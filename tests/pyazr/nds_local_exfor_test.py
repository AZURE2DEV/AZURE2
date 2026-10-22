#!/usr/bin/env python3
"""fetch_exfor and reference read a local copy of the EXFOR master files when the Web-API fails.

The IAEA EXFOR services are intermittently behind a Cloudflare challenge (on the Notre Dame CRC cluster
on 2026-10-22 x4get and the dataexplorer API both answered every request with one).  The EXFOR library
itself is public as https://github.com/IAEA-NDS/exfor_master (one exforall/<3 chars>/<entry>.x4 file per
entry), so ``fetch_exfor``/``reference`` fall back to a clone of it ($AZURE2_EXFOR_MASTER or
~/.azure2/exfor_master), and ``source="local"`` reads only that.

Offline: tests/pyazr/data/exfor_master holds one synthetic entry, T0001, written for this test in the
EXFOR exchange format with every feature the reader has to handle:
  002  a 7-column table wrapped onto two records per row; an angle in the subentry's COMMON and a
       per-cent error in subentry 001's COMMON, both repeated on every row; Fortran exponents without
       E ("1.23-4"), blanks inside a number (" 2.5   -04"), a blank field; a REACTION code split over
       two lines and followed by free text with parentheses;
  003  two pointers, the second a combination of reactions that wraps onto a line starting with "(".
Checks: columns, units and values; reaction, first author and year; the pointer rules (none given ->
error naming them, ".1" and the bare "1" suffix); the AZURE2 columns; reference(); a 403, a URLError and
a challenge page each trigger the fallback with a warning; source="web" re-raises; no local copy ->
re-raise; plus=1 and bad IDs are refused.

The reader was also checked against the IAEA dataexplorer API's own parse of the same entries
(13C+a/10-2-26_n16O_exfor/raw in the R-matrix archive, 439 entries saved 2026-10-02): 1,596 of 1,596
data tables and 1,663 of 1,663 reaction codes identical.

Run from anywhere:  python3 tests/pyazr/nds_local_exfor_test.py
"""
import importlib.util
import io
import math
import os
import sys
import urllib.error
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
os.environ["AZURE2_EXFOR_MASTER"] = os.path.join(HERE, "data", "exfor_master")
# load pyazr/nds.py by path: it needs numpy only, not the compiled engine
spec = importlib.util.spec_from_file_location("pyazr_nds_standalone", os.path.join(ROOT, "pyazr", "nds.py"))
N = importlib.util.module_from_spec(spec)
spec.loader.exec_module(N)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def same(a, b):
    return len(a) == len(b) and all((math.isnan(x) and math.isnan(y)) or abs(x - y) <= 1e-12 * max(1.0, abs(y))
                                    for x, y in zip(a, b))


def raises(fn, exc, text=""):
    try:
        fn()
    except exc as e:
        return text in str(e)
    return False


print("1. the local reader")
check("the copy is found through $AZURE2_EXFOR_MASTER", N.exfor_master_path() == os.environ["AZURE2_EXFOR_MASTER"],
      N.exfor_master_path())
d = N.fetch_exfor("T0001002", source="local")
cols = [(c.name, c.unit) for c in d.columns]
check("columns and units: the 7-column table, then the COMMON values (subentry, then 001)",
      cols == [("EN", "KEV"), ("EN-ERR", "KEV"), ("DATA", "B*KEV"), ("ERR-S", "B*KEV"), ("MISC1", "NO-DIM"),
               ("MISC2", "NO-DIM"), ("MISC3", "NO-DIM"), ("ANG-CM", "ADEG"), ("ERR-1", "PER-CENT")], cols)
check("Fortran exponents and blanks inside a number", same(list(d.arrays["DATA"]), [1.23e-4, 2.5e-4, 3.0e-4]),
      list(d.arrays["DATA"]))
check("a blank field is NaN", same(list(d.arrays["ERR-S"]), [1.0e-5, 2.0e-5, math.nan]), list(d.arrays["ERR-S"]))
check("the wrapped seventh column", same(list(d.arrays["MISC3"]), [3.0, 6.0, 9.0]), list(d.arrays["MISC3"]))
check("COMMON values repeat on every row", same(list(d.arrays["ANG-CM"]), [90.0] * 3)
      and same(list(d.arrays["ERR-1"]), [3.0] * 3))
check("the two-line REACTION code, free text dropped", d.reaction == "6-C-12(A,G)8-O-16,,SIG,,SFC", d.reaction)
check("first author and year from subentry 001", (d.author, d.year) == ("A.Author", 2026), (d.author, d.year))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    rows = N._azr_columns(d)
check("AZURE2 columns: three finite rows at 90 degrees, the blank error given the default",
      rows.shape == (3, 4) and all(math.isfinite(v) for v in rows.ravel()) and all(rows[:, 1] == 90.0), rows)

print("\n2. pointers")
check("no pointer given -> error naming them", raises(lambda: N.fetch_exfor("T0001003", source="local"),
                                                      ValueError, "pointers 1, 2"))
p1 = N.fetch_exfor("T0001003.1", source="local")
check("pointer 1: its reaction, its DATA column and subentry 001's COMMON",
      p1.reaction == "6-C-12(P,G)7-N-13,,SIG" and same(list(p1.arrays["DATA"]), [5.0, 7.0])
      and [c.name for c in p1.columns] == ["EN", "DATA", "ERR-1"], (p1.reaction, [c.name for c in p1.columns]))
p2 = N.fetch_exfor("T00010032", source="local")
check("pointer 2 (bare suffix): the combination that wraps onto a '(' line",
      p2.reaction == "(6-C-12(P,X)1-H-1,,SIG)+(6-C-12(P,X)1-H-2,,SIG)" and same(list(p2.arrays["DATA"]), [6.0, 8.0]),
      p2.reaction)

print("\n3. reference")
r = N.reference("T0001002", source="local")
check("reference from the local BIB", (r.entry, r.x4ref, r.year) == ("T0001", "J,PR/C,99,012345,2026", 2026),
      (r.entry, r.x4ref, r.year))

print("\n4. the fallback")
CHALLENGE = "<!DOCTYPE html><html><head><title>Just a moment...</title></head><body></body></html>"
real_get = N._http_get


def serve(answer):
    def fake(url, params=None, timeout=60):
        if isinstance(answer, Exception):
            raise answer
        return answer
    N._http_get = fake


for label, answer in (("HTTP 403", urllib.error.HTTPError("u", 403, "Forbidden", {}, io.BytesIO(b""))),
                      ("URLError", urllib.error.URLError("timed out")),
                      ("challenge page", CHALLENGE)):
    serve(answer)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        got = N.fetch_exfor("T0001002")
    check(f"{label}: fetch_exfor reads the local copy, with a warning",
          got.reaction == d.reaction and any("local EXFOR master" in str(x.message) for x in w), [str(x.message) for x in w])
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        check(f"{label}: reference too", N.reference("T0001").x4ref == "J,PR/C,99,012345,2026")
serve(CHALLENGE)
check("source='web' re-raises", raises(lambda: N.fetch_exfor("T0001002", source="web"), N._WebUnavailable, "Just a moment"))
real_path = N.exfor_master_path
N.exfor_master_path = lambda: None
check("no local copy: the web error is raised", raises(lambda: N.fetch_exfor("T0001002"), N._WebUnavailable))
N.exfor_master_path = real_path
N._http_get = real_get

print("\n5. refusals")
check("plus=1 is not available locally", raises(lambda: N.fetch_exfor("T0001002", plus=1, source="local"),
                                                ValueError, "original units only"))
check("a missing subentry", raises(lambda: N.fetch_exfor("T0001009", source="local"), ValueError, "no subentry 009"))
check("a missing entry", raises(lambda: N.fetch_exfor("T9999002", source="local"), ValueError, "not in the local copy"))
check("a malformed ID", raises(lambda: N.fetch_exfor("bad", source="local"), ValueError, "not an EXFOR dataset ID"))
check("a bad source", raises(lambda: N.fetch_exfor("T0001002", source="ftp"), ValueError, "source must be"))

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): " + ", ".join(failures))
    sys.exit(1)
print("all checks passed")
