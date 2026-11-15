#!/usr/bin/env python3
"""pyazr.endf.check_resonances (docs/upgrade/endf.md section 13.14) on the committed ENDF/B-VIII.1 3He + 4He extract
(tests/pyazr/data/endf/h-002_He_004_extract.endf): the design's three findings (section 6.2) and nothing else.

(a) reconstruction: MF=2 reconstructed by the engine reproduces MF=3/MT=600 to 1.2e-6 above 7.2 MeV; the two first
    points are off (a RECONSTRUCTION warning naming 2 points), the first 37 orders of magnitude above MF=2's
    (a THRESHOLD warning);
(b) masses: the 6Li mass of the file 21.7 keV below AME2020 (MASS), the file's Q of p + 6Li within 0.1 keV of
    AME2020's (a Q note), 22.6 keV from what the file's own masses give (Q_MASSES);
(c) text: the MF=1 text names LTP=2 while the extract has no MF=6/MT=2 (TEXT; on the full library file the record
    holds LTP=1, endf.md 13.14);
(d) check(path, resonances=True) puts the same findings among the own findings (what the page's Check screen lists).
Exits 77 (skipped) without the compiled engine module.  Run from anywhere.
"""
import os
import sys
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    import pyazr  # noqa: F401
    from pyazr import _azure2  # noqa: F401
    import pyazr.endf as E
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)

warnings.simplefilter("ignore")
EXTRACT = os.path.join(HERE, "data", "endf", "h-002_He_004_extract.endf")
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:600]}"))
    if not ok:
        failures.append(name)


r = E.check_resonances(EXTRACT)
print(r)
codes = sorted((f.severity, f.code) for f in r.findings)
check("the findings: RECONSTRUCTION, THRESHOLD, MASS, Q_MASSES, TEXT warnings and the Q note, nothing else",
      codes == sorted([("warning", "RECONSTRUCTION"), ("warning", "THRESHOLD"), ("warning", "MASS"),
                       ("note", "Q"), ("warning", "Q_MASSES"), ("warning", "TEXT")]), codes)
rec = r.parts["reconstruction"][0]
check("(a) MF=3/MT=600: 97 energies, 2 off, median below 1e-6", rec["mt"] == 600 and rec["points"] == 97 and
      rec["off"] == 2 and rec["median"] < 1e-6, rec)
th = next(f for f in r.findings if f.code == "THRESHOLD")
check("(a) the first point 37 orders of magnitude above MF=2", "37 orders" in th.message and "7.05e+06 eV" in
      th.message, th.message)
li6 = next(x for x in r.parts["masses"] if x.get("nuclide") == "6Li")
q = next(x for x in r.parts["masses"] if x.get("mt") == 600)
check("(b) 6Li 21.7 keV below AME2020; Q of MT=600 within 0.1 keV of AME2020's, 22.6 keV from the file's masses'",
      abs(li6["difference_keV"] + 21.7) < 0.1 and abs(q["q_file_eV"] - q["q_ame2020_eV"]) < 100 and
      abs((q["q_file_eV"] - q["q_file_masses_eV"]) / 1e3 + 22.6) < 0.1, (li6, q))
text = next(f for f in r.findings if f.code == "TEXT")
check("(c) the text's LTP=2 against a file without MF=6/MT=2 (listed in its directory)",
      "LTP=2" in text.message and "not in the file" in text.message and "directory lists it" in text.message,
      text.message)
result = E.check(EXTRACT, resonances=True)
own = [(f.code, f.source) for f in result.own.findings if f.source == "pyazr.endf resonances"]
check("(d) check(..., resonances=True): the six findings among the own findings",
      sorted(c for c, _ in own) == sorted(c for _, c in codes), own)
print(f"{len(failures)} failure(s)" + (": " + ", ".join(failures) if failures else ""))
sys.exit(1 if failures else 0)
