#!/usr/bin/env python3
"""pyazr.endf, the checking codes (docs/upgrade/endf.md section 9): finding them, their decks, their reports read
into findings, a crash as a finding, a rewritten file as a candidate whose difference counts the digits lost and that
replaces nothing until accepted.

The real codes are not part of the repository: tests/pyazr/endf_standin/standin.py stands in for them (linked here
under their names; its report lines are those of the real programs, see its text).  When AZURE2_ENDF_UTILITIES and
AZURE2_PREPRO name directories of the built codes, the last part runs the real ones on the ENDF/B-VIII.1 extract and
a hand-written file (skipped otherwise).
Exits 77 (skipped) without the compiled engine module.
Run from anywhere:  python3 tests/pyazr/endf_check_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
try:
    import pyazr
    from pyazr import _azure2  # noqa: F401
    import pyazr.endf as E
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)
from endf_testfiles import Writer, material

STANDIN = os.path.join(HERE, "endf_standin", "standin.py")
EXTRACT = os.path.join(HERE, "data", "endf", "h-002_He_004_extract.endf")
REAL = {k: os.environ.get(k) for k in ("AZURE2_ENDF_UTILITIES", "AZURE2_PREPRO")}
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:600]}"))
    if not ok:
        failures.append(name)


TMP = tempfile.mkdtemp(prefix="endf_check_")
for k in REAL:
    os.environ.pop(k, None)
BIN = os.path.join(TMP, "bin")
os.makedirs(BIN)
for name in ("checkr", "fizcon", "PSYCHE", "stanef", "fixup", "dictin", "endf2c"):
    os.symlink(STANDIN, os.path.join(BIN, name))
os.symlink(STANDIN, os.path.join(BIN, "notachecker"))
x3 = [1e6, 2e6, 3e6, 4e6]
sections = {(3, 800): Writer(form="widest").cont(6012.0, 11.89691, 0, 0, 0, 0).tab1(
    -1.0e6, -1.0e6, 0, 0, (4,), (2,), x3, [0.0, 0.123456789, 0.234567891, 0.345678912])}
hand = material(os.path.join(TMP, "hand.endf"), sections, form="widest")
original = open(hand).read()

print("finding the codes (section 9.2 of the design)")
found = E.checkers(BIN)
check("checkers(dir): the known names, any case, nothing else, in the table's order",
      [c.name for c in found] == ["checkr", "fizcon", "psyche", "stanef", "fixup", "dictin", "endf2c"], found)
check("... with their kind and family", {(c.name, c.kind, c.family) for c in found} >= {
    ("checkr", "check", "utility"), ("stanef", "rewrite", "utility"), ("dictin", "rewrite", "prepro")})
os.environ["AZURE2_ENDF_UTILITIES"] = BIN
check("checkers(): $AZURE2_ENDF_UTILITIES", [c.name for c in E.checkers()][:3] == ["checkr", "fizcon", "psyche"])
os.environ["PATH"] = os.environ.get("PATH", "")  # PATH is searched last
check("checkers([]): an empty search finds nothing", E.checkers([]) == [])

print("running the checking codes: clean, findings, crash")
os.environ["ENDF_STANDIN_ROLE"] = "clean"
r = E.check(hand)
check("check(): the own check first, named as its own, and the three codes", r.own.path == hand and
      [x.program for x in r.reports] == ["checkr", "fizcon", "psyche"] and all(
          f.source == "pyazr.endf" for f in r.own.findings), r)
deck = open(os.path.join(r.reports[0].workdir, "deck.txt")).read().splitlines()
check("the deck: input file, message file, Y, DONE (checkr.f 255-270, 806-861)", deck == [
    "hand.endf", "checkr.out", "Y", "DONE"], deck)
check("each code in a directory of its own, on a copy", len({x.workdir for x in r.reports}) == 3 and
      open(hand).read() == original)
check("a clean report: no finding but the unsequenced note; version and exit read", [
    (f.severity, f.code) for f in r.reports[0].findings] == [("note", "CHECKR")] and
      r.reports[0].version == "8.24" and r.reports[0].exit_status == 0, r.reports[0])
os.environ["ENDF_STANDIN_ROLE"] = "findings"
r = E.check(hand, codes=["checkr", "psyche"])
ck = r.reports[0]
errors = [(f.mat, f.mf, f.mt, f.record) for f in ck.findings if f.severity == "error"]
check("CHECKR's blocks: MAT, MF, MT, the record number, the text", errors == [
    (9999, 1, 451, 5), (9999, 6, 2, 189), (9999, 6, 2, 189)] and ck.totals == (3, 0), ck)
check("... and the meaning of a known message", "Appendix G" in ck.findings[2].meaning, ck.findings[2])
check("the report's count (3 errors, 'SECTION CANNOT BE CHECKED' one of them) equals what was read", not any(
    "were read" in f.message for f in ck.findings), ck)
ps = r.reports[1]
check("PSYCHE: the comment under FILE 2 SECTION 151 a warning, the heading not a finding", [
    (f.severity, f.mf, f.mt) for f in ps.findings] == [("warning", 2, 151)], ps)
os.environ["ENDF_STANDIN_ROLE"] = "crash"
r = E.check(hand, codes=["checkr"])
check("a crash is a finding of the run, not an exception", r.reports[0].exit_status == 2 and any(
    "stopped early" in f.message for f in r.reports[0].findings), r.reports[0])
for codes, part in ((["stanef"], "not a checking code"), (["inter"], "not found")):
    try:
        E.check(hand, codes=codes)
        check(f"refused: {part}", False, "nothing raised")
    except pyazr.FormatRefused as e:
        check(f"refused with three parts: {part}", part in str(e) and e.allowed and e.fix, e)

print("rewriting codes: candidates (section 9.4 of the design)")
os.environ["ENDF_STANDIN_ROLE"] = "lossy"
c = E.standardize(hand)
deck = open(os.path.join(os.path.dirname(c.path), "deck.txt")).read().splitlines()
check("STANEF's deck: input, output, blank label, blank options, DONE (stanef.f 145-180, 536-612)",
      deck == ["hand.endf", "candidate.endf", "", "", "DONE"], deck)
check("the original is not touched", open(hand).read() == original)
d = {(x.mf, x.mt): x for x in c.diff()}
ch = d[(3, 800)]
check("diff: MF=3/MT=800 numbers changed, the digits lost counted (widest 9-10 digits -> 7)", ch.kind == "numbers"
      and ch.numbers_changed == 3 and ch.digits_lost >= 2 and 1e-9 < ch.max_relative_change < 1e-6, ch)
check("diff: MF=1/MT=451 bookkeeping only", d[(1, 451)].kind in ("unchanged", "bookkeeping"), d[(1, 451)])
check("summary names the section and the digits lost", "MF=3 MT=800: numbers" in c.summary(), c.summary())
c2 = E.standardize(hand, mf2=True)
check("STANEF with MF=2: the options line 0,Y,Y,Y", open(os.path.join(os.path.dirname(c2.path), "deck.txt"))
      .read().splitlines()[3] == "0,Y,Y,Y")
written = c.accept()
check("accept(): a new file beside the original", written == os.path.join(TMP, "hand.stanef.endf") and
      os.path.isfile(written) and open(hand).read() == original, written)
try:
    c.accept(to=hand)
    check("accept(to=original) without overwrite is refused", False, "nothing raised")
except pyazr.FormatRefused as e:
    check("accept(to=original) without overwrite is refused (three parts)", e.fix and "overwrite=True" in e.fix, e)
check("accept(to=original, overwrite=True) replaces it", c.accept(to=hand, overwrite=True) == hand and
      open(hand).read() != original)
open(hand, "w").write(original)
os.environ["ENDF_STANDIN_ROLE"] = "clean"
c = E.standardize(hand, program="dictin")
check("PREPRO: ENDFB.IN, DICTIN.INP in the working directory, ENDFB.OUT read back", open(os.path.join(
    os.path.dirname(c.path), "deck.txt")).read() == "ENDFB.IN\nENDFB.OUT\n" and all(
        x.kind == "unchanged" for x in c.diff()), c.summary())
os.environ["ENDF_STANDIN_ROLE"] = "crash"
try:
    E.standardize(hand, program="fixup")
    check("a rewriting code that writes nothing is refused", False, "nothing raised")
except pyazr.FormatRefused as e:
    check("a rewriting code that writes nothing is refused (three parts)", "wrote no file" in str(e) and e.fix, e)

if REAL["AZURE2_ENDF_UTILITIES"] and REAL["AZURE2_PREPRO"]:
    print("the real codes (AZURE2_ENDF_UTILITIES, AZURE2_PREPRO)")
    real = [REAL["AZURE2_ENDF_UTILITIES"], REAL["AZURE2_PREPRO"]]
    r = E.check(EXTRACT, search=real)
    by = {x.program: x for x in r.reports}
    ck = by["checkr"]
    check("CHECKR on the extract: exit 0, LRP=2 and the three sections the directory lists, all read",
          ck.exit_status == 0 and ck.totals == (4, 0) and sum(f.severity == "error" for f in ck.findings) == 4
          and sum("IS MISSING" in f.message for f in ck.findings) == 3, ck)
    check("the own check finds the same: LRP=2 and the three sections", sum(
        f.code == "DIRECTORY" for f in r.own.findings) == 3 and any(f.code == "LRP" for f in r.own.findings))
    check("FIZCON on the extract: one error read of one counted", by["fizcon"].totals == (1, 0) and sum(
        f.severity == "error" for f in by["fizcon"].findings) == 1, by["fizcon"])
    clean = material(os.path.join(TMP, "real_hand.endf"), sections, form="widest", sequence=False)
    r = E.check(clean, search=real)
    check("CHECKR on a hand-written file: no error", not [f for f in r.reports[0].findings if f.severity == "error"],
          r.reports[0])
    c = E.standardize(clean, search=real)
    ch = {(x.mf, x.mt): x for x in c.diff()}[(3, 800)]
    check("STANEF on a widest-form file: numbers rounded to 7 digits", ch.kind == "numbers" and ch.digits_lost >= 2
          and ch.max_relative_change < 1e-6, c.summary())
else:
    print("the real codes: skipped (AZURE2_ENDF_UTILITIES and AZURE2_PREPRO not set)")

if os.environ.get("ENDF_TEST_KEEP"):  # keep the files written, for a run of the real codes on them
    print(f"files kept in {TMP}")
else:
    shutil.rmtree(TMP, ignore_errors=True)
print(f"{len(failures)} failure(s)" + (": " + ", ".join(failures) if failures else ""))
sys.exit(1 if failures else 0)
