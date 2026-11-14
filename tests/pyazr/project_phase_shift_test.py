#!/usr/bin/env python3
"""A Project edit that leaves a phase-shift data set or extrapolation without its partial wave is refused.

The engine checks at a session's set-up that each phase shift's (J, l, s) is a channel of its entrance pair in an
active level (src/WorkingCopy.cpp); Project makes the same check after every edit, so that the edit is refused (with
its three parts, the project unchanged) rather than the next session.  On a copy of tests/upgrade/T3 (data set
aa_phase_d: J = 2, l = 2; extrapolations x5, x6, x7: J = 0, 2, 4):

  (a) both 2+ levels made inactive: refused, naming aa_phase_d; the project unchanged; the same edit written past
      Project's check is refused by the engine (a session of the data);
  (b) the only 4+ level made inactive: refused, naming extrapolation x7; written past the check, a session of the
      data sets up and a session of the extrapolations is refused by the engine;
  (c) one of the two 2+ levels removed: accepted, and both sessions set up (the other 2+ level carries the wave);
  (d) an edit of a project already without the wave (written past the check) that does not touch it is accepted.

Run from anywhere:  python3 tests/pyazr/project_phase_shift_test.py
"""
import logging
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import pyazr  # noqa: E402
    from pyazr import nuclear_part as npart  # noqa: E402
except Exception as err:
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

warnings.simplefilter("ignore")
logging.getLogger("pyazr.engine").setLevel(logging.CRITICAL)


def inactive_past_the_check(path, labels, out):
    """``path`` with the levels ``labels`` made inactive by the engine's edit itself (as a hand edit of the file
    would be), not through Project's checks; saved as ``out``."""
    p = pyazr.Project.open(path)
    for label in labels:
        npart._replace(p._doc, npart._engine_call("level", npart._engine_module().set_level_edits, p._doc,
                                                  npart._level_address(label), energy=None, fixed=None,
                                                  active=False))
    return p.save(out)


def sets_up(path, extrapolations):
    try:
        pyazr.Session(path, extrapolations=extrapolations).close()
        return True
    except pyazr.Error:
        return False


tmp = tempfile.mkdtemp(prefix="project_phase_shift_")
try:
    shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T3"), os.path.join(tmp, "T3"),
                    ignore=shutil.ignore_patterns("output", "checks"))
    T3 = os.path.join(tmp, "T3", "T3.azr")
    twos = [lv.label for lv in pyazr.Project.open(T3).levels() if lv.jpi == "2+"]
    check("T3 has two 2+ levels and one 4+", len(twos) == 2 and
          [lv.label for lv in pyazr.Project.open(T3).levels() if lv.jpi == "4+"] == ["4+ 11.35"], twos)

    print("(a) both 2+ levels inactive")
    p = pyazr.Project.open(T3)
    text = p.text
    p.set_level(twos[0], active=False)
    try:
        p.set_level(twos[1], active=False)
        check("refused", False, "accepted")
    except pyazr.FormatError as e:
        parts = (e.what, e.allowed, e.fix)
        check("refused with what, allowed, fix in the message", all(parts) and all(x in str(e) for x in parts), parts)
        check("it names the data set and the wave", "dataset aa_phase_d" in e.what and "J = 2, l = 2" in e.what,
              e.what)
    check("the project is as before the refused edit", "2+ 3.03" not in [lv.label for lv in p.levels() if lv.active]
          and twos[1] in [lv.label for lv in p.levels() if lv.active], [lv.label for lv in p.levels() if lv.active])
    path = inactive_past_the_check(T3, twos, os.path.join(tmp, "T3", "a.azr"))
    check("the engine refuses that project too (a session of the data)", not sets_up(path, False), path)

    print("(b) the 4+ level inactive")
    p = pyazr.Project.open(T3)
    try:
        p.set_level("4+ 11.35", active=False)
        check("refused", False, "accepted")
    except pyazr.FormatError as e:
        check("refused, naming extrapolation x7", "extrapolation x7" in e.what and "J = 4, l = 4" in e.what, e.what)
    check("the project is unchanged", p.text == text, "")
    path = inactive_past_the_check(T3, ["4+ 11.35"], os.path.join(tmp, "T3", "b.azr"))
    check("the engine: a session of the data sets up, one of the extrapolations is refused",
          sets_up(path, False) and not sets_up(path, True), path)

    print("(c) one 2+ level removed")
    p = pyazr.Project.open(T3)
    p.remove_level(twos[1])
    path = p.save(os.path.join(tmp, "T3", "c.azr"))
    check("accepted, and both sessions set up", sets_up(path, False) and sets_up(path, True), path)

    print("(d) an edit elsewhere of a project already without the wave")
    path = inactive_past_the_check(T3, twos, os.path.join(tmp, "T3", "d.azr"))
    p = pyazr.Project.open(path)
    try:
        p.set_level("1+ 17.64", energy=17.65)
        check("accepted (the problem was there before)", True)
    except pyazr.FormatError as e:
        check("accepted (the problem was there before)", False, str(e))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"{'FAILED: ' + str(len(failures)) if failures else 'all ok'}")
sys.exit(1 if failures else 0)
