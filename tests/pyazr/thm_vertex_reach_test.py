#!/usr/bin/env python3
"""session.thm_vertex says where the spectator window does not reach E.

Beyond the kinematic reach of a ``ps`` window (no accepted spectator
direction at E) the engine uses the nodes of the nearest data point -- by
design, for the folding sub-points next to the data (EPoint.cpp,
ThmSpectatorWindow::Nodes).  thm_vertex returned those rows without a word,
so a vertex scan out there showed identical rows that read as physics (the
7Li cross-check, 1.88-4 MeV).  Now, on tests/18O_p_a_thm, 2H(18O,a15N)n at
54 MeV with a Hulthen window 0-40 MeV/c:

  1. ``reached`` is True at the data energies, False at 3 and 5 MeV; a
     UserWarning names the unreached energies; ``strict=True`` raises
     ValueError instead; the unreached rows are those of the last data point
     (the documented fallback).
  2. Energies all reached: no warning.
  3. A delta experiment (no window): ``reached`` all True, no warning.
  4. vertexModel=dw (Coulomb waves, spectatorAngles=cm:0-180): ``reached``
     False outside the DW vertex grid, with a warning saying so.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/thm_vertex_reach_test.py
"""
import glob
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "18O_p_a_thm")
AZR = "18O_p_a_thm.azr"
KIN = dict(beam="18O", target="d", spectator="n", Ebeam=54)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
try:
    import numpy as np
    from pyazr import azure2
    from pyazr.azrfile import AzrModel
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def fresh_copy(dst):
    shutil.copytree(SOURCE, dst)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(dst, junk), ignore_errors=True)
        os.makedirs(os.path.join(dst, junk))
    for f in glob.glob(os.path.join(dst, "run.log")):
        os.remove(f)


def vertex(s, energies, **kw):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r = s.thm_vertex("A", np.asarray(energies, float), **kw)
    return r, [str(x.message) for x in w if issubclass(x.category, UserWarning)]


with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "w")
    fresh_copy(work)
    src = os.path.join(work, AZR)

    print("1-2. Hulthen window 0-40 MeV/c")
    m = AzrModel.from_file(src)
    m.set_thm_experiment("A", [1, 2], ps="hulthen:0-40", psNodes=8, **KIN)
    m.write(os.path.join(work, "win.azr"))
    with azure2(os.path.join(work, "win.azr"), cwd=work) as s:
        data_e = np.sort(np.concatenate([np.asarray(e, float) for e in s.calculate_energies()]))
        E = [0.55, 0.7, 0.9, 3.0, 5.0]
        r, w = vertex(s, E)
        reached = np.asarray(r["reached"], bool)
        check("reached at the data energies, not at 3 and 5 MeV",
              reached.tolist() == [True, True, True, False, False], reached)
        check("a UserWarning names them", len(w) == 1 and "2 of 5" in w[0]
              and "3, 5" in w[0], w)
        last = s.thm_vertex("A", np.array([data_e[-1]]))
        same = all(np.array_equal(r["p_s"][k], last["p_s"][0]) for k in (3, 4))
        check("the unreached rows are the last data point's (the fallback)", same,
              (r["p_s"][3][:3], last["p_s"][0][:3]))
        try:
            s.thm_vertex("A", np.array(E), strict=True)
            check("strict=True raises", False, "no error")
        except ValueError as err:
            check("strict=True raises", "beyond the reach" in str(err), err)
        r, w = vertex(s, E[:3])
        check("all reached: no warning", np.all(r["reached"]) and not w, (r["reached"], w))

    print("3. delta (no window)")
    m = AzrModel.from_file(src)
    m.set_thm_experiment("A", [1, 2], **KIN)
    m.write(os.path.join(work, "delta.azr"))
    with azure2(os.path.join(work, "delta.azr"), cwd=work) as s:
        r, w = vertex(s, [0.55, 3.0, 5.0])
        check("reached everywhere, no warning", np.all(r["reached"]) and not w, (r["reached"], w))

    print("4. vertexModel=dw")
    m = AzrModel.from_file(src)
    m.set_thm_experiment("A", [1, 2], distortion="coulomb", vertexModel="dw",
                         spectatorAngles="cm:0-180", **KIN)
    m.write(os.path.join(work, "dw.azr"))
    with azure2(os.path.join(work, "dw.azr"), cwd=work) as s:
        r, w = vertex(s, [0.55, 0.9, 3.0])
        check("False outside the DW vertex grid only",
              np.asarray(r["reached"], bool).tolist() == [True, True, False], r["reached"])
        check("with a warning naming the grid", len(w) == 1 and "DW vertex grid" in w[0], w)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all thm_vertex reach checks passed")
