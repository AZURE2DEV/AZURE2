#!/usr/bin/env python3
"""pyazr.tabulate on a project with THM experiments.

tabulate() replaces the data segments of the project by its own (on-shell)
tabulation segments.  The <thm> block keys experiment[...] segments= and
weight[k] by data-segment line, so left in place they named the tabulation
segments, and the engine refused the file ("segment 1 is not a THM segment").

On a copy of tests/18O_p_a_thm with its two THM segments made one experiment
and a weight on segment 1: the tabulation runs, gives finite positive cross
sections at the requested points, and equals the tabulation of the same
project without the experiment and the weight.

Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/tabulate_thm_test.py
"""
import glob
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "18O_p_a_thm")
AZR = "18O_p_a_thm.azr"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
try:
    import numpy                                               # noqa: F401
    from pyazr.azrfile import AzrModel
    from pyazr.tabulate import tabulate
    import pyazr._azure2                                       # noqa: F401
except Exception as err:                                       # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def fresh_copy(dst):
    shutil.copytree(SOURCE, dst)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(dst, junk), ignore_errors=True)
        os.makedirs(os.path.join(dst, junk))
    for f in glob.glob(os.path.join(dst, "run.log")):
        os.remove(f)


with tempfile.TemporaryDirectory() as tmp:
    plain = os.path.join(tmp, "plain")
    thm = os.path.join(tmp, "thm")
    fresh_copy(plain)
    fresh_copy(thm)
    with open(os.path.join(thm, "w.dat"), "w") as fh:
        fh.write("0.4 1.0\n1.0 2.0\n")
    m = AzrModel.from_file(os.path.join(thm, AZR))
    m.set_thm_experiment("E", [1, 2], background="linear")
    m.set_thm_weight(1, "w.dat")
    m.write(os.path.join(thm, AZR))

    kw = dict(rel_tol=5e-2, min_points=9, max_points=9, max_rounds=0)
    try:
        a = tabulate(os.path.join(thm, AZR), [(1, 2)], 0.55, 0.9, **kw)[(1, 2)]
        check("a project with experiment[...] and weight[k] tabulates", True)
    except Exception as err:
        check("a project with experiment[...] and weight[k] tabulates", False,
              f"{type(err).__name__}: {err}")
        a = None
    b = tabulate(os.path.join(plain, AZR), [(1, 2)], 0.55, 0.9, **kw)[(1, 2)]
    if a is not None:
        check("9 finite positive cross sections",
              len(a.sigma) == 9 and all(math.isfinite(s) and s > 0 for s in a.sigma), list(a.sigma))
        check("the same as without the THM keys (the tabulation is on shell)",
              list(a.e_cm) == list(b.e_cm) and list(a.sigma) == list(b.sigma))

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all tabulate checks passed")
