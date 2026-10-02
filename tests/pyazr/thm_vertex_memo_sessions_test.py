#!/usr/bin/env python3
"""vertex=constant: the boundary B_c = S_c(E_lowest) follows the channel radius
across pyazr sessions and after set_channel_radius on a live session.

THMMatrixFunc keeps a per-thread memo of S_c at the lowest level of each J
group.  It used to be keyed by the PPair address, l and the level energy: the
memo outlives a CNuc, so a later session in the same process could get the
address back from the allocator and read S_c of another radius, and a radius
changed in place was never seen.  The model-averaging driver fits its radius
variants one after the other in one process; on the 18O(p,a) narrow-resonance
example a 6.1 fm variant read the 4.1 fm S_c (chi2 35.99 instead of 208.6 for
the same parameters).  Now the memo is keyed by the physical quantities.

On tests/18O_p_a_thm (two interfering 1/2+ levels, vertex=constant):
  1. chi2 of a 4.1 fm and a 6.1 fm copy, each computed in a fresh process,
     equal the in-process values in the orders 4.1, 6.1, 4.1 and 6.1, 4.1
     (rel 1e-10);
  2. set_channel_radius(1, 6.1) on a live 4.1 fm session gives the fresh
     6.1 fm chi2 (rel 1e-10), and back.

Needs the compiled engine; skips cleanly without it.
Run from anywhere:  python3 tests/pyazr/thm_vertex_memo_sessions_test.py
"""
import glob
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "18O_p_a_thm")
AZR = "18O_p_a_thm.azr"
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


_spec = importlib.util.spec_from_file_location("azrfile", os.path.join(ROOT, "pyazr", "azrfile.py"))
_azrfile = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_azrfile)
AzrModel = _azrfile.AzrModel

try:
    import numpy as np                                          # noqa: F401
    from pyazr import azure2
except Exception as err:                                       # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)

CHI2 = ("import sys, numpy as np; sys.path.insert(0, %r); from pyazr import azure2\n"
        "with azure2(sys.argv[1], cwd=sys.argv[2]) as m:\n"
        "    print('CHI2=' + repr(float(np.sum(np.asarray(m.residuals(np.asarray(m.params_rwa, float))) ** 2))))\n") % ROOT


def chi2(m):
    import numpy as np
    return float(np.sum(np.asarray(m.residuals(np.asarray(m.params_rwa, float))) ** 2))


with tempfile.TemporaryDirectory() as tmp:
    proj = {}
    for r in (4.1, 6.1):
        d = os.path.join(tmp, f"r{r}")
        shutil.copytree(SOURCE, d)
        for junk in ("output", "checks"):
            shutil.rmtree(os.path.join(d, junk), ignore_errors=True)
            os.makedirs(os.path.join(d, junk))
        for f in glob.glob(os.path.join(d, "run.log")):
            os.remove(f)
        m = AzrModel.from_file(os.path.join(d, AZR))
        m.set_channel_radius(1, r)
        m.write(os.path.join(d, AZR))
        proj[r] = d

    fresh = {}
    for r, d in proj.items():
        out = subprocess.run([sys.executable, "-c", CHI2, os.path.join(d, AZR), d],
                             capture_output=True, text=True, env=dict(os.environ))
        fresh[r] = float([ln for ln in out.stdout.splitlines() if ln.startswith("CHI2=")][-1][5:])
    print(f"1. fresh processes: chi2(4.1 fm) = {fresh[4.1]:.10g}, chi2(6.1 fm) = {fresh[6.1]:.10g}")
    check("the two radii differ (the test can see the memo)", rel(fresh[4.1], fresh[6.1]) > 1e-3,
          f"{fresh}")
    for order in ((4.1, 6.1, 4.1), (6.1, 4.1)):
        for i, r in enumerate(order):
            with azure2(os.path.join(proj[r], AZR), cwd=proj[r]) as m:
                c = chi2(m)
            check(f"order {order}, session {i + 1} ({r} fm) == fresh", rel(c, fresh[r]) < 1e-10,
                  f"{c!r} vs {fresh[r]!r}")

    print("2. set_channel_radius on a live session")
    with azure2(os.path.join(proj[4.1], AZR), cwd=proj[4.1]) as m:
        c0 = chi2(m)
        m.set_channel_radius(1, 6.1)
        c1 = chi2(m)
        m.set_channel_radius(1, 4.1)
        c2 = chi2(m)
    check("4.1 fm session == fresh", rel(c0, fresh[4.1]) < 1e-10, f"{c0!r}")
    check("-> 6.1 fm in place == fresh 6.1 fm", rel(c1, fresh[6.1]) < 1e-10, f"{c1!r} vs {fresh[6.1]!r}")
    check("-> back to 4.1 fm == fresh 4.1 fm", rel(c2, fresh[4.1]) < 1e-10, f"{c2!r}")

print("FAILED: " + ", ".join(failures) if failures else "all passed")
sys.exit(1 if failures else 0)
