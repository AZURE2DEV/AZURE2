#!/usr/bin/env python3
"""Pole guard (rmfit.poleguard, BoundsPolicy.pole_guard; 2026-11-27, the evaluator's choice).  Fixture
fixtures_p3/relvector_p12c/ngcov.azr (12C+p, format 3).  Checks:
  (a) narrow_poles finds every narrow pole of the model ON a level of its J^pi: the guard reports nothing;
  (b) polish with the guard on equals polish without it when the guard never fires (same objective, same vector);
  (c) a guard that flags every trial point rejects every step: polish returns its start, and says so in the log."""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if not _PYAZR3:
    print("SKIP: a pyazr 3 test (Session.jpi_block)")
    sys.exit(77)
from _fixtures import fixture  # noqa: E402
import numpy as np
from rmfit.engine import SingleSession
from rmfit.local import BoundsPolicy, polish
from rmfit import poleguard
failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)
def main():
    src = fixture("relvector_p12c")
    tmp = tempfile.mkdtemp(prefix="p3poleguard_")
    try:
        os.symlink(os.path.join(src, "data"), os.path.join(tmp, "data"))
        for d in ("output", "checks"):
            os.makedirs(os.path.join(tmp, d))
        shutil.copy(os.path.join(src, "ngcov.azr"), os.path.join(tmp, "ngcov.azr"))
        azr = os.path.join(tmp, "ngcov.azr")
        cfg = {"step": 0.002, "gamma_max": 0.01}
        with SingleSession(azr, cwd=tmp, flags=dict(use_brune=True), name="poleguard") as ev:
            x0 = np.array(ev.x0, float)
            g = poleguard.PoleGuard(ev, cfg, log=print)
            allp = [(g.groups[jg], poleguard.narrow_poles(g.sess, jg, g.e_min, g.e_max, g.vector(x0), step=0.002, gmax=0.01))
                    for jg in sorted(g.groups)]
            print("  narrow poles (Gamma < 10 keV) by J^pi:", [(j, [(round(e, 4), round(w * 1e3, 3)) for e, w in p]) for j, p in allp if p])
            bad = g(x0)
            check("(a) no spurious pole on the fixture's physical model", bad == [], poleguard.PoleGuard.text(bad))
            check("(a) the guard maps every level parameter", len(g.map) >= sum(1 for k in ev.keys if k.kind in ("energy", "width")))
            g.close()
            f0 = polish(ev, x0, policy=BoundsPolicy(), max_nfev=6)
            lines = []
            f1 = polish(ev, x0, policy=BoundsPolicy(pole_guard=cfg), max_nfev=6, log=lines.append)
            check("(b) guard on = guard off when it never fires", abs(f0.objective - f1.objective) < 1e-9 * max(1, abs(f0.objective))
                  and np.allclose(f0.x, f1.x), (f0.objective, f1.objective))
            check("(b) the log has the guard's start and end lines", any("pole guard at the start: none" in l for l in lines)
                  and any("pole guard at the end: none" in l for l in lines), lines)
            orig = poleguard.PoleGuard.__call__
            def flag_all(self, x):
                return [] if np.array_equal(np.asarray(x, float), x0) else [("test", 5.0, 1e-6, 9.0)]
            poleguard.PoleGuard.__call__ = flag_all
            try:
                lines = []
                f2 = polish(ev, x0, policy=BoundsPolicy(pole_guard=cfg), max_nfev=6, log=lines.append)
            finally:
                poleguard.PoleGuard.__call__ = orig
            check("(c) every step rejected: polish keeps the start", np.allclose(f2.x, x0) and abs(f2.objective - f2.start_objective) < 1e-9,
                  (f2.objective, f2.start_objective, f2.status))
            check("(c) the rejections are logged", any("pole guard: trial step rejected" in l for l in lines), lines[-3:])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("ALL PASSED" if not failures else f"FAILED: {failures}")
    return 1 if failures else 0
if __name__ == "__main__":
    sys.exit(main())
