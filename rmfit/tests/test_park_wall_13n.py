#!/usr/bin/env python3
"""Park's J > 0 wall as residual rows (engine.ParkWall, 2026-10-24) on the 13N model.

In Park mode (flags use_park) the engine adds sum (J/1e-3)^2 over levels with J = 1 - sum gamma^2 dS/dE < 0 to its own
scalar chi-squared; rmfit's objective is built from residual rows and used to miss it.  Checks:
  * Brune mode has no wall rows; Park mode one row per level, labels matching park_norms;
  * at the fitted 13N point (all J > 0) the rows are zero and the objective equals the engine's chi2 + penalties;
  * a level's widths scaled past the wall: J follows 1 - f^2 (1 - J0) (Park amplitudes enter J quadratically),
    the row is J/1e-3, the objective equals the engine's chi2 (which contains the wall) + penalties, and the score
    reports the wall in ``penalty``;
  * the wall's Jacobian rows match central differences of the rows;
  * a Sharded evaluator gives the same objective and wall;
  * a polish of that level's widths walks back to J >= 0 (to the wall's tolerance).
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3 (its pyazr 3 version is test_p3_park_wall_13n.py)")
    sys.exit(77)
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.engine import SingleSession, Sharded
from rmfit.keys import take
from rmfit.local import Bounds, polish
from rmfit.model import single_variant
from pyazr import AzrModel

failures = []
# the 13N fixture: from AZURE2_ROOT, or ~/AZURE2 when AZURE2_ROOT is a frozen build copy without tests/
FIX = next(d for d in (AZURE2_ROOT, os.path.expanduser("~/AZURE2")) if os.path.isdir(os.path.join(d, "tests", "13N")))
PARK = dict(use_park=True, use_brune=True)


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_park_")
    shutil.copytree(os.path.join(FIX, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(FIX, "tests", "13N", "13N.azr"))
    mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))

    with SingleSession(single_variant(base, "brune", tmp), cwd=tmp, name="brune") as b:
        check("Brune mode: no wall rows", b.park is None)
        nb = b.residuals(b.x0)[0].size
        check("Brune mode: rows = data + penalties", nb == b.npts + b.penalties.cols.size, f"{nb}")

    ev = SingleSession(single_variant(base, "park", tmp), cwd=tmp, flags=PARK, name="park")
    try:
        pw = ev.park
        check("Park mode: a wall", pw is not None)
        x0 = np.array(ev.x0, float)
        J0 = ev.park_norms(x0)
        L = len(pw.labels)
        print(f"  {L} levels, J at the fitted point: min {J0.min():.4g}")
        check("one label per park_norms entry", J0.size == L, f"{J0.size} vs {L}")
        check("every free R-matrix column belongs to a labelled level",
              sorted(c for cs in pw.cols for c in cs) == list(range(ev.n_rmatrix)),
              f"{sum(len(c) for c in pw.cols)} of {ev.n_rmatrix}")
        r0, _ = ev.residuals(x0)
        check("rows = data + penalties + one per level", r0.size == ev.npts + ev.penalties.cols.size + L, f"{r0.size}")
        check("all J > 0 at the fit: wall rows zero", np.all(J0 > 0) and np.all(r0[-L:] == 0), f"{J0}")
        eng0 = float(ev.worker.call("eval", x0)) + ev.penalties.value(x0)
        o0 = ev.objective(x0)
        check("objective = engine chi2 + penalties (no wall)", np.isclose(o0, eng0, rtol=1e-9), f"{o0} vs {eng0}")

        # the level with the most free width amplitudes, scaled past the wall
        wcols = [[c for c in cs if ev.keys[c].kind == "width"] for cs in pw.cols]
        l = int(np.argmax([len(c) for c in wcols]))
        f = np.sqrt(1.5 / (1.0 - J0[l]))
        x1 = x0.copy()
        x1[wcols[l]] *= f
        J1 = ev.park_norms(x1)
        print(f"  level {pw.labels[l]}: {len(wcols[l])} widths x{f:.3f}, J {J0[l]:.4f} -> {J1[l]:.4f}")
        check("J is quadratic in the Park amplitudes", np.isclose(J1[l], 1 - f * f * (1 - J0[l]), rtol=1e-6, atol=1e-9),
              f"{J1[l]} vs {1 - f * f * (1 - J0[l])}")
        check("the other levels' J unchanged", np.allclose(np.delete(J1, l), np.delete(J0, l), rtol=1e-12))
        r1, Jm = ev.residuals(x1, jac=True)
        check("its wall row is J/1e-3", np.isclose(r1[-L + l], J1[l] / 1e-3, rtol=1e-12), f"{r1[-L + l]}")
        wall = float(np.sum((J1[J1 < 0] / 1e-3) ** 2))
        eng1 = float(ev.worker.call("eval", x1)) + ev.penalties.value(x1)   # the engine's chi2 contains the wall
        o1 = ev.objective(x1)
        check("objective = engine chi2 (with its wall) + penalties", np.isclose(o1, eng1, rtol=1e-9),
              f"{o1} vs {eng1} (wall {wall})")
        s1 = ev.score(x1)
        check("score: wall counted in penalty", np.isclose(s1.penalty, ev.penalties.value(x1) + wall, rtol=1e-9),
              f"{s1.penalty} vs {ev.penalties.value(x1) + wall}")
        check("score objective = objective", np.isclose(s1.objective, o1, rtol=1e-12), f"{s1.objective} vs {o1}")

        # Jacobian rows vs central differences of the rows
        worst = 0.0
        for c in pw.cols[l]:
            h = 1e-6 if ev.keys[c].kind == "energy" else 1e-5 * max(abs(x1[c]), 1e-3)
            xp, xm = x1.copy(), x1.copy()
            xp[c] += h
            xm[c] -= h
            fd = (ev.residuals(xp)[0][-L + l] - ev.residuals(xm)[0][-L + l]) / (2 * h)
            an = Jm[-L + l, c]
            worst = max(worst, abs(an - fd) / max(abs(fd), 1e-6))
        check("wall Jacobian rows = finite differences (rel 1e-4)", worst < 1e-4, f"worst {worst:.2e}")
        check("wall Jacobian zero off the level", np.count_nonzero(np.delete(Jm[-L + l], pw.cols[l])) == 0)

        sh = Sharded(base, os.path.join(tmp, "shards"), nshards=2, flags=PARK)
        try:
            xs = take(x1, ev.keys, sh.keys, sh.x0)
            os_ = sh.objective(xs)
            check("Sharded objective = single (wall once, not per shard)", np.isclose(os_, o1, rtol=1e-9), f"{os_} vs {o1}")
            rs, _ = sh.residuals(xs)
            check("Sharded wall rows = single", np.allclose(rs[-L:], r1[-L:], rtol=1e-12))
        finally:
            sh.close()

        free = np.zeros(x1.size, bool)
        free[wcols[l]] = True
        bnd = Bounds(np.full(x1.size, -np.inf), np.full(x1.size, np.inf), np.zeros(x1.size, bool))
        res = polish(ev, x1, free=free, bounds=bnd, max_nfev=40)
        J2 = ev.park_norms(res.x)
        print(f"  polish {res.start_objective:.1f} -> {res.objective:.1f} ({res.status}); J {J1[l]:.4f} -> {J2[l]:.5f}")
        r2, _ = ev.residuals(res.x)
        check("a polish walks the level back across the wall (J >= 0 to its tolerance, wall share ~0)",
              J2[l] > -2e-3 and ev.park_value(r2) < 1.0 and res.objective < 0.1 * o1,
              f"J {J2[l]}, wall {ev.park_value(r2)}, objective {res.objective} vs start {o1}")
    finally:
        ev.close()
        shutil.rmtree(tmp, ignore_errors=True)
    print("FAILED: " + ", ".join(failures) if failures else "all passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
