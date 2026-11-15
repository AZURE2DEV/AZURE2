#!/usr/bin/env python3
"""BoundsPolicy.total_width_floor_windows (2026-10-12, the evaluator: "the floor should apply to the total width, not just the n0
width.  A state could have a large alpha width instead and a small neutron width and still be fine ... as long as at least one
of them is large enough to satisfy the floor").  A level whose seed energy lies in a window keeps
sum_c Gamma_W,c theta_c^2 >= gamma_min through a one-sided penalty row in polish; single amplitudes stay free (no sign lock).
Stubbed evaluator, no engine.  Run:  python3 rmfit/tests/test_total_width_floor.py"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.keys import ParamKey  # noqa: E402
from rmfit.local import BoundsPolicy, default_bounds, polish, floor_rows  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


GW_AMP = {0: 2.0, 1: 1.0, 2: 0.5, 3: 2.0}            # gamma_W (amplitude units)
GW_WID = {0: 4.0e6, 1: 1.0e6, 2: 0.0, 3: 4.0e6}       # Gamma_W (eV); channel 2 closed at the level's energy


class Pen:
    def value(self, x):
        return 0.0


class Stub:
    """New 3/2- at 10.30 with n0 (pair 1), a0 (pair 2) and a closed n1 (pair 3); an old 3/2- at 10.27 (pair 1).
    Data residuals: each open channel of the new level is pulled toward a target partial width (eV), with weight 1 per
    0.5 % of 100 keV, so the data alone decide the partition and the floor decides only the total."""
    keys = [ParamKey("width", level="3/2-#20", pair=1, L=1, S=0.5),
            ParamKey("width", level="3/2-#20", pair=2, L=2, S=0.5),
            ParamKey("width", level="3/2-#20", pair=3, L=1, S=0.5),
            ParamKey("width", level="3/2-#8", pair=1, L=1, S=0.5),
            ParamKey("energy", level="3/2-#20")]
    level_energy = {"3/2-#20": 10.30, "3/2-#8": 10.27}
    n_rmatrix = 5
    azr = None
    npts = 2
    prior = None
    ls_loss = None
    penalties = Pen()

    def __init__(self, target_n0, target_a0):
        self.t = np.array([target_n0, target_a0], float)

    def wigner_amplitudes(self):
        return dict(GW_AMP)

    def wigner_widths(self, x):
        return dict(GW_WID)

    def _gam(self, x):
        return np.array([GW_WID[0] * (x[0] / GW_AMP[0]) ** 2, GW_WID[1] * (x[1] / GW_AMP[1]) ** 2])

    def residuals(self, x, jac=False):
        x = np.asarray(x, float)
        r = (self._gam(x) - self.t) / 500.0
        J = None
        if jac:
            J = np.zeros((2, x.size))
            J[0, 0] = 2 * GW_WID[0] * x[0] / GW_AMP[0] ** 2 / 500.0
            J[1, 1] = 2 * GW_WID[1] * x[1] / GW_AMP[1] ** 2 / 500.0
        return r, J

    def rows_objective(self, r):
        return float(np.sum(np.asarray(r) ** 2))

    def rows_chi2(self, r):
        return float(np.sum(np.asarray(r)[:2] ** 2))


def amp(ch, gamma_ev):
    return GW_AMP[ch] * np.sqrt(gamma_ev / GW_WID[ch])


def main():
    base = dict(theta2_cap=3.0, bg_energy=11.5, e_window=0.3, threshold_margin_mev=0.0)
    pol = BoundsPolicy(**base, total_width_floor_windows=[[10.2995, 10.3005, 1.0e5]])
    check("total_width_floor_for: inside / outside", pol.total_width_floor_for(10.30) == 1.0e5 and pol.total_width_floor_for(10.27) is None)
    x = np.array([amp(0, 1.5e5), amp(1, 3.0e4), 0.1, 0.3, 10.30])
    ev = Stub(0.0, 0.0)
    b = default_bounds(ev, x, pol)
    check("one floor row, for the new level only", len(b.floors) == 1 and b.floors[0].level == "3/2-#20", [f.level for f in b.floors])
    f = b.floors[0]
    check("closed channel (Gamma_W = 0) not counted", sorted(f.cols.tolist()) == [0, 1], f.cols)
    check("total width = 150 + 30 keV", np.isclose(f.gamma(x), 1.8e5), f.gamma(x))
    check("no sign lock: amplitude bounds stay symmetric", np.allclose(b.lb[:2], -b.ub[:2]), (b.lb[:2], b.ub[:2]))
    check("above the floor: residual 0, Jacobian row 0", f.residual(x) == 0.0 and not np.any(f.jac_row(x, 5)))
    xb = x.copy(); xb[0] = amp(0, 4.0e4)                                  # 40 + 30 = 70 keV < 100
    r, J = floor_rows(b.floors, xb, 5, jac=True)
    check("below the floor: residual = (100 - 70) keV / (0.01 % of 100 keV)", np.isclose(r[0], 3.0e4 / (1e-4 * 1.0e5)), r)
    fd = np.zeros(5)
    for c in (0, 1):
        h = 1e-7; xp = xb.copy(); xp[c] += h; xm = xb.copy(); xm[c] -= h
        fd[c] = (floor_rows(b.floors, xp, 5)[0][0] - floor_rows(b.floors, xm, 5)[0][0]) / (2 * h)
    check("Jacobian row = finite difference", np.allclose(J[0], fd, rtol=1e-5, atol=1e-6), (J[0], fd))
    xoff = x.copy(); xoff[:3] = 0.0
    check("level switched off (all amplitudes 0, the level-add control): no floor", not default_bounds(ev, xoff, pol).floors)

    # polish: the data want n0 -> 0 and a0 -> 300 keV.  The total floor is satisfied by the alpha width: n0 goes to ~0.
    ev = Stub(0.0, 3.0e5)
    fr = polish(ev, x, free=np.array([True, True, False, False, False]), bounds=default_bounds(ev, x, pol), policy=pol,
                max_nfev=200)
    g = ev._gam(fr.x)
    check("alpha-dominated level passes the total floor: Gamma_n0 -> ~0, Gamma_a0 -> 300 keV, floor not binding",
          g[0] < 1e3 and abs(g[1] - 3.0e5) < 1e3 and not default_bounds(ev, x, pol).floors[0].binding(fr.x), g)
    # the OLD per-channel n0 floor forces Gamma_n0 >= 100 keV on the same data (what the evaluator flagged)
    old = BoundsPolicy(**base, width_floor_windows=[[10.2995, 10.3005, 1.0e5, 1]])
    fo = polish(ev, x, free=np.array([True, True, False, False, False]), bounds=default_bounds(ev, x, old), policy=old,
                max_nfev=200)
    check("old n0-only floor on the same data holds Gamma_n0 at 100 keV (the defect)", abs(ev._gam(fo.x)[0] - 1.0e5) < 2e3,
          ev._gam(fo.x))
    # the data want the level gone (both widths -> 0): the floor holds the total at ~100 keV, split freely, and binds
    ev = Stub(0.0, 0.0)
    bb = default_bounds(ev, x, pol)
    fr = polish(ev, x, free=np.array([True, True, False, False, False]), bounds=bb, policy=pol, max_nfev=300)
    gt = bb.floors[0].gamma(fr.x)
    check("unwanted level: total held within 0.1 % of the floor and reported binding", 0.999e5 < gt < 1.0e5 * 1.001
          and bb.floors[0].binding(fr.x), gt)
    check("polish objective includes the floor row (small, positive) and chi2 excludes it",
          fr.objective > fr.chi2 and fr.objective - fr.chi2 < 0.01 * fr.chi2, (fr.objective, fr.chi2))
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all total-width-floor checks passed")


if __name__ == "__main__":
    main()
