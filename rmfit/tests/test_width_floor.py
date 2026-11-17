#!/usr/bin/env python3
"""BoundsPolicy.width_floor_windows (2026-10-05, the evaluator: test a BROAD level at Ex 9.4 without letting it collapse): a level
whose seed energy lies in a window keeps Gamma_c >= gamma_min_ev in the window's pair, on the side of zero its amplitude is on
(one-sided sign lock), never above 0.95 of the cap.  Stubbed evaluator, no engine.  Run:  python3 rmfit/tests/test_width_floor.py"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.keys import ParamKey  # noqa: E402
from rmfit.local import BoundsPolicy, default_bounds  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


class Stub:
    """New 3/2- at 9.40 with n0 (pair 1, Gamma_W 4 MeV) and a0 (pair 2, Gamma_W 1 MeV); a second new level with a
    small n0 Wigner width (pair 1, Gamma_W 50 keV); an old 3/2- at 9.30 (pair 1)."""
    keys = [ParamKey("width", level="3/2-#9", pair=1, L=1, S=0.5),
            ParamKey("width", level="3/2-#9", pair=2, L=2, S=0.5),
            ParamKey("width", level="9/2+#4", pair=1, L=4, S=0.5),
            ParamKey("width", level="3/2-#3", pair=1, L=1, S=0.5),
            ParamKey("energy", level="3/2-#9")]
    level_energy = {"3/2-#9": 9.40, "9/2+#4": 9.40, "3/2-#3": 9.30}
    n_rmatrix = 5
    azr = None

    def wigner_amplitudes(self):
        return {0: 2.0, 1: 1.0, 2: 0.5, 3: 2.0}

    def wigner_widths(self, x):
        return {0: 4.0e6, 1: 1.0e6, 2: 5.0e4, 3: 4.0e6}


def main():
    ev = Stub()
    base = dict(theta2_cap=1.0, bg_energy=11.5, e_window=0.3, threshold_margin_mev=0.0)
    x = np.array([0.4, -0.2, -0.1, 0.3, 9.40])
    b0 = default_bounds(ev, x, BoundsPolicy(**base))
    check("no window: symmetric bounds", np.allclose(b0.lb[:4], -b0.ub[:4]))
    win = [9.3995, 9.4005, 1.0e5, 1]
    pol = BoundsPolicy(**base, width_floor_windows=[win])
    check("width_floor_for: pair match / other pair / outside", pol.width_floor_for(9.40, 1) == 1.0e5
          and pol.width_floor_for(9.40, 2) is None and pol.width_floor_for(9.30, 1) is None)
    b = default_bounds(ev, x, pol)
    amin = np.sqrt(1.0e5 / 4.0e6) * 2.0
    check("n0 of the new level: lb = sqrt(100 keV / 4 MeV) gamma_W, positive side", np.isclose(b.lb[0], amin) and np.isclose(b.ub[0], 2.0), (b.lb[0], b.ub[0]))
    check("a0 of the new level untouched (other pair)", np.isclose(b.lb[1], -1.0) and np.isclose(b.ub[1], 1.0), (b.lb[1], b.ub[1]))
    check("negative amplitude: floor on the negative side, capped at 0.95 of the cap when the floor exceeds it",
          np.isclose(b.ub[2], -0.95 * 0.5) and np.isclose(b.lb[2], -0.5), (b.lb[2], b.ub[2]))
    check("old level outside the window untouched", np.isclose(b.lb[3], -2.0), b.lb[3])
    lock = {ev.keys[0]: -1}
    b = default_bounds(ev, x, pol, sign_lock=lock)
    check("explicit sign lock wins the side, floor sets the magnitude", np.isclose(b.ub[0], -amin) and np.isclose(b.lb[0], -2.0), (b.lb[0], b.ub[0]))
    xoff = x.copy(); xoff[0] = 0.0
    b = default_bounds(ev, xoff, pol)
    check("a zeroed amplitude (level switched off, the level-add control) gets NO floor and stays at 0 after clip",
          np.isclose(b.lb[0], -2.0) and np.isclose(b.ub[0], 2.0) and b.clip(xoff)[0] == 0.0, (b.lb[0], b.ub[0]))
    pol_any = BoundsPolicy(**base, width_floor_windows=[[9.3995, 9.4005, 1.0e5]])
    check("window without a pair applies to every channel", pol_any.width_floor_for(9.40, 2) == 1.0e5)
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all width-floor checks passed")


if __name__ == "__main__":
    main()
