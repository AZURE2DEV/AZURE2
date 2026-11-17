#!/usr/bin/env python3
"""BoundsPolicy.width_bound_windows (2026-09-28): every free partial width of a level whose seed energy
lies in a window is bounded to Gamma_c <= gamma_max_ev, i.e. theta_c^2 <= gamma_max_ev / Gamma_W,c, on
top of the theta^2 cap.  Stubbed evaluator, no engine.  Run:  python3 rmfit/tests/test_width_bound.py"""
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
    """Two levels: the new 3/2+ at 13.226 with three channels (Wigner widths 1 MeV, 20 keV, 200 keV)
    and an old 5/2- at 12.0 with one channel (1 MeV)."""
    keys = [ParamKey("width", level="3/2+#9", pair=2, L=1, S=0.5),
            ParamKey("width", level="3/2+#9", pair=4, L=2, S=0.5),
            ParamKey("width", level="3/2+#9", pair=1, L=1, S=0.5),
            ParamKey("width", level="5/2-#2", pair=2, L=3, S=0.5),
            ParamKey("energy", level="3/2+#9")]
    level_energy = {"3/2+#9": 13.226, "5/2-#2": 12.0}
    n_rmatrix = 5
    azr = None

    def wigner_amplitudes(self):
        return {0: 2.0, 1: 0.3, 2: 1.0, 3: 2.0}

    def wigner_widths(self, x):
        return {0: 1.0e6, 1: 2.0e4, 2: 2.0e5, 3: 1.0e6}


def main():
    ev = Stub()
    x = np.array([0.1, 0.05, 0.1, 0.1, 13.226])
    pol = BoundsPolicy(theta2_cap=0.5, bg_energy=13.3, e_window=0.3, threshold_margin_mev=0.0)
    b0 = default_bounds(ev, x, pol)
    check("no window: theta2 cap 0.5 on every channel", np.allclose(b0.ub[:4], np.sqrt(0.5) * np.array([2.0, 0.3, 1.0, 2.0])), b0.ub[:4])

    pol = BoundsPolicy(theta2_cap=0.5, bg_energy=13.3, e_window=0.3, threshold_margin_mev=0.0,
                       width_bound_windows=[BoundsPolicy.seed_window(13.226, 4.0e4)[:2] + [4.0e4]])
    check("width_bound_for inside / outside", pol.width_bound_for(13.226) == 4.0e4 and pol.width_bound_for(12.0) is None)
    b = default_bounds(ev, x, pol)
    check("1 MeV-Wigner channel bounded to 40 keV: theta2 0.04", np.isclose(b.ub[0], np.sqrt(0.04) * 2.0), b.ub[0])
    check("20 keV-Wigner channel: bound (2) above the cap, cap 0.5 stays", np.isclose(b.ub[1], np.sqrt(0.5) * 0.3), b.ub[1])
    check("200 keV-Wigner channel bounded to theta2 0.2", np.isclose(b.ub[2], np.sqrt(0.2) * 1.0), b.ub[2])
    check("other level untouched", np.isclose(b.ub[3], np.sqrt(0.5) * 2.0), b.ub[3])
    check("bounds symmetric", np.allclose(b.lb[:4], -b.ub[:4]))

    pol3 = BoundsPolicy(theta2_cap=0.5, bg_energy=13.3, e_window=0.3, threshold_margin_mev=0.0,
                        theta2_cap_windows=[BoundsPolicy.seed_window(13.226, 3.0)],
                        width_bound_windows=[BoundsPolicy.seed_window(13.226, 4.0e4)[:2] + [4.0e4]])
    b = default_bounds(ev, x, pol3)
    check("with cap 3 the width bound still wins on the 1 MeV channel", np.isclose(b.ub[0], np.sqrt(0.04) * 2.0), b.ub[0])
    check("with cap 3 the 20 keV channel gets min(3, 2) = 2", np.isclose(b.ub[1], np.sqrt(2.0) * 0.3), b.ub[1])

    # 2026-10-05: a 4th element restricts the bound to one pair
    polp = BoundsPolicy(theta2_cap=0.5, bg_energy=13.3, e_window=0.3, threshold_margin_mev=0.0,
                        width_bound_windows=[BoundsPolicy.seed_window(13.226, 4.0e4)[:2] + [4.0e4, 2]])
    check("pair-specific window: width_bound_for(pair 2) / (pair 4) / (no pair)", polp.width_bound_for(13.226, 2) == 4.0e4
          and polp.width_bound_for(13.226, 4) is None and polp.width_bound_for(13.226) == 4.0e4)
    b = default_bounds(ev, x, polp)
    check("pair-specific window bounds only the pair-2 channel", np.isclose(b.ub[0], np.sqrt(0.04) * 2.0)
          and np.isclose(b.ub[1], np.sqrt(0.5) * 0.3) and np.isclose(b.ub[2], np.sqrt(0.5) * 1.0), b.ub[:3])
    polp2 = BoundsPolicy(theta2_cap=0.5, bg_energy=13.3, e_window=0.3, threshold_margin_mev=0.0,
                         width_bound_windows=[BoundsPolicy.seed_window(13.226, 4.0e4)[:2] + [4.0e4, 4],
                                              BoundsPolicy.seed_window(13.226, 1.0e5)[:2] + [1.0e5]])
    b = default_bounds(ev, x, polp2)
    check("a non-matching pair entry falls through to the next window", np.isclose(b.ub[0], np.sqrt(0.1) * 2.0)
          and np.isclose(b.ub[1], np.sqrt(0.5) * 0.3), b.ub[:2])

    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all width-bound checks passed")


if __name__ == "__main__":
    main()
