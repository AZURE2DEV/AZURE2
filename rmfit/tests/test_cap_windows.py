#!/usr/bin/env python3
"""theta^2 cap by energy window (the evaluator 2026-09-25: "don't relax the reduced width limit for all of the levels in the fit,
just those near the n1 threshold").  Pure checks of BoundsPolicy.cap_for, then default_bounds on the 13N fixture."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
import numpy as np
from rmfit.local import BoundsPolicy, default_bounds
ok = True
def check(name, cond, detail=""):
    global ok
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)
def main():
    p = BoundsPolicy(theta2_cap=1.0, bg_theta2_cap=3.0, bg_energy=11.5, theta2_cap_windows=[[9.95, 10.60, 10.0]])
    check("inside the window -> window cap", p.cap_for(10.2625) == 10.0)
    check("window edges are inclusive", p.cap_for(9.95) == 10.0 and p.cap_for(10.60) == 10.0)
    check("below the window -> physical cap", p.cap_for(9.9) == 1.0)
    check("above the window, below bg_energy -> physical cap", p.cap_for(10.61) == 1.0)
    check("background pole -> bg cap", p.cap_for(13.0) == 3.0)
    check("default (no windows) unchanged", BoundsPolicy().cap_for(10.2) == 1.0 and BoundsPolicy().cap_for(20.0) == 3.0)
    from rmfit.engine import SingleSession
    azr = os.path.join(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")), "tests", "13N", "13N.azr")
    with SingleSession(azr, cwd=os.path.dirname(azr), flags=dict(use_long_wavelength=False, ignore_externals=True, use_brune=True), name="capwin") as ev:
        x = np.array(ev.x0, float); gW = ev.wigner_amplitudes()
        E = sorted({round(v, 4) for v in ev.level_energy.values() if v < 10})
        tgt = E[len(E) // 2]; lo, hi = tgt - 1e-3, tgt + 1e-3
        b0 = default_bounds(ev, x, BoundsPolicy(bg_energy=15.0))
        b1 = default_bounds(ev, x, BoundsPolicy(bg_energy=15.0, theta2_cap_windows=[[lo, hi, 10.0]]))
        inside = [j for j, k in enumerate(ev.keys) if k.kind == "width" and j in gW and gW[j] > 0 and lo <= ev.level_energy.get(k.level, -1) <= hi]
        outside = [j for j, k in enumerate(ev.keys) if k.kind == "width" and j in gW and gW[j] > 0 and j not in inside]
        check(f"13N: widths of the level at {tgt} MeV get sqrt(10) x the bound", len(inside) > 0 and all(abs(b1.ub[j] / b0.ub[j] - np.sqrt(10)) < 1e-9 for j in inside), inside)
        check("13N: every other width bound unchanged", all(abs(b1.ub[j] - b0.ub[j]) < 1e-12 and abs(b1.lb[j] - b0.lb[j]) < 1e-12 for j in outside))
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1
if __name__ == "__main__":
    sys.exit(main())
