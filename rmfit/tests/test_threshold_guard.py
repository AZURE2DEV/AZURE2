#!/usr/bin/env python3
"""Free level energies keep a margin from particle thresholds (2026-09-23): pure test of keep_off_thresholds
plus particle_thresholds on the 13N fixture."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.local import keep_off_thresholds, particle_thresholds, BoundsPolicy
ok = True
def check(name, cond, detail=""):
    global ok
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)
def main():
    T = (10.2726,); m = 0.003
    check("a level 0.4 keV above the n2 threshold gets lb = T + margin", keep_off_thresholds(9.973, 10.573, 10.2730, T, m) == (10.2756, 10.573))
    check("a level below the threshold gets ub = T - margin", keep_off_thresholds(9.9, 10.5, 10.2500, T, m) == (9.9, 10.2696))
    check("a window that does not touch the band is untouched", keep_off_thresholds(9.0, 9.5, 9.2, T, m) == (9.0, 9.5))
    check("margin 0 disables the guard", keep_off_thresholds(9.973, 10.573, 10.2730, T, 0.0) == (9.973, 10.573))
    check("several thresholds: the nearest on each side bounds the window",
          keep_off_thresholds(9.4, 10.3, 10.1, (9.4481, 10.0432, 10.1925, 10.273), m) == (10.0462, 10.1895))
    thr = particle_thresholds(os.path.join(os.path.expanduser("~"), "AZURE2", "tests", "13N", "13N.azr"))
    check("13N fixture: one particle threshold, the proton one, photons excluded", len(thr) == 1 and abs(thr[0] - 1.94351) < 1e-4, thr)
    check("BoundsPolicy default margin is 3 keV", abs(BoundsPolicy().threshold_margin_mev - 0.003) < 1e-12)
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1
if __name__ == "__main__":
    sys.exit(main())
