#!/usr/bin/env python3
"""levelscan.scoring_segments: which segments a candidate level is judged on locally.
Regression for 2026-09-18: CM-frame differential data (isDiff 4, e.g. secondary gamma rays
such as (a,a1 gamma) and (a,n2 gamma)) were left out because only the literal observable
"differential" was accepted."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from types import SimpleNamespace as NS
from rmfit.levelscan import scoring_segments


class Segs:
    def __init__(self, obs): self.obs = obs
    def by_key(self, k): return NS(observable=self.obs[k])


def main():
    obs = {1: "differential", 2: "differential-cm", 3: "angle-integrated", 4: "total-capture", 5: "differential", 6: "analyzing-power"}
    ranges = {1: (10.7, 11.1), 2: (10.3, 11.2), 3: (7.0, 11.2), 4: (10.0, 11.0), 5: (8.0, 9.0), 6: (10.8, 11.0)}
    ev = NS(segments=Segs(obs))
    ok = True
    def check(name, cond, detail=""):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)
    got = scoring_segments(ev, 10.9, half_width=0.25, ex_ranges=ranges)
    check("lab-differential segment near the level is used", 1 in got, got)
    check("CM-differential segment near the level is used", 2 in got, got)
    check("angle-integrated and total segments are not", 3 not in got and 4 not in got, got)
    check("a differential segment far from the level is not", 5 not in got, got)
    allk = scoring_segments(ev, 10.9, half_width=0.25, differential_only=False, ex_ranges=ranges)
    check("differential_only=False takes everything in range", allk == [1, 2, 3, 4, 6], allk)
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
