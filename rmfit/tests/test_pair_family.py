#!/usr/bin/env python3
"""The two-channel ('pair') template family (2026-09-21): a level carried by two channels."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
import numpy as np
from rmfit.levelscan import _patterns, grid_coverage, partition_fractions
ok = True
def check(name, cond, detail=""):
    global ok
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)
def main():
    # carrying channels of the 9/2+ dummy at 11.07 (screen order) and their Wigner widths (eV)
    order = ["n0 L4", "a0 L5", "n1 L4", "a1 L4", "n2 L3 S2.5", "n2 L1", "n2 L3 S3.5", "a2 L3", "n3 L2", "a3 L2"]
    gw = [7.91e5, 1.96e5, 172, 536, 7.56e3, 2.35e6, 7.56e3, 71, 5.2, 88.9]
    # the group's priors (the other 9/2+ levels): none has the 44/56 n0/n2 balance
    pr = lambda d: dict(label="p", fractions=[d.get(c, 0.0) for c in order], signs=[-1 if d.get(c, 0.0) < 0 else 1 for c in order])
    priors = [pr({"n0 L4": 1.0}), pr({"n0 L4": 0.915, "a0 L5": 0.085}), pr({"n0 L4": -0.009, "a0 L5": -0.047, "a1 L4": -0.003, "n2 L1": -0.941}),
              pr({"n0 L4": -0.027, "a0 L5": 0.050, "n2 L1": 0.923}), pr({"n0 L4": -0.222, "a0 L5": 0.233, "n2 L1": -0.544})]
    for p in priors: s = sum(abs(f) for f in p["fractions"]); p["fractions"] = [abs(f) / s for f in p["fractions"]]
    old = _patterns(10, 32, priors=priors, pair_top=0)      # the grid the 09-21 screen used
    new = _patterns(10, 32, priors=priors, weights=gw)
    true = [-87860, -151.3, 0, 52.19, 0, -111200, 0, 0, 0, 0]
    totals = [199210 * s for s in (0.5, 1.0, 2.0)]
    c_old, c_new = grid_coverage(true, old, totals), grid_coverage(true, new, totals)
    check("the old grid had 74 shapes (as the 09-21 screen logged)", len(old) == 74, len(old))
    check("the old grid could not represent the level (distance > 0.1)", c_old["distance"] > 0.1, c_old["distance"])
    check("the new grid contains a pair template within 0.03 of it", c_new["distance"] < 0.03 and c_new["desc"].startswith("pair"), (c_new["distance"], c_new["desc"]))
    check("... with the right relative sign (both negative = same sign)", c_new["signs"][0] == c_new["signs"][5], c_new["signs"])
    check("... using n0 L4 and n2 L1", c_new["fractions"][0] > 0.2 and c_new["fractions"][5] > 0.4, c_new["fractions"])
    check("the pair family is bounded: at most 60 extra shapes for 5 candidate channels", 74 < len(new) <= 74 + 60, len(new))
    check("no pair family for a single channel", len(_patterns(1, 32, weights=[1.0])) == len(_patterns(1, 32)))
    nw = [d for _, _, d in _patterns(10, 32, priors=priors) if d.startswith("pair")]
    check("without weights the pair channels are the first pair_top plus the priors' channels", nw and all(d.split()[1] in ("0+1", "0+2", "0+3", "0+5", "1+2", "1+3", "1+5", "2+3", "2+5", "3+5") for d in nw), sorted(set(d.split()[1] for d in nw)))
    n2 = _patterns(3, 8, weights=[1.0, 2.0, 3.0]); pairs = [d for _, _, d in n2 if d.startswith("pair")]
    check("3 channels, no priors: 3 pairs x 3 splits x 2 signs = 18 pair shapes", len(pairs) == 18, len(pairs))
    check("pair fractions sum to one", all(abs(sum(f) - 1) < 1e-9 for f, _, d in new if d.startswith("pair")))
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1
if __name__ == "__main__":
    sys.exit(main())
