#!/usr/bin/env python3
"""diagnose.dataset_table on an evaluator whose segments are split into pieces (as Sharded
does): chi2 must be counted once per base segment, points once per piece.  Regression for
the 2026-09-18 bug (28 segments x 2 pieces, reported at exactly twice its
chi2 -- and therefore at twice its weight in compare()'s redistribution ratio)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from types import SimpleNamespace as NS
from rmfit.diagnose import dataset_table, compare


class Segs:
    def __init__(self, files): self.files = files
    def by_key(self, k): return NS(data_file=self.files[k])


def main():
    files = {1: "data/A.dat", 2: "data/A.dat", 3: "data/B.dat"}
    # segment 1 split in two pieces, segment 2 whole, segment 3 split in three
    ev = NS(segments=Segs(files), row_ranges=[(1, (0, 10)), (2, (10, 30)), (3, (30, 35)), (1, (35, 50)), (3, (50, 60)), (3, (60, 70))])
    before = NS(segment_chi2={1: 100.0, 2: 50.0, 3: 30.0}, objective=180.0)
    after = NS(segment_chi2={1: 60.0, 2: 55.0, 3: 30.0}, objective=145.0)
    t = {d["name"]: d for d in dataset_table(ev, before)}
    ok = True
    def check(name, cond, detail=""):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)
    check("split segment counted once (A = 150)", abs(t["A"]["chi2"] - 150.0) < 1e-12, t["A"])
    check("three-piece segment counted once (B = 30)", abs(t["B"]["chi2"] - 30.0) < 1e-12, t["B"])
    check("points summed over pieces (A = 45, B = 25)", t["A"]["n"] == 45 and t["B"]["n"] == 25, (t["A"]["n"], t["B"]["n"]))
    check("dataset chi2 sums to the total", abs(sum(d["chi2"] for d in t.values()) - 180.0) < 1e-12)
    c = compare(ev, before, after)
    # compare() nets per DATASET: A = segment 1 (-40) + segment 2 (+5) = -35, B unchanged
    check("compare: A nets to -35, no losses, total -35", abs(c["gains"] - 35) < 1e-12 and c["losses"] == 0 and abs(c["total_delta"] + 35) < 1e-12, c)
    check("compare: the dataset deltas sum to the objective change", abs(sum(r["delta"] for r in c["rows"]) - c["total_delta"]) < 1e-12, c)
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
