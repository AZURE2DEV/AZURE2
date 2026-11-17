#!/usr/bin/env python3
"""diagnose.window_delta (2026-10-12, the evaluator: "give the difference between the chi2 with and without the new resonance
based on something like +/- 3 times the width of the resonance to either side").  Stub evaluator, no engine.
Run:  python3 rmfit/tests/test_window_delta.py"""
import os
import sys
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.diagnose import window_delta, format_window_delta  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


class Segs:
    def __init__(self, files):
        self.f = files

    def by_key(self, key):
        return SimpleNamespace(data_file=self.f[key])


class Stub:
    """Two segments: A (key 1) on Ex 9.0-11.0, B (key 2) on 10.0-10.6, plus one row with no Ex.  Residual r_i = x[0] at
    rows near 10.3 (the 'level' helps there when x[0] = 0) and x[1] everywhere (a global offset, the noise term)."""
    def __init__(self):
        self.ex = np.r_[np.linspace(9.0, 11.0, 201), np.linspace(10.0, 10.6, 61), np.nan]
        self.npts = self.ex.size
        self.row_ranges = [(1, (0, 201)), (2, (201, 262)), (1, (262, 263))]
        self.segments = Segs({1: "data/setA.dat", 2: "data/setB.dat"})

    def row_excitation(self):
        return self.ex

    def residuals(self, x, jac=False):
        near = np.nan_to_num(np.abs(self.ex - 10.3), nan=9.0) <= 0.15
        r = np.where(near, x[0], 0.0) + x[1]
        return np.r_[r, 7.0], None          # one penalty row after the data rows: must be ignored

    def _plain_data_rows(self, r):
        return np.asarray(r, float)[:self.npts]


def main():
    ev = Stub()
    x_without = np.array([2.0, 0.0]); x_with = np.array([0.0, 0.0])
    res = window_delta(ev, x_with, x_without, 10.3, 0.05)            # Gamma 50 keV: windows +-50, 100, 150 keV
    near = int((np.nan_to_num(np.abs(ev.ex - 10.3), nan=9.0) <= 0.15).sum())
    check("total delta = -4 per helped row, penalty row ignored", np.isclose(res["total_delta"], -4.0 * near),
          (res["total_delta"], near))
    w3 = res["windows"][-1]
    check("+-3 Gamma window holds every helped row: outside = 0", np.isclose(w3["outside"], 0.0) and np.isclose(w3["delta"], res["total_delta"]), w3)
    w1 = res["windows"][0]
    check("+-1 Gamma window: smaller gain, the rest outside", w1["delta"] > w3["delta"] and np.isclose(w1["delta"] + w1["outside"], res["total_delta"]), w1)
    check("window bounds and point counts", np.isclose(w1["lo"], 10.25) and np.isclose(w1["hi"], 10.35) and w1["n"] < w3["n"], (w1, w3))
    names = {r["name"]: r for r in res["datasets"]}
    check("per data set names from the data files, both sets present", set(names) == {"setA", "setB"}, names)
    check("per data set deltas add up to the window delta", np.isclose(sum(r["delta"] for r in res["datasets"]), w3["delta"]))
    # a global offset (noise) lands mostly OUTSIDE the window
    res2 = window_delta(ev, np.array([0.0, 0.1]), np.array([2.0, 0.0]), 10.3, 0.05)
    check("a global offset shows up outside the window", res2["windows"][-1]["outside"] > 0, res2["windows"][-1])
    check("min_half_mev floors a tiny width", window_delta(ev, x_with, x_without, 10.3, 1e-6, min_half_mev=0.2)["windows"][0]["n"] > 0)
    lines = format_window_delta(res)
    check("formatter: header + 3 windows + table", len(lines) >= 5 and "+-3 Gamma" in lines[3], lines[:4])
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all window-delta checks passed")


if __name__ == "__main__":
    main()
