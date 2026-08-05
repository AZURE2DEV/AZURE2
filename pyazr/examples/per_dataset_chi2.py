#!/usr/bin/env python3
"""Per-dataset chi-squared, which the API does not report directly.

`calculate_chi2_rwa` gives one number per *segment*, but a fit is usually judged
per *experiment*, and an experiment often owns several segments -- one per angle,
or one per final state.  The residuals carry the information: they come back in
segment order, so slicing them by segment length and summing the squares gives
the per-segment chi-squared, which can then be grouped by data file.

The same slicing is how you find which experiment is fighting the fit.

Usage:
    python3 per_dataset_chi2.py <file.azr> [--params output/param.sav]
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

from pyazr import azure2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--params", type=Path, default=None,
                    help="a param.sav to load instead of the .azr values")
    args = ap.parse_args()

    with azure2(str(args.azr)) as azr:
        x = np.asarray(azr.params_rwa, float)
        if args.params and args.params.exists():
            saved = np.loadtxt(args.params, usecols=1)
            if saved.shape == x.shape:
                x = saved

        r, _ = azr.residual_jacobian(x)
        r = np.asarray(r, float)

        lengths = [len(azr.energies[i]) for i in range(azr.nsegments)]
        edges = np.cumsum([0] + lengths)
        per_segment = np.array([np.sum(r[edges[i]:edges[i + 1]] ** 2)
                                for i in range(azr.nsegments)])

        by_file = defaultdict(lambda: [0.0, 0])
        rows = []
        for i in range(azr.nsegments):
            d = azr.datasets[i]
            rows.append((d.name, lengths[i], per_segment[i]))
            by_file[d.name][0] += per_segment[i]
            by_file[d.name][1] += lengths[i]

        total, n = per_segment.sum(), sum(lengths)
        print(f"{args.azr.name}: chi2 = {total:.1f} over {n} points "
              f"(chi2/N = {total / n:.2f})\n")
        print(f"{'segment':>7}  {'dataset':28s} {'N':>6} {'chi2':>12} {'chi2/N':>9}")
        for i, (name, ln, c) in enumerate(rows, 1):
            print(f"{i:7d}  {name:28s} {ln:6d} {c:12.1f} {c/max(ln,1):9.2f}")

        print(f"\ngrouped by dataset, worst first:")
        for name, (c, ln) in sorted(by_file.items(), key=lambda kv: -kv[1][0] / max(kv[1][1], 1)):
            print(f"   {name:28s} {ln:6d} {c:12.1f} {c/max(ln,1):9.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
