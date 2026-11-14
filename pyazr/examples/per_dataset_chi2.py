#!/usr/bin/env python3
"""Per-data-set chi-squared, and which data set is fighting the fit.

``Session.evaluate(x)`` runs one evaluation; each of its ``data[i]`` carries the data set's chi-squared (``chi2``,
the row of the program's ``chiSquared.out``), its number of points and its normalization prior (``norm_chi2``).
Where a data covariance joins several data sets only the group's sum has a meaning: ``chi2`` splits it as
``chiSquared.out`` does, and the sum of the squares of the data set's rows of the objective
(``objective().residuals(x)[d.rows]``, whitened rows) splits it otherwise; both add up to the group's chi-squared.

Usage (from the repository's root):

    python3 pyazr/examples/per_dataset_chi2.py tests/upgrade/T2/T2.azr [--params output/param.sav]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np

from pyazr import Session


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path, help="the format-3 project (.azr)")
    ap.add_argument("--params", type=Path, default=None,
                    help="a param.sav to load instead of the .azr values")
    args = ap.parse_args()

    with Session(str(args.azr), parameter_file=args.params) as azr:
        x = np.asarray(azr.parameters.x0, float)

        ev = azr.evaluate(x)
        rows = azr.objective().residuals(x)
        total, n = ev.chi2, sum(d.n for d in ev.data)
        print(f"{args.azr.name}: chi2 = {total:.1f} over {n} points "
              f"(chi2/N = {total / n:.2f})\n")
        print(f"{'data set':28s} {'N':>6} {'chi2':>12} {'chi2/N':>9} {'rows^2':>12} {'norm prior':>10}")
        for d in sorted(ev.data, key=lambda d: -d.chi2 / max(d.n, 1)):   # worst first
            whitened = float(np.sum(rows[d.rows] ** 2))
            print(f"{d.name:28s} {d.n:6d} {d.chi2:12.1f} {d.chi2 / max(d.n, 1):9.2f} {whitened:12.1f} "
                  f"{d.norm_chi2:10.3g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
