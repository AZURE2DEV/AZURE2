#!/usr/bin/env python3
"""Scan the channel radius and watch what chi-squared does about it.

The channel radius is a boundary, not an observable, so the fit ought not to
care much where it is put.  In practice it does, because the external region is
assumed free of nuclear forces and it is not: a_c sets the hard-sphere phase and
the lower limit of the external-capture integral, so a scan of chi2(a_c) is the
standard way of choosing it.

The scan is a loop over *files*: the radius lives in the <levels> block, and a
session reads its model once at startup, so each radius is a new .azr and a new
process.  `AzrModel.set_channel_radius` rewrites every channel line of that pair.

The external-capture caches are keyed on the grid, not on the radius, so
`output/intEC.dat` and `intEC.extrap` must go with each change -- keeping them
reuses integrals computed from the previous radius, which is wrong and silent.

Usage:
    python3 channel_radius_scan.py <file.azr> --pair 1 --radii 3.5 4.0 4.5 5.0
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

from pyazr import AzrModel

SCAN = """
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from pyazr import azure2
with azure2(sys.argv[1]) as m:
    x = np.asarray(m.params_rwa, float)
    print(json.dumps({"chi2": float(np.sum(m.calculate_chi2_rwa(x))),
                      "n": int(sum(len(m.energies[i])
                                   for i in range(m.nsegments)))}))
"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--pair", type=int, default=1)
    ap.add_argument("--radii", type=float, nargs="+",
                    default=[3.5, 4.0, 4.5, 5.0, 5.5, 6.0])
    args = ap.parse_args()

    here = args.azr.parent
    runner = here / "_scan_one.py"
    runner.write_text(SCAN)
    out = here / "output"

    mdl = AzrModel.from_file(args.azr)
    print(f"channel radii in {args.azr.name}: {mdl.channel_radii()}")
    print(f"\n{'a_c (fm)':>9} {'chi2':>12} {'chi2/N':>9}")

    rows = []
    for a in args.radii:
        m2 = AzrModel.from_file(args.azr)
        m2.set_channel_radius(args.pair, a)
        path = m2.write(here / f"_radius_{a:g}.azr")
        for cache in ("intEC.dat", "intEC.extrap"):
            f = out / cache
            if f.exists():
                f.unlink()
        # A separate process per radius: one AZURE2 session per interpreter.
        p = subprocess.run([sys.executable, str(runner), str(path)],
                           cwd=here, capture_output=True, text=True)
        line = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ""
        try:
            import json
            r = json.loads(line)
        except Exception:
            print(f"{a:9.2f}   failed: {p.stderr.strip().splitlines()[-1:]}")
            continue
        rows.append((a, r["chi2"], r["chi2"] / max(r["n"], 1)))
        print(f"{a:9.2f} {r['chi2']:12.1f} {r['chi2']/max(r['n'],1):9.3f}")
        Path(path).unlink()

    runner.unlink()
    if rows:
        best = min(rows, key=lambda t: t[1])
        span = max(r[1] for r in rows) - min(r[1] for r in rows)
        print(f"\nminimum at a_c = {best[0]:.2f} fm, chi2 = {best[1]:.1f}")
        print(f"spread over the scan: {span:.1f} in chi2 "
              f"({span / max(best[1], 1) * 100:.1f}% of the minimum)")
        print("A large spread means a_c is acting as a fit parameter; see the "
              "hybrid\nexternal-region option for the alternative.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
