#!/usr/bin/env python3
"""The thermonuclear reaction rate N_A<sigma v> of a project, its band and its Gamow window (pyazr.rate).

The engine computes the rate as the program's ``AZURE2 rate`` command does (src/ReactionRate.cpp: sigma on a grid
with a cluster of points around every level, a spline of ln sigma against ln E, the integral at each temperature);
nothing is integrated here.  This script prints the rate on the rate tables' temperature grid, the Gamow window of
each temperature, and, with a fit's parameter file and covariance, the 1-sigma band by linear propagation; it writes
the table as CSV with its header.  For a rate per posterior sample, ``rate.rate_band(session, draws=..., names=...)``
gives the 16th, 50th and 84th percentiles over the samples.

Usage (from the repository's root):

    python3 pyazr/examples/reaction_rate.py tests/upgrade/T1/T1.azr --entrance 1 --exit 2 \\
        [--params output/param.sav --covariance output/covariance.dat] [--csv rate.csv]
"""

import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "2")

from pyazr import rate, read_covariance_dat  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr")
    ap.add_argument("--entrance", type=int, default=1)
    ap.add_argument("--exit", dest="exit_pair", type=int, default=2, help="-1: the total radiative capture")
    ap.add_argument("--grid", default="table", choices=sorted(rate.TEMPERATURE_GRIDS))
    ap.add_argument("--params", default=None, help="a parameter file (param.sav) of a fit")
    ap.add_argument("--covariance", default=None, help="the fit's covariance.dat (with its #columns line)")
    ap.add_argument("--csv", default=None, help="write the table here")
    args = ap.parse_args()

    with rate.RateSession(args.azr, args.entrance, args.exit_pair, parameter_file=args.params) as s:
        r = s.rate(grid=args.grid)
        band = None
        if args.covariance:
            C, names, _ = read_covariance_dat(args.covariance)
            band = rate.rate_band(s, covariance=C, names=names, grid=args.grid)
    g = rate.gamow_window(args.azr, args.entrance, r.T9)
    print(f"N_A<sigma v> of {r.reaction} ({r.units['rate']}), parameters: {r.source}")
    print(f"{'T9':>7} {'rate':>13} {'E0 (MeV)':>9} {'window (MeV)':>17}" + ("  1-sigma band" if band else ""))
    for i, t in enumerate(r.T9):
        line = f"{t:7.3g} {r.rate[i]:13.5e} {g['E0'][i]:9.4f} {g['low'][i]:8.4f}-{g['high'][i]:<8.4f}"
        if band:
            line += f"  {band.low[i]:.4e} - {band.high[i]:.4e}"
        print(line)
    if args.csv:
        rate.write_rate_table(args.csv, r, band=band)
        print(f"written: {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
