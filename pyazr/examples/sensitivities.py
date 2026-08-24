#!/usr/bin/env python3
"""Which parameters does a dataset actually constrain?

The sensitivity matrix d(observable)/d(parameter) answers that, and it is what
an uncertainty band is built from: sigma^2 = g^T C g.  Read on its own it says
something a fit cannot -- a parameter with a large uncertainty may still be
irrelevant here, and a well-determined one may be carrying the whole segment.

Raw derivatives are not comparable across parameters, because a level energy
(MeV) and a reduced-width amplitude (MeV^1/2) do not share units.  The
dimensionless form does:

    elasticity  =  d ln(sigma) / d ln(p)  =  (p / sigma) * d(sigma)/d(p)

"a 1% change in this parameter moves the cross section by this many percent",
comparable across every column and every segment.

Two ways to get the matrix:

    analytic    one reverse-mode adjoint per point -- about two forward
                evaluations for the whole matrix, whatever the parameter count.
    central     finite differences, 2 forward passes per column.  Slower by a
                factor of the column count, and the right way to check the
                analytic path.  --nprocs spreads the columns over worker
                processes, each with its own engine.

Usage:
    python3 sensitivities.py <file.azr> [--segment 1] [--params output/param.sav]
    python3 sensitivities.py <file.azr> --check              # analytic vs findiff
    python3 sensitivities.py <file.azr> --method central --nprocs 8
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

from pyazr import azure2
from pyazr.bands import best_fit_params, rmatrix_columns, sensitivities


def label(param):
    """A short, readable name for one R-matrix parameter."""
    if param.kind == "energy":
        return f"E({param.jpi})"
    channel = f"L={param.L} S={param.S:g} pair {param.pair}"
    if param.radiation_type in ("E", "M"):
        channel = f"{param.radiation_type}{param.L} pair {param.pair}"
    return f"g({param.jpi}, {channel})"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--segment", type=int, default=1,
                    help="1-based data segment (default 1)")
    ap.add_argument("--params", type=Path, default=None,
                    help="a param.sav / param.par to use instead of the .azr values")
    ap.add_argument("--method", default="analytic",
                    choices=("analytic", "central", "forward"))
    ap.add_argument("--step", type=float, default=1e-4,
                    help="finite-difference step, relative to each parameter")
    ap.add_argument("--nprocs", type=int, default=1,
                    help="worker processes for the finite-difference columns")
    ap.add_argument("--top", type=int, default=12, help="rows to print")
    ap.add_argument("--check", action="store_true",
                    help="also finite-difference, and compare with the analytic matrix")
    args = ap.parse_args()

    with azure2(str(args.azr)) as azr:
        x = np.asarray(azr.params_rwa, float)
        if args.params:
            x = np.asarray(best_fit_params(azr, str(args.params)), float)

        seg = args.segment - 1
        if not 0 <= seg < azr.nsegments:
            print(f"{args.azr.name} has {azr.nsegments} active segments; "
                  f"no {args.segment}.", file=sys.stderr)
            return 1

        # The columns are the free R-matrix parameters -- level energies and
        # reduced widths.  Normalizations and energy shifts are excluded
        # because no calculated observable depends on them.
        cols = rmatrix_columns(azr)
        free = sorted((p for p in azr.parameters if not p.fixed),
                      key=lambda p: p.free_index)
        names = [label(free[j]) for j in cols]

        t0 = time.time()
        G = sensitivities(azr, [seg], params=x, method=args.method,
                          step=args.step, cols=cols, nprocs=args.nprocs)[seg]
        elapsed = time.time() - t0

        sigma = np.asarray(azr.calculate_rwa(x)[seg], float)
        energy = np.asarray(azr.calculate_energies(x)[seg], float)

        how = args.method
        if args.method != "analytic" and args.nprocs > 1:
            how += f", {args.nprocs} processes"
        print(f"{args.azr.name}  segment {args.segment}: {len(energy)} points, "
              f"{len(cols)} free R-matrix parameters")
        print(f"sensitivity matrix {G.shape} in {elapsed:.2f} s ({how})\n")

        # d ln(sigma)/d ln(p): dimensionless, so columns are comparable.
        p = x[cols]
        with np.errstate(divide="ignore", invalid="ignore"):
            elasticity = G * p[None, :] / sigma[:, None]
        elasticity = np.nan_to_num(elasticity, nan=0.0, posinf=0.0, neginf=0.0)

        peak = np.abs(elasticity).max(axis=0)
        where = np.abs(elasticity).argmax(axis=0)
        order = np.argsort(-peak)

        print(f"{'parameter':<34} {'value':>12} {'|dlnS/dlnp|':>12} "
              f"{'at E (MeV)':>11}")
        for j in order[:args.top]:
            print(f"{names[j]:<34} {p[j]:12.4g} {peak[j]:12.3g} "
                  f"{energy[where[j]]:11.4f}")
        if len(order) > args.top:
            print(f"... {len(order) - args.top} more, all below "
                  f"{peak[order[args.top]]:.3g}")

        dead = [names[j] for j in order if peak[j] < 1e-6]
        if dead:
            print(f"\nthis segment says nothing about: {', '.join(dead)}")

        if args.check:
            t0 = time.time()
            F = sensitivities(azr, [seg], params=x, method="central",
                              step=args.step, cols=cols, nprocs=args.nprocs)[seg]
            t_fd = time.time() - t0
            # Scale each column by its own size, but never by less than a
            # sliver of the whole matrix: a column the observable does not
            # depend on is all round-off, and dividing round-off by round-off
            # manufactures a disagreement out of numbers that are both zero.
            floor = 1e-8 * np.abs(F).max()
            scale = np.maximum(np.abs(F).max(axis=0), floor)
            err = np.abs(G - F).max(axis=0) / scale
            print(f"\nanalytic vs central differences ({t_fd:.2f} s):")
            print(f"  worst column: {names[int(err.argmax())]}  "
                  f"rel. error {err.max():.2e}")
            print(f"  median over columns: {np.median(err):.2e}")
            if err.max() > 1e-3:
                print("  NOTE: larger than finite-difference noise should be; "
                      "try --step 1e-3 before believing either one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
