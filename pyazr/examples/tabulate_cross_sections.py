#!/usr/bin/env python3
"""Tabulate a cross section on an adaptive energy grid, for a table another code reads.

``pyazr.tabulate(project, pairs=[(entrance, exit)], e_min=, e_max=, rel_tol=)`` places knots so that log-linear
interpolation between them reproduces the engine's sigma(E) to the relative tolerance asked for: the knots pile up
across resonances and thin out on smooth stretches.  The object it returns interpolates (``tab(E_cm)``, in b) and
writes a three-column table (E_cm, E_lab in MeV, sigma in b) with ``tab.save(path)``.  Central values only: a band
comes from ``pyazr.uncertainty_bands`` (``uncertainty_band.py``).

``tabulate`` works in a scratch directory (``cwd=``; here a new temporary one, removed at the end) and writes the
table only where ``--out`` names it.

Usage (from the repository's root):

    python3 pyazr/examples/tabulate_cross_sections.py tests/upgrade/T1/T1.azr [--entrance 1 --exit 2] \\
        [--emin 0.1 --emax 2.0] [--tolerance 5e-3] [--out sigma.dat]
"""
import argparse
import os
import shutil
import sys
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np  # noqa: E402

from pyazr import tabulate  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr)")
    ap.add_argument("--entrance", type=int, default=1, help="the entrance pair's key")
    ap.add_argument("--exit", dest="exit_pair", type=int, default=2, help="the exit pair's key")
    ap.add_argument("--emin", type=float, default=0.1, help="lowest c.m. energy (MeV)")
    ap.add_argument("--emax", type=float, default=2.0, help="highest c.m. energy (MeV)")
    ap.add_argument("--tolerance", type=float, default=5e-3, help="relative interpolation tolerance")
    ap.add_argument("--out", default=None, help="write the table (E_cm, E_lab in MeV, sigma in b) here")
    args = ap.parse_args()

    scratch = tempfile.mkdtemp(prefix="tabulate.")
    try:
        pair = (args.entrance, args.exit_pair)
        tab = tabulate(args.project, pairs=[pair], e_min=args.emin, e_max=args.emax, rel_tol=args.tolerance,
                       cwd=scratch)[pair]
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    print(tab)
    spacing = np.diff(tab.e_cm)
    print(f"{len(tab.e_cm)} knots; spacing {1e3 * spacing.min():.2f} to {1e3 * spacing.max():.1f} keV")
    peak = int(np.argmax(tab.sigma))
    print(f"largest sigma {tab.sigma[peak]:.5g} b at E_cm = {tab.e_cm[peak]:.4f} MeV")
    e = 0.5 * (args.emin + args.emax)
    print(f"sigma({e:.3f} MeV) = {float(tab(e)):.5g} b (interpolated)")
    if args.out:
        tab.save(args.out)
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
