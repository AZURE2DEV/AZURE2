#!/usr/bin/env python3
"""The 1-sigma band of an S factor from a fit's covariance, on an extrapolation of the project and on a new grid.

The band on any quantity the model calculates is linear propagation of the fit's parameter covariance C,

    sigma_y^2 = g^T C g,        g_j = d y / d p_j,

with g from the engine's analytic sensitivities (``Session.model_gradients``; one reverse-mode pass for all
parameters).  ``pyazr.uncertainty_bands(session, keys=, params=, covariance=, quantity=)`` puts the two together on
the extrapolations of a session made for them, and ``method="central"`` checks it by central differences.

The fit's vector and covariance come from ``Session.fit_with_engine`` (MIGRAD and HESSE, as ``fit_engine.py``), or
from ``--params`` and ``--covariance`` (``param.sav`` and ``covariance.dat`` of ``AZURE2 fit``).  The extrapolations
asked for are put in a trimmed copy of the project (``pyazr.trimmed_model``: only those grids, so the band costs what
was asked); the copy and its outputs go to ``--out`` (default: a new temporary directory, removed at the end), and
its session reads the data files from the project's directory (``Session(..., cwd=)``).

Usage (from the repository's root):

    python3 pyazr/examples/uncertainty_band.py tests/upgrade/T1/T1.azr --key 1 \\
        [--grid "range 0.01 0.6 MeV cm points 30"] [--params output/param.sav --covariance output/covariance.dat] \\
        [--out DIR]
"""
import argparse
import os
import shutil
import sys
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np  # noqa: E402

from pyazr import Project, Session, read_covariance_dat, trimmed_model, uncertainty_bands  # noqa: E402


def band_of(project, here, work, x, C, **which):
    """The S-factor band (MeV b) on the trimmed extrapolations, analytic and by central differences."""
    tmp, outputs = trimmed_model(project, cwd=work, workdir=os.path.join(work, "output"), **which)
    p = Project.open(tmp)
    p.set_config(output_dir=outputs + "/")          # the engine writes there, not beside the project
    p.save()
    with Session(tmp, cwd=here, extrapolations=True) as m:
        b = uncertainty_bands(m, keys=[1], params=x, covariance=C, quantity="sfactor")[1]
        fd = uncertainty_bands(m, keys=[1], params=x, covariance=C, quantity="sfactor", method="central")[1]
    return b, fd


def show(title, b, fd):
    print(f"\n{title}: {len(b.energy)} points, E_cm = {b.energy.min():.4f} to {b.energy.max():.4f} MeV")
    print(f"{'E_cm (MeV)':>11} {'S (keV b)':>12} {'1 sigma (keV b)':>16} {'relative':>9}")
    for k in np.linspace(0, len(b.energy) - 1, 5).astype(int):
        print(f"{b.energy[k]:11.4f} {1e3 * b.value[k]:12.5g} {1e3 * b.sigma[k]:16.4g} {b.relative[k]:9.2%}")
    print(f"central differences against the analytic band: largest relative difference of sigma "
          f"{np.abs(fd.sigma / b.sigma - 1).max():.1e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr) with an extrapolation record")
    ap.add_argument("--key", type=int, default=1, help="the extrapolation, 1-based in file order")
    ap.add_argument("--grid", default="range 0.01 0.6 MeV cm points 30",
                    help="a new grid of the extrapolation's reaction (an extrapolation record's energies)")
    ap.add_argument("--params", default=None, help="a fit's param.sav (with --covariance)")
    ap.add_argument("--covariance", default=None, help="the fit's covariance.dat (with its #columns line)")
    ap.add_argument("--out", default=None, help="keep the trimmed copies and their outputs here")
    args = ap.parse_args()
    here = os.path.dirname(os.path.abspath(args.project))

    with Session(args.project, parameter_file=args.params) as s:
        if args.covariance:
            C, names, _ = read_covariance_dat(args.covariance)
            x = s.parameters.x0
            if list(names) != list(s.parameters.free.names):
                sys.exit(f"{args.covariance}: its columns are not the free parameters of {args.project}")
            print(f"the fit of {args.params}: objective {s.evaluate().objective:.6f}")
        else:
            r = s.fit_with_engine()
            x, C = r.x, r.covariance
            print(f"MIGRAD and HESSE: {r.status}, objective {r.objective:.6f}, {len(x)} free parameters")
    grids = [r for r in Project.open(args.project).records().values() if r.record_kind == "extrapolation"]
    record = grids[args.key - 1]

    work = os.path.abspath(args.out) if args.out else tempfile.mkdtemp(prefix="uncertainty_band.")
    os.makedirs(work, exist_ok=True)
    try:
        b, fd = band_of(args.project, here, work, x, C, keys=[args.key])
        show(f"extrapolation {args.key} ({record['reaction']})", b, fd)
        grid = dict(reaction=record["reaction"], quantity="sigma", unit="b", energies=args.grid)
        b, fd = band_of(args.project, here, work, x, C, grids=[grid])
        show(f"the grid {args.grid!r}", b, fd)
    finally:
        if not args.out:
            shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
