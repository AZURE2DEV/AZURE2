#!/usr/bin/env python3
"""Extrapolate an S factor below the measured range, and read S(0) off it.

The astrophysical S factor removes the Coulomb barrier from the cross section,

    S(E) = sigma(E) E exp(2 pi eta),

which is what makes a low-energy extrapolation meaningful at all: sigma falls
by decades where S is nearly flat.  The extrapolation itself is an ordinary
calculation on a `<segmentsTest>` grid, so the only work is defining the grid
and switching the session into extrapolation mode.

The grid is an ``extrapolation`` record of an edited copy of the project, given in c.m. energies of the entrance
pair (the engine takes them to the lab); ``evaluate(x).data[i].energies`` are c.m. as always.  The copy is written to
``--out`` (default: a new temporary directory, printed); its session reads the project's files from the project's
directory (``Session(..., cwd=)``).

Usage (from the repository's root):

    python3 pyazr/examples/sfactor_extrapolation.py tests/upgrade/T1/T1.azr --entrance 1 --exit 2 \\
        [--emin 0.01] [--emax 2.0] [--params output/param.sav] [--plot s.pdf] [--out DIR]
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np

from pyazr import Project, Session


def entrance_pair(mdl: Project, pair: int):
    """The particle pair ``pair`` of the project (masses in u, charges)."""
    for q in mdl.pairs():
        if q.key == pair and q.kind == 0:
            return q
    raise SystemExit(f"no particle pair with key {pair}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path, help="the format-3 project (.azr)")
    ap.add_argument("--entrance", type=int, default=1)
    ap.add_argument("--exit", dest="exit_pair", type=int, default=2)
    ap.add_argument("--emin", type=float, default=0.01, help="c.m. MeV")
    ap.add_argument("--emax", type=float, default=2.0, help="c.m. MeV")
    ap.add_argument("--npts", type=int, default=400)
    ap.add_argument("--params", type=Path, default=None, help="a parameter file (param.sav) of a fit")
    ap.add_argument("--plot", type=Path, default=None, help="draw S(E) into this file (needs matplotlib)")
    ap.add_argument("--out", default=None, help="the directory of the edited copy (default: a new temporary one)")
    args = ap.parse_args()

    mdl = Project.open(args.azr)
    q = entrance_pair(mdl, args.entrance)
    beam, target = q.light_mass, q.heavy_mass
    cm2lab = target / (beam + target)
    print(f"entrance pair {args.entrance}: m1 = {beam:.5f}, m2 = {target:.5f} u; "
          f"E_cm = {cm2lab:.5f} E_lab")

    step = (args.emax - args.emin) / max(args.npts - 1, 1)
    for name, r in mdl.records().items():
        if r.record_kind == "extrapolation":
            mdl.remove(name)
    # sigma on a c.m. grid of the entrance pair (the engine takes it to lab energies)
    mdl.add_extrapolation("grid", reaction=mdl.reaction(args.entrance, args.exit_pair), quantity="sigma", unit="b",
                          energies=f"range {args.emin!r} {args.emax!r} MeV cm step {step!r}")
    out_dir = args.out or tempfile.mkdtemp(prefix="sfactor_extrapolation.")
    os.makedirs(out_dir, exist_ok=True)
    out = mdl.save(os.path.join(out_dir, args.azr.name))
    print(f"the project with the grid: {out}")

    # a session for the extrapolations, with the fit's parameter file (fixed values included) when given
    with Session(out, extrapolations=True, parameter_file=args.params, cwd=str(args.azr.parent)) as azr:
        x = np.asarray(azr.parameters.x0, float)

        grid = azr.evaluate(x).data[0]          # the first extrapolation
        e = grid.energies
        s = grid.model * grid.extra["sfactor_conversion"]   # MeV b

    o = np.argsort(e)
    e, s = e[o], s[o]
    print(f"{len(e)} points, {e[0]:.4f}-{e[-1]:.4f} MeV c.m.")
    print(f"S({e[0]:.3f} MeV) = {s[0] * 1e3:.4g} keV b")
    print(f"peak S = {s.max() * 1e3:.4g} keV b at {e[np.argmax(s)]:.4f} MeV")

    # S(0) by a quadratic through the three lowest points -- the usual
    # convention, and honest only if the lowest point is well below any
    # resonance.
    if len(e) >= 3:
        c = np.polyfit(e[:3], s[:3], 2)
        print(f"S(0) = {np.polyval(c, 0.0) * 1e3:.4g} keV b "
              f"(quadratic through the three lowest points)")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.0, 3.0))
        ax.plot(e, s * 1e3, lw=1.2)
        ax.set_yscale("log")
        ax.set_xlabel(r"$E_{\rm c.m.}$ (MeV)")
        ax.set_ylabel(r"$S$ (keV b)")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=200)
        print(f"wrote {args.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
