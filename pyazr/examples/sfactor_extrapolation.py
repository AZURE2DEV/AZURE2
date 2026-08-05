#!/usr/bin/env python3
"""Extrapolate an S factor below the measured range, and read S(0) off it.

The astrophysical S factor removes the Coulomb barrier from the cross section,

    S(E) = sigma(E) E exp(2 pi eta),

which is what makes a low-energy extrapolation meaningful at all: sigma falls
by decades where S is nearly flat.  The extrapolation itself is an ordinary
calculation on a `<segmentsTest>` grid, so the only work is defining the grid
and switching the session into extrapolation mode.

Two things to keep straight, both of which silently give a wrong answer:

  * `add_extrapolation` takes LAB energies, while `calculate_energies` returns
    c.m.  Divide the c.m. range by m_target / (m_beam + m_target).
  * `output/intEC.extrap` caches external-capture integrals for the grid it was
    built on and is reused without checking.  A changed grid means a stale
    cache; delete it.

Usage:
    python3 sfactor_extrapolation.py <file.azr> --entrance 1 --exit 2 \\
        --emin 0.01 --emax 2.0 [--params output/param.sav] [--plot s.pdf]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

from pyazr import AzrModel, azure2


def pair_masses(mdl: AzrModel, pair: int) -> tuple[float, float]:
    """(beam, target) masses of a particle pair, from any channel that uses it.

    The masses are carried on the channel lines of the <levels> block rather
    than in a pair table of their own, so the first channel with this pair
    number answers the question.
    """
    for lvl in mdl.levels:
        for ch in lvl.channels:
            if ch.pair == pair and not ch.is_photon:
                return float(ch.M1), float(ch.M2)
    raise SystemExit(f"no particle channel found for pair {pair}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--entrance", type=int, default=1)
    ap.add_argument("--exit", dest="exit_pair", type=int, default=2)
    ap.add_argument("--emin", type=float, default=0.01, help="c.m. MeV")
    ap.add_argument("--emax", type=float, default=2.0, help="c.m. MeV")
    ap.add_argument("--npts", type=int, default=400)
    ap.add_argument("--params", type=Path, default=None)
    ap.add_argument("--plot", type=Path, default=None)
    ap.add_argument("--keep", action="store_true",
                    help="keep the generated .azr instead of removing it")
    args = ap.parse_args()

    mdl = AzrModel.from_file(args.azr)
    beam, target = pair_masses(mdl, args.entrance)
    cm2lab = target / (beam + target)
    print(f"entrance pair {args.entrance}: m1 = {beam:.5f}, m2 = {target:.5f} u; "
          f"E_cm = {cm2lab:.5f} E_lab")

    step = (args.emax - args.emin) / max(args.npts - 1, 1)
    mdl.clear_extrapolations()
    mdl.add_extrapolation(entrance=args.entrance, exit=args.exit_pair,
                          e_min=args.emin / cm2lab, e_max=args.emax / cm2lab,
                          e_step=step / cm2lab, observable="angle-integrated")
    out = args.azr.with_name("_sfactor.azr")
    mdl.write(out)

    stale = args.azr.parent / "output" / "intEC.extrap"
    if stale.exists():
        stale.unlink()
        print(f"removed stale {stale}")

    with azure2(str(out)) as azr:
        x = np.asarray(azr.params_rwa, float)
        if args.params and args.params.exists():
            saved = np.loadtxt(args.params, usecols=1)
            if saved.shape == x.shape:
                x = saved

        azr.extrap_mode()
        e = np.asarray(azr.calculate_energies(x)[0], float)
        s = np.asarray(azr.calculate_sfactor_rwa(x)[0], float)

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
    # The grid has to live beside the original, because AZURE2 resolves data
    # paths relative to the .azr; do not leave it there.
    if not args.keep:
        out.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
