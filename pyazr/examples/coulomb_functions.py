#!/usr/bin/env python3
"""Plot the Coulomb wave functions, penetrability and hard-sphere phase.

These are the quantities that describe the external region: everything outside
the channel radius, where the interaction is assumed to be purely Coulomb.  They
are what the penetrability factors in the level matrix are built from, what sets
the hard-sphere phase shift, and what the external-capture integrals integrate.

AZURE2 has always computed them; `azr.coulomb_functions` makes them askable.
The values follow the run's own configuration, so the same call returns the
accurate Coulomb routine's answer, GSL's (`--gsl-coul`), or the Numerov solution
through a nuclear potential (the hybrid model) -- which is the point of being
able to plot them.

Usage:
    python3 coulomb_functions.py <file.azr> [--pair 1] [--L 0 1 2] [--radius 0]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("this example needs matplotlib: pip install 'pyazr[examples]'")

from pyazr import azure2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--pair", type=int, default=1, help="particle-pair key")
    ap.add_argument("--L", type=int, nargs="*", default=[0, 1, 2])
    ap.add_argument("--radius", type=float, default=0.0,
                    help="fm; 0 uses the pair's own channel radius")
    ap.add_argument("--emin", type=float, default=0.01)
    ap.add_argument("--emax", type=float, default=5.0)
    ap.add_argument("--npoints", type=int, default=400)
    ap.add_argument("--out", type=Path, default=Path("coulomb_functions.png"))
    args = ap.parse_args()

    energies = np.linspace(args.emin, args.emax, args.npoints)

    with azure2(str(args.azr)) as azr:
        pair = azr.pairs[args.pair - 1]
        radius = args.radius or pair.channel_radius
        curves = {L: azr.coulomb_functions(args.pair, energies, L=L,
                                           radius=args.radius)
                  for L in args.L}
        stats = azr.cache_stats()

    fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)

    for L, c in curves.items():
        e = c["energy"]
        axes[0, 0].plot(e, c["F"], label=rf"$\ell={L}$")
        # G diverges as E -> 0 (it is the irregular solution), so log scale.
        axes[0, 1].semilogy(e, np.maximum(c["G"], 1e-300), label=rf"$\ell={L}$")
        axes[1, 0].semilogy(e, np.maximum(c["P"], 1e-300), label=rf"$\ell={L}$")
        axes[1, 1].plot(e, np.degrees(c["delta_hs"]), label=rf"$\ell={L}$")

    axes[0, 0].set_ylabel(r"$F_\ell(\rho)$")
    axes[0, 1].set_ylabel(r"$G_\ell(\rho)$")
    axes[0, 1].set_ylim(bottom=1e-1)
    axes[1, 0].set_ylabel(r"penetrability $P_\ell$")
    axes[1, 1].set_ylabel(r"hard-sphere phase $\delta_\ell^{\rm hs}$ (deg)")
    for ax in axes.ravel():
        ax.set_xlabel(r"$E_{\rm c.m.}$ (MeV)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(f"{args.azr.name}, pair {args.pair}, a = {radius:g} fm")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")

    print(f"\nCoulomb cache after {args.npoints * len(args.L)} evaluations:")
    print(f"  {stats['queries']} queries, {stats['hits']} hits "
          f"({100 * stats['hit_rate']:.1f}%), {stats['entries']} entries "
          f"over {stats['keys']} keys")
    return 0


if __name__ == "__main__":
    sys.exit(main())
