#!/usr/bin/env python3
"""The Coulomb wave functions, penetrability, shift function and hard-sphere phase of a channel, from the engine.

These describe the external region: everything outside the channel radius, where the interaction is taken to be
Coulomb only.  The penetrability P_l and shift S_l enter the level matrix, the hard-sphere phase the elastic
amplitude, F_l and G_l the external-capture integrals.  ``Session.channel_functions(pair, E_cm, l=, radius=)``
returns them as the session computes them -- the Coulomb table, the exact routine outside it or the
Numerov solution through a nuclear potential (the hybrid model) -- with c.m. energies in MeV and the phase in rad.

The script prints a table at a few energies and, with ``--plot FILE``, draws F, G, P and the phase (matplotlib).

Usage (from the repository's root):

    python3 pyazr/examples/channel_functions.py tests/upgrade/T1/T1.azr [--pair 1] [--L 0 1 2] [--plot cf.png]
"""
import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np  # noqa: E402

from pyazr import Session  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr)")
    ap.add_argument("--pair", type=int, default=1, help="the particle pair's key")
    ap.add_argument("--L", type=int, nargs="*", default=[0, 1, 2], help="orbital angular momenta")
    ap.add_argument("--radius", type=float, default=0.0, help="channel radius in fm; 0: the pair's own")
    ap.add_argument("--emin", type=float, default=0.05, help="lowest c.m. energy (MeV)")
    ap.add_argument("--emax", type=float, default=3.0, help="highest c.m. energy (MeV)")
    ap.add_argument("--npoints", type=int, default=200)
    ap.add_argument("--plot", default=None, help="draw the four functions into this file (needs matplotlib)")
    args = ap.parse_args()

    energies = np.linspace(args.emin, args.emax, args.npoints)
    with Session(args.project) as s:
        pair = next(q for q in s.pairs if q.key == args.pair)
        curves = {L: s.channel_functions(args.pair, energies, l=L, radius=args.radius) for L in args.L}
    radius = args.radius or pair.channel_radius
    print(f"{os.path.basename(args.project)}, pair {args.pair}, a = {radius:g} fm")
    print(f"{'l':>2} {'E_cm (MeV)':>11} {'F_l':>12} {'G_l':>12} {'P_l':>12} {'S_l':>12} {'phi_hs (deg)':>13}")
    for L, c in curves.items():
        for k in np.linspace(0, len(energies) - 1, 4).astype(int):
            print(f"{L:2d} {c['e_cm'][k]:11.4f} {c['F'][k]:12.5g} {c['G'][k]:12.5g} {c['P'][k]:12.5g} "
                  f"{c['S'][k]:12.5g} {np.degrees(c['delta_hs'][k]):13.5g}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
        for L, c in curves.items():
            e = c["e_cm"]
            axes[0, 0].plot(e, c["F"], label=rf"$\ell={L}$")
            axes[0, 1].semilogy(e, np.maximum(c["G"], 1e-300), label=rf"$\ell={L}$")  # G diverges as E -> 0
            axes[1, 0].semilogy(e, np.maximum(c["P"], 1e-300), label=rf"$\ell={L}$")
            axes[1, 1].plot(e, np.degrees(c["delta_hs"]), label=rf"$\ell={L}$")
        for ax, name in zip(axes.ravel(), (r"$F_\ell(\rho)$", r"$G_\ell(\rho)$", r"penetrability $P_\ell$",
                                           r"hard-sphere phase (deg)")):
            ax.set_ylabel(name)
            ax.set_xlabel(r"$E_{\rm c.m.}$ (MeV)")
            ax.legend(fontsize=8)
        fig.suptitle(f"{os.path.basename(args.project)}, pair {args.pair}, a = {radius:g} fm")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"wrote {args.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
