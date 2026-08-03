"""Angular distributions: Legendre coefficients at energies you choose.

AZURE2 expresses an angular distribution as

    W(theta) = sum_k a_k P_k(cos theta)

and computes the a_k only for segments declared as angular distributions in the
.azr file. Two ways to get at them:

  * ``azure2.calculate_angular_dists`` -- for the grids a model already
    declares. No temporary files, works on a live instance.

  * ``pyazr.angular_distribution`` -- for arbitrary energies. It writes a
    temporary model requesting exactly those energies, evaluates it, and cleans
    up.

Usage:
    python angular_distribution.py <model.azr> [--entrance 1] [--exit 2]
                                   [--order 4] [--emin 0.05] [--emax 1.0] [-n 10]

The a_0 coefficient is the normalisation, so a distribution with a_0 = 1 and
everything else zero is isotropic -- which is what a single s-wave resonance
gives. Anisotropy shows up as the higher orders growing away from zero.
"""

import argparse
import os
import sys

import numpy as np

from pyazr import angular_distribution, azure2


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", help="the .azr model")
    ap.add_argument("--entrance", type=int, default=1, help="entrance pair key")
    ap.add_argument("--exit", type=int, default=2, dest="exit_pair",
                    help="exit pair key")
    ap.add_argument("--order", type=int, default=4,
                    help="highest Legendre order (default 4)")
    ap.add_argument("--emin", type=float, default=0.05,
                    help="lowest LAB energy, MeV (default 0.05)")
    ap.add_argument("--emax", type=float, default=1.0,
                    help="highest LAB energy, MeV (default 1.0)")
    ap.add_argument("-n", "--num", type=int, default=8,
                    help="number of energies (default 8)")
    ap.add_argument("--plot", action="store_true",
                    help="plot the coefficients against energy (needs matplotlib)")
    args = ap.parse_args()

    if not os.path.isfile(args.azr):
        sys.exit(f"no such file: {args.azr}")

    # --- what the model already declares ------------------------------------
    # Cheap, and it tells you whether the model has angular-distribution
    # segments of its own before you go making new ones.
    with azure2(args.azr) as m:
        x = np.asarray(m.params_rwa, float)
        declared = m.calculate_angular_dists_rwa(x)
    with_coeffs = sum(1 for seg in declared for row in seg if row.size)
    print(f"{args.azr}: {len(declared)} data segments, "
          f"{with_coeffs} points carrying angular-distribution coefficients")
    if with_coeffs == 0:
        print("  (none declared as angular distributions -- expected for most models)")

    # --- coefficients at energies of our choosing ---------------------------
    energies = np.linspace(args.emin, args.emax, args.num)
    print(f"\nRequesting order {args.order} coefficients for pair "
          f"{args.entrance} -> {args.exit_pair} at {args.num} energies "
          f"({args.emin} to {args.emax} MeV lab):\n")

    e_cm, coeffs = angular_distribution(
        args.azr, energies, entrance=args.entrance, exit=args.exit_pair,
        order=args.order)

    header = "  E_lab      E_cm   " + "".join(f"     a_{k}" for k in range(args.order + 1))
    print(header)
    print("  " + "-" * (len(header) - 2))
    for e_lab, ecm, row in zip(energies, e_cm, coeffs):
        cells = "".join(f"{v:8.4f}" for v in row)
        print(f"  {e_lab:6.3f} {ecm:9.5f}  {cells}")

    # a_0 is the normalisation; the rest measure how far from isotropic it is.
    anisotropy = np.nanmax(np.abs(coeffs[:, 1:]), axis=1)
    print(f"\n  largest |a_k>0| over the range: {np.nanmax(anisotropy):.4g}")
    if np.nanmax(anisotropy) < 1e-6:
        print("  -> isotropic throughout: a single s-wave resonance dominates,")
        print("     or no higher partial wave contributes at these energies.")

    if args.plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for k in range(coeffs.shape[1]):
            ax.plot(e_cm, coeffs[:, k], marker="o", label=f"$a_{k}$")
        ax.set_xlabel("c.m. energy (MeV)")
        ax.set_ylabel("Legendre coefficient")
        ax.set_title(f"Angular distribution, pair {args.entrance} "
                     f"$\\rightarrow$ {args.exit_pair}")
        ax.legend()
        ax.grid(alpha=.3)
        fig.tight_layout()
        out = "angular_distribution.png"
        fig.savefig(out, dpi=150)
        print(f"\n  wrote {out}")


if __name__ == "__main__":
    main()
