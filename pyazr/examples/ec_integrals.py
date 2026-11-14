#!/usr/bin/env python3
"""The external-capture radial integrals of a capture model, pathway by pathway, from the engine.

External capture is the non-resonant part of a radiative-capture cross section: the particle radiates from the
external region straight into a bound state.  Its amplitude is a radial integral from the channel radius outward,

    R^EC = C * int_{a_c}^{inf} dr r^L [F_l cos(phi_hs) + G_l sin(phi_hs)] W_{-eta,lf+1/2}(2 k_f r),

one per (l_i, l_f, s_i, s_f, multipolarity, E/M) pathway the compound nucleus allows.
``Session.ec_integrals(pair, E_cm)`` returns them as dicts (quantum numbers, the complex ``value`` per energy and its
``unit``, fm^(L+1); the final state is not part of the engine's answer: pathways to two final states with the same
quantum numbers are told apart only by their position in the list, printed as #k).  The script prints the largest
pathways at a few energies, the time of two identical requests, and with ``--plot FILE`` draws |R^EC| and its phase
(matplotlib).

Usage (from the repository's root):

    python3 pyazr/examples/ec_integrals.py tests/upgrade/T1/T1.azr [--pair 1] [--top 4] [--plot ec.png]
"""
import argparse
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np  # noqa: E402

from pyazr import Session  # noqa: E402


def label(p):
    return (f"l_i={p['li']} l_f={p['lf']} s_i={p['si']:g} s_f={p['sf']:g} {p['radiation']}{p['multipolarity']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="the format-3 project (.azr) with external capture")
    ap.add_argument("--pair", type=int, default=1, help="the entrance particle pair's key")
    ap.add_argument("--emin", type=float, default=0.05, help="lowest c.m. energy (MeV)")
    ap.add_argument("--emax", type=float, default=3.0, help="highest c.m. energy (MeV)")
    ap.add_argument("--npoints", type=int, default=60)
    ap.add_argument("--top", type=int, default=4, help="the N pathways with the largest integrals")
    ap.add_argument("--plot", default=None, help="draw the N largest into this file (needs matplotlib)")
    args = ap.parse_args()

    energies = np.linspace(args.emin, args.emax, args.npoints)
    with Session(args.project) as s:
        t0 = time.perf_counter()
        paths = s.ec_integrals(args.pair, energies)
        t1 = time.perf_counter()
        s.ec_integrals(args.pair, energies)            # the same request again
        t2 = time.perf_counter()
    if not paths:
        sys.exit(f"no external-capture pathway from pair {args.pair} in {args.project}")

    print(f"{len(paths)} external-capture pathways from pair {args.pair}, {args.npoints} energies; "
          f"first request {t1 - t0:.2f} s, second {t2 - t1:.2f} s")
    order = np.argsort([-np.max(np.abs(p["value"])) for p in paths])[:args.top]
    cols = np.linspace(0, len(energies) - 1, 4).astype(int)
    print(f"\n{'pathway':<36} {'unit':<6}" + "".join(f"  |R| at {energies[k]:.2f} MeV" for k in cols))
    for k in order:
        p = paths[k]
        print(f"#{k:<3d}{label(p):<32} {p['unit']:<6}" + "".join(f"  {abs(p['value'][c]):18.5g}" for c in cols))

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        for k in order:
            p = paths[k]
            ax1.semilogy(p["energy"], np.abs(p["value"]), label=f"#{k} " + label(p))
            ax2.plot(p["energy"], np.degrees(np.angle(p["value"])))
        ax1.set_ylabel(r"$|R^{\rm EC}|$")
        ax2.set_ylabel(r"$\arg R^{\rm EC}$ (deg)")
        for ax in (ax1, ax2):
            ax.set_xlabel(r"$E_{\rm c.m.}$ (MeV)")
        ax1.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"wrote {args.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
