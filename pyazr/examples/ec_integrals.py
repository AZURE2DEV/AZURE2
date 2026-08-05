#!/usr/bin/env python3
"""Plot the external-capture radial integrals, and show what caching them buys.

External capture is the non-resonant part of a radiative-capture cross section:
the particle radiates from the external region straight into a bound state,
without forming a compound resonance.  Its amplitude is a radial integral from
the channel radius outward,

    R^EC = C * int_{a_c}^{inf} dr r^L [F_l cos(delta_hs) + G_l sin(delta_hs)]
                                       W_{-eta,lf+1/2}(2 k_f r),

over Coulomb functions and a Whittaker function.  There is one of these for
every (l_i, l_f, s_i, s_f, L, E/M) pathway the compound nucleus allows -- for a
typical capture model, tens of them -- and each is an adaptive quadrature over
expensive special functions.  They are the single most costly thing in a capture
calculation, and they are why AZURE3 caches Coulomb functions at all.

This script asks for them over an energy grid, plots them, and then asks a
second time to show the cache doing its job.

Usage:
    python3 ec_integrals.py <file.azr> [--pair 1] [--top 8]
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("this example needs matplotlib: pip install 'pyazr[examples]'")

from pyazr import azure2


def label(p) -> str:
    return (rf"$\ell_i={p['li']},\ \ell_f={p['lf']}$, "
            rf"$s_i={p['si']:g},\ s_f={p['sf']:g}$, "
            rf"{p['radiation']}{p['multipolarity']}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--pair", type=int, default=1, help="entrance particle-pair key")
    ap.add_argument("--emin", type=float, default=0.05)
    ap.add_argument("--emax", type=float, default=3.0)
    ap.add_argument("--npoints", type=int, default=120)
    ap.add_argument("--top", type=int, default=8,
                    help="plot the N pathways with the largest integrals")
    ap.add_argument("--out", type=Path, default=Path("ec_integrals.png"))
    args = ap.parse_args()

    energies = np.linspace(args.emin, args.emax, args.npoints)

    with azure2(str(args.azr)) as azr:
        before = azr.cache_stats()

        t0 = time.perf_counter()
        paths = azr.ec_integrals(args.pair, energies)
        t_cold = time.perf_counter() - t0
        mid = azr.cache_stats()

        # Exactly the same request again.  Every Coulomb function it needs is
        # now memoized at exactly the energies and radii it asks for.
        t0 = time.perf_counter()
        azr.ec_integrals(args.pair, energies)
        t_warm = time.perf_counter() - t0
        after = azr.cache_stats()

    if not paths:
        sys.exit(f"no external-capture pathways from pair {args.pair} "
                 f"in {args.azr.name}")

    print(f"{len(paths)} external-capture pathways, {args.npoints} energies")
    print(f"  first pass  {t_cold:6.2f} s")
    print(f"  second pass {t_warm:6.2f} s   ({t_cold / max(t_warm, 1e-9):.1f}x faster)")
    print(f"\nCoulomb cache")
    for name, s in (("before", before), ("after first pass", mid),
                    ("after second", after)):
        print(f"  {name:17s} {s['queries']:9d} queries  {s['hits']:9d} hits  "
              f"({100 * s['hit_rate']:5.1f}%)  {s['entries']:7d} entries  "
              f"{s['disabled_keys']} keys given up on")

    order = np.argsort([-np.max(np.abs(p["value"])) for p in paths])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    for k in order[:args.top]:
        p = paths[k]
        ax1.semilogy(p["energy"], np.abs(p["value"]), label=label(p))
        ax2.plot(p["energy"], np.degrees(np.angle(p["value"])))

    ax1.set_xlabel(r"$E_{\rm c.m.}$ (MeV)")
    ax1.set_ylabel(r"$|R^{\rm EC}|$")
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=7)
    ax2.set_xlabel(r"$E_{\rm c.m.}$ (MeV)")
    ax2.set_ylabel(r"$\arg R^{\rm EC}$ (deg)")
    ax2.grid(alpha=0.3)
    fig.suptitle(f"{args.azr.name}: external-capture integrals from pair {args.pair}")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
