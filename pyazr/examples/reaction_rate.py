#!/usr/bin/env python3
"""A thermonuclear reaction rate from an extrapolated cross section.

AZURE2's own mode 5 writes `reactionrates.out`; this does the same integral in
Python, which is the useful form when the rate has to be recomputed for every
sample of a posterior, or split by level, or differentiated.

    N_A <sigma v> = N_A (8 / pi mu)^(1/2) (kT)^(-3/2)
                    int_0^inf sigma(E) E exp(-E / kT) dE

with E and kT in the centre of mass.  Two practical points:

  * The integrand is the Gamow window -- a narrow product of a rising
    penetrability and a falling Boltzmann factor -- so the grid must resolve it
    at the lowest temperature of interest, not just at the highest.  The script
    reports the window it actually integrated over and warns if the grid ends
    inside it.
  * A resonance narrower than the grid spacing is simply missed.  The default
    grid here is uniform; for a real evaluation, refine it around every
    resonance in the range, or use mode 5, which does that itself.

Usage:
    python3 reaction_rate.py <file.azr> --entrance 1 --exit 2 \\
        --emin 0.005 --emax 3.0 --npts 3000 [--params output/param.sav]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

from pyazr import Project, Session

NA = 6.02214076e23            # 1/mol
KB = 8.617333262e-2           # MeV / GK
AMU = 931.49410242            # MeV
C_CM_S = 2.99792458e10        # cm/s
BARN_CM2 = 1e-24

T9 = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0,
               3.0, 5.0, 10.0])


def entrance_pair(mdl: Project, pair: int):
    """The particle pair ``pair`` of the project (masses in u, charges)."""
    for q in mdl.pairs():
        if q.key == pair and q.kind == 0:
            return q
    raise SystemExit(f"no particle pair with key {pair}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--entrance", type=int, default=1)
    ap.add_argument("--exit", dest="exit_pair", type=int, default=2)
    ap.add_argument("--emin", type=float, default=0.005, help="c.m. MeV")
    ap.add_argument("--emax", type=float, default=3.0, help="c.m. MeV")
    ap.add_argument("--npts", type=int, default=3000)
    ap.add_argument("--params", type=Path, default=None)
    ap.add_argument("--keep", action="store_true",
                    help="keep the generated .azr instead of removing it")
    args = ap.parse_args()

    mdl = Project.open(args.azr)
    q = entrance_pair(mdl, args.entrance)
    m1, m2 = q.light_mass, q.heavy_mass
    mu = m1 * m2 / (m1 + m2)          # amu
    cm2lab = m2 / (m1 + m2)

    step = (args.emax - args.emin) / max(args.npts - 1, 1)
    for name, r in mdl.records().items():
        if r.record_kind == "extrapolation":
            mdl.remove(name)
    # sigma on a c.m. grid of the entrance pair (the engine takes it to lab energies)
    mdl.add_extrapolation("grid", reaction=mdl.reaction(args.entrance, args.exit_pair), quantity="sigma", unit="b",
                          energies=f"range {args.emin!r} {args.emax!r} MeV cm step {step!r}")
    out = args.azr.with_name("_rate.azr")
    mdl.save(out)

    stale = args.azr.parent / "output" / "intEC.extrap"
    if stale.exists():
        stale.unlink()

    with Session(str(out)) as azr:
        x = np.asarray(azr.params_rwa, float)
        if args.params:
            # param.sav has one line per parameter, fixed ones included, while
            # params_rwa holds only the free ones -- so a length comparison
            # never matches and would silently leave x at the .azr's values.
            # best_fit_params does the mapping (and the extrapolation-mode
            # case) properly.
            from pyazr.bands import best_fit_params
            x = np.asarray(best_fit_params(azr, str(args.params)), float)
        azr.extrap_mode()
        e = np.asarray(azr.calculate_energies(x)[0], float)
        sig = np.asarray(azr.calculate_rwa(x)[0], float)

    o = np.argsort(e)
    e, sig = e[o], sig[o]
    print(f"mu = {mu:.5f} u; sigma on {len(e)} points, "
          f"{e[0]:.4f}-{e[-1]:.4f} MeV, spacing {np.diff(e).mean()*1e3:.2f} keV")

    v_pref = np.sqrt(8.0 / (np.pi * mu * AMU)) * C_CM_S   # cm/s, times sqrt(MeV)
    print(f"\n{'T9':>7} {'kT (MeV)':>10} {'NA<sv>':>13} {'E0 (MeV)':>10} "
          f"{'window':>18}")
    for t9 in T9:
        kt = KB * t9
        w = sig * e * np.exp(-e / kt) * BARN_CM2
        rate = NA * v_pref * kt ** -1.5 * np.trapezoid(w, e)

        # Gamow peak and its 1/e width, for reporting the window integrated.
        eta_c = 0.1574854 * q.light[0] * q.heavy[0] * np.sqrt(mu)
        e0 = (0.5 * eta_c * np.pi * kt) ** (2.0 / 3.0) if eta_c else np.nan
        dw = 4.0 * np.sqrt(e0 * kt / 3.0) if eta_c else np.nan
        flag = "  <- grid ends inside" if e[-1] < e0 + dw / 2 else ""
        print(f"{t9:7.2f} {kt:10.4f} {rate:13.4e} {e0:10.4f} "
              f"{e0-dw/2:7.3f}-{e0+dw/2:.3f}{flag}")

    print("\nUnits: cm^3 mol^-1 s^-1.  Compare against AZURE2 mode 5, which "
          "refines\nthe grid around each resonance; a rate that disagrees at "
          "low T9 usually\nmeans this uniform grid missed a narrow one.")
    # The grid has to live beside the original, because AZURE2 resolves data
    # paths relative to the .azr; do not leave it there.
    if not args.keep:
        out.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
