#!/usr/bin/env python3
"""Decompose a cross section into level, interference and external contributions.

"Which resonance does this dataset constrain" is a question about a
decomposition, and the pieces are not independent: the total is not a sum of
level contributions, because levels of the same J^pi interfere.  Four
evaluations at the *fitted* parameters, with no refitting, separate them:

    full          everything
    off(L)        L's reduced widths zeroed -- everything except L
    only(L)       every *other* level zeroed -- L on its own background
    bare          every free width zeroed -- the non-resonant background

from which

    full - off    = what L does, resonance plus its interference
    only - bare   = the bare resonance on the background
    difference    = the interference alone

The interference is block-diagonal in J^pi -- only levels of the same spin and
parity interfere -- which is a check on any decomposition you compute this way.

Usage:
    python3 decompose_cross_section.py <file.azr> [--segment 1]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

from pyazr import azure2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--segment", type=int, default=1, help="1-based")
    ap.add_argument("--params", type=Path, default=None)
    args = ap.parse_args()

    with azure2(str(args.azr)) as azr:
        x = np.asarray(azr.params_rwa, float)
        if args.params:
            # param.sav has one line per parameter, fixed ones included, while
            # params_rwa holds only the free ones -- so a length comparison
            # never matches and would silently leave x at the .azr's values.
            # best_fit_params does the mapping (and the extrapolation-mode
            # case) properly.
            from pyazr.bands import best_fit_params
            x = np.asarray(best_fit_params(azr, str(args.params)), float)

        i = args.segment - 1
        energies = np.asarray(azr.calculate_energies(x)[i], float)
        full = np.asarray(azr.calculate_rwa(x)[i], float)

        # The background: every free width off, leaving hard-sphere scattering
        # and external capture.
        bare_x = x.copy()
        for p in azr.parameters.widths:
            if p.free_index is not None:
                bare_x[p.free_index] = 0.0
        bare = np.asarray(azr.calculate_rwa(bare_x)[i], float)

        print(f"{args.azr.name}, segment {args.segment}: "
              f"{len(energies)} points, {energies.min():.3f}-{energies.max():.3f} MeV")
        print(f"\n{'level':>16} {'max |full-off|':>15} {'max |only-bare|':>16} "
              f"{'max |interf.|':>14}")

        for key in azr.physical_levels():
            # Select by the LevelKey itself.  A (J^pi, energy) pair is
            # ambiguous when a level's energy is fixed and reported as None,
            # which is exactly the case for bound states and background poles.
            off = np.asarray(azr.calculate_rwa(
                azr.without_level(x, level=key))[i], float)
            only = np.asarray(azr.calculate_rwa(
                azr.only_level(x, level=key))[i], float)
            contrib = full - off
            bare_res = only - bare
            interference = contrib - bare_res
            print(f"{str(key):>16} {np.max(np.abs(contrib)):15.3e} "
                  f"{np.max(np.abs(bare_res)):16.3e} "
                  f"{np.max(np.abs(interference)):14.3e}")

        print("\nThe third column is interference and should be non-zero only "
              "where\nanother level of the same J^pi overlaps.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
