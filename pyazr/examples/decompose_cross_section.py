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

from pyazr import Session


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--segment", type=int, default=1, help="1-based")
    ap.add_argument("--params", type=Path, default=None)
    args = ap.parse_args()

    with Session(str(args.azr), parameter_file=args.params) as azr:
        x = np.asarray(azr.parameters.x0, float)

        i = args.segment - 1
        energies = np.asarray(azr.evaluate(x).data[i].energies, float)
        full = np.asarray(azr.evaluate(x).data[i].model, float)

        # The background: every free width off, leaving hard-sphere scattering
        # and external capture.
        bare_x = x.copy()
        P = azr.parameters
        bare_x = P.put(bare_x, P.of_kind("width").free, 0.0)
        bare = np.asarray(azr.evaluate(bare_x).data[i].model, float)

        print(f"{args.azr.name}, segment {args.segment}: "
              f"{len(energies)} points, {energies.min():.3f}-{energies.max():.3f} MeV")
        print(f"\n{'level':>16} {'max |full-off|':>15} {'max |only-bare|':>16} "
              f"{'max |interf.|':>14}")

        widths = P.of_kind("width").free
        for e in P.of_kind("energy"):
            # A level by its label (level_label): a (J^pi, energy) pair is
            # ambiguous when a level's energy is fixed and reported as None,
            # which is exactly the case for bound states and background poles.
            key = e.level_label
            off = np.asarray(azr.evaluate(
                P.put(x, P.of_level(key).of_kind("width").free, 0.0)).data[i].model, float)
            only = np.asarray(azr.evaluate(
                P.put(x, [w for w in widths if w.level_label != key], 0.0)).data[i].model, float)
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
