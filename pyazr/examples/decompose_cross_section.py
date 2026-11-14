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

In an angle-integrated reaction cross section the interference is block-diagonal in J^pi -- only levels of the
same spin and parity interfere -- which is a check on a decomposition computed this way; in a differential cross
section, and in elastic scattering with its Coulomb amplitude, every amplitude interferes.

The same evaluations probe what one level does, with no file written: ``--level LABEL`` prints full, off and only
at every fifth point of the segment.  To keep a level switched off in a project, edit the file
(``Project.set_level(label, active=False)``, as ``edit_model.py`` does).

Usage (from the repository's root):

    python3 pyazr/examples/decompose_cross_section.py tests/upgrade/T1/T1.azr [--segment 1] [--level "1/2+ 2.365"]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np

from pyazr import Session


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path, help="the format-3 project (.azr)")
    ap.add_argument("--segment", type=int, default=1, help="the data set, 1-based among the active ones")
    ap.add_argument("--params", type=Path, default=None, help="a parameter file (param.sav) to take")
    ap.add_argument("--level", default=None, help="a level's label (Project.levels()) to print point by point")
    args = ap.parse_args()

    with Session(str(args.azr), parameter_file=args.params) as azr:
        x = np.asarray(azr.parameters.x0, float)

        i = args.segment - 1
        data = azr.evaluate(x).data[i]
        energies = np.asarray(data.energies, float)
        full = np.asarray(data.model, float)

        # The background: every free width off, leaving hard-sphere scattering
        # and external capture.
        bare_x = x.copy()
        P = azr.parameters
        bare_x = P.put(bare_x, P.of_kind("width").free, 0.0)
        bare = np.asarray(azr.evaluate(bare_x).data[i].model, float)

        print(f"{args.azr.name}, segment {args.segment} ({data.name}, {data.quantity} in {data.unit}): "
              f"{len(energies)} points, {energies.min():.3f}-{energies.max():.3f} MeV (c.m.)")
        print(f"\n{'level':>16} {'max |full-off|':>15} {'max |only-bare|':>16} "
              f"{'max |interf.|':>14}   ({data.unit})")

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

        print("\nThe third column is the interference: in an angle-integrated reaction cross section only levels of "
              "the same J^pi\ninterfere; in a differential or elastic one every amplitude does (Coulomb and "
              "hard-sphere scattering too).")

        if args.level:
            key = args.level
            off = azr.evaluate(P.put(x, P.of_level(key).of_kind("width").free, 0.0)).data[i].model
            only = azr.evaluate(P.put(x, [w for w in widths if w.level_label != key], 0.0)).data[i].model
            print(f"\nthe level {key}, {data.unit}:")
            print(f"{'E_cm (MeV)':>11} {'full':>12} {'off':>12} {'only':>12}")
            for k in range(0, len(energies), 5):
                print(f"{energies[k]:11.4f} {full[k]:12.5g} {off[k]:12.5g} {only[k]:12.5g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
