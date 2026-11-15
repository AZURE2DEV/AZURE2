#!/usr/bin/env python3
"""The model of every data set on a fine grid over its range, and what a line through the points misses.

``Session.evaluate(x).data[i].model`` is the model at the data set's points; ``pyazr.model_curve.model_curves(s)``
gives it between them as well (``Session.model_curve(name)`` for one): the engine's evaluation of extrapolation
records made from each data set's record (reaction, quantity, target, beam, detector, ...) on a grid of ``points``
lab energies between the first and the last point, seeded with the level energies and refined where the curve bends
(``grid="refined"``, the default), and at the points themselves, one branch per fixed angle of an excitation
function or per fixed energy of an angular distribution.  The script prints, per data set, the largest model value
at the points and on the curve (a resonance between two points shows as a larger maximum on the curve), and writes
each curve as a table with ``--out``.  A data set the curve cannot be made for (a beam-profile window per point,
points that are not one function of one variable) prints its refusal: what, what is allowed, the way out.

Usage (from the repository's root):

    python3 pyazr/examples/model_curve.py tests/upgrade/T2/T2.azr [--points 400] [--out DIR]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np

from pyazr import Session
from pyazr.model_curve import model_curves


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path, help="the format-3 project (.azr)")
    ap.add_argument("--points", type=int, default=200, help="grid values per branch between the first and last point")
    ap.add_argument("--out", type=Path, default=None, help="a directory for one table per curve")
    args = ap.parse_args()

    with Session(str(args.azr)) as s:
        x = np.asarray(s.parameters.x0, float)
        ev = s.evaluate(x)
        curves = model_curves(s, x=x, points=args.points)
        print(f"{'data set':18s} {'unit':>10} {'points':>6} {'grid':>6} {'max at the points':>18} "
              f"{'max on the curve':>17} {'at E_cm/MeV':>12}")
        for d in ev.data:
            c = curves[d.name]
            if isinstance(c, Exception):
                print(f"{d.name:18s} no curve: {c.what}; allowed: {c.allowed}; {c.fix}")
                continue
            k = int(np.argmax(c.model))
            print(f"{d.name:18s} {c.unit:>10} {d.n:6d} {c.n:6d} {np.max(d.model):18.6g} {c.model[k]:17.6g} "
                  f"{c.energies[k]:12.6g}")
            if args.out:
                args.out.mkdir(parents=True, exist_ok=True)
                (args.out / f"{d.name}_curve.txt").write_text(c.table() + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
