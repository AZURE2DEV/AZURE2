#!/usr/bin/env python3
"""Add/remove resonances and data sets, recalc EC integrals, save a new file.

``pyazr.Project`` edits the two things a model is made of -- the level scheme
(``<levels>``: levels by label or by J^pi and E_x, channels by pair, l and s) and
the data (``dataset`` and ``extrapolation`` records) -- and writes a **new**
project, leaving the original untouched.  This example walks through the full
edit loop:

  1. remove a resonance, add a resonance,
  2. remove a data set, add a data set,
  3. **recalculate the external-capture integrals** -- required whenever the
     data change, because AZURE2 caches them in ``output/intEC.dat`` /
     ``output/intEC.extrap`` keyed on the *grid*, not on the segment selection,
     and silently reuses the stale file otherwise,
  4. save the edited file and run a fresh session from it.

The result is deliberately a wreck (one level gone, one added with tiny widths,
one data set dropped, a new one included) -- the point is the *mechanics* of
editing, not a physical model.

Usage:
    python3 edit_model.py <file.azr> [--out edited.azr]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

from pyazr import Project, Session


def recalc_ec_caches(outdir: Path):
    """Drop AZURE2's cached external-capture integrals.

    They are keyed on the integration grid, so a change to the segment
    selection (or the channel radius) leaves a stale file that AZURE2 would
    otherwise read back silently.  Deleting it forces a recomputation on the
    next session start; they cost time but are always correct.
    """
    for name in ("intEC.dat", "intEC.extrap"):
        f = outdir / name
        if f.exists():
            f.unlink()
            print(f"    deleted stale cache {f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azr", type=Path)
    ap.add_argument("--out", type=Path, default=None,
                    help="where to write the edited model "
                         "(default: <name>_edited.azr next to the input)")
    args = ap.parse_args()

    here = args.azr.parent
    out = args.out or here / f"{args.azr.stem}_edited.azr"

    # -- 1. resonances -------------------------------------------------------
    model = Project.open(args.azr)
    print(f"{args.azr.name}: {len(model.levels())} levels")

    model.remove_level(("1/2+", 20))                          # a background pole
    print("\nremoved level: 1/2+ at 20 MeV")

    model.add_level(1.5, +1, 4.1,                             # a new 3/2+ resonance (E_x, MeV)
                    channels=[dict(pair=1, l=2, s=0.5, width=1000.0),
                              dict(pair=2, l=1, s=0.5, width=0.1)])
    print(f"added level:   {model.levels()[-1].label}")

    # -- 2. data ---------------------------------------------------------------
    print("\ndata sets in the original file:")
    for name, r in model.records().items():
        if r.record_kind == "dataset":
            print(f"  {name:<24} {r['quantity']:<14} {r['reaction']}")

    model.remove("artemov")
    print("\nremoved the data set artemov")

    # add a data set for a data file present but not yet used (roughton.dat):
    # the cross section of 12C(p,g)13N, every capture pair summed.
    model.add_dataset("roughton", reaction=model.reaction(1, -1), quantity="sigma", unit="b",
                      file="data/roughton.dat", columns="E_lab[MeV] _ value error",
                      energy_window="0.30 2.30 MeV", norm="1 fixed")
    print("added the data set roughton (total capture, E=0.3-2.3 MeV)")

    # -- 3. recalculate the EC integrals ---------------------------------------
    # data changed => the cached integrals belong to the old grids.  Delete
    # them so the fresh session below recomputes rather than reuses them.
    print(f"\nrecalculating external-capture integrals (data changed):")
    recalc_ec_caches(here / "output")

    # -- 4. save + run ----------------------------------------------------------
    model.save(out)
    print(f"\nwrote edited model to {out}")

    with Session(str(out), cwd=str(here)) as azr:
        x = np.asarray(azr.parameters.x0, float)
        chi2 = azr.evaluate(x).chi2
        npts = sum(len(e) for e in azr.energies)
        print(f"fresh session on the edited file: {azr.nsegments} segments, "
              f"{npts} data points, chi2 = {chi2:.1f}")
        # EC integrals were recomputed on startup (cache was deleted); asking
        # again here would hit the running cache.  To force a recomputation in
        # a *live* session instead of a fresh one:
        #   azr.recalculate_external_capture()

    # cleanup the variant file
    Path(out).unlink(missing_ok=True)
    print(f"(removed {out.name})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
