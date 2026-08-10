#!/usr/bin/env python3
"""Add/remove resonances and data segments, recalc EC integrals, save a new file.

``AzrModel`` edits the two things a model is made of -- the level scheme
(``<levels>``) and the data (``<segmentsData>`` / ``<segmentsTest>``) -- and
writes a **new** ``.azr``, leaving the original untouched.  This example walks
through the full edit loop:

  1. remove a resonance, add a resonance (and add a channel to one),
  2. remove a data segment, add a data segment,
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

from pyazr import AzrModel, azure2


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
    model = AzrModel.from_file(args.azr)
    n_before = len(model.levels)
    print(f"{args.azr.name}: {n_before} levels, "
          f"{len(model.find())} levels found")

    removed = model.remove_level(jpi="1/2+", energy=20)      # a background pole
    print(f"\nremoved level: {[str(r) for r in removed]}")

    added = model.add_level(J=1.5, parity=+1, energy=4.1,    # a new 3/2+ resonance
                            channels=[dict(pair=1, L=2, S=0.5, gamma=1000.0),
                                      dict(pair=2, L=1, S=0.5, gamma=0.1)])
    print(f"added level:   {added}")

    # -- 2. data ---------------------------------------------------------------
    print(f"\ndata segments in the original file:")
    from pyazr import datasets
    for s in datasets.SegmentSet.from_file(args.azr):
        print(f"  #{s.key:<2} {s.name:<24} {s.observable:<18} "
              f"E={s.energy_min:g}-{s.energy_max:g}")

    # remove every segment that reads data/artemov.dat
    n_removed = model.remove_data_segments("artemov.dat")
    print(f"\nremoved {n_removed} data segment(s) (artemov)")

    # add a segment for a data file present but not yet used (roughton.dat).
    # total capture: entrance pair 1 -> exit pair 2 (summed).
    model.add_data_segment("data/roughton.dat", entrance=1, exit=2,
                           observable="total-capture",
                           energy_min=0.30, energy_max=2.30,
                           angle_min=0.0, angle_max=180.0,
                           norm=1.0, vary_norm=False, norm_error=5.0)
    print("added 1 data segment (roughton, total capture, E=0.3-2.3)")

    # -- 3. recalculate the EC integrals ---------------------------------------
    # data changed => the cached integrals belong to the old grids.  Delete
    # them so the fresh session below recomputes rather than reuses them.
    print(f"\nrecalculating external-capture integrals (data changed):")
    recalc_ec_caches(here / "output")

    # -- 4. save + run ----------------------------------------------------------
    model.write(out)
    print(f"\nwrote edited model to {out}")

    with azure2(str(out), cwd=str(here)) as azr:
        x = np.asarray(azr.params_rwa, float)
        chi2 = float(np.sum(azr.calculate_chi2_rwa(x)))
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
