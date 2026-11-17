"""Cheap regression tests for the two gaps that let a level-search gate fail twice.

1. template_variant writing into a SUBDIRECTORY must link the base model's
   data/ beside the variant -- the existing 13N test writes its variant beside
   the base model, where the link already exists, so it never covered this.
2. grid_coverage must say, as pure geometry, whether the template grid can
   represent a given level at all.  The first gate run spent 77 minutes
   discovering something this answers in milliseconds.

No fitting, no engine session: runs in a second.
"""
import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(1, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if not _PYAZR3:
    print("SKIP: a pyazr 3 test; this process runs pyazr 2 (the pyazr 2 version is test_levelscan_grid_13n.py)")
    sys.exit(77)

from rmfit.levelscan import (_patterns, grid_coverage, partition_fractions,
                             pattern_distance, template_variant)

AZURE2_ROOT = os.environ.get("AZURE2_ROOT") or os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
BASE = os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr")
fails = []


def check(cond, msg):
    print(("  ok    " if cond else "  FAIL  ") + msg)
    if not cond:
        fails.append(msg)


def test_variant_links_data_in_a_subdirectory():
    print("variant written into a subdirectory")
    tmp = tempfile.mkdtemp(prefix="lsgrid_")
    try:
        wd = os.path.join(tmp, "scan_X")
        var = template_variant(BASE, os.path.join(wd, "templates.azr"), ex=9.0,
                               jpis=["3/2-"], output_dir=os.path.join(wd, "output"))
        check(os.path.exists(var["path"]), "the variant file is written")
        link = os.path.join(wd, "data")
        check(os.path.islink(link), "data/ is linked beside the variant")
        check(os.path.exists(link), "the data link resolves")
        check(os.path.samefile(link, os.path.join(os.path.dirname(BASE), "data")),
              "it points at the base model's data/")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_variant_records_partition_priors():
    print("priors recorded from the group's existing levels")
    tmp = tempfile.mkdtemp(prefix="lsgrid_")
    try:
        var = template_variant(BASE, os.path.join(tmp, "sub", "t.azr"), ex=9.0,
                               output_dir=os.path.join(tmp, "sub", "out"))
        withp = [j for j, i in var["levels"].items() if i.get("priors")]
        check(len(withp) >= 3, f"at least three J^pi groups carry priors (got {len(withp)})")
        info = var["levels"][withp[0]]
        pr = info["priors"][0]
        n_part = sum(1 for ph in info["photon"] if not ph)
        check(len(pr["fractions"]) == n_part, "a prior spans the particle channels")
        check(abs(sum(abs(f) for f in pr["fractions"]) - 1.0) < 1e-9, "its fractions are normalised")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_pattern_distance_symmetries():
    print("pattern_distance")
    f = [0.5, 0.3, 0.2]
    check(pattern_distance(f, [1, 1, -1], f, [1, 1, -1]) < 1e-12, "identical pattern has distance 0")
    check(pattern_distance(f, [1, 1, -1], f, [-1, -1, 1]) < 1e-12,
          "a global sign flip is free (it is an exact symmetry of a level)")
    d_rel = pattern_distance(f, [1, 1, -1], f, [1, -1, -1])
    check(d_rel > 0.1, "a relative sign change is not free")
    d_shape = pattern_distance([1.0, 0.0, 0.0], [1, 1, 1], [0.0, 0.0, 1.0], [1, 1, 1])
    check(d_shape > 0.9, "orthogonal shapes are far apart")


def test_grid_coverage_sees_a_width_out_of_reach():
    print("grid_coverage")
    # a level whose partition is nearly even but whose total is far above the grid
    widths = [70000.0, 38000.0, -1900.0, 80000.0, 66000.0]
    pats = _patterns(len(widths), 8)
    near = grid_coverage(widths, pats, [9e3, 30e3, 90e3])
    far = grid_coverage(widths, pats, [s * 30e3 for s in (0.3, 0.6, 1.2, 2.5, 5, 10)])
    print(f"        old totals penalty {near['penalty']:.3f} (ratio {near['width_ratio']:.2f});"
          f" wide totals penalty {far['penalty']:.3f} (ratio {far['width_ratio']:.2f})")
    check(near["penalty"] > far["penalty"], "a grid that reaches the true width scores better")
    check(far["width_ratio"] > 0.5 and far["width_ratio"] < 2.0,
          "the wide grid lands within a factor two of the true total width")
    check(near["distance"] < 0.3, "the shape itself was never the problem for this level")
    f, s = partition_fractions(widths)
    check(abs(float(np.sum(f)) - 1.0) < 1e-12 and list(s) == [1, 1, -1, 1, 1],
          "partition_fractions splits magnitude and sign")


if __name__ == "__main__":
    test_variant_links_data_in_a_subdirectory()
    test_variant_records_partition_priors()
    test_pattern_distance_symmetries()
    test_grid_coverage_sees_a_width_out_of_reach()
    print(("\nFAILED: " + "; ".join(fails)) if fails else "\nall grid/variant checks passed")
    raise SystemExit(1 if fails else 0)
