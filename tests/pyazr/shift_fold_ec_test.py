#!/usr/bin/env python3
"""A folded capture data set with external capture follows its energy shift beyond the reach.

When a shifted folded point's shift leaves the reach its nodes were planned for, the hand-over plans the mesh again
and re-places its sharing group (include/ShiftedFold.h, question 2).  A segment with external capture is re-placed
too: its new sites take their capture amplitudes at their energies from the tables (Site::GetECAmplitudeWithShift),
as a shifted site always has.  Before, such a segment threw in the middle of a fit.

Model: tests/13N with only skowronski_fels active (12C(p,g)13N angle-integrated, external capture on), a Gaussian beam
spread of 2 keV on it (quadrature points 20) and a free constant shift without a prior (the reach is then 3 x the
support's width, about 55 keV c.m.; +0.2 MeV lab is 3.4 reaches).
In one session: chi^2 at the file value, at +0.2 MeV (beyond the reach: the mesh is planned again), and at the file
value again.  Checked: finite chi^2 beyond the reach, and the file value's chi^2 again on return to 1e-12.  The value
at +0.2 MeV is printed for the comparison with the engine before (which moved every sub-point with the shift).

Run from anywhere:  python3 tests/pyazr/shift_fold_ec_test.py
"""
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import numpy as np                         # noqa: E402
    from pyazr import azure2, Project          # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

SRC = os.path.join(ROOT, "tests", "13N")
FAR = float(os.environ.get("SHIFT_FOLD_EC_FAR", "0.2"))  # MeV lab

work = tempfile.mkdtemp(prefix="pyazr_shift_ec_")
try:
    shutil.copytree(os.path.join(SRC, "data"), os.path.join(work, "data"))
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    p = Project.open(os.path.join(SRC, "13N.azr"))
    for name, r in p.records().items():
        if r.record_kind == "dataset" and name != "skowronski_fels":
            p.set(name, active="no")
    p.add_beam("beam_fels", spread="0.002 MeV lab")
    p.set("skowronski_fels", shift="offset 0 MeV fit", beam="beam_fels", quadrature="points 20")
    a = azure2(p.save(os.path.join(work, "ec.azr")), cwd=work)
    x0 = np.array(a.params_rwa, float)
    free = [q for q in a.parameters if not q.fixed]
    k = [i for i, q in enumerate(free) if q.kind == "shift"][0]
    c0 = a.segment_chi2(x0.copy())[0]
    x = x0.copy()
    x[k] = x0[k] + FAR
    cfar = a.segment_chi2(x)[0]
    cback = a.segment_chi2(x0.copy())[0]
    print(f"    chi2 at the file value {c0:.12g}, at +{FAR} MeV {cfar:.12g}, back {cback:.12g}")
    check("finite chi2 beyond the reach", math.isfinite(cfar), cfar)
    check("the shift moves chi2", abs(cfar - c0) > 1e-6 * abs(c0), (c0, cfar))
    check("the file value's chi2 again on return (1e-12)", abs(cback - c0) <= 1e-12 * abs(c0), (c0, cback))
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
