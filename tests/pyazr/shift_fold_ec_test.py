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

Then the shift moved from 0 to +0.2, -0.2, +0.05 and +0.02 MeV against the same shift written in the file, fixed (no
reach: planned for exactly that shift) and free (planned around it): the three chi^2 agree to 1e-12 relative.  They
are equal to all digits today, because the planner's nodes are a function of the anchors (the levels) alone and a
re-planned support covers the same nodes as one planned at the shift from the start; a re-plan that depended on where
the shift came from would break this (the director's measurement of the shift kernel's review, 2026-11-12).

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
    from pyazr import Session, Project          # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
    sys.exit(77)

SRC = os.path.join(ROOT, "tests", "13N")
FAR = float(os.environ.get("SHIFT_FOLD_EC_FAR", "0.2"))  # MeV lab

work = tempfile.mkdtemp(prefix="pyazr_shift_ec_")


def session(shift):
    """A session of the model with skowronski_fels's shift record `shift` (format 3), in its own folder."""
    folder = tempfile.mkdtemp(dir=work)
    shutil.copytree(os.path.join(SRC, "data"), os.path.join(folder, "data"))
    os.makedirs(os.path.join(folder, "output"))
    os.makedirs(os.path.join(folder, "checks"))
    p = Project.open(os.path.join(SRC, "13N.azr"))
    for name, r in p.records().items():
        if r.record_kind == "dataset" and name != "skowronski_fels":
            p.set(name, active="no")
    p.add_beam("beam_fels", spread="0.002 MeV lab")
    p.set("skowronski_fels", shift=shift, beam="beam_fels", quadrature="points 20")
    return Session(p.save(os.path.join(folder, "ec.azr")), cwd=folder)


def shift_index(a):
    """The free parameter vector's index of the shift."""
    free = [q for q in a.parameters if not q.fixed]
    return [i for i, q in enumerate(free) if q.kind == "shift"][0]


try:
    with session("offset 0 MeV fit") as a:
        x0 = np.array(a.params_rwa, float)
        k = shift_index(a)
        c0 = np.array([d.chi2 for d in a.evaluate(x0.copy()).data])[0]
        x = x0.copy()
        x[k] = x0[k] + FAR
        cfar = np.array([d.chi2 for d in a.evaluate(x).data])[0]
        cback = np.array([d.chi2 for d in a.evaluate(x0.copy()).data])[0]
    print(f"    chi2 at the file value {c0:.12g}, at +{FAR} MeV {cfar:.12g}, back {cback:.12g}")
    check("finite chi2 beyond the reach", math.isfinite(cfar), cfar)
    check("the shift moves chi2", abs(cfar - c0) > 1e-6 * abs(c0), (c0, cfar))
    check("the file value's chi2 again on return (1e-12)", abs(cback - c0) <= 1e-12 * abs(c0), (c0, cback))

    for far in (0.2, -0.2, 0.05, 0.02):
        with session("offset 0 MeV fit") as a:
            x = np.array(a.params_rwa, float)
            x[shift_index(a)] += far
            moved = np.array([d.chi2 for d in a.evaluate(x).data])[0]
        with session(f"offset {far} MeV fixed") as a:
            fixed = np.array([d.chi2 for d in a.evaluate(np.array(a.params_rwa, float)).data])[0]
        with session(f"offset {far} MeV fit") as a:
            free = np.array([d.chi2 for d in a.evaluate(np.array(a.params_rwa, float)).data])[0]
        print(f"    shift {far:+.2f} MeV: moved from 0 {moved:.12g}, fixed in the file {fixed:.12g}, free in the file "
              f"{free:.12g}")
        check(f"shift {far:+.2f} MeV: moved from 0 = written in the file, fixed and free (1e-12)",
              math.isfinite(fixed) and abs(moved - fixed) <= 1e-12 * abs(fixed)
              and abs(free - fixed) <= 1e-12 * abs(fixed), (moved, fixed, free))
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
