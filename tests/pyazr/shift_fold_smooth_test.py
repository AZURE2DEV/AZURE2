#!/usr/bin/env python3
"""chi^2 of a folded data set is smooth in its energy shift: the kernel moves on the fixed nodes.

A folded point with an energy shift is folded on the nodes of the shared mesh, which do not move with the shift; the
kernel moves with the point and is clipped at its moved support (include/ShiftedFold.h, energy_grid.md section 2.5).
The clip makes the fold a smooth function of the shift: the integral of the linear interpolant up to a moving edge has
a continuous first derivative, and only its second derivative changes where an edge crosses a node.  A fold that cut
its range at whole nodes instead (the support's edge rounded to the enclosing node) would step at every crossing.

Model: tests/shift_sqrt_subpoints/models/const.azr (14N(p,g) through the 259.7 keV resonance, a 3.2e19 at/cm2 layer)
with the shift of its first data set varied (prior 0.5 keV, so a reach of 1.4 keV c.m.), in one session: chi^2 of
that data set at 21 shifts across 0.2 keV lab around its file value.  The back edge of the layer lies in the
resonance's core lattice (nodes 20 eV apart), so the sweep crosses about 10 nodes there.  Checked:
  (a) no re-plan inside the reach: every evaluation at the file value gives the same chi^2 to 1e-13;
  (c) a shift 5 keV away (beyond the reach: the mesh is planned again) and back gives the file value's chi^2 again;
  (b) smooth: the third differences of the 21 values stay below 5 % of the largest second difference.  Measured:
      second differences -0.79 .. -0.69, third differences 0.002 .. 0.009 (the trend of the curvature and its small
      changes where an edge crosses a node, every other step here).  A step at a node crossing has the size of the
      content of one 20 eV interval, hundreds in chi^2: a range cut at whole nodes (the kernel not clipped at the moved
      support) gives third differences of that size.

Run from anywhere:  python3 tests/pyazr/shift_fold_smooth_test.py
"""
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

MODELS = os.path.join(ROOT, "tests", "shift_sqrt_subpoints", "models")

work = tempfile.mkdtemp(prefix="pyazr_shift_smooth_")
try:
    shutil.copytree(os.path.join(MODELS, "data"), os.path.join(work, "data"))
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    p = Project.open(os.path.join(MODELS, "const.azr"))
    before = p["yield_1"]["shift"]
    p.set("yield_1", shift="offset 0.01043168251 MeV fit prior 0 0.0005")
    check("the shift of yield_1 is varied now", p["yield_1"]["shift"] != before, before)
    a = azure2(p.save(os.path.join(work, "smooth.azr")), cwd=work)
    x0 = np.array(a.params_rwa, float)
    free = [q for q in a.parameters if not q.fixed]
    k = [i for i, q in enumerate(free) if q.kind == "shift"][0]
    step = 1.0e-5  # MeV lab: 21 values over 0.2 keV
    shifts = x0[k] + step * np.arange(-10, 11)
    chi = []
    for v in shifts:
        x = x0.copy()
        x[k] = v
        chi.append(a.segment_chi2(x)[0])
    x = x0.copy()
    again = a.segment_chi2(x)[0]
    chi = np.array(chi)
    print("    chi2(shift): " + ", ".join(f"{c:.9g}" for c in chi))
    print("    d2: " + ", ".join(f"{c:.4g}" for c in np.diff(chi, 2)))
    print("    d3: " + ", ".join(f"{c:.3g}" for c in np.diff(chi, 3)))
    check("(a) the file value again, after the scan: the same chi2 (no re-plan inside the reach)",
          abs(again - chi[10]) <= 1e-13 * abs(chi[10]), (again, chi[10]))
    # (c) beyond the reach: the mesh is planned again for the moved support, and planned back on return
    far = x0.copy()
    far[k] = x0[k] + 0.005
    c_far = a.segment_chi2(far)[0]
    back = a.segment_chi2(x0.copy())[0]
    check("(c) 5 keV beyond the reach and back: the file value's chi2 again (the plan is a function of the shift)",
          abs(back - chi[10]) <= 1e-13 * abs(chi[10]) and np.isfinite(c_far), (back, chi[10], c_far))
    d2 = np.diff(chi, 2)
    d3 = np.diff(chi, 3)
    scale = float(np.max(np.abs(d2)))
    ratio = float(np.max(np.abs(d3))) / scale if scale > 0 else float("inf")
    print(f"    largest |second difference| {scale:.3e}, largest |third difference| {np.max(np.abs(d3)):.3e} "
          f"(ratio {ratio:.2e})")
    check("the scan moves chi2 (the shift is active)", abs(chi[-1] - chi[0]) > 1e-6 * abs(chi[10]), (chi[0], chi[-1]))
    check("(b) no step: third differences below 5 % of the second differences", ratio < 0.05, ratio)
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
