#!/usr/bin/env python3
"""Identity check (ARCHITECTURE.md section 10): the analytic gradient of chi^2 equals finite differences.

For the free parameters x of T1 and T3 (the engine's reduced-width vector in the Brune basis: level energies E_x in
MeV, reduced width amplitudes and ANC-derived amplitudes in MeV^1/2, beta amplitudes), the gradient of the data chi^2
from azure2.chi2_and_grad (reverse-mode adjoint of the engine) must equal the Richardson-extrapolated central
difference of azure2.calculate_chi2_rwa,

    D(h) = [chi^2(x + h e_i) - chi^2(x - h e_i)] / (2 h),   d chi^2 / d x_i ~ [4 D(h/2) - D(h)] / 3,
    h = 1e-6 max(|x_i|, 1e-2),

to 1e-5 of the largest gradient component.  The plain central difference is not enough: for the energy of the
3/2- level of T1 (sigma_E = 5e-5 MeV) its truncation error, h^2/6 times the third derivative of chi^2, is ~20
against a largest gradient component of 9e4 at the truth (2.4e-4); the extrapolation removes the h^2 term (observed
1.2e-8).  The check is made at two points: the truth, where the gradient comes from the noise of the data (largest
component 9e4 for T1, 6e3 for T3), and the start model, where it is 50 times larger.  A missing derivative
term (the Brune shift-function derivative, an external-capture amplitude, the identical-particle symmetrization, a
scattering-length constraint) shows up as a component off by order 1.

Needs the compiled engine module (pyazr); skips without it.

    python3 tests/pyazr/upgrade_gradient_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2
    sys.path.insert(0, os.path.join(ROOT, "tests", "upgrade", "lib"))
    import upgrade_project as up
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)   # ctest SKIP_RETURN_CODE (tests/pyazr/CMakeLists.txt)


def gradient_check(project, azr_name):
    work = tempfile.mkdtemp(prefix=f"upgrade_grad_{project}.")
    try:
        up.project_copy(project, work)
        with azure2(os.path.join(work, azr_name + ".azr")) as m:
            x = np.asarray(m.params_rwa, float)
            chi2, grad = m.chi2_and_grad(x)
            chi2_ref = float(np.sum(m.calculate_chi2_rwa(x)))
            fd = np.empty_like(x)

            def central(i, h):
                xp, xm = x.copy(), x.copy()
                xp[i] += h
                xm[i] -= h
                return (np.sum(m.calculate_chi2_rwa(xp)) - np.sum(m.calculate_chi2_rwa(xm))) / (2 * h)

            for i in range(x.size):
                h = 1e-6 * max(abs(x[i]), 1e-2)
                fd[i] = (4.0 * central(i, h / 2) - central(i, h)) / 3.0
    finally:
        shutil.rmtree(work, ignore_errors=True)
    label = f"{project} at {'the truth' if azr_name == project else 'the start model'}"
    check(f"{label}: chi2 from chi2_and_grad = calculate_chi2_rwa", abs(chi2 / chi2_ref - 1) < 1e-10,
          f"{chi2:.6f} vs {chi2_ref:.6f}")
    scale = np.max(np.abs(fd))
    worst = np.max(np.abs(grad - fd)) / scale
    i = int(np.argmax(np.abs(grad - fd)))
    check(f"{label}: analytic gradient = central differences ({x.size} parameters)", worst < 1e-5,
          f"largest |difference| / max|grad| = {worst:.1e} (parameter {i}); max|grad| = {scale:.3g}")


for project in ("T1", "T3"):
    gradient_check(project, project)
    gradient_check(project, project + "_start")

print("FAIL" if failures else "ok")
sys.exit(1 if failures else 0)
