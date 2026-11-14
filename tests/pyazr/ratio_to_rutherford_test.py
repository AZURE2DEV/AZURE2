#!/usr/bin/env python3
"""Ratio to Rutherford through pyazr (data_kinds.md section 4.3, item 4 of section 5.2), against the independent
numbers of tests/ratio_to_rutherford/truth.py (mpmath hard-sphere scattering in the Coulomb field):

1. without shift or fold, the chi^2 of ratio data (models/aC_data.azr) equals, to 1e-12, the chi^2 of the same data
   multiplied by the Rutherford cross section of truth.py and given as dsigma_domega; both equal truth.py's chi^2;
2. with a fixed 5 keV lab energy shift the two forms differ, each as truth.py computes it (the ratio form
   re-evaluates sigma_R at the shifted energy, the multiplied data keep sigma_R(E));
3. with a Gaussian beam spread (0.2 MeV lab) the model is the ratio of the two folds <sigma>/<sigma_R>, not the fold
   of the ratio <sigma/sigma_R>, and the numerator data set (dsigma_domega with the same beam) over the ratio gives
   truth.py's <sigma_R>;
4. the binding's rutherford_cross_section and mott_cross_section against truth.py's |f_C|^2 and
   |f_C(theta) + f_C(pi - theta)|^2;
5. the analytic gradient of the chi^2 (levels, norm; with and without the beam) against central differences
   (models/aC_story.azr, the worked example of project_format3.md section 9.2);
6. the worked example's lab-angle data and the same data at c.m. angles converted by an independent formula give the
   same chi^2 (the ratio is the same in both frames: the value is not converted at load).

Needs the compiled module and mpmath; skips cleanly without them.

Run from anywhere:  python3 tests/pyazr/ratio_to_rutherford_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
CASE = os.path.join(ROOT, "tests", "ratio_to_rutherford")
sys.path.insert(0, ROOT)
sys.path.insert(0, CASE)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


try:
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    import mpmath  # noqa: F401  (truth.py)
    import truth
    from pyazr import Session, _azure2
except Exception as err:
    print(f"skip: the engine or mpmath is not importable ({type(err).__name__}: {err})")
    sys.exit(77)

DATA = os.path.join(CASE, "models", "data", "aC_rth_cm.dat")


def project(work, name, text):
    d = os.path.join(work, name)
    os.makedirs(os.path.join(d, "output"))
    os.makedirs(os.path.join(d, "checks"))
    shutil.copytree(os.path.join(CASE, "models", "data"), os.path.join(d, "data"))
    with open(os.path.join(d, "m.azr"), "w") as f:
        f.write(text)
    return d


def chi2_of(d):
    with Session(os.path.join(d, "m.azr"), cwd=d) as m:
        x = np.asarray(m.parameters.x0, float)
        return m.evaluate(x).chi2


def values_of(d):
    with Session(os.path.join(d, "m.azr"), cwd=d) as m:
        return np.concatenate([np.asarray(v, float) for v in [d.model for d in m.evaluate(np.asarray(m.parameters.x0, float)).data]])


base = open(os.path.join(CASE, "models", "aC_data.azr")).read()
SHIFT = "    norm           1 fixed\n    shift          offset 0.005 MeV fixed\n"
BEAM = "    norm           1 fixed\n    beam           B\n    quadrature     points 200\n"


def as_dsdo(text):
    return (text.replace("quantity       ratio_to_rutherford", "quantity       dsigma_domega")
                .replace("unit           1\n", "unit           b/sr\n")
                .replace("data/aC_rth_cm.dat", "data/aC_dsdo_cm.dat"))


with tempfile.TemporaryDirectory() as work:
    print("1. ratio data against the same data multiplied by sigma_R (no shift, no fold)")
    r0 = project(work, "r0", base)
    truth.multiply(DATA, os.path.join(r0, "data", "aC_dsdo_cm.dat"))
    m0 = project(work, "m0", as_dsdo(base))
    shutil.copy(os.path.join(r0, "data", "aC_dsdo_cm.dat"), os.path.join(m0, "data"))
    cr, cm = chi2_of(r0), chi2_of(m0)
    tr0 = truth.chi2_values(DATA, "0")
    check("chi2 ratio = chi2 multiplied to 1e-12", rel(cr, cm) <= 1e-12, f"{cr!r} {cm!r} ({rel(cr, cm):.2e})")
    check("chi2 ratio = truth.py to 1e-9", rel(cr, float(tr0[0])) <= 1e-9, f"{cr!r} {float(tr0[0])!r}")

    print("2. a fixed 5 keV lab energy shift")
    r5 = project(work, "r5", base.replace("    norm           1 fixed\n", SHIFT))
    m5 = project(work, "m5", as_dsdo(base).replace("    norm           1 fixed\n", SHIFT))
    shutil.copy(os.path.join(r0, "data", "aC_dsdo_cm.dat"), os.path.join(m5, "data"))
    cr5, cm5 = chi2_of(r5), chi2_of(m5)
    tr5, tm5 = (float(v) for v in truth.chi2_values(DATA, "0.005"))
    check("ratio form = truth.py (sigma_R at the shifted energy) to 1e-9", rel(cr5, tr5) <= 1e-9, f"{cr5!r} {tr5!r}")
    check("multiplied form = truth.py (sigma_R frozen at E) to 1e-9", rel(cm5, tm5) <= 1e-9, f"{cm5!r} {tm5!r}")
    check("the two forms differ by truth.py's amount (%.6g)" % (tm5 - tr5),
          abs((cm5 - cr5) - (tm5 - tr5)) <= 1e-8 * max(cm5, cr5) and abs(tm5 - tr5) > 1e-3 * tr5,
          f"engine {cm5 - cr5!r}")

    print("3. a Gaussian beam spread: the ratio of the folds")
    beam = base.replace("    norm           1 fixed\n", BEAM).replace("<targets>\n", "<targets>\nbeam B\n    spread         0.2 MeV lab\n")
    rb = project(work, "rb", beam)
    nb = project(work, "nb", as_dsdo(beam))
    shutil.copy(os.path.join(r0, "data", "aC_dsdo_cm.dat"), os.path.join(nb, "data"))
    rows = [3, 6, 24, 27, 55]  # where the two differ by 1.8e-4 .. 7e-3 (row 0, 20 deg at 2 MeV: 3e-5)
    ratio, numerator = values_of(rb), values_of(nb)
    for k, (e_lab, th, rof, fr, ref) in zip(rows, truth.fold(DATA, "0.2", rows)):
        rof, fr, ref = float(rof), float(fr), float(ref)
        check(f"E_lab {float(e_lab):g} MeV, {float(th):g} deg: model = <sigma>/<sigma_R> to 5e-6 (fold of the ratio: "
              f"{(fr - rof) / rof:+.2e})", rel(ratio[k], rof) <= 5e-6 and abs(fr - rof) > 20 * 5e-6 * rof,
              f"{ratio[k]!r} {rof!r} {fr!r}")
        # each fold alone carries the kernel's quadrature error (sigma linear between sub-points: 1.4e-5 at 2 MeV),
        # which the ratio of two folds on the same sub-points cancels to the 5e-6 above
        check("  the numerator's fold over the model = truth.py's <sigma_R> to 5e-5",
              rel(numerator[k] / ratio[k], ref) <= 5e-5, f"{numerator[k] / ratio[k]!r} {ref!r}")

    print("4. the binding's Coulomb cross sections")
    e, th = np.array([1.5, 2.25, 3.75]), np.array([20.0, 90.0, 160.0])
    import pyazr.physics
    ru = pyazr.physics.rutherford_cross_section(2, 6, 4.00260325, 12.0, e, th)
    mo = pyazr.physics.mott_cross_section(2, 4.00260325, 0, e, th)
    check("rutherford_cross_section = |f_C|^2 (alpha + 12C) to 1e-13",
          all(rel(ru["value"][i], float(truth.coulomb("aC", e[i], th[i]))) <= 1e-13 for i in range(3)), str(ru["value"]))
    check("mott_cross_section = |f_C(theta) + f_C(pi - theta)|^2 (alpha + alpha) to 1e-13",
          all(rel(mo["value"][i], float(truth.coulomb("aa", e[i], th[i]))) <= 1e-13 for i in range(3)), str(mo["value"]))
    check("d/dE of Rutherford = -2 sigma_R / E", all(rel(ru["d_dE"][i], -2 * ru["value"][i] / e[i]) <= 1e-13 for i in range(3)))

    print("5. the gradient against central differences (worked example, section 9.2)")
    story = open(os.path.join(CASE, "models", "aC_story.azr")).read()
    with_beam = story.replace("    norm           1 fit prior 1 3%\n",
                              "    norm           1 fit prior 1 3%\n    beam           B\n    quadrature     points 60\n").replace(
        "<targets>\n", "<targets>\nbeam B\n    spread         0.05 MeV lab\n")
    for name, text in (("plain", story), ("beam", with_beam)):
        d = project(work, "g_" + name, text)
        with Session(os.path.join(d, "m.azr"), cwd=d) as m:
            x = np.asarray(m.parameters.x0, float)
            c, g = m.objective().value_and_grad(x)
            worst = 0.0
            for i in range(len(x)):
                h = 1e-6 * max(abs(x[i]), 1e-3)
                xp, xm = x.copy(), x.copy()
                xp[i] += h
                xm[i] -= h
                fd = (m.objective()(xp) - m.objective()(xm)) / (2 * h)
                worst = max(worst, abs(g[i] - fd) / max(np.max(np.abs(g)), 1e-300))
            check(f"{name}: {len(x)} free parameters, analytic = central differences to 1e-6 of the largest",
                  len(x) >= 7 and worst <= 1e-6, f"worst {worst:.2e}, chi2 {c!r}")

    print("6. the same ratio data at lab and at c.m. angles (the value is not converted)")
    # elastic kinematics, independent of the engine: theta_cm = theta_lab + asin((m1/m2) sin theta_lab)
    import math
    rows = [l.split() for l in open(os.path.join(CASE, "models", "data", "aC_rth.dat")) if l.split() and not l.startswith("#")]
    lab = project(work, "lab", story)
    cm = project(work, "cm", story.replace("theta_lab[deg]", "theta_cm[deg]").replace("data/aC_rth.dat", "data/aC_rth_cm_angles.dat"))
    with open(os.path.join(cm, "data", "aC_rth_cm_angles.dat"), "w") as f:
        for e_lab, th, r, dr in rows:
            t = math.radians(float(th))
            f.write("%s %.17g %s %s\n" % (e_lab, math.degrees(t + math.asin(4.00260325 / 12.0 * math.sin(t))), r, dr))
    cl, cc = chi2_of(lab), chi2_of(cm)
    check("chi2 with lab angles = chi2 with the c.m. angles of an independent conversion, to 1e-9", rel(cl, cc) <= 1e-9,
          f"{cl!r} {cc!r} ({rel(cl, cc):.2e})")

if failures:
    print(f"\n{len(failures)} FAILED")
    sys.exit(1)
print("\nall passed")
