#!/usr/bin/env python3
"""Two levels of one J^pi at one energy: the limits of the shift-function divided differences.

The Brune level matrix holds (S_a (E - E_b) - S_b (E - E_a))/(E_a - E_b) and Park's overlap (S_a - S_b)/(E_a - E_b)
for every pair of levels of a J-group.  At E_a = E_b these are 0/0: before the fix the chi-squared of a pair exactly
on top of each other came back finite but wrong (here 74,726 instead of 47,983; the matrix is not a number and the
J-group's contribution is lost) and the gradient NaN.  Below kLevelCoincidenceGap (1e-9 MeV) the engine now takes the limits, S + S' (E - E_a) and S', and in
the analytic derivatives their limits in S''.

Model: tests/identical_pp_res (p+p, a 2- group with two channels) with a second 2- level, free energy, amplitudes not
parallel to the first level's.  The energy of that level is put ON the file's 1.2 MeV level through the parameter
vector, in the Brune and Park formalisms:

  * the chi-squared is finite;
  * it is continuous across the switch: gap 0, 5e-10 (limit) and 2e-9 MeV (divided difference) agree;
  * it joins the ordinary divided differences at 1e-7 MeV;
  * `chi2_and_grad` and `residual_jacobian` at gap 0 agree with central differences.

A file cannot hold the pair at the identical energy (the reader takes them as one level, tests/level_ids), so the
standard parametrization, whose transformation is done from the file, is checked at a gap of 5e-10 MeV against 1e-7.

Run from anywhere:  python3 tests/pyazr/coincident_levels_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "4")
    import numpy as np
    from pyazr import azure2, AzrModel
except Exception as err:
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)

NAME = "identical_pp_res"
E1, E2 = 1.2, 1.21


def stage(tmp, tag, energy=E2):
    work = os.path.join(tmp, tag)
    shutil.copytree(os.path.join(ROOT, "tests", NAME), work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
        os.makedirs(os.path.join(work, junk))
    azr = os.path.join(work, NAME + ".azr")
    m = AzrModel.from_file(azr)
    lv = [l for l in m.levels if l.active and l.jpi == "2-" and abs(float(l.energy) - E1) < 1e-9][0]
    ch = [dict(pair=c.pair, L=c.L, S=c.S, gamma=g, fixed=False)
          for c, g in zip([c for c in lv.channels if c.active], (30000.0, -40000.0))]
    assert len(ch) == 2, ch
    m.add_level(lv.J, lv.parity, energy, ch, level_fixed=False)
    m.set_output_dir(os.path.join(work, "output"))
    return m.write(azr)


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(a)), np.max(np.abs(b)), 1e-300))


def central(f, x, h):
    r0 = np.atleast_1d(f(x))
    J = np.zeros((r0.size, x.size))
    for i in range(x.size):
        z = x.copy(); z[i] += h[i]; fp = np.atleast_1d(f(z))
        z[i] += h[i]; fpp = np.atleast_1d(f(z))
        z[i] = x[i] - h[i]; fm = np.atleast_1d(f(z))
        z[i] -= h[i]; fmm = np.atleast_1d(f(z))
        J[:, i] = (8.0 * (fp - fm) - (fpp - fmm)) / (12.0 * h[i])
    return J


with tempfile.TemporaryDirectory() as tmp:
    for label, kw in (("Brune", dict(use_brune=True)), ("Park", dict(use_park=True))):
        print(label, flush=True)
        with azure2(stage(tmp, label), **kw) as m:
            x0 = np.asarray(m.params_rwa, float)
            free = [p.free_index for p in m.parameters.energies if p.free_index is not None and abs(x0[p.free_index] - E2) < 1e-9]
            assert len(free) == 1, free
            j = free[0]

            def at(gap):
                z = x0.copy(); z[j] = E1 + gap
                return z

            def chi2(z):
                return float(np.sum(m.calculate_chi2_rwa(z)))

            c0 = chi2(at(0.0))
            check(f"{label}: chi2 at gap 0 is finite ({c0:.8g})", bool(np.isfinite(c0)))
            if not np.isfinite(c0):
                continue
            c1, c2, c3, c4 = chi2(at(5e-10)), chi2(at(2e-9)), chi2(at(1e-7)), chi2(at(-5e-10))
            check(f"{label}: 5e-10 MeV (limit) = gap 0", abs(c1 - c0) < 1e-6 * c0, f"{c1} {c0}")
            check(f"{label}: continuous across the switch (5e-10 / 2e-9 MeV)", abs(c2 - c1) < 1e-6 * c0, f"{c2} {c1}")
            check(f"{label}: joins the divided differences at 1e-7 MeV", abs(c3 - c0) < 1e-5 * c0, f"{c3} {c0}")
            check(f"{label}: symmetric in the sign of the gap", abs(c4 - c1) < 1e-5 * c0, f"{c4} {c1}")

            z = at(0.0)
            h = 1e-4 * (np.abs(z) + 1.0)
            ref_val, ref_g = m.chi2_and_grad(at(0.01))
            ref_fd = central(chi2, at(0.01), h).ravel()
            ref = rel(ref_g, ref_fd)
            val, g = m.chi2_and_grad(z)
            gfd = central(chi2, z, h).ravel()
            d = rel(g, gfd)
            check(f"{label}: gradient at gap 0 = central differences (rel {d:.1e}; {ref:.1e} at 10 keV)", d < max(1e-4, 3.0 * ref),
                  "\n analytic " + np.array2string(g, precision=4) + "\n numeric  " + np.array2string(gfd, precision=4))
            r, Jac = m.residual_jacobian(z)
            Jfd = central(lambda y: np.asarray(m.residual_jacobian(y)[0], float), z, 0.1 * h)
            d = rel(Jac, Jfd)
            check(f"{label}: Jacobian at gap 0 = central differences (rel {d:.1e})", d < max(2e-4, 3.0 * ref))

    print("Standard parametrization", flush=True)
    out = {}
    for tag, gap in (("near", 5e-10), ("far", 1e-7)):
        with azure2(stage(tmp, "std_" + tag, E1 + gap), use_brune=False) as m:
            out[tag] = float(np.sum(m.calculate_chi2_rwa(np.asarray(m.params_rwa, float))))
    check(f"standard: finite at 5e-10 MeV ({out['near']:.8g})", bool(np.isfinite(out["near"])))
    check("standard: limit (5e-10 MeV) joins the divided differences (1e-7 MeV)",
          abs(out["near"] - out["far"]) < 1e-4 * abs(out["far"]), str(out))

print("ALL PASSED" if not failures else f"FAILED: {failures}")
sys.exit(1 if failures else 0)
