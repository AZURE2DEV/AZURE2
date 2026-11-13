#!/usr/bin/env python3
"""The fitted sample thickness of transmission / self-shielded-yield data sets through pyazr, on the format-3 models
of tests/transmission_yield (make_models.py, n12c_models, on the tree's n12C).

On grad_fit (eight transmission / yield data sets on sample targets whose areal_density is fitted: 1, 2 without a
beam, 3, 4 under a Gaussian spread with a 0.01 atoms/b prior, 5, 6 one target under doppler alone, 7, 8 under doppler
+ the Gaussian; every width free): the type-10 parameter records (one per target: data sets 5 and 6 share one),
penalties()["thickness"] and objective() at moved thicknesses, and the thickness entries of chi2_and_grad and the
residual_jacobian columns against central differences.  (The record grammar -- areal_density on a sample target, its
fit form and refusals -- is tests/transmission_yield/check.sh's.)

Needs the compiled module and skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/thickness_test.py
"""
import importlib.util
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SUITE = os.path.join(ROOT, "tests", "transmission_yield")
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({detail})"))
    if not ok:
        failures.append(name)


tmp = tempfile.mkdtemp(prefix="pyazr_thickness_")
try:
    spec = importlib.util.spec_from_file_location("ty_models", os.path.join(SUITE, "make_models.py"))
    mm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mm)
    mm.n12c_models(os.path.join(ROOT, "tests", "reaction_rate", "models", "n12C.azr"), tmp)

    # --- engine ----------------------------------------------------------------------------------
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import numpy as np
        from pyazr import Session
        from pyazr import _azure2  # noqa: F401
    except Exception as err:
        print(f"skip engine part: {type(err).__name__}: {err}")
    else:
        work = os.path.join(tmp, "run_grad_fit")
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        shutil.copy(os.path.join(tmp, "grad_fit.azr"), work)
        os.symlink(os.path.join(tmp, "data"), os.path.join(work, "data"))
        with Session(os.path.join(work, "grad_fit.azr"), cwd=work) as s:
            tp = s.parameters.thicknesses
            print("    " + "; ".join(f"{p.name} (data set {p.segment_key}, prior {p.prior_sigma}, fixed {p.fixed})"
                                      for p in tp))
            free = [p for p in tp if not p.fixed]
            keys = {p.segment_key: p for p in tp}
            ok = (all(p.prior_nominal == 0.15 for p in tp) and
                  [keys[k].prior_sigma for k in (1, 2, 3, 4)] == [0, 0, 0.01, 0.01] and len(free) == 7)
            check("type-10 records: priors on 3, 4; 5 and 6 one parameter (7 free in all)", ok, list(tp))
            x = np.asarray(s.params_rwa, float)
            pen0 = s.penalties(x)["thickness"]
            xm = x.copy()
            for p in free:
                xm[p.free_index] = {3: 0.16, 4: 0.17}.get(p.segment_key, 0.15 + 0.004 * p.segment_key)
            pen = s.penalties(xm)["thickness"]
            obj = s.objective(xm) - float(np.sum(s.calculate_chi2_rwa(xm)))
            check("penalties()['thickness'] = 1 + 4 at data set 3 = 0.16, data set 4 = 0.17; 0 at the start",
                  abs(pen - 5.0) < 1e-9 and pen0 == 0.0, (pen, pen0))
            check("objective() carries it", abs(obj - 5.0) < 1e-6, obj)
            c0 = float(np.sum(s.calculate_chi2_rwa(x)))
            c1 = float(np.sum(s.calculate_chi2_rwa(xm)))
            check("moving the thicknesses moves the data chi2", abs(c1 - c0) > 1e-6 * abs(c0), (c0, c1))
            for name, xx in (("start", x), ("moved", xm)):
                c, g = s.chi2_and_grad(xx)
                g = np.asarray(g, float)
                worst = 0.0
                for p in free:
                    h = 1e-6
                    xp, xn = xx.copy(), xx.copy()
                    xp[p.free_index] += h
                    xn[p.free_index] -= h
                    fd = (float(np.sum(s.calculate_chi2_rwa(xp))) - float(np.sum(s.calculate_chi2_rwa(xn)))) / (2 * h)
                    worst = max(worst, abs(g[p.free_index] - fd) / max(abs(fd), 1e-300))
                check(f"chi2_and_grad d chi2/d thickness = central differences ({name}, 7 of them, rel 1e-6)",
                      worst < 1e-6, worst)
            r0, J = s.residual_jacobian(xm)
            J = np.asarray(J, float)
            worst = 0.0
            for p in free:
                h = 1e-6
                xp, xn = xm.copy(), xm.copy()
                xp[p.free_index] += h
                xn[p.free_index] -= h
                rp, _ = s.residual_jacobian(xp)
                rn, _ = s.residual_jacobian(xn)
                fd = (np.asarray(rp, float) - np.asarray(rn, float)) / (2 * h)
                worst = max(worst, float(np.max(np.abs(J[:, p.free_index] - fd)) / np.max(np.abs(fd))))
            check("residual_jacobian columns of the thicknesses = central differences (1e-6 of the largest)",
                  worst < 1e-6, worst)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: pyazr sample thickness")
