#!/usr/bin/env python3
"""Free-gas Doppler broadening through pyazr, on the format-3 models of tests/doppler_fgm (make_models.py, on the
tree's n12C): a Doppler-only target changes sigma_t at the narrowed n + 12C resonance by tens of percent and a
transmission gives exp(-N sigma_D) (the API path assembles transformed segments through Segment::NeedsAssembly, which
must include the Doppler step).

The fitted temperature (``doppler free_gas T K fit ...`` on a target, parameter <target>.doppler), on grad_fit
(targets D0 (data sets 1, 2) and D1 (5, 6) "fit grid_max 900", D2 (3, 4) "fit prior 300 20"): the type-9 parameter
records, evaluate().terms["Doppler priors"] and the objective at moved temperatures, and the doppler entries of the
objective's gradient and Jacobian columns against central finite differences, at the start and at moved temperatures.  (The record
grammar -- the doppler field, its fit form and refusals -- is tests/doppler_fgm/check.sh's.)

Needs the compiled module and skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/doppler_test.py
"""
import importlib.util
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SUITE = os.path.join(ROOT, "tests", "doppler_fgm")
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({detail})"))
    if not ok:
        failures.append(name)


tmp = tempfile.mkdtemp(prefix="pyazr_doppler_")
try:
    spec = importlib.util.spec_from_file_location("doppler_models", os.path.join(SUITE, "make_models.py"))
    mm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mm)
    mm.main(os.path.join(ROOT, "tests", "reaction_rate", "models", "n12C.azr"), tmp)

    # --- engine ----------------------------------------------------------------------------------
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import numpy as np
        from pyazr import Session
        from pyazr import _azure2  # noqa: F401
    except Exception as err:
        print(f"skip engine part: {type(err).__name__}: {err}")
    else:
        def evaluate(name):
            work = os.path.join(tmp, "run_" + name)
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))
            shutil.copy(os.path.join(tmp, name + ".azr"), work)
            os.symlink(os.path.join(tmp, "data"), os.path.join(work, "data"))
            with Session(os.path.join(work, name + ".azr"), cwd=work) as s:
                return [np.asarray(v, float) for v in [d.model for d in s.evaluate(np.asarray(s.parameters.x0, float)).data]]
        def session(name):
            work = os.path.join(tmp, "run_" + name)
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))
            shutil.copy(os.path.join(tmp, name + ".azr"), work)
            os.symlink(os.path.join(tmp, "data"), os.path.join(work, "data"))
            return Session(os.path.join(work, name + ".azr"), cwd=work)

        bare, dop, dopT = evaluate("bare"), evaluate("dop"), evaluate("dopT")
        rel = float(np.max(np.abs(dop[0] / bare[0] - 1.0)))
        check("a Doppler-only line changes sigma_t by more than 5 % somewhere", rel > 0.05, rel)
        want = np.exp(-mm.N * dop[0])
        err = float(np.max(np.abs(dopT[0] - want) / want))
        check("transmission with doppler alone = exp(-N sigma_D) of the Doppler-only sigma_t (1e-9)", err < 1e-9, err)

        # --- fitted temperature ------------------------------------------------------------------
        with session("grad_fit") as s:
            dp = s.parameters.of_kind("doppler")
            by = {p.name: p for p in dp}
            ok = (sorted(by) == ["D0.doppler", "D1.doppler", "D2.doppler"] and all(not p.fixed for p in dp) and
                  all(p.prior_centre == 300.0 for p in dp) and
                  [by[n].prior_sigma for n in ("D0.doppler", "D1.doppler", "D2.doppler")] == [0.0, 0.0, 20.0] and
                  abs(by["D0.doppler"].layout_max - 900.0) < 1e-9 and abs(by["D1.doppler"].layout_max - 900.0) < 1e-9
                  and abs(by["D2.doppler"].layout_max - 363.0) < 1e-9)
            check("three type-9 records: D0, D1 (no prior, grids for 900 K), D2 (prior 20 K, grids for 363 K)", ok, dp)
            x = np.asarray(s.parameters.x0, float)
            pen0 = s.evaluate(x).terms["Doppler priors"]
            xm = x.copy()
            xm[by["D2.doppler"].free_index] = 320.0
            xm[by["D0.doppler"].free_index] = 450.0
            xm[by["D1.doppler"].free_index] = 450.0
            pen = s.evaluate(xm).terms["Doppler priors"]
            obj = s.objective()(xm) - s.evaluate(xm).chi2
            check("terms['Doppler priors'] = ((320 - 300)/20)^2 = 1 at D2 = 320 (D0, D1 have no prior), 0 at 300",
                  abs(pen - 1.0) < 1e-12 and pen0 == 0.0, (pen, pen0))
            check("objective() carries it", abs(obj - 1.0) < 1e-9, obj)
            c_lo = s.evaluate(x).chi2
            c_mv = s.evaluate(xm).chi2
            check("moving the temperatures moves the data chi2", abs(c_mv - c_lo) > 1e-6 * abs(c_lo), (c_lo, c_mv))
            # the objective's gradient: the doppler entries against central differences of the objective
            for name, xx in (("start", x), ("moved", xm)):
                c, g = s.objective().value_and_grad(xx)
                g = np.asarray(g, float)
                for p in dp:
                    h = 1e-3
                    xp, xn = xx.copy(), xx.copy()
                    xp[p.free_index] += h
                    xn[p.free_index] -= h
                    fd = (s.objective()(xp) - s.objective()(xn)) / (2 * h)
                    rel = abs(g[p.free_index] - fd) / max(abs(fd), 1e-300)
                    check(f"objective gradient d f/d {p.name} = central difference ({name}, rel 1e-6)", rel < 1e-6,
                          (g[p.free_index], fd, rel))
            # objective().residuals_and_jacobian: the doppler columns against central differences of the residuals
            r0, J = s.objective().residuals_and_jacobian(xm)
            J = np.asarray(J, float)
            for p in dp:
                h = 1e-3
                xp, xn = xm.copy(), xm.copy()
                xp[p.free_index] += h
                xn[p.free_index] -= h
                rp, _ = s.objective().residuals_and_jacobian(xp)
                rn, _ = s.objective().residuals_and_jacobian(xn)
                fd = (np.asarray(rp, float) - np.asarray(rn, float)) / (2 * h)
                err = float(np.max(np.abs(J[:, p.free_index] - fd)) / np.max(np.abs(fd)))
                check(f"objective Jacobian column {p.name} = central differences (1e-6 of the largest)", err < 1e-6, err)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: pyazr doppler")
