#!/usr/bin/env python3
"""Free-gas Doppler broadening through pyazr: AzrModel.add_target_effect(doppler=T) writes the ``doppler`` keyword
and refuses T <= 0 (pure Python); through the engine, a Doppler-only line changes sigma_t at the narrowed n + 12C
resonance of tests/doppler_fgm by tens of percent and a transmission line gives exp(-N sigma_D) (the API path
assembles transformed segments through ESegment::NeedsAssembly, which must include the Doppler step).

The fitted temperature (``doppler T vary sigma [Tmax]``, parameter doppler_<n>): AzrModel.dopplers / set_doppler /
add_target_effect(doppler=(T, vary, sigma, tmax)) round trips and refusals (pure Python); through the engine, on
tests/doppler_fgm's grad_fit model (doppler_1: "300 1 0 900" on a Doppler-only + transmission/yield line, doppler_2:
"300 1 20" on a Doppler + Gaussian transmission/yield line), the type-9 parameter records, penalties()["doppler"]
and objective() at a moved T, and the doppler_n entries of chi2_and_grad and the residual_jacobian columns against
central finite differences, at the start and at moved temperatures.

The engine part needs the compiled module and skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/doppler_test.py
"""
import importlib.util
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SUITE = os.path.join(ROOT, "tests", "doppler_fgm")
sys.path.insert(0, ROOT)

try:
    from pyazr.azrfile import AzrModel  # noqa: E402
except ImportError as e:
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(0)

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

    # --- AzrModel (no engine) ---------------------------------------------------------------------
    b = AzrModel.from_file(os.path.join(tmp, "bare.azr"))
    b.add_target_effect("1,2", n_points=400, doppler=293.6)
    path = os.path.join(tmp, "added.azr")
    b.write(path)
    lines = open(path).read().split("<targetInt>")[1].split("</targetInt>")[0].strip().splitlines()
    check("add_target_effect(doppler=293.6) writes one line ending in 'doppler 293.6'",
          len(lines) == 1 and lines[0].split()[-2:] == ["doppler", "293.6"], lines)
    for bad in (0, -5):
        try:
            AzrModel.from_file(os.path.join(tmp, "bare.azr")).add_target_effect("1", doppler=bad)
            check(f"add_target_effect(doppler={bad}) refused", False)
        except ValueError:
            check(f"add_target_effect(doppler={bad}) refused", True)

    def target_lines(model):
        path = os.path.join(tmp, "rt.azr")
        model.write(path)
        return open(path).read().split("<targetInt>")[1].split("</targetInt>")[0].strip().splitlines()

    f = AzrModel.from_file(os.path.join(tmp, "bare.azr")).add_target_effect("1,2", n_points=400,
                                                                            doppler=(300, True, 20, 900))
    check("add_target_effect(doppler=(300, 1, 20, 900)) writes the fit form",
          target_lines(f)[0].split()[-5:] == ["doppler", "300", "1", "20", "900"], target_lines(f))
    d = f.dopplers
    check("dopplers reads it back (fit, vary, sigma 20, tmax 900, max 900)",
          d == {1: dict(temperature=300.0, fit=True, vary=True, sigma=20.0, tmax=900.0, max=900.0)}, d)
    f.set_doppler(1, tmax=0)
    d = f.dopplers[1]
    check("set_doppler(tmax=0) removes Tmax: the layout is max(T + 3 sigma, 1.21 T) = 363 K",
          target_lines(f)[0].split()[-4:] == ["doppler", "300", "1", "20"] and d["max"] == 363.0, (target_lines(f), d))
    f.set_doppler(1, T=310, vary=False, sigma=0)
    d = f.dopplers[1]
    check("set_doppler(T=310, vary=False, sigma=0): 'doppler 310 0 0', layout 2 T",
          target_lines(f)[0].split()[-4:] == ["doppler", "310", "0", "0"] and d["max"] == 620.0 and not d["vary"],
          (target_lines(f), d))
    g = AzrModel.from_file(os.path.join(tmp, "dop.azr"))
    g.set_doppler(1, T=350)
    check("set_doppler(T=350) on a plain line keeps it plain",
          target_lines(g)[0].split()[-2:] == ["doppler", "350"] and not g.dopplers[1]["fit"], target_lines(g))
    g.set_doppler(1, vary=True, sigma=25, tmax=800)
    written = target_lines(g)
    check("set_doppler(vary, sigma, tmax) turns it into the fit form, round trip through the file",
          AzrModel.from_file(os.path.join(tmp, "rt.azr")).dopplers == g.dopplers and
          written[0].split()[-5:] == ["doppler", "350", "1", "25", "800"], written)
    for what, kw in (("a negative sigma", dict(sigma=-1)), ("Tmax not above T", dict(tmax=300)),
                     ("T <= 0", dict(T=0))):
        try:
            AzrModel.from_file(os.path.join(tmp, "dop.azr")).set_doppler(1, vary=True, **kw)
            check(f"set_doppler refuses {what}", False)
        except ValueError:
            check(f"set_doppler refuses {what}", True)
    try:
        AzrModel.from_file(os.path.join(tmp, "bare.azr")).add_target_effect("1", n_points=10).set_doppler(1, vary=True)
        check("set_doppler on a line without a doppler block needs T", False)
    except ValueError:
        check("set_doppler on a line without a doppler block needs T", True)

    # --- engine ----------------------------------------------------------------------------------
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import numpy as np
        from pyazr import azure2
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
            with azure2(os.path.join(work, name + ".azr"), cwd=work) as s:
                return [np.asarray(v, float) for v in s.calculate_rwa(np.asarray(s.params_rwa, float))]
        def session(name):
            work = os.path.join(tmp, "run_" + name)
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))
            shutil.copy(os.path.join(tmp, name + ".azr"), work)
            os.symlink(os.path.join(tmp, "data"), os.path.join(work, "data"))
            return azure2(os.path.join(work, name + ".azr"), cwd=work)

        bare, dop, dopT = evaluate("bare"), evaluate("dop"), evaluate("dopT")
        rel = float(np.max(np.abs(dop[0] / bare[0] - 1.0)))
        check("a Doppler-only line changes sigma_t by more than 5 % somewhere", rel > 0.05, rel)
        want = np.exp(-mm.N * dop[0])
        err = float(np.max(np.abs(dopT[0] - want) / want))
        check("transmission with doppler alone = exp(-N sigma_D) of the Doppler-only sigma_t (1e-9)", err < 1e-9, err)

        # --- fitted temperature ------------------------------------------------------------------
        with session("grad_fit") as s:
            dp = s.parameters.dopplers
            ok = (len(dp) == 2 and [p.effect for p in dp] == [1, 2] and all(not p.fixed for p in dp) and
                  [p.prior_nominal for p in dp] == [300.0, 300.0] and [p.prior_sigma for p in dp] == [0.0, 20.0] and
                  abs(dp[0].doppler_max - 900.0) < 1e-9 and abs(dp[1].doppler_max - 363.0) < 1e-9 and
                  [p.name for p in dp] == ["doppler_1", "doppler_2"])
            check("two type-9 records: doppler_1 (no prior, layout 900 K), doppler_2 (prior 20 K, layout 363 K)", ok, dp)
            x = np.asarray(s.params_rwa, float)
            i1, i2 = dp[0].free_index, dp[1].free_index
            pen0 = s.penalties(x)["doppler"]
            xm = x.copy()
            xm[i2] = 320.0
            xm[i1] = 450.0
            pen = s.penalties(xm)["doppler"]
            obj = s.objective(xm) - float(np.sum(s.calculate_chi2_rwa(xm)))
            check("penalties()['doppler'] = ((320 - 300)/20)^2 = 1 at doppler_2 = 320 (doppler_1 has no prior), 0 at 300",
                  abs(pen - 1.0) < 1e-12 and pen0 == 0.0, (pen, pen0))
            check("objective() carries it", abs(obj - 1.0) < 1e-9, obj)
            c_lo = float(np.sum(s.calculate_chi2_rwa(x)))
            c_mv = float(np.sum(s.calculate_chi2_rwa(xm)))
            check("moving the temperatures moves the data chi2", abs(c_mv - c_lo) > 1e-6 * abs(c_lo), (c_lo, c_mv))
            # chi2_and_grad: the doppler entries against central differences of the data chi2
            for name, xx in (("start", x), ("moved", xm)):
                c, g = s.chi2_and_grad(xx)
                g = np.asarray(g, float)
                for p in dp:
                    h = 1e-3
                    xp, xn = xx.copy(), xx.copy()
                    xp[p.free_index] += h
                    xn[p.free_index] -= h
                    fd = (float(np.sum(s.calculate_chi2_rwa(xp))) - float(np.sum(s.calculate_chi2_rwa(xn)))) / (2 * h)
                    rel = abs(g[p.free_index] - fd) / max(abs(fd), 1e-300)
                    check(f"chi2_and_grad d chi2/d {p.name} = central difference ({name}, rel 1e-6)", rel < 1e-6,
                          (g[p.free_index], fd, rel))
            # residual_jacobian: the doppler columns against central differences of the residuals
            r0, J = s.residual_jacobian(xm)
            J = np.asarray(J, float)
            for p in dp:
                h = 1e-3
                xp, xn = xm.copy(), xm.copy()
                xp[p.free_index] += h
                xn[p.free_index] -= h
                rp, _ = s.residual_jacobian(xp)
                rn, _ = s.residual_jacobian(xn)
                fd = (np.asarray(rp, float) - np.asarray(rn, float)) / (2 * h)
                err = float(np.max(np.abs(J[:, p.free_index] - fd)) / np.max(np.abs(fd)))
                check(f"residual_jacobian column {p.name} = central differences (1e-6 of the largest)", err < 1e-6, err)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s)")
    sys.exit(1)
print("PASS: pyazr doppler")
