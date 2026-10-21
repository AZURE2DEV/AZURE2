#!/usr/bin/env python3
"""Electron screening through pyazr: the AzrModel helpers for the ``screening`` block of a <targetInt>
line (pure Python), the ``screening_<n>`` parameter record, its prior in penalties() / objective(),
and the U_e entries of chi2_and_grad and residual_jacobian against central finite differences.
Model: tests/electron_screening/grad.azr (12C(p,g), LUNA segments, U_e free with a 200 eV prior).

The engine part needs the compiled module and skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/screening_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "electron_screening")
DATA = os.path.join(ROOT, "tests", "13N", "data")
sys.path.insert(0, ROOT)

try:
    from pyazr.azrfile import AzrModel  # noqa: E402
except ImportError as e:  # as the other pyazr tests: no package dependencies (e.g. mpmath) -> nothing to test
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(0)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def close(a, b, tol):
    return abs(a - b) <= tol * max(abs(b), 1e-300)


tmp = tempfile.mkdtemp(prefix="pyazr_screening_")
try:
    # --- AzrModel: the block (no engine) ----------------------------------------------------------
    m = AzrModel.from_file(os.path.join(SRC, "scr.azr"))
    check("scr.azr: one screening line, effect 1 (13N's own lines are inactive)",
          m.screenings == {1: (1000.0, False, 0.0)}, m.screenings)
    m.set_screening(1, ue=900.0, vary=True, sigma=100.0)
    path = os.path.join(tmp, "set.azr")
    m.write(path)
    check("set_screening round trip", AzrModel.from_file(path).screenings == {1: (900.0, True, 100.0)},
          AzrModel.from_file(path).screenings)
    m.set_screening(1, sigma=0.0)
    check("set_screening keeps the values not given", m.screenings == {1: (900.0, True, 0.0)}, m.screenings)
    b = AzrModel.from_file(os.path.join(SRC, "bare.azr"))
    check("bare.azr: no screening", b.screenings == {}, b.screenings)
    b.add_target_effect("6", gaussian_sigma=0.001)
    b.add_target_effect("7", screening=(500.0, True, 50.0))
    check("add_target_effect(screening=...): numbered among the active lines",
          b.screenings == {2: (500.0, True, 50.0)}, b.screenings)
    b.set_screening(1, ue=700.0)
    check("set_screening adds a block to a line without one",
          b.screenings == {1: (700.0, False, 0.0), 2: (500.0, True, 50.0)}, b.screenings)
    try:
        b.set_screening(5, ue=1.0)
        check("set_screening refuses a missing line", False)
    except ValueError:
        check("set_screening refuses a missing line", True)

    # --- engine ----------------------------------------------------------------------------------
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import numpy as np
        from pyazr import azure2
        from pyazr import _azure2  # noqa: F401  the compiled engine
    except Exception as err:  # engine not built
        print(f"skip engine part: {type(err).__name__}: {err}")
        session = None
    else:
        work = os.path.join(tmp, "grad")
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        shutil.copy(os.path.join(SRC, "grad.azr"), work)
        os.symlink(DATA, os.path.join(work, "data"))
        session = azure2(os.path.join(work, "grad.azr"))

    if session is not None:
        scr = list(session.parameters.screenings)
        check("one screening parameter", len(scr) == 1, scr)
        p = scr[0]
        check("screening_1: free, U_e 1000 eV, prior 1000 +- 200, effect 1",
              p.name == "screening_1" and not p.fixed and close(p.value, 1000.0, 1e-12)
              and p.effect == 1 and close(p.prior_nominal, 1000.0, 1e-12) and close(p.prior_sigma, 200.0, 1e-12), p)

        x = np.asarray(session.params_rwa, float).copy()
        k = p.free_index
        x[k] = 1100.0
        pen = session.penalties(x)["screening"]
        check("penalties()['screening'] = ((1100 - 1000)/200)^2", close(pen, 0.25, 1e-12), pen)
        others = session.penalties(x)
        rest = (float(np.sum(others["norm"])) + float(np.sum(others["shift"])) +
                float(np.sum(others["shift_sqrt"])) + float(np.sum(others["shift_fn"])))
        obj = session.objective(x)
        data = float(np.sum(session.calculate_chi2_rwa(x)))
        check("objective = data chi2 + norm/shift terms + screening prior", close(obj, data + rest + 0.25, 1e-12),
              (obj, data, rest))
        x0 = np.asarray(session.params_rwa, float)
        check("U_e moves the data chi2", abs(data - float(np.sum(session.calculate_chi2_rwa(x0)))) > 1e-6)

        value, grad = session.chi2_and_grad(x)
        h = 1e-3 * max(1.0, abs(x[k]))
        e = np.zeros_like(x)
        e[k] = h
        fd = (float(np.sum(session.calculate_chi2_rwa(x + e))) - float(np.sum(session.calculate_chi2_rwa(x - e)))) / (2 * h)
        check(f"chi2_and_grad: d chi2/d U_e = finite difference ({grad[k]:.6g} vs {fd:.6g})",
              abs(grad[k] - fd) <= 1e-6 * max(1.0, abs(fd)), (grad[k], fd))

        r, J = session.residual_jacobian(x)
        r, J = np.asarray(r, float), np.asarray(J, float)
        rp = np.asarray(session.residual_jacobian(x + e)[0], float)
        rm = np.asarray(session.residual_jacobian(x - e)[0], float)
        col = (rp - rm) / (2 * h)
        worst = float(np.max(np.abs(J[:, k] - col) / np.maximum(1e-12, np.max(np.abs(col)))))
        check(f"residual_jacobian: U_e column = finite differences ({worst:.1e})", worst < 1e-6)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
print("PASS")
