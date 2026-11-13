#!/usr/bin/env python3
"""Electron screening through pyazr: the ``screening`` field of a format-3 target record (Project: read, set,
added, refused), the ``<target>.screening`` parameter record, its prior in evaluate().terms / objective(),
and the U_e entries of objective().value_and_grad and objective().residuals_and_jacobian against central finite differences.
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
    from pyazr import Project, Nuisance, FormatError  # noqa: E402
except ImportError as e:  # as the other pyazr tests: no package dependencies (e.g. mpmath) -> nothing to test
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(77)

failures = []
engine_skipped = False


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def close(a, b, tol):
    return abs(a - b) <= tol * max(abs(b), 1e-300)


tmp = tempfile.mkdtemp(prefix="pyazr_screening_")
try:
    # --- the target records (no engine session) ---------------------------------------------------
    def screenings(p):
        return {n: r["screening"] for n, r in p.records().items() if r.record_kind == "target" and "screening" in r}

    m = Project.open(os.path.join(SRC, "scr.azr"))
    check("scr.azr: one screening target, named by the two LUNA data sets",
          screenings(m) == {"tgt_targetInt3": "1000 eV fixed"} and
          [n for n, r in m.records().items() if r.get("target") == "tgt_targetInt3"]
          == ["skowronski_luna_hpge", "skowronski_luna_bgo"], screenings(m))
    m.set("tgt_targetInt3", screening=Nuisance(900, "fit", unit="eV", prior=(900, 100)))
    path = os.path.join(tmp, "set.azr")
    m.save(path)
    check("set screening round trip", screenings(Project.open(path)) == {"tgt_targetInt3": "900 eV fit prior 900 100"},
          screenings(Project.open(path)))
    b = Project.open(os.path.join(SRC, "bare.azr"))
    check("bare.azr: no screening", screenings(b) == {}, screenings(b))
    b.add("target", "scr7", kind="none", screening=Nuisance(500, "fit", unit="eV", prior=(500, 50)))
    b.set("skowronski_luna_hpge", target="scr7")
    check("a new target with screening, named by a data set",
          screenings(b) == {"scr7": "500 eV fit prior 500 50"} and b["skowronski_luna_hpge"]["target"] == "scr7",
          screenings(b))
    try:
        b.set("scr7", screening="700 MeV/sr fixed")
        check("a screening of the wrong unit is refused", False)
    except FormatError as err:
        check("a screening of the wrong unit is refused", "E_UNKNOWN_UNIT" in err.codes and
              screenings(b)["scr7"] == "500 eV fit prior 500 50", err.codes)

    # --- engine ----------------------------------------------------------------------------------
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import numpy as np
        from pyazr import Session
        from pyazr import _azure2  # noqa: F401  the compiled engine
    except Exception as err:  # engine not built
        print(f"skip engine part: {type(err).__name__}: {err}")
        session = None
        engine_skipped = True
    else:
        work = os.path.join(tmp, "grad")
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        shutil.copy(os.path.join(SRC, "grad.azr"), work)
        os.symlink(DATA, os.path.join(work, "data"))
        session = Session(os.path.join(work, "grad.azr"))

    if session is not None:
        scr = list(session.parameters.screenings)
        check("one screening parameter", len(scr) == 1, scr)
        p = scr[0]
        check("tgt_targetInt3.screening: free, U_e 1000 eV, prior 1000 +- 200, effect 1",
              p.name == "tgt_targetInt3.screening" and not p.fixed and close(p.value, 1000.0, 1e-12)
              and p.effect == 1 and close(p.prior_nominal, 1000.0, 1e-12) and close(p.prior_sigma, 200.0, 1e-12), p)

        x = np.asarray(session.params_rwa, float).copy()
        k = p.free_index
        x[k] = 1100.0
        pen = session.evaluate(x).terms["screening priors"]
        check("terms['screening priors'] = ((1100 - 1000)/200)^2", close(pen, 0.25, 1e-12), pen)
        others = session.evaluate(x).data
        rest = (float(np.sum([d.norm_chi2 for d in others])) + float(np.sum([d.shift_constant_chi2 for d in others])) +
                float(np.sum([d.shift_sqrt_chi2 for d in others])) + float(np.sum([d.shift_function_chi2 for d in others])))
        obj = session.objective()(x)
        data = session.evaluate(x).chi2
        check("objective = data chi2 + norm/shift terms + screening prior", close(obj, data + rest + 0.25, 1e-12),
              (obj, data, rest))
        x0 = np.asarray(session.params_rwa, float)
        check("U_e moves the data chi2", abs(data - session.evaluate(x0).chi2) > 1e-6)

        value, grad = session.objective().value_and_grad(x)
        h = 1e-3 * max(1.0, abs(x[k]))
        e = np.zeros_like(x)
        e[k] = h
        fd = (session.objective()(x + e) - session.objective()(x - e)) / (2 * h)
        check(f"objective gradient: d f/d U_e = finite difference ({grad[k]:.6g} vs {fd:.6g})",
              abs(grad[k] - fd) <= 1e-6 * max(1.0, abs(fd)), (grad[k], fd))

        r, J = session.objective().residuals_and_jacobian(x)
        r, J = np.asarray(r, float), np.asarray(J, float)
        rp = np.asarray(session.objective().residuals_and_jacobian(x + e)[0], float)
        rm = np.asarray(session.objective().residuals_and_jacobian(x - e)[0], float)
        col = (rp - rm) / (2 * h)
        worst = float(np.max(np.abs(J[:, k] - col) / np.maximum(1e-12, np.max(np.abs(col)))))
        check(f"objective().residuals_and_jacobian: U_e column = finite differences ({worst:.1e})", worst < 1e-6)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
if engine_skipped:
    # thm's convention (tests/pyazr/CMakeLists.txt): ctest lists the test as Skipped.
    print("PASS (pure-Python part; engine part skipped)")
    sys.exit(77)
print("PASS")
