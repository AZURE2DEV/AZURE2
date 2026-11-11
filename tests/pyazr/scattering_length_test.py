#!/usr/bin/env python3
"""Scattering lengths through pyazr: the constraints as `datum` records of <modelData> (Project: read, added,
removed),
azure2.scattering_lengths against the reference values of tests/scattering_length, the constraint term
in calculate_chi2_rwa / penalties(), chi2_and_grad and the residual_jacobian rows against central
finite differences.  Model: tests/scattering_length/be7p_fit.azr (7Be+p, every kind of constraint; in format 3
the constraints are scattering_length data of <modelData>).

The engine part needs the compiled module and skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/scattering_length_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "scattering_length")
sys.path.insert(0, ROOT)

failures = []
engine_skipped = False


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def close(a, b, tol):
    return abs(a - b) <= tol * max(abs(b), 1e-300)


# --- the constraints as records (no engine session) ----------------------------------------------
try:
    from pyazr import Project  # noqa: E402
except ImportError as e:  # as the other pyazr tests: no package dependencies (e.g. mpmath) -> nothing to test
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(77)


def constraints(p):
    return [r for r in p.records().values() if r.record_kind == "datum" and r["quantity"] == "scattering_length"]


tmp = tempfile.mkdtemp(prefix="pyazr_scatlen_")
try:
    m = Project.open(os.path.join(SRC, "be7p_fit.azr"))
    lines = constraints(m)
    check("4 scattering_length data in <modelData>", len(lines) == 4, lines)
    check("coherent datum: s coherent, 5.0 fm", lines[2]["s"] == "coherent" and lines[2]["value"] == "5.0 fm",
          lines[2].fields)
    check("complex datum", lines[3]["imag"] == "-0.3 fm" and lines[3]["imag_error"] == "0.05 fm", lines[3].fields)
    m.add_datum("a_new", quantity="scattering_length", pair="p+7Be", s=2, value=(-5.0, "fm"), error=(0.2, "fm"))
    out = os.path.join(tmp, "rt.azr")
    m.save(out)
    back = constraints(Project.open(out))
    check("round trip: one added, at the end", len(back) == 5 and back[-1].fields ==
          {"quantity": "scattering_length", "pair": "p+7Be", "s": "2", "value": "-5.0 fm", "error": "0.2 fm"},
          [r.fields for r in back])
    plain = Project.open(os.path.join(SRC, "be7p.azr"))
    check("no <modelData>: no constraint", constraints(plain) == [] and "<modelData>" not in plain.text)
    plain.add_datum("a01", quantity="scattering_length", pair="p+7Be", s=1, value=(17.0, "fm"), error=(0.5, "fm"))
    plain.save(os.path.join(tmp, "new.azr"))
    check("a new <modelData> section holds it", len(constraints(Project.open(os.path.join(tmp, "new.azr")))) == 1)
    plain.remove("a01")
    plain.save(os.path.join(tmp, "none.azr"))
    check("removing the record removes the constraint",
          constraints(Project.open(os.path.join(tmp, "none.azr"))) == [])

    # --- engine ----------------------------------------------------------------------------------
    try:
        os.environ.setdefault("OMP_NUM_THREADS", "2")
        import numpy as np
        from pyazr import azure2
        from pyazr import _azure2  # noqa: F401  the compiled engine
    except Exception as err:  # engine not built
        print(f"skip engine part: {type(err).__name__}: {err}")
        session = None
        engine_skipped = True
    else:
        work = os.path.join(tmp, "be7p")
        shutil.copytree(SRC, work, ignore=shutil.ignore_patterns("output", "checks", "run.log"))
        for d in ("output", "checks"):
            os.makedirs(os.path.join(work, d), exist_ok=True)
        # be7p.azr has no data segments: an extrapolation-mode session
        session = azure2(os.path.join(work, "be7p.azr"), data_mode=False)

    if session is not None:
        values = {(v.pair, v.jpi): v for v in session.scattering_lengths()}
        check("a01 = 17.06776 fm", close(values[(1, "1-")].a0.real, 17.06776438, 1e-5), values[(1, "1-")])
        check("a02 = -5.01840 fm", close(values[(1, "2-")].a0.real, -5.018395912, 1e-5), values[(1, "2-")])
        check("7Be*+p 0-: hard sphere", not values[(2, "0-")].has_levels, values[(2, "0-")])
        check("7Be*+p 1-: complex", values[(2, "1-")].a0.imag < -1e-3, values[(2, "1-")])

        fit = azure2(os.path.join(work, "be7p_fit.azr"))
        x = np.asarray(fit.params_rwa, float)
        pen = fit.penalties(x)["scattering_length"]
        # CODATA 2018 constants since 2026-11-05 (were 48.3471155 and 275.3872387 with the CODATA 2006 ones)
        check("penalties()['scattering_length'] = 48.3472 (the CLI's)", close(pen, 48.3471861, 1e-6), pen)
        total = float(np.sum(fit.calculate_chi2_rwa(x)))
        check("calculate_chi2_rwa includes it (275.3878)", close(total, 275.3878019, 1e-7), total)
        check("objective = calculate_chi2_rwa (no norm/shift terms here)", close(fit.objective(x), total, 1e-12))

        value, grad = fit.chi2_and_grad(x)
        h = 1e-6 * np.maximum(1.0, np.abs(x))
        fd = np.array([(float(np.sum(fit.calculate_chi2_rwa(x + h[i] * e))) -
                        float(np.sum(fit.calculate_chi2_rwa(x - h[i] * e)))) / (2 * h[i])
                       for i, e in enumerate(np.eye(x.size))])
        worst = float(np.max(np.abs(grad - fd) / np.maximum(1.0, np.abs(fd))))
        check(f"chi2_and_grad vs finite differences ({worst:.1e})", worst < 1e-5, (grad, fd))

        r, J = fit.residual_jacobian(x)
        r, J = np.asarray(r, float), np.asarray(J, float)
        n_data = 34
        check("residual_jacobian: 34 data rows + 5 constraint rows", r.size == n_data + 5, r.size)
        check("sum of the constraint rows squared = the penalty", close(float(np.sum(r[n_data:] ** 2)), pen, 1e-8))

        def rows(z):
            return np.asarray(fit.residual_jacobian(z)[0], float)[n_data:]
        fdJ = np.array([(rows(x + h[i] * e) - rows(x - h[i] * e)) / (2 * h[i]) for i, e in enumerate(np.eye(x.size))]).T
        worstJ = float(np.max(np.abs(J[n_data:] - fdJ) / np.maximum(1.0, np.abs(fdJ))))
        check(f"constraint Jacobian rows vs finite differences ({worstJ:.1e})", worstJ < 1e-5)
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
