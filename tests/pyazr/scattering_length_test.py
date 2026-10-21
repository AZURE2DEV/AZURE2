#!/usr/bin/env python3
"""Scattering lengths through pyazr: the <scatteringLengths> block helpers of AzrModel (pure Python),
azure2.scattering_lengths against the reference values of tests/scattering_length, the constraint term
in calculate_chi2_rwa / penalties(), chi2_and_grad and the residual_jacobian rows against central
finite differences.  Model: tests/scattering_length/be7p_fit.azr (7Be+p, every kind of constraint line).

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


# --- AzrModel: the block (no engine) ------------------------------------------------------------
try:
    from pyazr.azrfile import AzrModel, ScatteringLengthLine  # noqa: E402
except ImportError as e:  # as the other pyazr tests: no package dependencies (e.g. mpmath) -> nothing to test
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(77)

tmp = tempfile.mkdtemp(prefix="pyazr_scatlen_")
try:
    m = AzrModel.from_file(os.path.join(SRC, "be7p_fit.azr"))
    lines = m.scattering_lengths
    check("block read: 5 lines", len(lines) == 5, lines)
    check("coherent line s = -1", lines[2].s == -1 and lines[2].value == 5.0, lines[2])
    check("complex line", lines[3].im_value == -0.3 and lines[3].im_sigma == 0.05, lines[3])
    check("inactive line", not lines[4].active, lines[4])
    m.add_scattering_length(1, 2, -5.0, 0.2)
    m.set_scattering_lengths([e for e in m.scattering_lengths if e.active])
    out = os.path.join(tmp, "rt.azr")
    m.write(out)
    back = AzrModel.from_file(out).scattering_lengths
    check("round trip: inactive dropped, one added", len(back) == 5 and back[-1] == ScatteringLengthLine(1, 2.0, -5.0, 0.2),
          back)
    text = open(out).read()
    check("block kept where it was (after </segmentsData>, before <segmentsTest>)",
          text.index("</segmentsData>") < text.index("<scatteringLengths>") < text.index("<segmentsTest>"))
    plain = AzrModel.from_file(os.path.join(SRC, "be7p.azr"))
    check("no block: empty list", plain.scattering_lengths == [])
    plain.add_scattering_length(1, 1, 17.0, 0.5)
    plain.write(os.path.join(tmp, "new.azr"))
    t2 = open(os.path.join(tmp, "new.azr")).read()
    check("new block after </segmentsTest>",
          t2.index("</segmentsTest>") < t2.index("<scatteringLengths>") < t2.index("<targetInt>"))
    plain.set_scattering_lengths([])
    plain.write(os.path.join(tmp, "none.azr"))
    check("empty list removes the block", "<scatteringLengths>" not in open(os.path.join(tmp, "none.azr")).read())

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
        check("penalties()['scattering_length'] = 48.3471 (the CLI's)", close(pen, 48.3471155, 1e-6), pen)
        total = float(np.sum(fit.calculate_chi2_rwa(x)))
        check("calculate_chi2_rwa includes it (275.3872)", close(total, 275.3872387, 1e-7), total)
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
