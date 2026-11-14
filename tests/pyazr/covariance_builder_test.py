#!/usr/bin/env python3
"""The data covariance builder (pyazr/covariance.py).

  1. each kind of part on three rows, against matrices worked by hand (no engine);
  2. the checks: a scale part alone is refused (the engine needs every fitted row's own variance), a correlation
     matrix that is not positive semi-definite is refused naming its part, a remainder larger than the column is
     refused, a wrong field is refused with the kind's fields;
  3. the file: the text of a small case line by line, the parts read back from it;
  4. with the engine: the rows are the engine's (61 rows of p12c_60.dat, the window of the data set), and the
     engine's chi-squared of the data set with the written file equals v^T C^-1 v, v = model - data from a session
     without it and C formed here (own errors, a 5 % model scale, a 3 % part over a correlation length of 0.2 MeV);
     the engine replaces the error column by sqrt(C_ii) (the session's errors);
  5. the example pyazr/examples/covariance_builder.py runs and gives that chi-squared.

Run from anywhere:  python3 tests/pyazr/covariance_builder_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import numpy as np
    from pyazr import covariance as cv
except Exception as err:
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)


def rows3():
    return cv.Rows(files=[("d.dat", 3)], dataset=np.array(["s", "s", "s"], dtype=object), record=np.array([1, 2, 3]),
                   energy=np.array([1.0, 1.5, 3.0]), angle=np.zeros(3), value=np.array([2.0, 4.0, 5.0]),
                   error=np.array([0.2, 0.3, 0.5]), fitted=np.array([True, True, True]), unit="b", norms={})


def dense(parts, rows=None):
    rows = rows or rows3()
    parts = [cv._normal_part(p) for p in parts]
    terms = cv._terms(rows, parts, ".")
    return sum((t.dense(rows) for t in terms), np.zeros((3, 3))), terms


print("1. the parts against matrices worked by hand")
C, _ = dense([{"kind": "own"}])
check("own (column): diag(0.04, 0.09, 0.25)", np.allclose(C, np.diag([0.04, 0.09, 0.25]), rtol=0, atol=1e-15), C)
C, _ = dense([{"kind": "own", "source": "percent", "percent": 10}])
check("own (10 % of the value): diag(0.04, 0.16, 0.25)", np.allclose(C, np.diag([0.04, 0.16, 0.25]), atol=1e-15), C)
C, t = dense([{"kind": "scale", "percent": 10, "basis": "data"}])
hand = np.array([[0.04, 0.08, 0.10], [0.08, 0.16, 0.20], [0.10, 0.20, 0.25]])
check("scale 10 % (data): u u^T, u = (0.2, 0.4, 0.5); a vector term", np.allclose(C, hand, atol=1e-15)
      and t[0].kind == "vector", (C, t[0].kind))
C, t = dense([{"kind": "scale", "percent": 10}])
check("scale 10 % (model): the same with the data standing in; a relvector term of fractions 0.1",
      np.allclose(C, hand, atol=1e-15) and t[0].kind == "relvector" and np.allclose(t[0].values, 0.1), t[0].values)
C, t = dense([{"kind": "range", "percent": 10, "e_min": 1.2, "e_max": 3.5}])
hand_r = np.array([[0, 0, 0], [0, 0.16, 0.20], [0, 0.20, 0.25]])
check("range 10 % on 1.2-3.5 MeV: rows 2 and 3 only", np.allclose(C, hand_r, atol=1e-15)
      and np.allclose(t[0].values, [0, 0.1, 0.1]), C)
C, t = dense([{"kind": "length", "percent": 10, "length": 1.0}])
hand_l = np.array([[0.04, 0.2 * 0.4 * 0.60653065971263342, 0.2 * 0.5 * 0.13533528323661270],
                   [0.2 * 0.4 * 0.60653065971263342, 0.16, 0.4 * 0.5 * 0.22313016014842982],
                   [0.2 * 0.5 * 0.13533528323661270, 0.4 * 0.5 * 0.22313016014842982, 0.25]])
check("length 10 %, 1 MeV, exponential: s_i s_j exp(-|dE|) (C12 = 0.0485225, C13 = 0.0135335, C23 = 0.0446260)",
      np.allclose(C, hand_l, rtol=1e-15, atol=0) and t[0].kind == "block", C)
C, _ = dense([{"kind": "length", "percent": 10, "length": 1.0, "shape": "gaussian", "e_max": 2.0}])
g = 0.2 * 0.4 * np.exp(-0.125)
check("length, gaussian, up to 2 MeV: a 2 x 2 block, C12 = 0.08 exp(-1/8)",
      np.allclose(C, [[0.04, g, 0], [g, 0.16, 0], [0, 0, 0]], rtol=1e-15, atol=0), C)
C, _ = dense([{"kind": "scale", "percent": 5, "basis": "data"}, {"kind": "own", "source": "remainder"}])
check("own (remainder) beside a 5 % scale: diagonal 0.04-0.01, 0.09-0.04, 0.25-0.0625; C_ii = the column squared",
      np.allclose(np.diag(C), [0.04, 0.09, 0.25], atol=1e-15)
      and np.allclose(C - np.outer([0.1, 0.2, 0.25], [0.1, 0.2, 0.25]), np.diag([0.03, 0.05, 0.1875]), atol=1e-15), C)
work = tempfile.mkdtemp(prefix="pyazr_covbuild_")
rho = np.array([[1, 0.5, 0], [0.5, 1, -0.2], [0, -0.2, 1]])
np.savetxt(os.path.join(work, "rho.txt"), rho)
C, _ = dense([{"kind": "correlation", "file": os.path.join(work, "rho.txt")}])
s = np.array([0.2, 0.3, 0.5])
check("correlation (error column): s rho s", np.allclose(C, s[:, None] * rho * s[None, :], rtol=1e-15, atol=0), C)

print("\n2. the checks")


class FakeProject:
    directory = work


def built(parts, rows=None):
    """build() on rows3() without the engine (rows_of replaced)."""
    keep = cv.rows_of
    cv.rows_of = lambda project, datasets: rows or rows3()
    try:
        return cv.build(FakeProject(), ["s"], parts)
    finally:
        cv.rows_of = keep


c = built([{"kind": "scale", "percent": 5}])
check("a scale part alone is refused: no own variance; the part named", c.refusal is not None
      and "own" in c.refusal and c.refusal_parts == [0], (c.refusal, c.refusal_parts))
np.savetxt(os.path.join(work, "bad.txt"), np.array([[1, 0.9, -0.9], [0.9, 1, 0.9], [-0.9, 0.9, 1]]))
c = built([{"kind": "own"}, {"kind": "correlation", "file": "bad.txt", "percent": 10}])
check("a correlation that is not positive semi-definite is refused, naming part 2, though own errors hide it",
      c.refusal is not None and c.refusal_parts == [1] and "2 (correlation)" in c.refusal, (c.refusal, c.refusal_parts))
c = built([{"kind": "own"}, {"kind": "scale", "percent": 5}])
check("own + 5 % scale: positive definite, no refusal, 3 rows off the column by more than 1 %",
      c.refusal is None and c.checks["positive_definite"] and c.checks["own_variance"]
      and c.checks["sigma_mismatch"] == 3, c.checks)
for parts, words in (([{"kind": "scale", "percent": 30, "basis": "data"}, {"kind": "own", "source": "remainder"}],
                      "exceed the error column"),
                     ([{"kind": "scale", "size": 3}], "no field size"),
                     ([{"kind": "scael", "percent": 3}], "unknown kind"),
                     ([{"kind": "length", "percent": 3}], "needs length"),
                     ([{"kind": "own", "source": "percent"}], "needs percent")):
    try:
        built(parts)
        check(f"refused: {words}", False, "no error")
    except ValueError as e:
        check(f"refused: {words}", words in str(e), str(e))
try:
    c.write(os.path.join(work, "x.txt")) if False else built([{"kind": "scale", "percent": 5}]).write(
        os.path.join(work, "x.txt"))
    check("a refused covariance is not written", False)
except ValueError as e:
    check("a refused covariance is not written", not os.path.exists(os.path.join(work, "x.txt")), str(e))

print("\n3. the file")
c = built([{"kind": "own"}, {"kind": "range", "percent": 10, "e_min": 1.2, "e_max": 3.5},
           {"kind": "length", "percent": 10, "length": 1.0, "e_max": 2.0}])
lines = c.text().splitlines()
expect = ["file d.dat 3", "# part 1 (own)", "diagonal", "0.040000000000000008 0.089999999999999997 0.25",
          "# part 2 (range)", "relvector", "0 0.10000000000000001 0.10000000000000001", "# part 3 (length)", "block 2",
          "1 2", "0.040000000000000008 0.048522452777010679", "0.048522452777010679 0.16000000000000003"]
check("the lines after the parts' line", lines[2:] == expect, "\n".join(lines))
path = c.write(os.path.join(work, "c.txt"))
back = cv.read_parts(path)
check("the parts read back from the #pyazr.covariance line", back == {"datasets": ["s"], "parts": c.parts}, back)
with open(os.path.join(work, "hand.txt"), "w") as fh:
    fh.write("file d.dat 3\ndiagonal 1 1 1\n")
check("a file written by hand has no parts", cv.read_parts(os.path.join(work, "hand.txt")) is None)

print("\n4. with the engine")
try:
    from pyazr import Project, Session
    have_engine = True
except Exception as e:  # pragma: no cover
    print(f"  skip  engine not built ({e})")
    have_engine = False

if have_engine:
    d = os.path.join(work, "p")
    os.makedirs(os.path.join(d, "output"))
    os.makedirs(os.path.join(d, "checks"))
    shutil.copytree(os.path.join(ROOT, "tests", "angle_mapping", "data"), os.path.join(d, "data"))
    shutil.copy(os.path.join(ROOT, "tests", "data_covariance", "plain.azr"), os.path.join(d, "t.azr"))
    p = Project.open(os.path.join(d, "t.azr"))
    rows = cv.rows_of(p, ["p12c_60"])
    raw = np.loadtxt(os.path.join(d, "data", "p12c_60.dat"), comments="#", ndmin=2)
    check("the rows are the file's 61 data rows, values and errors as written",
          rows.n == 61 and rows.files == [("data/p12c_60.dat", 61)] and np.array_equal(rows.value, raw[:, 2])
          and np.array_equal(rows.error, raw[:, 3]) and np.array_equal(rows.record, np.arange(1, 62)),
          (rows.n, rows.files))
    parts = [{"kind": "own"}, {"kind": "scale", "percent": 5}, {"kind": "length", "percent": 3, "length": 0.2}]
    c = cv.build(p, ["p12c_60"], parts)
    check("p12c_60 with own, 5 % scale, 3 % over 0.2 MeV: no refusal; a note on the fitted 5 % norm prior",
          c.refusal is None and any("counted twice" in n for n in c.notes), (c.refusal, c.notes))
    with Session(os.path.join(d, "t.azr")) as s:
        x0 = np.asarray(s.parameters.x0, float)
        ev = s.evaluate(x0)
        k = [i for i, e in enumerate(ev.data) if e.name == "p12c_60"][0]
        r = np.asarray(ev.data[k].residuals, float)
        e_lab = np.asarray(ev.data[k].extra["lab_energy"], float) if "lab_energy" in ev.data[k].extra else None
    fit = c.fitted
    check("the session's points are the fitted rows", len(r) == fit.size, (len(r), fit.size))
    if e_lab is not None:
        check("... in file order (lab energies)", np.allclose(e_lab, rows.energy[fit], rtol=1e-12), e_lab[:3])
    v = r * rows.error[fit]                      # model - data in the file's units (n cancels: M = D^-1 C D^-1)
    m = rows.value[fit] + v                      # the model the relvector term freezes at the first evaluation
    terms = c.terms
    Cm = np.zeros((fit.size, fit.size))
    for t in terms:
        if t.kind == "relvector":
            u = t.values[fit] * m
            Cm += np.outer(u, u)
        else:
            Cm += t.dense(rows)[np.ix_(fit, fit)]
    formula = float(v @ np.linalg.solve(Cm, v))
    c.write(os.path.join(d, "cov_built.txt"))
    p.set("p12c_60", covariance="cov_built.txt")
    p.save(os.path.join(d, "t_cov.azr"))
    with Session(os.path.join(d, "t_cov.azr")) as s:
        ev2 = s.evaluate(np.asarray(s.parameters.x0, float))
        k2 = [i for i, e in enumerate(ev2.data) if e.name == "p12c_60"][0]
        engine = float(ev2.data[k2].chi2)
        errs = np.asarray(ev2.data[k2].errors, float) / np.asarray(ev.data[k].errors, float) * rows.error[fit]
    check(f"the engine's chi2 with the file = v^T C^-1 v ({engine:.12g} against {formula:.12g})",
          abs(engine - formula) < 1e-9 * formula, (engine, formula))
    check("... and differs from the chi2 without it", abs(engine - float(ev.data[k].chi2)) > 1e-3 * engine)
    check("the engine's error column is sqrt(C_ii) with the data standing in (the builder's total error)",
          np.allclose(errs, c.total_error[fit], rtol=1e-12), (errs[:3], c.total_error[fit][:3]))

    print("\n5. the example pyazr/examples/covariance_builder.py")
    import subprocess
    out = subprocess.run([sys.executable, os.path.join(ROOT, "pyazr", "examples", "covariance_builder.py"),
                          os.path.join(d, "t.azr"), "p12c_60", "--out", os.path.join(d, "ex.azr")],
                         capture_output=True, text=True, timeout=600, cwd=d,
                         env=dict(os.environ, PYTHONPATH=ROOT + os.pathsep + os.environ.get("PYTHONPATH", "")))
    with_line = [ln for ln in out.stdout.splitlines() if ln.startswith("chi2 of p12c_60")]
    check("it runs, writes the file and the project, and prints both chi2",
          out.returncode == 0 and os.path.isfile(os.path.join(d, "cov_p12c_60.txt"))
          and os.path.isfile(os.path.join(d, "ex.azr")) and len(with_line) == 1, out.stdout[-600:] + out.stderr[-600:])
    check("... with the engine's chi2 of the test above", f"{engine:.10g} with" in with_line[0] if with_line else False,
          with_line)

shutil.rmtree(work, ignore_errors=True)
print(f"\n{len(failures)} failure(s)" if failures else "\nall ok")
sys.exit(1 if failures else 0)
