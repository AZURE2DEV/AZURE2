#!/usr/bin/env python3
"""Data covariance matrices as pyazr sees them.

A format-3 data set may name ``covariance <file>`` (include/DataCovariance.h, tests/data_covariance):
the data sets naming one file share one data covariance matrix and their chi2 is r^T M^-1 r.
pyazr must read and write the field and,
with the engine, give a chi2, gradient and whitened residual Jacobian consistent with it:

  1. the field is read, written on a new data set (Project.add) and read back;
  2. a diagonal covariance diag(e_i^2) gives the chi2 of no covariance;
  3. dense covariance: sum of the whitened data rows squared = the data chi2, the
     analytic gradient of the objective = central differences of the objective, and
     = 2 J^T w from the whitened residual Jacobian;
  4. relvector terms: the model is frozen at the first evaluation; refreeze_covariance()
     makes the next evaluation freeze it again (at that point chi2 changes).

Run from anywhere:  python3 tests/pyazr/data_covariance_test.py
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


try:                                           # pyazr itself needs numpy
    import numpy as np                         # noqa: E402
    from pyazr import Project                  # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)

PROJECT = os.path.join(ROOT, "tests", "data_covariance")

print("1. the covariance field is read and written")
dense = Project.open(os.path.join(PROJECT, "dense.azr"))
sets = [r for r in dense.records().values() if r.record_kind == "dataset"]
with_cov = [r.record_name for r in sets if r.get("covariance") == "cov_dense.txt"]
check("eight data sets name cov_dense.txt", len(with_cov) == 8, with_cov)
shifted = [r for r in sets if r.get("shift") == "G"]
check("... beside a shared energy scale, which is still read", len(shifted) == 2 and all(r.get("covariance") for r in shifted),
      [(r.record_name, r.get("shift"), r.get("covariance")) for r in shifted])
rt = tempfile.mkdtemp(prefix="pyazr_covariance_rt_")  # the project beside its data (tests/angle_mapping/data)
shutil.copy(os.path.join(PROJECT, "plain.azr"), rt)
os.symlink(os.path.join(ROOT, "tests", "angle_mapping", "data"), os.path.join(rt, "data"))
with open(os.path.join(rt, "cov.txt"), "w") as fh:  # the validator reads the file: 61 points of p12c_60
    fh.write("file data/p12c_60.dat 61\ndiagonal\n" + " ".join(["1e-2"] * 61) + "\n")
m = Project.open(os.path.join(rt, "plain.azr"))
spec = {k: v for k, v in m[[n for n, r in m.records().items() if r.record_kind == "dataset"][0]].fields.items()
        if k not in ("norm", "shift", "covariance")}
m.add("normalization", "N", value="1 fixed")
m.add("dataset", "extra", norm="N", covariance="cov.txt", **spec)
check("the new data set holds the fields", {"    norm           N", "    covariance     cov.txt"}
      <= set(m["extra"].text.split("\n")), m["extra"].text)
out = os.path.join(rt, "rt.azr")
m.save(out)
back = Project.open(out)["extra"]
check("... and reads it back", back["covariance"] == "cov.txt" and back["norm"] == "N", back.fields)
shutil.rmtree(rt, ignore_errors=True)

print("\n2-4. with the engine")
try:
    from pyazr import Session
    have_engine = True
except Exception as e:                              # pragma: no cover
    print(f"  skip  engine not built ({e})")
    have_engine = False


def project(work, model):
    d = os.path.join(work, model)
    os.makedirs(os.path.join(d, "output"))
    os.makedirs(os.path.join(d, "checks"))
    shutil.copytree(os.path.join(ROOT, "tests", "angle_mapping", "data"), os.path.join(d, "data"))
    for f in os.listdir(PROJECT):
        if f.startswith("cov_"):
            shutil.copy(os.path.join(PROJECT, f), d)
    shutil.copy(os.path.join(PROJECT, model + ".azr"), os.path.join(d, "t.azr"))
    return os.path.join(d, "t.azr"), d


if have_engine:
    work = tempfile.mkdtemp(prefix="pyazr_covariance_")
    try:
        chi = {}
        for model in ("plain", "diag"):
            azr, d = project(work, model)
            with Session(azr, cwd=d) as mdl:
                x0 = np.asarray(mdl.parameters.x0, float)
                chi[model] = (mdl.evaluate(x0).chi2, np.asarray(np.array([d.chi2 for d in mdl.evaluate(x0).data]), float))
        check("diag(e^2) covariance: chi2 = no covariance",
              abs(chi["diag"][0] - chi["plain"][0]) < 1e-9 * chi["plain"][0], (chi["diag"][0], chi["plain"][0]))
        check("... segment by segment", float(np.max(np.abs(chi["diag"][1] - chi["plain"][1]))) < 1e-9 * chi["plain"][0])

        azr, d = project(work, "dense")
        with Session(azr, cwd=d) as mdl:
            x0 = np.asarray(mdl.parameters.x0, float)
            ev0 = mdl.evaluate(x0)
            c0 = ev0.chi2
            w, J = mdl.objective().residuals_and_jacobian(x0)
            w = np.asarray(w, float)
            wd = w[:ev0.data[-1].rows.stop]  # the data rows of the objective's rows
            check("dense: sum of whitened residuals^2 = chi2", abs(float(wd @ wd) - c0) < 1e-9 * c0, (float(wd @ wd), c0))
            check("... and data[i].chi2 adds up to it", abs(float(np.sum([d.chi2 for d in ev0.data])) - c0) < 1e-9 * c0)
            check("... and differs from the plain chi2", abs(c0 - chi["plain"][0]) > 0.1 * chi["plain"][0], (c0, chi["plain"][0]))
            out = mdl.objective().value_and_grad(x0)
            cg, g = float(out[0]), np.asarray(out[1], float)
            f0 = mdl.objective()(x0)
            check("value_and_grad's value = the objective", abs(cg - f0) < 1e-9 * f0, (cg, f0))
            jtw = 2.0 * np.asarray(J, float).T @ w
            scale = float(np.max(np.abs(g)))
            check("gradient = 2 J^T w (whitened Jacobian)", float(np.max(np.abs(jtw - g))) < 1e-7 * scale,
                  float(np.max(np.abs(jtw - g))) / scale)
            # central differences at two steps, the closer one counted: the J = 3/2- level is
            # narrow, and its energy needs a 1e-9 relative step (as in tests/angle_mapping)
            diff = np.zeros_like(x0)
            for i in range(len(x0)):
                best = np.inf
                for step in (1e-7, 1e-9):
                    h = step * max(abs(x0[i]), 1e-2)
                    xp, xm = x0.copy(), x0.copy()
                    xp[i] += h
                    xm[i] -= h
                    fd = (mdl.objective()(xp) - mdl.objective()(xm)) / (2 * h)
                    best = min(best, abs(fd - g[i]))
                diff[i] = best
            worst = float(np.max(diff)) / scale
            check(f"gradient = central differences for all {len(x0)} free parameters", worst < 1e-5, worst)

        azr, d = project(work, "rel")
        with Session(azr, cwd=d) as mdl:
            x0 = np.asarray(mdl.parameters.x0, float)
            c0 = mdl.evaluate(x0).chi2            # freezes m at x0
            g0 = np.asarray(mdl.objective().value_and_grad(x0)[1], float)
            x1 = x0 - 1e-3 * g0 / max(float(np.max(np.abs(g0))), 1e-30) * np.maximum(np.abs(x0), 1e-2)
            c1 = mdl.evaluate(x1).chi2            # m still frozen at x0
            c1b = mdl.evaluate(x1).chi2
            check("rel: chi2 at a new point with the model frozen is reproducible", c1 == c1b, (c1, c1b))
            n = mdl.refreeze_covariance()
            check("refreeze_covariance() reports one relvector group", n == 1, n)
            c1r = mdl.evaluate(x1).chi2           # freezes m at x1
            check("... after it the model is frozen at the new point (chi2 changes)", abs(c1r - c1) > 1e-9 * c1,
                  (c1, c1r))
            mdl.refreeze_covariance()
            c0r = mdl.evaluate(x0).chi2           # back at x0
            check("... and refreezing at the start point gives the start chi2 back", abs(c0r - c0) < 1e-10 * c0, (c0r, c0))
    finally:
        shutil.rmtree(work, ignore_errors=True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): " + ", ".join(failures))
    sys.exit(1)
print("all checks passed")
