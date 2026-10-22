#!/usr/bin/env python3
"""Data covariance matrices as pyazr sees them.

A <segmentsData> line may end with ``covariance <file>`` (include/SegLine.h, include/
DataCovariance.h, tests/data_covariance): the segments naming one file share one data
covariance matrix and their chi2 is r^T M^-1 r.  pyazr must read and write the block and,
with the engine, give a chi2, gradient and whitened residual Jacobian consistent with it:

  1. the block is read, written (add_data_segment) and read back;
  2. a diagonal covariance diag(e_i^2) gives the chi2 of no covariance;
  3. dense covariance: sum of the whitened residuals squared = calculate_chi2_rwa, the
     analytic gradient (chi2_and_grad) = central differences of calculate_chi2_rwa, and
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
    from pyazr.datasets import SegmentSet      # noqa: E402
    from pyazr.azrfile import AzrModel         # noqa: E402
except Exception as err:                       # same skip as the other pyazr tests
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(0)

PROJECT = os.path.join(ROOT, "tests", "data_covariance")

print("1. the covariance block is read and written")
s = list(SegmentSet.from_file(os.path.join(PROJECT, "dense.azr")))
with_cov = [x.key for x in s if x.covariance == "cov_dense.txt"]
check("eight segments name cov_dense.txt", len(with_cov) == 8, with_cov)
shifted = [x for x in s if x.shift_group == "G"]
check("... next to a shiftgroup block, which is still read", len(shifted) == 2 and all(x.covariance for x in shifted),
      [(x.key, x.shift_group, x.covariance) for x in shifted])
m = AzrModel.from_file(os.path.join(PROJECT, "plain.azr"))
m.add_data_segment("data/p12c_90.dat", 1, 1, observable="differential", norm_group="N", covariance="cov.txt")
lines = [l for l in m.to_text().splitlines() if l.strip()]
last = lines[[i for i, l in enumerate(lines) if l.strip() == "</segmentsData>"][0] - 1]
check("add_data_segment writes the block", last.split()[-4:] == ["normgroup", "N", "covariance", "cov.txt"], last)
back = SegmentSet._parse(last.split(), 99)
check("... and reads it back", back.covariance == "cov.txt" and back.norm_group == "N", (back.covariance, back.norm_group))

print("\n2-4. with the engine")
try:
    from pyazr import azure2
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
            with azure2(azr, cwd=d) as mdl:
                x0 = np.asarray(mdl.params_rwa, float)
                chi[model] = (float(np.sum(mdl.calculate_chi2_rwa(x0))), np.asarray(mdl.segment_chi2(x0), float))
        check("diag(e^2) covariance: chi2 = no covariance",
              abs(chi["diag"][0] - chi["plain"][0]) < 1e-9 * chi["plain"][0], (chi["diag"][0], chi["plain"][0]))
        check("... segment by segment", float(np.max(np.abs(chi["diag"][1] - chi["plain"][1]))) < 1e-9 * chi["plain"][0])

        azr, d = project(work, "dense")
        with azure2(azr, cwd=d) as mdl:
            x0 = np.asarray(mdl.params_rwa, float)
            c0 = float(np.sum(mdl.calculate_chi2_rwa(x0)))
            w, J = mdl.residual_jacobian(x0)
            w = np.asarray(w, float)
            check("dense: sum of whitened residuals^2 = chi2", abs(float(w @ w) - c0) < 1e-9 * c0, (float(w @ w), c0))
            check("... and segment_chi2 adds up to it", abs(float(np.sum(mdl.segment_chi2(x0))) - c0) < 1e-9 * c0)
            check("... and differs from the plain chi2", abs(c0 - chi["plain"][0]) > 0.1 * chi["plain"][0], (c0, chi["plain"][0]))
            out = mdl.chi2_and_grad(x0)
            cg, g = float(out[0]), np.asarray(out[1], float)
            check("chi2_and_grad's chi2 = calculate_chi2_rwa", abs(cg - c0) < 1e-9 * c0, (cg, c0))
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
                    fd = (float(np.sum(mdl.calculate_chi2_rwa(xp))) - float(np.sum(mdl.calculate_chi2_rwa(xm)))) / (2 * h)
                    best = min(best, abs(fd - g[i]))
                diff[i] = best
            worst = float(np.max(diff)) / scale
            check(f"gradient = central differences for all {len(x0)} free parameters", worst < 1e-5, worst)

        azr, d = project(work, "rel")
        with azure2(azr, cwd=d) as mdl:
            x0 = np.asarray(mdl.params_rwa, float)
            c0 = float(np.sum(mdl.calculate_chi2_rwa(x0)))            # freezes m at x0
            g0 = np.asarray(mdl.chi2_and_grad(x0)[1], float)
            x1 = x0 - 1e-3 * g0 / max(float(np.max(np.abs(g0))), 1e-30) * np.maximum(np.abs(x0), 1e-2)
            c1 = float(np.sum(mdl.calculate_chi2_rwa(x1)))            # m still frozen at x0
            c1b = float(np.sum(mdl.calculate_chi2_rwa(x1)))
            check("rel: chi2 at a new point with the model frozen is reproducible", c1 == c1b, (c1, c1b))
            n = mdl.refreeze_covariance()
            check("refreeze_covariance() reports one relvector group", n == 1, n)
            c1r = float(np.sum(mdl.calculate_chi2_rwa(x1)))           # freezes m at x1
            check("... after it the model is frozen at the new point (chi2 changes)", abs(c1r - c1) > 1e-9 * c1,
                  (c1, c1r))
            mdl.refreeze_covariance()
            c0r = float(np.sum(mdl.calculate_chi2_rwa(x0)))           # back at x0
            check("... and refreezing at the start point gives the start chi2 back", abs(c0r - c0) < 1e-10 * c0, (c0r, c0))
    finally:
        shutil.rmtree(work, ignore_errors=True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): " + ", ".join(failures))
    sys.exit(1)
print("all checks passed")
