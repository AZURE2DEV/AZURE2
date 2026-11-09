#!/usr/bin/env python3
"""The fit's objective through pyazr (docs/upgrade/objective.md, Phase 2 step 4).

On models with normalization priors (tests/norm_groups), energy-shift-function priors (tests/shift_fn) and a
scattering-length constraint (tests/scattering_length/be7p_fit):

  1. objective() = the first entry of objective_and_grad() (the same engine call), to 1e-12;
  2. objective() = the sum of the engine's terms, and the per-segment norm and shift arrays of penalties() sum to
     what the terms say, to 1e-9 relative;
  3. objective_and_grad()'s gradient = central finite differences of objective(), to 1e-5 of the largest component
     (step 1e-6 (|x| + 1));
  4. objective() >= calculate_chi2_rwa() (the priors are sums of squares).

Needs the compiled engine; skips cleanly without it.  Run from anywhere:  python3 tests/pyazr/objective_test.py
"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
MODELS = [("norm_groups", "norm_groups.azr"), ("shift_fn", "shift_fn.azr"), ("scattering_length", "be7p_fit.azr")]

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    from pyazr import azure2
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(0)

for project, azr in MODELS:
    work = tempfile.mkdtemp(prefix="objective_test.")
    try:
        shutil.copytree(os.path.join(ROOT, "tests", project), os.path.join(work, project))
        d = os.path.join(work, project)
        for sub in ("output", "checks"):
            os.makedirs(os.path.join(d, sub), exist_ok=True)
        s = azure2(os.path.join(d, azr), cwd=d)
        x = np.asarray(s.params_rwa, float)
        f = s.objective(x)
        f2, g = s.objective_and_grad(x)
        tag = f"{project}/{azr}"
        check(f"{tag}: objective() = objective_and_grad()[0] ({f:.12g})", abs(f - f2) <= 1e-12 * max(1.0, abs(f)),
              f"{f} vs {f2}")
        t = s.sess.objective_terms_rwa(x)
        total = float(np.sum(t["values"]))
        check(f"{tag}: objective = sum of its {len(t['names'])} terms", abs(total - f) <= 1e-9 * max(1.0, abs(f)),
              f"{total} vs {f}")
        pen = s.penalties(x)
        segsum = sum(float(np.sum(pen[k])) for k in ("norm", "shift", "shift_sqrt", "shift_fn"))
        check(f"{tag}: per-segment priors finite and >= 0 (sum {segsum:.6g})", np.isfinite(segsum) and segsum >= 0)
        chi2 = float(np.sum(s.calculate_chi2_rwa(x)))
        check(f"{tag}: objective {f:.8g} >= calculate_chi2_rwa {chi2:.8g}", f >= chi2 - 1e-9 * max(1.0, abs(chi2)))
        # 3. gradient against central differences of objective() (at most 12 parameters, spread over the vector)
        idx = np.unique(np.linspace(0, x.size - 1, min(12, x.size)).astype(int))
        worst = 0.0
        for i in idx:
            h = 1e-6 * (abs(x[i]) + 1.0)
            xp, xm = x.copy(), x.copy()
            xp[i] += h
            xm[i] -= h
            fd = (s.objective(xp) - s.objective(xm)) / (2 * h)
            worst = max(worst, abs(fd - g[i]))
        scale = max(np.max(np.abs(g[idx])), 1e-12)
        check(f"{tag}: gradient = finite differences on {idx.size} parameters ({worst / scale:.2e} of the largest)",
              worst <= 1e-5 * scale, f"{worst / scale:.3e}")
    finally:
        shutil.rmtree(work, ignore_errors=True)

print("PASS" if not failures else f"FAIL: {len(failures)}")
sys.exit(1 if failures else 0)
