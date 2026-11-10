#!/usr/bin/env python3
"""The gradient of composite data sets whose terms read different parts of a cross section (data_model.md S11).

A SUM whose leading term is E1-only (isDiff 5) or E2-only (isDiff 6) reads, in the forward model, the E1 / E2 part of
its own point and the FULL cross section of each component (src/ModelPass.cpp; before S10, ComposeModelValue:
`compValue = componentPoint->GetFitCrossSection(results)`).  The analytic gradient must differentiate the same thing.
On 12C(p,g)13N (tests/13N levels; the 3/2- level decays by E1 and M1/E2, so sigma_E1, sigma_E2 and sigma differ), with
synthetic data at 0.4-1.7 MeV, the gradient of the objective (pyazr objective_and_grad) is compared with central
differences of the objective, for E1 and E2 SUMs, by the shared-component route and by the per-term route
(AZURE_NO_COMPONENT_MAPPING=1, a process of its own).  Tolerance: |analytic - FD| <= 1e-5 of the largest gradient
entry (the differences themselves agree to ~1e-8 here).

  PYAZR_ENGINE_DIR=<build>/pyazr-out python3 tests/pyazr/composite_gradient_test.py
Needs the compiled engine; skips (exit 77) without it.
"""
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
TOL = 1e-5

try:
    sys.path.insert(0, ROOT)
    import numpy as np
    import pyazr
    from pyazr import _azure2  # noqa: F401
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

ENERGIES = [0.40, 0.44, 0.46, 0.50, 1.40, 1.52, 1.56, 1.60, 1.70]


def make_project(work, isdiff):
    """tests/13N with one data set: a SUM of the E1 (isDiff 5) or E2 (6) part of p + 12C -> 13N + g and, as its
    component, the full cross section of the same pair."""
    with open(os.path.join(ROOT, "tests", "13N", "13N.azr")) as f:
        text = f.read()
    line = f"1 1 2 0.3 1.8 0 180 {isdiff} 1 0 0 0 0 0 data/sum.dat 1 0 1 1 2 -999\n"
    text = re.sub(r"<segmentsData>.*?</segmentsData>", "<segmentsData>\n" + line + "</segmentsData>", text, flags=re.S)
    os.makedirs(os.path.join(work, "data"), exist_ok=True)
    os.makedirs(os.path.join(work, "output"), exist_ok=True)
    os.makedirs(os.path.join(work, "checks"), exist_ok=True)
    with open(os.path.join(work, "data", "sum.dat"), "w") as f:
        for e in ENERGIES:
            f.write(f"{e} 0 1.0e-6 1.0e-7\n")
    with open(os.path.join(work, "p.azr"), "w") as f:
        f.write(text)


def check_one(work):
    """Analytic gradient against central differences; prints and returns the relative error."""
    os.chdir(work)
    m = pyazr.azure2("p.azr")
    x = np.asarray(m.params_rwa, float)
    f0, g = m.objective_and_grad(x)
    fd = np.zeros_like(g)
    for k in range(len(x)):
        h = 1e-6 * max(abs(x[k]), 1e-3)
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] -= h
        fd[k] = (m.objective(xp) - m.objective(xm)) / (2 * h)
    err = np.max(np.abs(g - fd)) / np.max(np.abs(fd))
    return err, f0, len(x)


if __name__ == "__main__" and len(sys.argv) == 2:  # one case, in its own process (the route is read at start)
    err, f0, n = check_one(sys.argv[1])
    print(f"{err:.3e} {f0:.6e} {n}")
    sys.exit(0)

failures = []
with tempfile.TemporaryDirectory() as tmp:
    for isdiff, part in ((5, "E1"), (6, "E2")):
        for route, env in (("shared-component route", {}), ("per-term route", {"AZURE_NO_COMPONENT_MAPPING": "1"})):
            work = os.path.join(tmp, f"{part}_{len(env)}")
            make_project(work, isdiff)
            r = subprocess.run([sys.executable, os.path.abspath(__file__), work], stdout=subprocess.PIPE,
                               stderr=subprocess.DEVNULL, env=dict(os.environ, OMP_NUM_THREADS="1", **env), text=True)
            last = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else ""
            try:
                err, f0, n = last.split()
                ok = float(err) <= TOL
            except ValueError:
                err, ok = f"no result (exit {r.returncode})", False
            print(f"  {'ok  ' if ok else 'FAIL'}  {part}-only SUM + full component, {route}: gradient vs central "
                  f"differences {err} of the largest (tolerance {TOL:g})")
            if not ok:
                failures.append(f"{part} {route}")
if failures:
    print(f"{len(failures)} FAILURE(S)")
    sys.exit(1)
print("all passed")
