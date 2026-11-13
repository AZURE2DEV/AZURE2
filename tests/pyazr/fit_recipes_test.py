#!/usr/bin/env python3
"""The three fit recipes of pyazr/examples on pyazr.Objective: scipy, iminuit, emcee.

- T2 (tests/upgrade): the program's own fit (MIGRAD with its analytic gradient, menu 2, --max-calls=200) gives the
  minimum; from the project's start, iminuit's MIGRAD (fit_iminuit.fit, then HESSE) reaches it to 1e-6 relative,
  and scipy's L-BFGS-B with the gradient (fit_scipy.minimize) and least_squares on the objective's rows with their
  Jacobian (fit_scipy.least_squares) reach it or go below it (the program's MIGRAD stops 9e-4 above the minimum the
  two scipy recipes agree on to 1e-8: 297.8603 against 297.8594), and agree with each other to 1e-8.
- emcee (fit_emcee.sample) on the 15O gamma_decay record with a branching under a prior centred below 0
  (the near-zero variant of tests/reference/sampler_objective_test.cpp, three free branchings): a chain of 8
  walkers x 60 steps stays inside [0, 1] (the support), is the same for the same seed, and the same in a pool of
  two worker processes (one session each) as in this process.

iminuit and emcee parts are skipped when the package is missing.  Slow (about 4 min at one thread).

Run: python3 tests/pyazr/fit_recipes_test.py   (AZURE2_BIN: the program)
"""
import importlib
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
TESTS = os.path.join(ROOT, "tests")
EXAMPLES = os.path.join(ROOT, "pyazr", "examples")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


NEAR_ZERO = [
    ("    total_width  0.0385027 eV fit\n", "    total_width  0.0385027 eV fixed prior 0.040 0.004\n"),
    ("    remainder    g+15O[5.183]\n", "    remainder    g+15O[5.183] prior 0.15 0.01\n"),
    ("    branching    g+15O[6.1763] 0.5843746023006179 fit\n",
     "    branching    g+15O[6.1763] 0.5843746023006179 fixed prior 0.6 0.01\n"),
    ("    branching    g+15O[6.8594] 0.0014024990455214828 fit\n",
     "    branching    g+15O[6.8594] 0.0014024990455214828 fit prior -0.01 0.00001\n"),
]


def copy_dir(src, dst, project, edits=()):
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("output", "checks"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(dst, d), exist_ok=True)
    path = os.path.join(dst, project)
    text = open(path).read()
    for old, new in edits:
        assert old in text, old
        text = text.replace(old, new)
    open(path, "w").write(text)
    return path


def program_minimum(azr):
    proc = subprocess.run([os.environ["AZURE2_BIN"], "--no-readline", "--max-calls=200", os.path.basename(azr)],
                          input="2\n\n\n\n7\n", text=True, cwd=os.path.dirname(azr), capture_output=True, timeout=1800)
    found = re.findall(r"^Total Chi-Squared: (\S+)", proc.stdout, re.M)
    return float(found[-1]) if found else None


def main():
    try:
        sys.path.insert(0, ROOT)
        sys.path.insert(0, EXAMPLES)
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        import numpy as np
        import pyazr
        pyazr.kinds()
        fit_scipy = importlib.import_module("fit_scipy")
    except Exception as err:  # engine not built
        print(f"skip: engine not available ({type(err).__name__}: {err})")
        return 77
    if not os.environ.get("AZURE2_BIN"):
        print("skip: AZURE2_BIN (the program) is not set")
        return 77
    with tempfile.TemporaryDirectory() as tmp:
        t2 = os.path.join(TESTS, "upgrade", "T2")
        best = program_minimum(copy_dir(t2, os.path.join(tmp, "cli"), "T2.azr"))
        check(f"the program's MIGRAD minimum of T2: {best}", best is not None)
        if best is None:
            return 1

        def near(v):
            return abs(v - best) <= 1e-6 * best

        def reached(v):
            return v <= best * (1 + 1e-6)

        with pyazr.Session(copy_dir(t2, os.path.join(tmp, "api"), "T2.azr")) as s:
            f = s.objective()
            try:
                fit_iminuit = importlib.import_module("fit_iminuit")
                importlib.import_module("iminuit")
            except ImportError:
                print("  skip  iminuit is not installed")
                fit_iminuit = None
            if fit_iminuit is not None:
                m = fit_iminuit.fit(f)
                check(f"iminuit MIGRAD + HESSE: {m.fval:.12g} in {m.nfcn} calls ({m.ngrad} gradients), valid",
                      near(m.fval) and m.valid and m.accurate, f"want {best}")
            r, x = fit_scipy.minimize(f)
            check(f"scipy L-BFGS-B with the gradient: {r.fun:.12g} in {r.nfev} evaluations", reached(f(x)),
                  f"want {best}; {r.message}")
            q = fit_scipy.least_squares(f)
            check(f"scipy least_squares on the rows and Jacobian: {f(q.x):.12g} in {q.nfev} evaluations, "
                  f"{q.njev} Jacobians", reached(f(q.x)), f"want {best}; {q.message}")
            check(f"the two scipy recipes agree ({f(x) - f(q.x):.2e})", abs(f(x) - f(q.x)) <= 1e-8 * best)

        try:
            importlib.import_module("emcee")
        except ImportError:
            print("  skip  emcee is not installed")
            return 1 if failures else 0
        fit_emcee = importlib.import_module("fit_emcee")
        record = copy_dir(os.path.join(TESTS, "gamma_decay", "models"), os.path.join(tmp, "record"), "record.azr",
                          NEAR_ZERO)
        chains = []
        for processes in (0, 0, 2):
            sampler, names = fit_emcee.sample(record, nwalkers=8, nsteps=60, spread=0.3, seed=7, processes=processes)
            chains.append(sampler.get_chain())
        chain = chains[0]
        check(f"emcee: {chain.shape[0]} steps x {chain.shape[1]} walkers of {names}, all inside [0, 1] "
              f"(smallest {chain.min():.3g})",
              chain.shape == (60, 8, len(names)) and chain.min() >= 0 and chain.max() <= 1)
        check("emcee: the same seed gives the same chain", np.array_equal(chains[0], chains[1]))
        check("emcee: a pool of two processes (a session each) gives the same chain",
              np.array_equal(chains[0], chains[2]))
        check("emcee: the walkers moved (acceptance > 0)", np.mean(sampler.acceptance_fraction) > 0)
    if failures:
        print(f"FAIL: {len(failures)} check(s): {failures}")
        return 1
    print("ok: the fit recipes")
    return 0


if __name__ == "__main__":  # the pool's workers import this file without running it
    sys.exit(main())
