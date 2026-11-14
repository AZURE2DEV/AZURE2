#!/usr/bin/env python3
"""Self-test of the azure2-calculate scripts: extrapolate.py, decompose.py.

On a copy of the synthetic project T1 (tests/upgrade/T1: p + 12C -> 13N, capture to two final states) in a temporary
directory:

- extrapolate.py on a c.m. grid at the energies of the capture data set gives the data-mode model of that data set (to
  1e-8: the same calculation; the grid passes through a c.m. to lab conversion), and its S column equals sigma E exp(2
  pi eta) computed here from the masses (CODATA 2018 constants; to 1e-6);
- decompose.py: the model with the 1/2+ level's widths zeroed equals that of a copy of the project without the level
  (``Project.remove_level``; to 1e-6, measured 3e-8), the level does something at its resonance (full - off), and the
  pairwise term of the 1/2+ and the 3/2- level, which have different J^pi, vanishes for the angle-integrated capture
  cross section (to 1e-9 of the cross section).

Needs pyazr with its engine (PYAZR_ENGINE_DIR) and AZURE2_COULOMB_TABLE as the pyazr tests have them.
PASS/FAIL lines; exit 0, 1 on a failure, 77 (skip) without the engine. About 10 s at one thread.
"""
import argparse
import math
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
T1 = REPO / "tests" / "upgrade" / "T1"
FAILED = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILED.append(name)


def run(script, *args, cwd):
    t = time.time()
    r = subprocess.run([sys.executable, str(HERE / script), *map(str, args)], cwd=cwd, capture_output=True, text=True)
    print(f"  ran {script} in {time.time() - t:.1f} s (exit {r.returncode})")
    if r.returncode:
        print(r.stdout[-2000:], r.stderr[-2000:])
    return r


def table(path):
    return np.loadtxt(path, comments="#", ndmin=2)


def capture_set(work):
    """(name, c.m. energies, model) of T1's ground-state capture data set (angle-integrated sigma)."""
    import pyazr
    with pyazr.Session(str(work / "T1.azr")) as s:
        ev = s.evaluate(s.parameters.x0)
        d = next(d for d in ev.data if d.quantity == "sigma" and d.n > 5)
        return d.name, np.asarray(d.energies, float), np.asarray(d.model, float)


def test_extrapolate(work):
    import pyazr
    name, e_cm, model = capture_set(work)
    reaction = pyazr.Project.open(work / "T1.azr")[name]["reaction"]
    grid = "list " + " ".join(repr(float(e)) for e in e_cm) + " MeV cm"
    r = run("extrapolate.py", "T1.azr", "--reaction", reaction, "--energies", grid, "--out", "x.txt", cwd=work)
    if r.returncode:
        check("extrapolate: runs", False)
        return
    t = table(work / "x.txt")
    rel = np.max(np.abs(t[:, 3] / model - 1))
    check("extrapolate: the grid at the data energies gives the data-mode model", rel < 1e-8 and
          np.allclose(t[:, 0], e_cm, rtol=1e-9), f"largest rel. difference {rel:.1e} over {len(e_cm)} points")
    pair = next(q for q in pyazr.Project.open(work / "T1.azr").pairs() if q.key == 1)
    alpha, amu = 1 / 137.035999084, 931.49410242
    mu = pair.light_mass * pair.heavy_mass / (pair.light_mass + pair.heavy_mass) * amu
    eta = pair.light[0] * pair.heavy[0] * alpha * np.sqrt(mu / (2 * t[:, 0]))
    s_ref = t[:, 3] * t[:, 0] * np.exp(2 * math.pi * eta)
    rel = np.max(np.abs(t[:, 4] / s_ref - 1))
    check("extrapolate: S = sigma E exp(2 pi eta)", rel < 1e-6, f"largest rel. difference {rel:.1e}")


def test_decompose(work):
    name, e_cm, model = capture_set(work)
    r = run("decompose.py", "T1.azr", "--level", "1/2+ 2.365", "--with", "3/2- 3.502", "--dataset", name,
            "--out", "d.txt", cwd=work)
    if r.returncode:
        check("decompose: runs", False)
        return
    t = table(work / "d.txt")
    full, off, pair_term = t[:, 2], t[:, 3], t[:, 6]
    check("decompose: full is the data-mode model", np.allclose(full, model, rtol=1e-12))
    check("decompose: the 1/2+ level carries the cross section at its resonance",
          np.max(np.abs(full - off) / full) > 0.5, f"max (full - off)/full {np.max(np.abs(full - off) / full):.3f}")
    import pyazr
    p = pyazr.Project.open(work / "T1.azr")
    p.remove_level(("1/2+", 2.365))
    p.save(work / "T1_no12.azr")
    with pyazr.Session(str(work / "T1_no12.azr")) as s:
        removed = next(d for d in s.evaluate(s.parameters.x0).data if d.name == name).model
    rel = np.max(np.abs(off / removed - 1))
    check("decompose: 'off' equals the project with the level removed", rel < 1e-6,
          f"largest rel. difference {rel:.1e}")
    worst = np.max(np.abs(pair_term)) / np.max(full)
    check("decompose: levels of different J^pi do not interfere in an angle-integrated sigma", worst < 1e-9,
          f"max |pair term| / max sigma = {worst:.1e}")


def main():
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    try:
        import pyazr
        pyazr.describe()
    except Exception as err:
        print(f"skip: the pyazr engine is not available ({type(err).__name__}: {err})")
        return 77
    t = time.time()
    with tempfile.TemporaryDirectory(prefix="azure2_calculate_selftest_") as tmp:
        work = Path(tmp)
        shutil.copytree(T1 / "data", work / "data")
        shutil.copy(T1 / "T1.azr", work / "T1.azr")
        for d in ("output", "checks"):
            (work / d).mkdir()
        for test in (test_extrapolate, test_decompose):
            try:
                test(work)
            except Exception as err:
                check(test.__name__, False, f"{type(err).__name__}: {err}")
    print(f"{'FAIL' if FAILED else 'PASS'} azure2-calculate selftest ({len(FAILED)} failed, {time.time() - t:.1f} s)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
