#!/usr/bin/env python3
"""Self-test of the azure2-plot scripts: plot_fit.py.

On a copy of the synthetic project T1 (tests/upgrade/T1) in a temporary directory, with the norm of two data sets set to
1.2 (fixed) and 0.9 (fitted), plot_fit.py writes one figure per data set, and the pulls it draws -- (model - datum
n)/(error n) from the evaluated arrays -- equal the engine's residuals of every data set (lab-frame differential data,
A_y, a phase shift and capture: the figure shows what the fit compares). It is run once more with ``--sfactor``, which
must not change the pulls.

Needs pyazr with its engine (PYAZR_ENGINE_DIR) and matplotlib. PASS/FAIL lines; exit 0, 1 on a failure, 77 (skip)
without the engine. About 10 s.
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

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
    print(f"  ran {script} {' '.join(map(str, args))} in {time.time() - t:.1f} s (exit {r.returncode})")
    if r.returncode:
        print(r.stdout[-2000:], r.stderr[-2000:])
    return r


def main():
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    try:
        import pyazr
        pyazr.describe()
    except Exception as err:
        print(f"skip: the pyazr engine is not available ({type(err).__name__}: {err})")
        return 77
    t = time.time()
    with tempfile.TemporaryDirectory(prefix="azure2_plot_selftest_") as tmp:
        work = Path(tmp)
        shutil.copytree(T1 / "data", work / "data")
        shutil.copy(T1 / "T1.azr", work / "T1.azr")
        for d in ("output", "checks"):
            (work / d).mkdir()
        try:
            p = pyazr.Project.open(work / "T1.azr")     # norms other than 1, so that the test sees how data are scaled
            p.set("el_lab090", norm="1.2 fixed")
            p.set("cap_g0", norm="0.9 fit prior 1 10%")
            p.save(work / "T1.azr")
            with pyazr.Session(str(work / "T1.azr")) as s:
                names = [d.name for d in s.evaluate(s.parameters.x0).data]
            r = run("plot_fit.py", "T1.azr", "--outdir", "fig", cwd=work)
            check("plot_fit: drawn pulls equal the engine's residuals for every data set", r.returncode == 0,
                  f"{len(names)} data sets")
            made = sorted(p.stem for p in (work / "fig").glob("*.png"))
            check("plot_fit: one figure per data set", made == sorted(names), f"{len(made)} figures")
            r = run("plot_fit.py", "T1.azr", "--outdir", "figS", "--sfactor", cwd=work)
            check("plot_fit --sfactor: pulls unchanged", r.returncode == 0)
        except Exception as err:
            check("plot_fit", False, f"{type(err).__name__}: {err}")
    print(f"{'FAIL' if FAILED else 'PASS'} azure2-plot selftest ({len(FAILED)} failed, {time.time() - t:.1f} s)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
