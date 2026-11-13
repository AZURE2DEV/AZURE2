#!/usr/bin/env python3
"""Self-test of the azure2-thm scripts: thm_inspect.py.

On a copy of the synthetic project T2 (tests/upgrade/T2: p + 15N -> 16O with a THM experiment of two HOES data sets
that share one profiled norm) in a temporary directory, thm_inspect.py reports the experiment LC09 with the norm n*
that every one of its data sets carries in the evaluation (``Evaluation.data[i].norm``, to 1e-10) and with the
experiment's chi2 equal to the sum of its data sets' chi2 (to 1e-8), and the plane-wave vertex reaches every data
energy.

Needs pyazr with its engine (PYAZR_ENGINE_DIR) and AZURE2_COULOMB_TABLE. PASS/FAIL lines; exit 0, 1 on a failure,
77 (skip) without the engine.
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
T2 = REPO / "tests" / "upgrade" / "T2"
FAILED = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {name}" + (f": {detail}" if detail else ""))
    if not ok:
        FAILED.append(name)


def main():
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    try:
        import pyazr
        pyazr.describe()
    except Exception as err:
        print(f"skip: the pyazr engine is not available ({type(err).__name__}: {err})")
        return 77
    sys.path.insert(0, str(HERE))
    import thm_inspect
    t = time.time()
    with tempfile.TemporaryDirectory(prefix="azure2_thm_selftest_") as tmp:
        work = Path(tmp)
        work = work / "T2"
        shutil.copytree(T2, work)
        for d in ("output", "checks"):
            (work / d).mkdir(exist_ok=True)
        try:
            t0 = time.time()
            r = subprocess.run([sys.executable, str(HERE / "thm_inspect.py"), "T2.azr"], cwd=work,
                               capture_output=True, text=True)
            print(f"  ran thm_inspect.py T2.azr in {time.time() - t0:.1f} s (exit {r.returncode})")
            print(r.stdout.rstrip())
            check("thm_inspect: runs on T2 without a flag", r.returncode == 0, r.stderr.strip()[-300:])
            with pyazr.Session(str(work / "T2.azr")) as s:
                rep = thm_inspect.inspect(s, s.parameters.x0)["LC09"]
            check("thm_inspect: every data set of LC09 carries n*", len(rep["datasets"]) >= 2 and
                  all(abs(n / rep["norm"] - 1) < 1e-10 for n in rep["dataset_norms"]),
                  f"n* {rep['norm']:.8g}, data sets {rep['dataset_norms']}")
            check("thm_inspect: the experiment's chi2 is its data sets' sum", abs(rep["chi2"] - rep["dataset_chi2"])
                  < 1e-8 * max(1.0, rep["chi2"]), f"{rep['chi2']:.8f} vs {rep['dataset_chi2']:.8f}")
            check("thm_inspect: the vertex reaches every data energy", rep.get("vertex_unreached") == 0,
                  f"{rep.get('vertex_model')} {rep.get('vertex_window')}")
        except Exception as err:
            check("thm_inspect", False, f"{type(err).__name__}: {err}")
    print(f"{'FAIL' if FAILED else 'PASS'} azure2-thm selftest ({len(FAILED)} failed, {time.time() - t:.1f} s)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
