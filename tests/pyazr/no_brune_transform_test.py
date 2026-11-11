#!/usr/bin/env python3
"""pyazr, standard formalism (use_brune=False) on a project with no standard-R-matrix equivalent: a clear error.

tests/7Li_p_ay is a Brune fit whose 0+ background poles have partial widths up to 135 MeV; the transformation of its
observed parameters to formal ones (CNuc::TransformIn, N a = E M a with M = 1 - sum gamma gamma' dS/dE) has no real
solution because M is not positive definite (Brune, PRC 66, 044611).  Before 2026-11-01 GSL's failure went unnoticed
and the session died at load with "RuntimeError: z is not finite in log_Gamma." (a crash, with a core dump, in the
CLI).  Now opening the session raises AZURE2Error after the engine printed which J^pi group has no standard
equivalent; with Brune (the default) the project opens and calculates as before.  The CLI side:
tests/no_brune_transform/check.sh.

Run from anywhere:  python3 tests/pyazr/no_brune_transform_test.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

WORKER = r"""
import sys
sys.path.insert(0, sys.argv[1])
import numpy as np
from pyazr import azure2
try:
    with azure2(sys.argv[2], cwd=sys.argv[3], use_brune=sys.argv[4] == "1") as a:
        v = a.calculate_rwa(np.array(a.params_rwa, float))
        print("RESULT opened", int(all(np.all(np.isfinite(np.asarray(x, float))) for x in v)))
except Exception as err:
    print("RESULT", type(err).__name__, str(err).replace(chr(10), " "))
"""


def main():
    try:
        from pyazr import azure2  # noqa: F401
    except Exception as err:  # same skip as the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 0
    failures = []

    def check(name, ok, detail=""):
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
        if not ok:
            failures.append(name)

    work = tempfile.mkdtemp(prefix="pyazr_no_brune_")
    try:
        src = os.path.join(ROOT, "tests", "7Li_p_ay")
        shutil.copytree(os.path.join(src, "data"), os.path.join(work, "data"))
        shutil.copy(os.path.join(src, "7Li_p_ay.azr"), work)
        os.makedirs(os.path.join(work, "output"))
        os.makedirs(os.path.join(work, "checks"))
        azr = os.path.join(work, "7Li_p_ay.azr")
        for brune in ("0", "1"):
            p = subprocess.run([sys.executable, "-c", WORKER, ROOT, azr, work, brune], capture_output=True, text=True,
                               timeout=600, env=dict(os.environ, OMP_NUM_THREADS="2"))
            out = p.stdout + p.stderr
            result = [l for l in p.stdout.splitlines() if l.startswith("RESULT")]
            result = result[-1] if result else f"no result (exit {p.returncode})"
            if brune == "0":
                check("use_brune=False: the process survives", p.returncode == 0, f"exit {p.returncode}")
                check("use_brune=False: opening the session raises AZURE2Error", result.startswith("RESULT AZURE2Error"),
                      result)
                check("use_brune=False: the engine says which J^pi has no standard equivalent",
                      "No standard R-matrix (constant boundary conditions) corresponds to the levels of J^pi = 0+"
                      in out, "message missing")
            else:
                check("use_brune=True: opens and calculates (finite)", result == "RESULT opened 1", result)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
