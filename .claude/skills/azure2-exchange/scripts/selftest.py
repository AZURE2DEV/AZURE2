#!/usr/bin/env python3
"""Self-test of the azure2-exchange scripts: basis_compare.py.

On a copy of the synthetic project T1 (tests/upgrade/T1: resonances, sub-threshold states with ANCs, capture with
external capture) the Brune, Park and standard bases give the same model for every data set (to 1e-6: one collision
matrix in three parametrizations), and the Park norms J of the levels lie in (0, 1].

Needs pyazr with its engine (PYAZR_ENGINE_DIR). PASS/FAIL lines; exit 0, 1 on a failure, 77 (skip) without the
engine. About 10 s.
"""
import argparse
import re
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


def main():
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    try:
        import pyazr
        pyazr.describe()
    except Exception as err:
        print(f"skip: the pyazr engine is not available ({type(err).__name__}: {err})")
        return 77
    t = time.time()
    with tempfile.TemporaryDirectory(prefix="azure2_exchange_selftest_") as tmp:
        work = Path(tmp)
        shutil.copytree(T1 / "data", work / "data")
        shutil.copy(T1 / "T1.azr", work / "T1.azr")
        for d in ("output", "checks"):
            (work / d).mkdir()
        t0 = time.time()
        r = subprocess.run([sys.executable, str(HERE / "basis_compare.py"), "T1.azr"], cwd=work, capture_output=True,
                           text=True)
        print(f"  ran basis_compare.py T1.azr in {time.time() - t0:.1f} s (exit {r.returncode})")
        print(r.stdout.rstrip())
        if r.returncode:
            print(r.stderr[-2000:])
        lines = {l.split()[0]: l for l in r.stdout.splitlines() if l and not l.startswith(" ")}
        for basis in ("park", "standard"):
            m = re.search(r"from Brune (\S+)", lines.get(basis, ""))
            check(f"basis_compare: {basis} gives Brune's model", m is not None and float(m.group(1)) < 1e-6,
                  m.group(0) if m else lines.get(basis, "missing"))
        m = re.search(r"Park norms J: (.*)", r.stdout)
        norms = [float(v) for v in m.group(1).split()] if m else []
        check("basis_compare: Park norms J in (0, 1]", norms and all(0 < j <= 1 + 1e-12 for j in norms),
              " ".join(f"{j:.4f}" for j in norms))
    print(f"{'FAIL' if FAILED else 'PASS'} azure2-exchange selftest ({len(FAILED)} failed, {time.time() - t:.1f} s)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
