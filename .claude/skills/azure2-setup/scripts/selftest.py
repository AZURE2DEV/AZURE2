#!/usr/bin/env python3
"""Self-test of the azure2-setup scripts: audit_project.py.

On copies of the synthetic project T1 (tests/upgrade/T1: p + 12C -> 13N, AME2020 masses) in a temporary directory:

- the project as shipped raises no flag (its S_p agrees with the masses to well below 1 keV, every data set keeps
  its rows);
- a copy whose S_p was raised by 10 keV (``Project.set_pair(1, separation_energy=...)``) is flagged, with the
  difference reported as +10.00 keV;
- a copy whose proton mass is the mass number 1 is flagged as a mass number.

Needs pyazr with its engine (PYAZR_ENGINE_DIR). PASS/FAIL lines; exit 0, 1 on a failure, 77 (skip) without the
engine. About 5 s.
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


def run(script, *args, cwd):
    t = time.time()
    r = subprocess.run([sys.executable, str(HERE / script), *map(str, args)], cwd=cwd, capture_output=True, text=True)
    print(f"  ran {script} {' '.join(map(str, args))} in {time.time() - t:.1f} s (exit {r.returncode})")
    if r.returncode not in (0, 1):
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
    with tempfile.TemporaryDirectory(prefix="azure2_setup_selftest_") as tmp:
        work = Path(tmp)
        shutil.copytree(T1 / "data", work / "data")
        shutil.copy(T1 / "T1.azr", work / "T1.azr")
        for d in ("output", "checks"):
            (work / d).mkdir()
        try:
            r = run("audit_project.py", "T1.azr", cwd=work)
            check("audit: T1 as shipped raises no flag", r.returncode == 0, r.stdout.strip().splitlines()[-1])

            p = pyazr.Project.open(work / "T1.azr")
            s_p = next(q for q in p.pairs() if q.key == 1).separation_energy
            p.set_pair(1, separation_energy=s_p + 0.010)
            p.save(work / "T1_s.azr")
            r = run("audit_project.py", "T1_s.azr", cwd=work)
            m = re.search(r"pair 1: S differs .* by ([+-]\d+\.\d+) keV", r.stdout)
            check("audit: S raised by 10 keV is flagged as +10.00 keV", r.returncode == 1 and m is not None
                  and abs(float(m.group(1)) - 10.0) < 0.01, m.group(0) if m else r.stdout[-300:])

            p = pyazr.Project.open(work / "T1.azr")
            p.set_pair(1, light_mass=1.0)
            p.save(work / "T1_m.azr")
            r = run("audit_project.py", "T1_m.azr", cwd=work)
            check("audit: a proton mass of 1 u is flagged as a mass number", r.returncode == 1
                  and "mass number" in r.stdout)
        except Exception as err:
            check("audit", False, f"{type(err).__name__}: {err}")
    print(f"{'FAIL' if FAILED else 'PASS'} azure2-setup selftest ({len(FAILED)} failed, {time.time() - t:.1f} s)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
