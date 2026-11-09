#!/usr/bin/env python3
"""pyazr.convert_azr: an old-format project (format 1) in the updated format (docs/upgrade/project.md section 4).

On tests/17O (with <lastRun>, an <mcmc> block of run settings and a <potential> block without a hybrid
potential) and tests/hybrid_potential (a hybrid potential for one of two pairs):
  1. the converted file is format 2 (_azure2.project_format) and lists its changes; converting it again changes
     nothing; an existing target is not overwritten;
  2. <lastRun> and <mcmc>'s run settings are gone, and <potential> remains only where a pair uses the hybrid
     potential;
  3. a session on the converted file gives the chi-squared of a session on the original, to 1e-12.
Needs the compiled engine; skips cleanly without it.  Run from anywhere:  python3 tests/pyazr/convert_azr_test.py
"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
CASES = [("17O", "17O.azr", False), ("hybrid_potential", "hybrid_potential.azr", True)]

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    from pyazr import azure2, convert_azr
    from pyazr import _azure2
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(0)

for project, azr, hybrid in CASES:
    work = tempfile.mkdtemp(prefix="convert_azr.")
    try:
        d = os.path.join(work, project)
        shutil.copytree(os.path.join(ROOT, "tests", project), d)
        for sub in ("output", "checks"):
            os.makedirs(os.path.join(d, sub), exist_ok=True)
        old = os.path.join(d, azr)
        new = os.path.join(d, "converted.azr")
        changes = convert_azr(old, new)
        text = open(new, "rb").read()
        check(f"{project}: format 2, {len(changes)} change(s) listed",
              _azure2.project_format(text) == 2 and len(changes) > 0, f"{changes}")
        check(f"{project}: converting again changes nothing", convert_azr(new, new + ".2") == [])
        try:
            convert_azr(old, new)
            check(f"{project}: an existing target is not overwritten", False)
        except FileExistsError:
            check(f"{project}: an existing target is not overwritten", True)
        check(f"{project}: <lastRun> gone", b"<lastRun>" not in text)
        check(f"{project}: <mcmc>'s run settings gone", b"nwalkers" not in text)
        check(f"{project}: <potential> {'kept (hybrid)' if hybrid else 'gone (no hybrid)'}",
              (b"<potential>" in text) == hybrid)
        a = azure2(old, cwd=d)
        b = azure2(new, cwd=d)
        x = np.asarray(a.params_rwa, float)
        ca, cb = float(np.sum(a.calculate_chi2_rwa(x))), float(np.sum(b.calculate_chi2_rwa(x)))
        check(f"{project}: the converted file's chi2 = the original's ({ca:.10g})",
              abs(ca - cb) <= 1e-12 * max(1.0, abs(ca)), f"{ca} vs {cb}")
    finally:
        shutil.rmtree(work, ignore_errors=True)

print("PASS" if not failures else f"FAIL: {len(failures)}")
sys.exit(1 if failures else 0)
