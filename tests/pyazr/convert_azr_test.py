#!/usr/bin/env python3
"""pyazr.convert_azr: an old-format project (format 1) in the updated formats (docs/upgrade/project.md section 4,
project_format3.md section 8), and pyazr's one refusal of an old file.

On the originals of tests/17O (with <lastRun>, an <mcmc> block of run settings and a <potential> block without a
hybrid potential) and tests/hybrid_potential (a hybrid potential for one of two pairs), kept for the converter in
tests/convert/fixtures:
  1. to=2: the converted file is format 2 (_azure2.project_format) and lists its changes; converting it again
     changes nothing; an existing target is not overwritten; <lastRun> and <mcmc>'s run settings are gone, and
     <potential> remains only where a pair uses the hybrid potential;
  2. to=3 (the default): a session on the conversion gives the chi-squared of a session on the tree's project
     (the committed conversion), to 1e-12;
  3. pyazr reads format 3 only: a session (azure2), Project.open and the functions that take a path (tabulate,
     angular_distribution, model_average's writer, levels_from_azr, merge_projects) refuse the original and its
     format-2 conversion with one FormatError, the engine's message naming `AZURE2 --convert ... --to 3`.
Needs the compiled engine; skips cleanly without it.  Run from anywhere:  python3 tests/pyazr/convert_azr_test.py
"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
FIXTURES = os.path.join(ROOT, "tests", "convert", "fixtures", "tests")  # the originals of converted projects
CASES = [("17O", "17O.azr", False), ("hybrid_potential", "hybrid_potential.azr", True)]

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import pyazr
    from pyazr import azure2, convert_azr, FormatError, Project
    from pyazr import _azure2
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(0)

for project, azr, hybrid in CASES:
    work = tempfile.mkdtemp(prefix="convert_azr.")
    try:
        d = os.path.join(work, project)
        shutil.copytree(os.path.join(ROOT, "tests", project), d)
        old = os.path.join(d, "original.azr")
        shutil.copy(os.path.join(FIXTURES, project, azr), old)
        for sub in ("output", "checks"):
            os.makedirs(os.path.join(d, sub), exist_ok=True)
        new = os.path.join(d, "converted.azr")
        changes = convert_azr(old, new, to=2)
        text = open(new, "rb").read()
        check(f"{project}: format 2, {len(changes)} change(s) listed",
              _azure2.project_format(text) == 2 and len(changes) > 0, f"{changes}")
        check(f"{project}: converting again changes nothing", convert_azr(new, new + ".2", to=2) == [])
        try:
            convert_azr(old, new, to=2)
            check(f"{project}: an existing target is not overwritten", False)
        except FileExistsError:
            check(f"{project}: an existing target is not overwritten", True)
        check(f"{project}: <lastRun> gone", b"<lastRun>" not in text)
        check(f"{project}: <mcmc>'s run settings gone", b"nwalkers" not in text)
        check(f"{project}: <potential> {'kept (hybrid)' if hybrid else 'gone (no hybrid)'}",
              (b"<potential>" in text) == hybrid)

        f3 = os.path.join(d, "converted3.azr")
        rep = convert_azr(old, f3)
        check(f"{project}: to=3 is the default", rep["ok"] and _azure2.project_format(open(f3, "rb").read()) == 3)
        a = azure2(os.path.join(d, azr), cwd=d)
        b = azure2(f3, cwd=d)
        x = np.asarray(a.params_rwa, float)
        ca, cb = float(np.sum(a.calculate_chi2_rwa(x))), float(np.sum(b.calculate_chi2_rwa(x)))
        check(f"{project}: the conversion's chi2 = the tree's project's ({ca:.10g})",
              abs(ca - cb) <= 1e-12 * max(1.0, abs(ca)), f"{ca} vs {cb}")

        for label, path in (("format 1", old), ("format 2", new)):
            calls = {"azure2": lambda: azure2(path, cwd=d), "Project.open": lambda: Project.open(path),
                     "tabulate": lambda: pyazr.tabulate(path, [(1, 1)], 1.0, 2.0),
                     "angular_distribution": lambda: pyazr.angular_distribution(path, [1.0]),
                     "write_averaged_azr": lambda: pyazr.write_averaged_azr(path, os.path.join(d, "avg.azr"), {}),
                     "levels_from_azr": lambda: pyazr.levels_from_azr(path),
                     "merge_projects": lambda: pyazr.merge_projects([path], os.path.join(d, "merged.azr"))}
            for name, call in calls.items():
                try:
                    call()
                    check(f"{project}: {name} refuses the {label} file", False, "no error")
                except FormatError as e:
                    check(f"{project}: {name} refuses the {label} file, naming the conversion",
                          "AZURE2 --convert" in str(e) and "--to 3" in str(e), str(e)[:300])
    finally:
        shutil.rmtree(work, ignore_errors=True)

print("PASS" if not failures else f"FAIL: {len(failures)}")
sys.exit(1 if failures else 0)
