#!/usr/bin/env python3
"""pyazr.convert_azr and pyazr.convert_params: an old-format project (format 1) in project format 3
(project_format3.md section 8) and a parameter file of it renamed (section 8.4), and pyazr's one refusal of an old file.

On the originals of tests/17O (with <lastRun>, an <mcmc> block of run settings and a <potential> block without a
hybrid potential) and tests/hybrid_potential (a hybrid potential for one of two pairs), kept for the converter in
tests/convert/fixtures:
  1. the converter writes format 3 only: to=2 is refused (ValueError); an existing target is not overwritten;
  2. a session on the conversion gives the chi-squared of a session on the tree's project (the committed
     conversion), to 1e-12;
  3. pyazr reads format 3 only: a session (azure2), Project.open and the functions that take a path (tabulate,
     angular_distribution, model_average's writer, merge_projects) refuse the original, the same
     file with format=2, the format-3 project with an old <segmentsData> appended, an empty file and a file that is
     not a project, each with one FormatError: the engine's message (io::NotFormat3Message), the conversion named
     where there is one;
  4. convert_params renames a parameter file naming every old parameter of the converter's map (the report's
     "parameters") to the conversion's names, every other byte kept; it refuses a format-3 project as the map's
     source.
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
    from pyazr import Session, convert_azr, FormatError, Project
    from pyazr import _azure2
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

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
        try:
            convert_azr(old, new, to=2)
            check(f"{project}: to=2 is refused (format 3 only)", False, "no error")
        except ValueError as e:
            check(f"{project}: to=2 is refused (format 3 only)", "format 3 only" in str(e), str(e))
        with open(new, "wb") as f:
            f.write(b"exists")
        try:
            convert_azr(old, new)
            check(f"{project}: an existing target is not overwritten", False)
        except FileExistsError:
            check(f"{project}: an existing target is not overwritten", open(new, "rb").read() == b"exists")
        # the original with format=2 in <config>, the format-3 project with an old section, an empty file, not a project
        text = open(old, "rb").read()
        f2 = os.path.join(d, "format2.azr")
        with open(f2, "wb") as f:
            f.write(text.replace(b"</config>", b"format=2\n</config>", 1))
        check(f"{project}: format=2 written", _azure2.project_format(open(f2, "rb").read()) == 2)
        f3old = os.path.join(d, "f3old.azr")
        with open(f3old, "wb") as f:
            f.write(open(os.path.join(d, azr), "rb").read() + b"<segmentsData>\n</segmentsData>\n")
        empty = os.path.join(d, "empty.azr")
        open(empty, "wb").close()
        notproject = os.path.join(d, "text.azr")
        with open(notproject, "w") as f:
            f.write("Not a project.\n")

        f3 = os.path.join(d, "converted3.azr")
        rep = convert_azr(old, f3, report=os.path.join(d, "report.json"))
        check(f"{project}: converted to format 3", rep["ok"] and _azure2.project_format(open(f3, "rb").read()) == 3)
        a = Session(os.path.join(d, azr), cwd=d)
        b = Session(f3, cwd=d)
        x = np.asarray(a.parameters.x0, float)
        ca, cb = a.evaluate(x).chi2, b.evaluate(x).chi2
        check(f"{project}: the conversion's chi2 = the tree's project's ({ca:.10g})",
              abs(ca - cb) <= 1e-12 * max(1.0, abs(ca)), f"{ca} vs {cb}")

        expected = {"format 1": f"is a format-1 project", "format 2": f"is a format-2 project",
                    "format 3 over an old section": "says format=3 but holds <segmentsData>",
                    "empty": "is not an AZURE2 project: it has no <config> section (the file is empty)",
                    "not a project": "is not an AZURE2 project: it has no <config> section."}
        for label, path in (("format 1", old), ("format 2", f2), ("format 3 over an old section", f3old),
                            ("empty", empty), ("not a project", notproject)):
            message = _azure2.not_format3_message(open(path, "rb").read(), path)
            calls = {"azure2": lambda: Session(path, cwd=d), "Project.open": lambda: Project.open(path),
                     "tabulate": lambda: pyazr.tabulate(path, [(1, 1)], 1.0, 2.0),
                     "angular_distribution": lambda: pyazr.angular_distribution(path, [1.0]),
                     "write_averaged_azr": lambda: pyazr.write_averaged_azr(path, os.path.join(d, "avg.azr"), {}),
                     "merge_projects": lambda: pyazr.merge_projects([path], os.path.join(d, "merged.azr"))}
            for name, call in calls.items():
                try:
                    call()
                    check(f"{project}: {name} refuses the {label} file", False, "no error")
                except FormatError as e:
                    check(f"{project}: {name} refuses the {label} file with the engine's message",
                          str(e) == message and expected[label] in str(e)
                          and ("AZURE2 --convert" in str(e)) == label.startswith("format"), str(e)[:300])

        # 4. a parameter file with every old name of the conversion's map (the report's "parameters"), renamed by
        # convert_params; the level parameter energy_1 keeps its name
        names = rep["parameters"]
        lines = [f"{n:>20}  1.0000000e+00  1.0000000e-02\n" for n in names]
        sav_old, sav_new = os.path.join(d, "old.sav"), os.path.join(d, "new.sav")
        with open(sav_old, "w") as f:
            f.write("    #parametrization  1.0000000e+00  0.0000000e+00\n" + "".join(lines) + "         energy_1  1.0  0.1\n")
        renamed = pyazr.convert_params(sav_old, sav_new, old)
        out = open(sav_new).read().splitlines()
        check(f"{project}: convert_params renames every name of the map ({len(names)})",
              len(names) > 0 and renamed == list(names.items()) and [l.split()[0] for l in out[1:-1]] == list(names.values()), f"{renamed}")
        check(f"{project}: ... and keeps every other byte",
              out[0] == "    #parametrization  1.0000000e+00  0.0000000e+00" and out[-1] == "         energy_1  1.0  0.1"
              and all(l.split()[1:] == ["1.0000000e+00", "1.0000000e-02"] for l in out[1:-1]))
        try:
            pyazr.convert_params(sav_old, os.path.join(d, "new2.sav"), f3)
            check(f"{project}: convert_params refuses a format-3 project as the map's source", False, "no error")
        except ValueError as e:
            check(f"{project}: convert_params refuses a format-3 project as the map's source", "format 3" in str(e), str(e))
    finally:
        shutil.rmtree(work, ignore_errors=True)

print("PASS" if not failures else f"FAIL: {len(failures)}")
sys.exit(1 if failures else 0)
