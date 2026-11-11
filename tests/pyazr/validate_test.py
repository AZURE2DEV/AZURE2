#!/usr/bin/env python3
"""AZURE2 --validate and pyazr.validate(path); pyazr.convert_azr(..., to=3) (docs/upgrade/project_format3.md 6.2, 8).

  1. every fixture of experiment/tests/data/ validates through the command line (exit 0, no diagnostic) and through
     pyazr.validate (an empty list);
  2. a file with several errors: every one is reported in one run, the same codes by pyazr and by --validate --json
     (parsed here), the text form 'file:line: record.field: ...', exit 2;
  3. a format-2 file is refused with the conversion command: exit 3, pyazr raises ValueError;
  4. pyazr.convert_azr(T2, to=3) writes the committed conversion convert/tests/fixtures/T2.f3.azr byte for byte and its
     prior-covariance copy, returns the report (ok, the parameter names), and the result validates.
The command-line parts need AZURE2_BIN, the pyazr parts the compiled engine module; each part without what it needs
is skipped (exit 77 when nothing ran).  Run from anywhere:  python3 tests/pyazr/validate_test.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "experiment", "tests", "data")
FIXTURES = os.path.join(ROOT, "convert", "tests", "fixtures")
BIN = os.environ.get("AZURE2_BIN")
ORIGINALS = os.path.join(ROOT, "tests", "convert", "fixtures", "tests")  # the originals of converted projects

failures = []
ran = 0


def check(name, ok, detail=""):
    global ran
    ran += 1
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    import pyazr
    from pyazr import _azure2  # noqa: F401
except ImportError as e:
    print(f"pyazr parts skipped: {e}")
    pyazr = None


def cli(*args):
    p = subprocess.run([BIN, *args], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"))
    return p.returncode, p.stdout, p.stderr


BROKEN = """<config>
format=3
</config>
<data>
dataset d1
    reaction       15N(p,a)12C
    quantity       sigma
    unit           mb/sr
    file           d1.dat
    columns        E_lab[MeV] value error
    norm           1 fixed
    taget          t1
dataset d2
    reaction       15N(p,a)12C
    quantity       sigma
    unit           b
    file           d2.dat
    columns        E_lab[MeV] value error
    norm           1 fixed
    target         nowhere
</data>
"""
BROKEN_CODES = {"E_UNIT_NOT_ALLOWED", "E_UNKNOWN_FIELD", "E_UNKNOWN_REFERENCE"}

fixtures = sorted(f for f in os.listdir(DATA) if f.endswith(".azr"))
with tempfile.TemporaryDirectory() as tmp:
    broken = os.path.join(tmp, "broken.azr")
    with open(broken, "w") as f:
        f.write(BROKEN)
    # the old T2 (its original, tests/convert/fixtures) in a copy of its folder, beside its data and covariance files
    old = os.path.join(shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T2"), os.path.join(tmp, "T2old")),
                       "T2.azr")
    shutil.copy(os.path.join(ORIGINALS, "upgrade", "T2", "T2.azr"), old)

    # 1. the fixtures
    for name in fixtures:
        path = os.path.join(DATA, name)
        if BIN:
            code, out, err = cli("--validate", path, "--json")
            check(f"--validate --json {name}: exit 0, []", code == 0 and json.loads(out) == [], out[:300])
            code, out, _ = cli("--validate", path)
            check(f"--validate {name}: exit 0, valid", code == 0 and "0 diagnostic(s), valid" in out, out[:300])
        if pyazr:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                check(f"pyazr.validate {name}: []", pyazr.validate(path) == [])

    # 2. several errors, all of them, in one run
    if BIN:
        code, out, err = cli("--validate", broken, "--json")
        diags = json.loads(out)
        check("--validate --json broken: exit 2", code == 2, f"exit {code}")
        check("--validate --json broken: every error", {d["code"] for d in diags} == BROKEN_CODES, str([d["code"] for d in diags]))
        check("--validate --json broken: file, line, record and field named",
              all(d["file"] == broken and d["line"] > 0 and d["record"] in ("d1", "d2") for d in diags))
        code, out, _ = cli("--validate", broken)
        lines = [l for l in out.splitlines() if l.startswith(broken + ":") and l[len(broken) + 1:].split(":")[0].isdigit()]
        check("--validate broken (text): one line per error, file:line: record.field: ...",
              code == 2 and len(lines) == len(BROKEN_CODES) and all(": d1." in l or ": d2." in l for l in lines), out)
        check("--validate: the note on the checks skipped", "check(s) skipped" in out)
    if pyazr:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            diags = pyazr.validate(broken)
        check("pyazr.validate broken: every error", {d["code"] for d in diags} == BROKEN_CODES, str(diags))
        if BIN:
            check("pyazr.validate and --validate --json: the same diagnostics",
                  diags == json.loads(cli("--validate", broken, "--json")[1]))

    # 3. a format-2 file
    if BIN:
        code, out, _ = cli("--validate", old)
        check("--validate a format-1/2 file: exit 3 with the conversion command",
              code == 3 and "AZURE2 --convert" in out and "<segmentsData> at line" in out, out)
    if pyazr:
        try:
            pyazr.validate(old)
            check("pyazr.validate a format-1/2 file raises", False)
        except ValueError as e:
            check("pyazr.validate a format-1/2 file raises", "convert" in str(e))

    # 4. pyazr.convert_azr(..., to=3)
    if pyazr:
        new = os.path.join(tmp, "T2.f3.azr")
        rep = pyazr.convert_azr(old, new, to=3, report=os.path.join(tmp, "T2.report.json"))
        with open(new, "rb") as f, open(os.path.join(FIXTURES, "T2.f3.azr"), "rb") as g:
            check("convert_azr(to=3): the committed conversion, byte for byte", f.read() == g.read())
        with open(os.path.join(tmp, "prior_cov.f3.txt"), "rb") as f, open(os.path.join(FIXTURES, "T2.prior_cov.f3.txt"), "rb") as g:
            check("convert_azr(to=3): the prior-covariance copy", f.read() == g.read())
        check("convert_azr(to=3): the report", rep["ok"] and rep["parameters"].get("segment_7_norm") == "s07_profile.norm",
              str(rep.get("parameters")))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # data files are read relative to the working directory, as the engine opens them: the project's folder
            here = os.getcwd()
            os.chdir(os.path.dirname(old))
            try:
                check("convert_azr(to=3): the result validates and builds", pyazr.validate(new) == [])
                # 5. pyazr.describe: what the engine makes of it, as AZURE2 --describe
                text = pyazr.describe(new)
                check("pyazr.describe: the pair keys and the resolved data file",
                      "entrance pair 1 (p+15N), exit pair 2 (a+12C)" in text and
                      os.path.join(os.path.dirname(old), "data", "s01_int_thin.dat") in text, text[:400])
                d = pyazr.describe(new, as_json=True)
                check("pyazr.describe(as_json=True): the data sets with their parameters",
                      d["datasets"][0]["name"] == "s01_int_thin" and
                      any(p["name"] == "s01_int_thin.norm" for p in d["datasets"][0]["parameters"]), str(d)[:300])
                if BIN:
                    code, out, _ = cli("--describe", new)
                    check("--describe and pyazr.describe: the same text", code == 0 and out == text, out[:300])
            finally:
                os.chdir(here)
        try:
            pyazr.describe(broken)
            check("pyazr.describe of a project with errors raises", False)
        except ValueError as e:
            check("pyazr.describe of a project with errors raises", "errors" in str(e))
        try:
            pyazr.convert_azr(new, os.path.join(tmp, "again.azr"), to=3)
            check("convert_azr(to=3) of a format-3 file is refused", False)
        except ValueError as e:
            check("convert_azr(to=3) of a format-3 file is refused", "already format 3" in str(e))

if ran == 0:
    print("SKIP: neither AZURE2_BIN nor the engine module")
    sys.exit(77)
print(f"{ran} checks, {len(failures)} failures")
sys.exit(1 if failures else 0)
