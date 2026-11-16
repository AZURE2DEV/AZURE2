#!/usr/bin/env python3
"""The recipe "write a project from Python" (docs/source/user_guide/project_from_python.rst, the azure2-setup skill's
references/project-from-python.md) is pyazr/examples/project_from_python.py; this runs it.

  1. the script runs from schema to calculation in a temporary directory, the refusal it shows is the engine's
     (E_FIELD_REQUIRED, charge_state);
  2. its chi^2 equals, bit for bit, the chi^2 of the committed project it rebuilds
     (tests/units_scales/models/yield_story.azr) in a session;
  3. the docs page includes the script's recipe block and the skill names the script.

Needs the compiled module.  Run from anywhere:  python3 tests/pyazr/project_recipe_test.py
"""
import contextlib
import io
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "pyazr", "examples"))
os.environ.setdefault("OMP_NUM_THREADS", "1")

import pyazr  # noqa: E402
import project_from_python as recipe  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


work = tempfile.mkdtemp(prefix="project_recipe.")
try:
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        p, chi2 = recipe.main(os.path.join(work, "recipe"))
    text = out.getvalue()
    print(text)
    check("the recipe runs; its refusal is the engine's (E_FIELD_REQUIRED of charge_state)",
          "refused: ['E_FIELD_REQUIRED'] ['charge_state']" in text and "dataset Y_curve" in text, text)
    committed = os.path.join(work, "committed")
    os.makedirs(os.path.join(committed, "output"))
    shutil.copytree(os.path.join(recipe.STORY, "data"), os.path.join(committed, "data"))
    shutil.copy(os.path.join(recipe.STORY, "yield_story.azr"), committed)
    with pyazr.Session(os.path.join(committed, "yield_story.azr")) as m:
        want = m.evaluate().chi2
    check(f"chi^2 of the recipe's project = the committed project's, bit for bit ({chi2!r})", chi2 == want, (chi2, want))
finally:
    shutil.rmtree(work, ignore_errors=True)

page = open(os.path.join(ROOT, "docs", "source", "user_guide", "project_from_python.rst"), encoding="utf-8").read()
check("the docs page includes the script's recipe block",
      "pyazr/examples/project_from_python.py" in page and ":start-after: # [recipe]" in page)
skill = open(os.path.join(ROOT, ".claude", "skills", "azure2-setup", "references", "project-from-python.md"), encoding="utf-8").read()
check("the skill's reference names the script", "pyazr/examples/project_from_python.py" in skill)
print(f"{len(failures)} failure(s)" if failures else "PASS")
sys.exit(1 if failures else 0)
