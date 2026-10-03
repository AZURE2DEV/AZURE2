#!/usr/bin/env python3
"""write_output_files reports the directory the engine wrote to.

A project names its output directory in ``<config>`` (second line, relative
to the run's working directory unless absolute).  write_output_files used to
return ``cwd/output`` whatever that line said, so a project writing to
``results/`` got a path to a directory that was empty or did not exist.
Checked on tests/hybrid_potential (three elastic segments, seconds):

  1. ``output_dir`` parses the line as the engine does (comment cut, trimmed;
     relative to cwd, or absolute).
  2. With a relative ``results/``: the files are there, nothing in
     ``output/``, and the returned path is ``cwd/results``.
  3. With an absolute directory elsewhere: likewise.
  4. The usual ``output/``: unchanged, ``cwd/output``.

Needs the compiled engine for 2-4; skips them cleanly.

Run from anywhere:  python3 tests/pyazr/output_dir_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "hybrid_potential")
AZR = "hybrid_potential.azr"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
try:
    from pyazr import azure2
    from pyazr.azrfile import AzrModel
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)


def fresh_copy(dst):
    os.makedirs(dst)
    shutil.copy(os.path.join(SOURCE, AZR), dst)
    shutil.copytree(os.path.join(SOURCE, "data"), os.path.join(dst, "data"))
    os.makedirs(os.path.join(dst, "checks"))
    os.makedirs(os.path.join(dst, "output"))


def outfiles(d):
    return sorted(f for f in os.listdir(d) if f.startswith("AZUREOut_") or f == "chiSquared.out") \
        if os.path.isdir(d) else []


with tempfile.TemporaryDirectory() as tmp:
    # 2. relative results/ --------------------------------------------------
    print("1-2. relative output directory")
    work = os.path.join(tmp, "rel")
    fresh_copy(work)
    os.makedirs(os.path.join(work, "results"))
    azr = os.path.join(work, AZR)
    text = open(azr).read().replace("output/", "   results/   ", 1)
    open(azr, "w").write(text)
    with azure2(azr, cwd=work) as m:
        check("output_dir is cwd/results", m.output_dir == os.path.join(work, "results"),
              m.output_dir)
        got = m.write_output_files()
    check("write_output_files returns it", got == os.path.join(work, "results"), got)
    check("the files are there", "chiSquared.out" in outfiles(got), outfiles(got))
    check("and not in output/", outfiles(os.path.join(work, "output")) == [],
          outfiles(os.path.join(work, "output")))

    # 3. absolute -------------------------------------------------------------
    print("3. absolute output directory")
    work = os.path.join(tmp, "abs")
    fresh_copy(work)
    elsewhere = os.path.join(tmp, "elsewhere")
    m0 = AzrModel.from_file(os.path.join(work, AZR))
    m0.set_output_dir(elsewhere)
    m0.write(os.path.join(work, AZR))
    with azure2(os.path.join(work, AZR), cwd=work) as m:
        got = m.write_output_files()
    check("write_output_files returns the absolute directory",
          os.path.normpath(got) == os.path.normpath(elsewhere), got)
    check("the files are there", "chiSquared.out" in outfiles(elsewhere), outfiles(elsewhere))
    check("and not in output/", outfiles(os.path.join(work, "output")) == [])

    # 4. the usual output/ ----------------------------------------------------
    print("4. output/")
    work = os.path.join(tmp, "usual")
    fresh_copy(work)
    with azure2(os.path.join(work, AZR), cwd=work) as m:
        got = m.write_output_files()
    check("cwd/output, as before", got == os.path.join(work, "output"), got)
    check("the files are there", "chiSquared.out" in outfiles(got))

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all output-directory checks passed")
