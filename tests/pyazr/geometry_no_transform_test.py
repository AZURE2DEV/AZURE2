#!/usr/bin/env python3
"""The energies and angles of an evaluation come with no observed-to-formal transformation.

Point energies and angles do not depend on the R-matrix parameters, but
the methods that returned them once took the physical path, so a caller handing them ``parameters.x0`` (an easy slip: every
other call takes it) had each amplitude read as a width or ANC.  With a
sub-threshold ANC channel that prints "Denominator less than zero while
transforming ... input 0.118102 fm^(-1/2)" (the 17O cross-check), which the
CLI never prints because it never transforms an amplitude vector.  Now, on
examples/o17_guardo2017_fit (1- level 7.6 keV below the n + 17O threshold,
ANC channel; an extrapolation record added here), in data and extrapolation
mode:

  1. control: the physical path fed parameters.x0 still prints the warning
     (so the check below can see it);
  2. evaluate(parameters.x0) and evaluate() print nothing, and the energies,
     excitation energies and angles of their data sets agree exactly (a
     physical vector is no input of the session's calls).

The engine prints on the process's stdout, so the session runs in a child
process.  Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/geometry_no_transform_test.py

This test reads the binding itself (`_azure2.Session`), so it uses the private handle `Session._engine`.
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "examples", "o17_guardo2017_fit")
AZR = "o17_guardo2017_fit.azr"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


sys.path.insert(0, ROOT)
try:
    import pyazr  # noqa: F401
    from pyazr import Session  # noqa: F401
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)

CHILD = r'''
import sys
sys.path.insert(0, sys.argv[1])
import numpy as np
from pyazr import Session
s = Session(sys.argv[2])
if sys.argv[4] == "extrap":
    s.extrap_mode()
if sys.argv[3] == "control":
    s._engine.update_segments(np.asarray(s.parameters.x0, float))
else:
    same = True
    a, c = s.evaluate(s.parameters.x0).data, s.evaluate().data
    for field in ("energies", "excitation_energies", "angles"):
        same &= all(np.array_equal(getattr(x, field), getattr(z, field)) for x, z in zip(a, c))
        same &= sum(len(getattr(x, field)) for x in a) > 0
    print("SAME", same)
'''

# One child per check: the engine writes through C stdio, buffered
# separately from Python's own stdout when piped, so the two cannot be
# interleaved reliably within one process.
with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "w")
    shutil.copytree(SOURCE, work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
        os.makedirs(os.path.join(work, junk))
    from pyazr import Project
    m = Project.open(os.path.join(work, AZR))
    m.add("extrapolation", "ex", reaction=m.reaction(1, 2), quantity="sigma", unit="b",
          energies="range 0.01 0.3 MeV lab step 0.01")
    m.save(os.path.join(work, "ex.azr"))
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "2"))

    def child(what, mode):
        p = subprocess.run([sys.executable, "-c", CHILD, ROOT, os.path.join(work, "ex.azr"),
                            what, mode], capture_output=True, text=True, timeout=600, env=env)
        if p.returncode != 0:
            return None
        return p.stdout

    for mode in ("data", "extrap"):
        print(f"{mode} mode")
        ctl = child("control", mode)
        geo = child("geometry", mode)
        check("control: the physical path fed parameters.x0 warns",
              ctl is not None and "Denominator less than zero" in ctl, (ctl or "")[-300:])
        check("energies, excitation energies, angles: no warning",
              geo is not None and "Denominator" not in geo, (geo or "")[-300:])
        check("rwa, physical and None agree exactly", geo is not None and "SAME True" in geo,
              (geo or "")[-300:])

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all geometry checks passed")
