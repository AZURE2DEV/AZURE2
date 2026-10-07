#!/usr/bin/env python3
"""calculate_energies and friends run no observed-to-formal transformation.

Point energies and angles do not depend on the R-matrix parameters, but
calculate_energies / calculate_excitation_energy / calculate_angles took the
physical path, so a caller handing them ``params_rwa`` (an easy slip: every
other call takes it) had each amplitude read as a width or ANC.  With a
sub-threshold ANC channel that prints "Denominator less than zero while
transforming ... input 0.118102 fm^(-1/2)" (the 17O cross-check), which the
CLI never prints because it never transforms an amplitude vector.  Now, on
examples/o17_guardo2017_fit (1- level 7.6 keV below the n + 17O threshold,
ANC channel; an extrapolation segment added here), in data and extrapolation
mode:

  1. control: the physical path fed params_rwa still prints the warning
     (so the check below can see it);
  2. calculate_energies, calculate_excitation_energy and calculate_angles
     with params_rwa, with the physical params and with None print nothing
     and agree exactly.

The engine prints on the process's stdout, so the session runs in a child
process.  Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/geometry_no_transform_test.py
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
    from pyazr import azure2  # noqa: F401
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)

CHILD = r'''
import sys
sys.path.insert(0, sys.argv[1])
import numpy as np
from pyazr import azure2
s = azure2(sys.argv[2])
if sys.argv[4] == "extrap":
    s.extrap_mode()
if sys.argv[3] == "control":
    s.sess.update_segments(np.asarray(s.params_rwa, float))
else:
    same = True
    for f in (s.calculate_energies, s.calculate_excitation_energy, s.calculate_angles):
        a, b, c = f(s.params_rwa), f(s.params), f()
        same &= all(np.array_equal(x, y) and np.array_equal(x, z) for x, y, z in zip(a, b, c))
        same &= sum(len(x) for x in a) > 0
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
    from pyazr.azrfile import AzrModel
    m = AzrModel.from_file(os.path.join(work, AZR))
    m.add_extrapolation(1, 2, 0.01, 0.3, 0.01)
    m.write(os.path.join(work, "ex.azr"))
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
        check("control: the physical path fed params_rwa warns",
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
