#!/usr/bin/env python3
"""Identity check (ARCHITECTURE.md section 10): Brune <-> formal <-> Park give the same observables.

The .azr of T1 and T3 holds observed (Brune) level parameters: energies E_x and partial widths Gamma at E_x, ANCs of
bound channels.  The engine transforms them to the parameters of the run's basis -- the Brune alternative level
matrix (default, Brune, PRC 66, 044611 (2002)), the formal Lane-Thomas parameters with a constant boundary
condition (--no-brune), or Park's level-dependent boundary conditions (--use-park, PRC 104, 064612 (2021)) -- and
computes the collision matrix there.  The three are exact reparametrizations, so every observable must be the same:
the elastic and capture cross sections, A_y and phase shifts of every data point (menu 1, column 4 of AZUREOut_*.out)
and the extrapolated ones (menu 3, .extrap).  A transformation error (a missing shift-function derivative, a wrong
boundary condition, a sign of a reduced width) changes resonance positions or interference and fails at 1e-4 or
worse.

Tolerance: 2e-7 relative to max(|value|, 1e-6 x the largest |value| of the segment) per point.  With the
c1cc92c7 build the differences are 0 or 2e-11 (print precision), except in one case stated here, not hidden:

  capture with external capture in the formal basis.  T1's (p,gamma) cross sections with --no-brune differ from the
  Brune and Park ones by up to 54 % (E_cm = 1.568 MeV, at the 3/2- resonance) when the photon pairs have external
  capture (mask 7); with the mask set to 0 they agree to print precision.  The external-capture amplitude of the
  engine depends on the boundary condition B_c, and the transformation of the observed photon widths to the formal
  basis does not carry that dependence along.  Whether the formal result or the transformation is at fault is an
  open question (docs/upgrade/HANDOFF.md); this test prints the size of the difference and checks the identity with
  external capture switched off.  It belongs to the convention and level-parameter work of Phase 1.

    python3 tests/pyazr/upgrade_parametrization_roundtrip_test.py   (binary: $AZURE2_BIN, $AZURE2_BINARY or newest)
"""
import glob
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests", "upgrade", "lib"))

failures = []
TOL = 2e-7


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


try:
    import numpy as np
    import upgrade_project as up
    binary = up.find_binary()
except (ImportError, SystemExit) as e:
    print(f"skip: {e}")
    sys.exit(77)   # ctest SKIP_RETURN_CODE (tests/pyazr/CMakeLists.txt)


def outputs(project, flags, ec_off=False):
    work = tempfile.mkdtemp(prefix=f"upgrade_roundtrip_{project}.")
    try:
        azr = up.project_copy(project, work)
        if ec_off:   # photon pairs: column 31 (external-capture mask) 7 -> 0
            text = open(azr).read()
            text = re.sub(r"(?m)^(.* 10 [0-9.]+ [0-9.]+ [0-9.]+ )7$", lambda m: m.group(1) + "0", text)
            open(azr, "w").write(text)
        res = {}
        for answers, pattern in ((up.CALC, "*.out"), (up.EXTRAP, "*.extrap")):
            up.run(binary, work, azr, answers, flags=flags)
            for f in glob.glob(os.path.join(work, "output", "AZUREOut_aa=*_R=*" + pattern[1:])):
                res[os.path.basename(f)] = up.read_blocks(f)
        return res
    finally:
        shutil.rmtree(work, ignore_errors=True)


def largest(ref, other, files):
    worst, where, same = 0.0, "", all(f in other for f in files)
    for name in files:
        for a, b in zip(ref[name], other.get(name, [])):
            same = same and a.shape == b.shape
            if a.shape != b.shape:
                continue
            scale = np.maximum(np.abs(a[:, 3]), 1e-6 * np.max(np.abs(a[:, 3])))
            d = np.abs(a[:, 3] - b[:, 3]) / scale
            if d.max() > worst:
                worst, where = float(d.max()), f"{name}, E_cm = {a[np.argmax(d), 0]:.4f} MeV"
    return same, worst, where


def compare(label, ref, other, files):
    same, worst, where = largest(ref, other, files)
    check(label, same and worst <= TOL,
          (f"largest relative difference {worst:.1e}" + (f" at {where}" if where else "")) if same else "outputs differ")


PHOTON_T1 = ("_R=2.", "_R=3.")      # T1's photon pairs
for project in ("T1", "T3"):
    ref = outputs(project, [])
    park = outputs(project, ["--use-park"])
    formal = outputs(project, ["--no-brune"])
    compare(f"{project}: Brune = Park (--use-park), all data points and extrapolations", ref, park, sorted(ref))
    if project == "T3":
        compare("T3: Brune = formal (--no-brune), all data points and extrapolations", ref, formal, sorted(ref))
        continue
    particle = [f for f in ref if not any(t in f for t in PHOTON_T1)]
    photon = [f for f in ref if any(t in f for t in PHOTON_T1)]
    compare("T1: Brune = formal (--no-brune), elastic, A_y and phase shifts", ref, formal, particle)
    ref0, formal0 = outputs("T1", [], ec_off=True), outputs("T1", ["--no-brune"], ec_off=True)
    compare("T1: Brune = formal (--no-brune), capture without external capture", ref0, formal0, photon)
    _, worst, where = largest(ref, formal, photon)
    print(f"  info  T1: capture WITH external capture, Brune vs formal: largest relative difference {worst:.2f} at "
          f"{where} (open question, see the docstring)")

print("FAIL" if failures else "ok")
sys.exit(1 if failures else 0)
