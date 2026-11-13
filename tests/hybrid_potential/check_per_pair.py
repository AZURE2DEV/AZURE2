#!/usr/bin/env python3
"""Behavioural checks on the per-pair hybrid nuclear potential.

run_tests.sh already pins this project's chi-squared, which proves the
<potentials> section reaches the calculation at all.  What it cannot show is that
the *per-pair* semantics are right, because one project produces one number.
That is what this script covers:

  1. The record of pair 1 (p + 12C) is pair 1's setting; there is no default.
  2. (The old <potential> block's default, inherited by a pair without a section of
     its own, is the converter's: it writes one record per particle pair, and
     convert_spec_equivalence compares the potential of every particle pair.)
  3. Switching a pair off returns exactly the no-hybrid answer.
  4. Setting a potential through the API equals loading the same potential
     from a file -- the two paths must not drift.
  5. The Coulomb memo is keyed on the potential, so changing one is seen.
  6. The parameter vector is re-derived by the change, and reusing a vector
     captured beforehand is what quietly gives a different answer.

The project and its variants are run in a temporary copy of this folder.
Needs the compiled engine (ctest pyazr_hybrid_per_pair; PYAZR_ENGINE_DIR for a
module built elsewhere):  python3 check_per_pair.py
"""
import atexit
import os
import shutil
import sys
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "2")
SOURCE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(SOURCE, "..", "..")))

try:
    import numpy as np
    from pyazr import _azure2  # noqa: F401  (the compiled engine)
    from pyazr import Session
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

HERE = tempfile.mkdtemp(prefix="check_per_pair.")
atexit.register(shutil.rmtree, HERE, ignore_errors=True)
os.makedirs(os.path.join(HERE, "output"))
shutil.copytree(os.path.join(SOURCE, "data"), os.path.join(HERE, "data"))
AZR = os.path.join(HERE, "hybrid_potential.azr")
shutil.copy(os.path.join(SOURCE, "hybrid_potential.azr"), AZR)


def committed_total():
    """Total-Chi-Squared of expected/chiSquared.out, the reference tests/run_tests.sh compares this project with."""
    with open(os.path.join(SOURCE, "expected", "chiSquared.out")) as f:
        for line in f:
            if line.startswith("Total-Chi-Squared:"):
                return float(line.split()[1])
    raise RuntimeError("no Total-Chi-Squared line in expected/chiSquared.out")


def committed_tolerance():
    """The project's relative tolerance on that total (the first number of ./tolerance, as tests/run_tests.sh reads
    it): the reference is printed to six digits and recorded with one toolchain."""
    with open(os.path.join(SOURCE, "tolerance")) as f:
        for line in f:
            if line.strip() and not line.lstrip().startswith("#"):
                return float(line.split()[0])
    raise RuntimeError("no number in tolerance")


# The project as committed: one record, pair 1 (p + 12C) Woods-Saxon V0 = 20 (the
# old file's pair section over a default of V0 = 40); pair 2 (the photon pair) has none.  Its chi-squared is the one of
# expected/chiSquared.out (check 1, within the project's tolerance); every other
# path that must give the same project is compared with the value computed here
# (FILE_CHI2, set in check 1) to round-off.
EXPECTED_CHI2 = committed_total()
EXPECTED_TOL = committed_tolerance()
TOL = 1e-6           # relative; the two paths should agree to round-off

failures = []


def check(name, got, want, tol=TOL):
    ok = abs(got - want) <= tol * max(abs(want), 1.0)
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}: {got:.6f}"
          f"{'' if ok else f' (expected {want:.6f})'}")
    if not ok:
        failures.append(name)


def check_true(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}{'' if ok else f' -- {detail}'}")
    if not ok:
        failures.append(name)


def variant(tmp, records):
    """Write a copy of the project with other records in <potentials>."""
    text = open(AZR).read()
    start = text.index("<potentials>")
    end = text.index("</potentials>") + len("</potentials>")
    path = os.path.join(HERE, tmp)
    open(path, "w").write(text[:start] + "<potentials>\n" + records + "</potentials>" + text[end:])
    return path


NO_HYBRID = ""      # no record: no pair has a potential

print("1. the project as committed (pair 1's record)")
with Session(AZR, cwd=HERE) as m:
    x = np.asarray(m.params_rwa, float)
    chi2_file = m.evaluate(x).chi2
    check("chi2 from the file equals expected/chiSquared.out", chi2_file, EXPECTED_CHI2, EXPECTED_TOL)
    FILE_CHI2 = chi2_file
    check_true("pair 1 carries its own setting", m.nuclear_potential(1).own)
    check_true("pair 1 uses its record's V0 = 20; there is no default",
               m.nuclear_potential(1).V0 == 20.0
               and not m.nuclear_potential(0).enabled,
               f"pair1 V0={m.nuclear_potential(1).V0}, "
               f"default {m.nuclear_potential(0)}")
    check_true("pair 2 is off", not m.nuclear_potential(2).enabled)

print("\n3. no hybrid at all")
path = variant("_no_hybrid.azr", NO_HYBRID)
with Session(path, cwd=HERE) as m:
    chi2_off = m.evaluate(np.asarray(m.params_rwa, float)).chi2
    print(f"        baseline chi2 = {chi2_off:.6f}")
    # 1e-4, not the 1e-3 this used to demand. The potential's effect on this
    # project's chi-squared is 2.8e-4 relative (2466827.58 against 2466142.31
    # with the distorted bound states; 2.9e-4 before them). The larger figure the old threshold assumed came from
    # an unconverged target integral, in which the grid itself responded to the
    # potential through the penetrability in its resonance-width estimate;
    # refining the grid removes that in the pre-fix code too (see ../tolerance).
    # Both values here come from the same process, so the comparison is
    # deterministic and 1e-4 still leaves a factor of three of margin.
    check_true("the potential changes the answer",
               abs(chi2_off - FILE_CHI2) > 1e-4 * FILE_CHI2,
               f"{chi2_off} vs {FILE_CHI2}")

print("\n4. switching pair 1 off at runtime returns the baseline")
with Session(AZR, cwd=HERE) as m:
    m.set_nuclear_potential(1, enabled=False)
    m.set_nuclear_potential(0, enabled=False)
    # Re-read: the potential changes the penetrabilities, which are what map
    # physical widths to reduced-width amplitudes, so the vector is re-derived.
    x = np.asarray(m.params_rwa, float)
    check("chi2 with everything off", m.evaluate(x).chi2, chi2_off)

print("\n5. setting a potential through the API equals loading it from a file")
with Session(variant("_no_hybrid.azr", NO_HYBRID), cwd=HERE) as m:
    m.set_nuclear_potential(1, type="WoodsSaxon", enabled=True, V0=20.0, R=3.6, a=0.6)
    x = np.asarray(m.params_rwa, float)
    check("chi2 built through the API", m.evaluate(x).chi2, FILE_CHI2)

print("\n6. the Coulomb memo is keyed on the potential, so a change is seen")
with Session(variant("_no_hybrid.azr", NO_HYBRID), cwd=HERE) as m:
    seen = []
    for v0 in (10.0, 20.0, 30.0):
        m.set_nuclear_potential(1, type="WoodsSaxon", enabled=True, V0=v0, R=3.6, a=0.6)
        seen.append(m.evaluate(np.asarray(m.params_rwa, float)).chi2)
    print("        chi2 at V0 = 10, 20, 30: " + ", ".join(f"{c:.1f}" for c in seen))
    check_true("every depth gives its own answer", len(set(seen)) == 3, str(seen))
    check("V0=20 matches the file", seen[1], FILE_CHI2)

print("\n7. a parameter vector captured before the change no longer fits the model")
with Session(variant("_no_hybrid.azr", NO_HYBRID), cwd=HERE) as m:
    stale = np.asarray(m.params_rwa, float)
    m.set_nuclear_potential(1, type="WoodsSaxon", enabled=True, V0=20.0, R=3.6, a=0.6)
    with_stale = m.evaluate(stale).chi2
    with_fresh = m.evaluate(np.asarray(m.params_rwa, float)).chi2
    check("re-read vector reproduces the file", with_fresh, FILE_CHI2)
    check_true("the stale vector does not", abs(with_stale - FILE_CHI2) > 1e-3 * FILE_CHI2,
               f"stale {with_stale} vs file {FILE_CHI2}")

for tmp in ("_no_hybrid.azr",):
    p = os.path.join(HERE, tmp)
    if os.path.exists(p):
        os.remove(p)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all checks passed")
