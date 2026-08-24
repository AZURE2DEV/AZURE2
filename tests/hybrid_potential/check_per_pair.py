#!/usr/bin/env python3
"""Behavioural checks on the per-pair hybrid nuclear potential.

run_tests.sh already pins this project's chi-squared, which proves the
<potential> block reaches the calculation at all.  What it cannot show is that
the *per-pair* semantics are right, because one project produces one number.
That is what this script covers:

  1. A pair's own setting beats the default.
  2. A pair without one follows the default.
  3. Switching a pair off returns exactly the no-hybrid answer.
  4. Setting a potential through the API equals loading the same potential
     from a file -- the two paths must not drift.
  5. The Coulomb memo is keyed on the potential, so changing one is seen.
  6. The parameter vector is re-derived by the change, and reusing a vector
     captured beforehand is what quietly gives a different answer.

Run from this directory:  python3 check_per_pair.py
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np

from pyazr import AzrModel, azure2

HERE = os.path.dirname(os.path.abspath(__file__))
AZR = os.path.join(HERE, "hybrid_potential.azr")

# The project as committed: default Woods-Saxon V0 = 40, pair 1 overriding it
# with V0 = 20, pair 2 (the photon pair) off.
FILE_CHI2 = 3659795.66031825
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


def variant(tmp, potential_block):
    """Write a copy of the project with a different <potential> block."""
    text = open(AZR).read()
    start = text.index("<potential>")
    end = text.index("</potential>") + len("</potential>")
    path = os.path.join(HERE, tmp)
    open(path, "w").write(text[:start] + potential_block + text[end:])
    return path


NO_HYBRID = """<potential>
useHybridPotential=0
useAdaptiveGrid=1
potentialType=0
V0=40
R=3.6
a=0.6
</potential>"""

DEFAULT_ONLY = """<potential>
useHybridPotential=1
useAdaptiveGrid=1
potentialType=0
V0=20
R=3.6
a=0.6
</potential>"""

print("1. the project as committed (pair 1 overrides the default)")
with azure2(AZR, cwd=HERE) as m:
    x = np.asarray(m.params_rwa, float)
    chi2_file = float(np.sum(m.calculate_chi2_rwa(x)))
    check("chi2 from the file", chi2_file, FILE_CHI2)
    check_true("pair 1 carries its own setting", m.nuclear_potential(1).own)
    check_true("pair 1 uses its own V0, not the default's",
               m.nuclear_potential(1).V0 == 20.0
               and m.nuclear_potential(0).V0 == 40.0,
               f"pair1 V0={m.nuclear_potential(1).V0}, "
               f"default V0={m.nuclear_potential(0).V0}")
    check_true("pair 2 is off", not m.nuclear_potential(2).enabled)

print("\n2. the same potential on the default alone gives the same answer")
print("   (pair 1 is the only particle pair, so inheriting V0=20 must match)")
path = variant("_default_only.azr", DEFAULT_ONLY)
with azure2(path, cwd=HERE) as m:
    chi2_default = float(np.sum(m.calculate_chi2_rwa(np.asarray(m.params_rwa, float))))
    check_true("pair 1 inherits, it has no setting of its own",
               not m.nuclear_potential(1).own)
    check("chi2 inheriting the default", chi2_default, FILE_CHI2)

print("\n3. no hybrid at all")
path = variant("_no_hybrid.azr", NO_HYBRID)
with azure2(path, cwd=HERE) as m:
    chi2_off = float(np.sum(m.calculate_chi2_rwa(np.asarray(m.params_rwa, float))))
    print(f"        baseline chi2 = {chi2_off:.6f}")
    check_true("the potential changes the answer",
               abs(chi2_off - FILE_CHI2) > 1e-3 * FILE_CHI2,
               f"{chi2_off} vs {FILE_CHI2}")

print("\n4. switching pair 1 off at runtime returns the baseline")
with azure2(AZR, cwd=HERE) as m:
    m.set_nuclear_potential(1, enabled=False)
    m.set_nuclear_potential(0, enabled=False)
    # Re-read: the potential changes the penetrabilities, which are what map
    # physical widths to reduced-width amplitudes, so the vector is re-derived.
    x = np.asarray(m.params_rwa, float)
    check("chi2 with everything off", float(np.sum(m.calculate_chi2_rwa(x))), chi2_off)

print("\n5. setting a potential through the API equals loading it from a file")
with azure2(variant("_no_hybrid.azr", NO_HYBRID), cwd=HERE) as m:
    m.set_nuclear_potential(0, type="WoodsSaxon", enabled=True, V0=40.0, R=3.6, a=0.6,
                            reinitialize=False)
    m.set_nuclear_potential(1, type="WoodsSaxon", enabled=True, V0=20.0, R=3.6, a=0.6)
    x = np.asarray(m.params_rwa, float)
    check("chi2 built through the API", float(np.sum(m.calculate_chi2_rwa(x))), FILE_CHI2)

print("\n6. the Coulomb memo is keyed on the potential, so a change is seen")
with azure2(variant("_no_hybrid.azr", NO_HYBRID), cwd=HERE) as m:
    seen = []
    for v0 in (10.0, 20.0, 30.0):
        m.set_nuclear_potential(1, type="WoodsSaxon", enabled=True, V0=v0, R=3.6, a=0.6)
        seen.append(float(np.sum(m.calculate_chi2_rwa(np.asarray(m.params_rwa, float)))))
    print("        chi2 at V0 = 10, 20, 30: " + ", ".join(f"{c:.1f}" for c in seen))
    check_true("every depth gives its own answer", len(set(seen)) == 3, str(seen))
    check("V0=20 matches the file", seen[1], FILE_CHI2)

print("\n7. a parameter vector captured before the change no longer fits the model")
with azure2(variant("_no_hybrid.azr", NO_HYBRID), cwd=HERE) as m:
    stale = np.asarray(m.params_rwa, float)
    m.set_nuclear_potential(1, type="WoodsSaxon", enabled=True, V0=20.0, R=3.6, a=0.6)
    with_stale = float(np.sum(m.calculate_chi2_rwa(stale)))
    with_fresh = float(np.sum(m.calculate_chi2_rwa(np.asarray(m.params_rwa, float))))
    check("re-read vector reproduces the file", with_fresh, FILE_CHI2)
    check_true("the stale vector does not", abs(with_stale - FILE_CHI2) > 1e-3 * FILE_CHI2,
               f"stale {with_stale} vs file {FILE_CHI2}")

for tmp in ("_default_only.azr", "_no_hybrid.azr"):
    p = os.path.join(HERE, tmp)
    if os.path.exists(p):
        os.remove(p)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all checks passed")
