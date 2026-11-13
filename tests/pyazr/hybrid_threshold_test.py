#!/usr/bin/env python3
"""The shift function of a hybrid pair is continuous through threshold (docs/upgrade/spin_orbit_potential.md
sections 1.4 item 1 and 2.4).

Below threshold S is a W~'/W~ of the decaying solution of the radial equation with the pair's potential (physics/
DecayingSolution, HybridBoundShift); above it rho H'/H of the outgoing one.  Builds before 2026-11-10 took S below
threshold from the undistorted Whittaker function, so S of tests/hybrid_potential (p + 12C, Woods-Saxon 20 / 3.6 /
0.6 MeV, fm, a_c = 3.6 fm) jumped at E = 0 by 1.07 (l = 0), 0.92 (l = 1), 0.75 (l = 2) (work/so-note/threshold.py).

Checks through `Session.coulomb_functions` (S at the channel radius, E <= 0 included), pair 1, l = 0-2:
  (1) |S(E) - S(0)| <= 1e-3 for |E| <= 1e-3 MeV, on both sides and in the closed band (eta > 100, E < 8e-5 MeV,
      where S = 2 S(0) - S(-E)); dS/dE is about 0.3-0.8 MeV^-1 there, so 1e-3 MeV moves S by < 1e-3;
  (2) S(-1.94351 MeV) (the 13N ground state) against the independent DOP853 integration
      (physics/tests/data/hybrid_bound.csv, make_hybrid_bound_reference.py) to 1e-6.
Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/hybrid_threshold_test.py
"""
import csv
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "hybrid_potential")
REF = os.path.join(ROOT, "physics", "tests", "data", "hybrid_bound.csv")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import Session
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

ref = {}
with open(REF) as f:
    for row in csv.DictReader(f):
        ref[(row["kind"], int(row["l"]), float(row["E_MeV"]))] = float(row["value"])

work = tempfile.mkdtemp(prefix="hybrid_threshold.")
try:
    os.makedirs(os.path.join(work, "output"))
    shutil.copytree(os.path.join(SRC, "data"), os.path.join(work, "data"))
    shutil.copy(os.path.join(SRC, "hybrid_potential.azr"), work)
    A = Session(os.path.join(work, "hybrid_potential.azr"), cwd=work)
    near = np.array([-1e-3, -1e-4, -1e-5, -1e-7, 0.0, 1e-7, 1e-5, 5e-5, 1e-4, 1e-3])
    for L in (0, 1, 2):
        S = A.coulomb_functions(1, near, L=L)["S"]
        s0 = S[near == 0.0][0]
        jump = float(np.max(np.abs(S - s0)))
        check(f"l = {L}: S continuous through threshold (max |S(E) - S(0)| {jump:.1e} for |E| <= 1 keV, S(0) {s0:.5f})",
              jump <= 1e-3, f"{jump:.3e} > 1e-3")
        gs = A.coulomb_functions(1, np.array([-1.94351]), L=L)["S"][0]
        d = abs(gs - ref[("S", L, -1.94351)])
        check(f"l = {L}: S(-1.94351 MeV) = {gs:.8f} against DOP853 to {d:.1e}", d <= 1e-6, f"{d:.3e} > 1e-6")
except Exception as e:  # pragma: no cover
    check("session", False, repr(e))
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("all checks passed")
