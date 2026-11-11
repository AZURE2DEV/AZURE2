#!/usr/bin/env python3
"""The hybrid nuclear potential belongs to the session, not to the process (audit 2026-10-25, item 30; 2026-10-27).

`NuclearPotentialManager` is a process singleton.  Up to dev 2363ada `Config::ReadPotentialBlock` of a later pyazr
session (or `set_potential` on it) replaced the settings an earlier live session was using: session A, built on a file
with V0 = 80 MeV, computed with B's V0 = 120 MeV from the moment B existed.  Now each session keeps its own potential
state and every engine call installs it (ConfigScope), so sessions cannot reconfigure each other.

Checks on tests/hybrid_potential (p + 12C, hybrid model on):
  (1) the shift function of session A (V0 = 80) is the same before and after session B (V0 = 120) is created and used;
  (2) A and B differ (teeth);
  (3) after B is deleted, A is still itself.
Needs the compiled engine; skips cleanly without it.

Run from anywhere:  python3 tests/pyazr/potential_sessions_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "hybrid_potential", "hybrid_potential.azr")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(0)

work = tempfile.mkdtemp(prefix="potential_sessions.")
try:
    os.makedirs(os.path.join(work, "output"))
    shutil.copytree(os.path.join(os.path.dirname(SRC), "data"), os.path.join(work, "data"))
    text = open(SRC).read()
    s = text.index("<potential>")
    e = text.index("</potential>") + len("</potential>")

    def variant(name, v0):
        path = os.path.join(work, name + ".azr")
        open(path, "w").write(text[:s] + "<potential>\nuseHybridPotential=1\nuseAdaptiveGrid=1\npotentialType=0\n"
                              f"V0={v0}\nR=3.6\na=0.6\n</potential>" + text[e:])
        return path

    energies = np.array([0.5, 1.0, 2.0, 3.0])
    A = azure2(variant("a80", 80.0), cwd=work)
    before = A.coulomb_functions(1, energies, L=1)["S"].copy()
    B = azure2(variant("b120", 120.0), cwd=work)
    b_values = B.coulomb_functions(1, energies, L=1)["S"].copy()
    after = A.coulomb_functions(1, energies, L=1)["S"].copy()
    check(f"the two potentials give different shift functions (max |dS| {np.max(np.abs(before - b_values)):.2e})",
          np.max(np.abs(before - b_values)) > 1e-4)
    check(f"session A is unchanged by session B (max |dS| {np.max(np.abs(after - before)):.2e})",
          np.max(np.abs(after - before)) < 1e-12, "A computed with B's potential")
    del B
    again = A.coulomb_functions(1, energies, L=1)["S"].copy()
    check("session A is unchanged after B is gone", np.max(np.abs(again - before)) < 1e-12)
    reported = A.nuclear_potential(0)
    v0 = reported["V0"] if isinstance(reported, dict) else getattr(reported, "V0", None)
    check(f"A still reports its own V0 = 80 ({reported})", v0 is not None and abs(v0 - 80.0) < 1e-9,
          "nuclear_potential reports the other session's setting")
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
