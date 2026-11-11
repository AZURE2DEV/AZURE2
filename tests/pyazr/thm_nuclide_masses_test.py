#!/usr/bin/env python3
"""One mass table (docs/upgrade/project.md section 5, step 5): the engine's built-in THM nuclei (ThmExperiment.cpp
kNuclides, which pyazr reads through _azure2.thm_nuclides) carry the nuclear masses of the AME2020 file that
pyazr.masses reads (pyazr/data/mass_1.mas20):

    M_nuc = M_atomic - Z m_e + B_e(Z),   B_e = 14.4381 Z^2.39 + 1.55468e-6 Z^5.35 eV  (Lunney et al., RMP 75 (2003)
                                                                                     1021, eq. A4)

to 0.6e-10 u (the table's 10 decimals), with m_e = physics Constants.h kElectronMass_MeV.  A nuclide added to the
engine's table, or a newer mass file, fails here until the two agree.  Needs the compiled engine; skips without it.
Run from anywhere:  python3 tests/pyazr/thm_nuclide_masses_test.py
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
try:
    from pyazr import masses, _azure2
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

consts = open(os.path.join(ROOT, "physics", "include", "azure2", "physics", "Constants.h")).read()
me = float(re.search(r"kElectronMass_MeV\s*=\s*([0-9.]+)", consts).group(1)) / masses.U_MEV   # u
failures = []
table = _azure2.thm_nuclides()
for name, Z, A, mass in table:
    entry = masses.lookup((Z, A))
    be = (14.4381 * Z ** 2.39 + 1.55468e-6 * Z ** 5.35) * 1e-6 / masses.U_MEV               # u
    expected = entry.atomic_mass - Z * me + be
    ok = abs(mass - expected) < 0.6e-10
    print(f"  {'ok  ' if ok else 'FAIL'}  {name:5s} {mass:.10f}  AME2020 {expected:.10f}  ({(mass - expected) * 1e10:+.2f}e-10 u)")
    if not ok:
        failures.append(name)
print(f"{len(table)} nuclei: " + ("PASS" if not failures else f"FAIL {failures}"))
sys.exit(1 if failures else 0)
