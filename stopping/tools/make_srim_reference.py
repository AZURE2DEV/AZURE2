#!/usr/bin/env python3
"""Reference values of the stopping module's tests from LunaSRIM (stopping/tests/data/srim_reference.csv).

LunaSRIM (github.com/skowrons94/LunaSRIM) is the Python re-implementation of SRIM-2013 that stopping/ is ported
from; these values are its own, computed with its own data files, so the C++ tests compare the port with the
original.  Run with any Python that has the package `lunasrim` (version 1.0.0 was used):

    python3 stopping/tools/make_srim_reference.py      # rewrites stopping/tests/data/srim_reference.csv

Rows (comma-separated, first field the kind; energies in keV, stopping in eV/(1e15 atoms/cm^2), areal densities in
atoms/cm^2, widths in keV):
  se,Z1,M1,Z2,phase,E,S_e            lunasrim.stopping.electronic_stopping (phase factor included, no compound corr.)
  sn,Z1,M1,Z2,M2,E,S_n               lunasrim.stopping.nuclear_stopping, M2 SRIM's atomic weight of Z2
  compound,formula,phase,C,Z1,M1,E,S_e,S_n   Bragg's rule with the compound correction C applied once, the sum
                                     lunasrim.straggling.energy_loss and SRIM's table mode use
  layer,formula,phase,Z1,M1,E,N,E_out,sigma_trim,sigma_bohr   lunasrim.straggling.energy_loss (RK4, 2000 steps)
phase is 0 (solid) or 1 (gas), as in LunaSRIM.
"""
import math
import os
import re
import sys

import lunasrim
from lunasrim import data_files, stopping, straggling

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "tests", "data", "srim_reference.csv")

# Projectiles: Z1, mass (u; AME2020 atomic masses, as stopping/src/Material.cpp's Proton() ... and the tests use)
PROJECTILES = [("p", 1, 1.00782503223), ("d", 1, 2.01410177812), ("3He", 2, 3.01602932197),
               ("a", 2, 4.00260325413), ("7Li", 3, 7.0160034366), ("12C", 6, 12.0),
               ("14N", 7, 14.00307400443), ("16O", 8, 15.99491461957)]
ELEMENTS = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "Ca", "Ti", "Fe", "Ni",
            "Cu", "Zr", "Ag", "Ta", "W", "Au", "Pb"]
COMPOUNDS = [("TiN", 0), ("Ta2O5", 0), ("TiH2", 0), ("LiF", 0), ("CaF2", 0), ("H2O", 0), ("CH4", 1), ("CO2", 1)]


def boundaries(z1):
    """SRIM's branch boundaries in keV/u for ion z1 (lunasrim/stopping.py)."""
    if z1 == 1:
        return [2.0, 1000.0, 1200.0]
    if z1 <= 3:
        return [2.0, 2000.0, 2400.0]
    return [2.0, 2000.0, 2400.0, 30000.0]


def energies_keV_u(z1):
    """Every branch: below 2 keV/u, the fit, the join, Bethe (and above 30 MeV/u for heavy ions); each boundary exactly
    and 1e-3 on each side."""
    grid = [0.3, 1.0, 10.0, 50.0, 150.0, 400.0, 800.0, 1100.0, 1600.0, 2200.0, 3000.0, 8000.0, 20000.0]
    if z1 >= 4:
        grid += [45000.0, 100000.0]
    for b in boundaries(z1):
        grid += [b * (1 - 1e-3), b, b * (1 + 1e-3)]
    return sorted(set(grid))


def formula(text):
    """[(Z, count)] of a chemical formula."""
    out = []
    for sym, count in re.findall(r"([A-Z][a-z]?)([0-9.]*)", text):
        out.append((data_files.symbol_to_z(sym), float(count) if count else 1.0))
    return out


def main():
    sc = data_files.load_scoef()
    weight = {z: sc[z].natural_mass for z in sc}
    rows = ["# written by stopping/tools/make_srim_reference.py with lunasrim " + lunasrim.__version__
            + "; do not edit (kinds: see the script)"]
    elements = [data_files.symbol_to_z(s) for s in ELEMENTS]
    for name, z1, m1 in PROJECTILES:
        es = energies_keV_u(z1)
        for z2 in elements:
            phases = (0, 1) if z2 in stopping.GAS_ELEMENTS else (0,)
            for phase in phases:
                for e in es:
                    ek = e * m1
                    rows.append(f"se,{z1},{m1!r},{z2},{phase},{ek!r},{stopping.electronic_stopping(z1, m1, z2, ek, phase)!r}")
            # nuclear: a grid in the reduced energy, with the change of form at 30 (1e-6 on each side: exactly at 30
            # rounding could take either branch)
            m2 = weight[z2]
            eps1 = stopping.reduced_energy(z1, m1, z2, m2, 1.0)
            for eps in (1e-4, 1e-3, 0.01, 0.1, 1.0, 5.0, 30.0 * (1 - 1e-6), 30.0 * (1 + 1e-6), 100.0, 1000.0):
                ek = eps / eps1
                rows.append(f"sn,{z1},{m1!r},{z2},{m2!r},{ek!r},{stopping.nuclear_stopping(z1, m1, z2, m2, ek)!r}")
    for text, phase in COMPOUNDS:
        parts = formula(text)
        tot = sum(c for _, c in parts)
        for corr in (1.0, 0.93, 1.08):
            for name, z1, m1 in PROJECTILES:
                for e in (1.0, 20.0, 300.0, 1500.0, 5000.0):
                    ek = e * m1
                    se = sum(c / tot * stopping.electronic_stopping(z1, m1, z, ek, phase, corr=corr) for z, c in parts)
                    sn = sum(c / tot * stopping.nuclear_stopping(z1, m1, z, weight[z], ek) for z, c in parts)
                    rows.append(f"compound,{text},{phase},{corr!r},{z1},{m1!r},{ek!r},{se!r},{sn!r}")
    # layers: (formula, phase, projectile, E keV, N atoms/cm^2)
    layers = [("Ta", 0, "p", 300.0, 1e17), ("Ta", 0, "p", 1000.0, 2e18), ("TiN", 0, "p", 200.0, 5e17),
              ("TiN", 0, "p", 2000.0, 5e18), ("C", 0, "a", 1000.0, 1e18), ("C", 0, "a", 5000.0, 1e19),
              ("C", 0, "12C", 12000.0, 5e17), ("C", 0, "12C", 30000.0, 3e18)]
    proj = {name: (z1, m1) for name, z1, m1 in PROJECTILES}
    for text, phase, pname, ek, n in layers:
        z1, m1 = proj[pname]
        parts = formula(text)
        els = [(z, c, weight[z]) for z, c in parts]
        r = straggling.energy_loss(z1, m1, ek, els, 5.0, n, unit="atoms/cm2", phase=phase, steps=2000)
        rows.append(f"layer,{text},{phase},{z1},{m1!r},{ek!r},{n!r},{r['energy_out']!r},{r['sigma']!r},{r['sigma_bohr']!r}")
    with open(OUT, "w") as f:
        f.write("\n".join(rows) + "\n")
    print(f"{len(rows) - 1} rows -> {os.path.normpath(OUT)}", file=sys.stderr)


if __name__ == "__main__":
    main()
