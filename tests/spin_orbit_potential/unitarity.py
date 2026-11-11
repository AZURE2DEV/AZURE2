#!/usr/bin/env python3
"""Unitarity with the spin-orbit term (spin_orbit_potential.md section 4.2, test 6): n + 16O with a Woods-Saxon and a
spin-orbit term on the neutron (format 3, converted from a model written by tests/upgrade/lib), s, p1/2, p3/2, d3/2,
d5/2 levels; elastic scattering the only open channel, so U_J = exp(2 i delta_J) and the total cross section from the
optical theorem equals the angle-integrated elastic one at every energy (296 energies, to 1e-9; the outputs carry 11
digits).  Also: the term moves sigma_el (it is applied).  Usage: unitarity.py AZURE2 WORKDIR"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "upgrade", "lib"))
import upgrade_project as up  # noqa: E402

binary, work = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
M16O, M17O = 15.99491461957, 16.99913175650
n16o = up.Pair(1, 0.5, 1, 0.0, 1, up.MASS["n"], M16O, 0, 8, up.separation(up.MASS["n"], M16O, M17O), R=4.0)
C = up.Channel
model = up.Project(
    "unitarity", [n16o],
    [up.Level(0.5, 1, 5.0, [C(n16o, 1, 0, 1.0e5)]), up.Level(1.5, -1, 5.5, [C(n16o, 1, 2, 5.0e4)]),
     up.Level(0.5, -1, 5.7, [C(n16o, 1, 2, 3.0e4)]), up.Level(2.5, 1, 6.0, [C(n16o, 1, 4, 2.0e4)]),
     up.Level(1.5, 1, 6.2, [C(n16o, 1, 4, 1.0e4)]), up.Level(0.5, 1, 15.0, [C(n16o, 1, 0, 5.0e6)])],
    data=[], tests=[up.TestSegment(1, 1, 0.05, 3.0, 0.01), up.TestSegment(1, 1, 0.05, 3.0, 0.01, extra="total 1")])
os.makedirs(work, exist_ok=True)
azr = model.write(work, binary)   # format 3 (upgrade_project converts its old text)
text = open(azr).read()
so = "    spin_orbit     woods_saxon depth 7 MeV radius 3.2 fm diffuseness 0.6 fm on light operator 2l.s\n"
record = ("<potentials>\npotential n16O\n    pair           n+16O\n    coulomb_radius 3.5 fm\n"
          "    central        woods_saxon depth 40 MeV radius 3.2 fm diffuseness 0.6 fm\n" + so + "</potentials>\n")
out = {}
for name, rec in (("so", record), ("central", record.replace(so, ""))):
    d = os.path.join(work, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "m.azr"), "w") as f:
        f.write(text.replace("<levels>\n", rec + "<levels>\n", 1))
    up.run(binary, d, os.path.join(d, "m.azr"), up.EXTRAP)
    out[name] = up.read_blocks(os.path.join(d, "output", "AZUREOut_aa=1_R=1.extrap"))
el, tot = out["so"]
worst = float(np.max(np.abs(tot[:, 3] / el[:, 3] - 1.0)))
moved = float(np.max(np.abs(el[:, 3] / out["central"][0][:, 3] - 1.0)))
print(f"{len(el)} energies: |sigma_tot / sigma_el - 1| <= {worst:.1e}; the term moves sigma_el by up to {moved:.1e}")
sys.exit(0 if len(el) == len(tot) == 296 and worst <= 1e-9 and moved > 1e-3 else 1)
