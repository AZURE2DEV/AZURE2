#!/usr/bin/env python3
"""The phase shifts of models/phase.azr (p + 12C with the spin-orbit term; levels without widths, so delta_{J l} is the
potential's: delta = -phi of the channel at a_c) against physics/tests/data/spin_orbit.csv (make_spin_orbit_reference.py,
an integration independent of the engine's).  Usage: phase_check.py AZUREOut_aa=1_R=1.extrap spin_orbit.csv
Blocks of the extrapolation file, in the order of the records: p1half (J = 1/2, l = 1), p3half (3/2, 1), s1half (1/2, 0);
column 1 E_cm (MeV), column 4 the phase shift (deg).  Bound 1e-6 rad; exit 1 above it."""
import csv
import math
import sys

blocks, cur = [], []
for line in open(sys.argv[1]):
    if line.strip():
        cur.append([float(x) for x in line.split()])
    elif cur:
        blocks.append(cur)
        cur = []
if cur:
    blocks.append(cur)
ref = {}
for row in csv.DictReader(open(sys.argv[2])):
    if row["kind"] == "phi":
        ref[(int(row["l"]), int(row["two_j"]), float(row["E_MeV"]))] = float(row["value"])
worst, n = 0.0, 0
for (l, two_j), block in zip(((1, 1), (1, 3), (0, 1)), blocks):
    for r in block:
        delta = math.radians(r[3])
        d = abs(delta + ref[(l, two_j, round(r[0], 9))])
        worst, n = max(worst, d), n + 1
        print(f"  l = {l} 2J = {two_j} E_cm = {r[0]:.3f} MeV: delta {delta:+.9f} rad, -phi {-ref[(l, two_j, round(r[0], 9))]:+.9f}, diff {d:.1e}")
print(f"{n} phase shifts, worst {worst:.1e} rad")
sys.exit(0 if n == 12 and worst <= 1e-6 else 1)
