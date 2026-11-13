#!/usr/bin/env python3
"""Pseudo-data of tests/gamma_decay: the capture cross sections of partial.azr (the 14N(p,g)15O model of
tests/thick_target_resonance, level 5 the 259.7 keV resonance with its <levels> radiative widths, the truth) at 18 lab
energies across the resonance, to five final states, with a 3 % error and a fixed deterministic scatter of
0.03 sin(1.7 k + 0.9 m) (k the point, m the data set), so that the chi^2 minimum is not 0.

    python3 make_data.py OUTPUT_DIR    # OUTPUT_DIR: output/ of a menu-1 run of partial.azr (any data values)
"""
import math
import sys

out = sys.argv[1]
sets = [("g0", 2), ("g6793", 3), ("g6176", 6), ("g5183", 4), ("g6859", 7)]  # data set, exit pair key
for m, (name, key) in enumerate(sets):
    rows = [l.split() for l in open(f"{out}/AZUREOut_aa=1_R={key}.out") if l.strip()]
    elab = [l.split()[0] for l in open(f"data/{name}.dat")]
    with open(f"data/{name}.dat", "w") as f:
        for k, (e, r) in enumerate(zip(elab, rows)):
            model = float(r[3])
            f.write(f"{e} 0 {model * (1 + 0.03 * math.sin(1.7 * k + 0.9 * m)):.8e} {0.03 * model:.8e}\n")
