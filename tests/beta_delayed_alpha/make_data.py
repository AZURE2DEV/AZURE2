#!/usr/bin/env python3
"""Regenerate data/bda.dat: the model's own beta-delayed alpha spectrum on a
0.6-3.0 MeV grid (50 keV steps), 5% uncertainties, with 5% Gaussian pseudo-noise
from a fixed seed so the chi-squared of the shipped levels is of order the
number of points (a value ~0 would make run_tests.sh's relative comparison
meaningless).  Any sign error on a beta amplitude shows up as a large one.

Run from this directory with the engine built:  python3 make_data.py
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
import numpy as np
from pyazr import azure2

azr = os.path.join(HERE, "16N.azr")
energies = [round(0.6 + 0.05 * i, 3) for i in range(49)]
with open(os.path.join(HERE, "data", "bda.dat"), "w") as f:
    for e in energies:
        f.write(f"{e:.3f}\t0\t1.0\t0.05\n")
with azure2(azr, cwd=HERE, use_brune=True, ignore_externals=True,
            use_long_wavelength=False) as m:
    y = np.asarray(m.calculate_rwa(m.params_rwa)[0], float)
assert y.size == len(energies), (y.size, len(energies))
noise = np.random.default_rng(16).normal(0.0, 0.05, y.size)
with open(os.path.join(HERE, "data", "bda.dat"), "w") as f:
    f.write("# 16N beta-delayed alpha spectrum computed from 16N.azr (make_data.py);\n"
            "# columns: E_cm(12C+a) [MeV], angle, yield x (1 + 5% noise, seed 16), 5% uncertainty\n")
    for e, v, n in zip(energies, y, noise):
        f.write(f"{e:.3f}\t0\t{v * (1.0 + n):.6e}\t{0.05 * abs(v):.6e}\n")
print(f"wrote {len(energies)} points, max yield {y.max():.4g}")
