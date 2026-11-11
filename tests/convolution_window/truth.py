#!/usr/bin/env python3
"""Full-Gaussian convolution of the bare 16O(p,g) excitation function at the conv.azr test points.

Reads models/points.txt (k, E_lab) and a bare mode-3 output (AZUREOut_aa=1_R=3.extrap: E_cm, E_x, angle, sigma, ...),
convolves sigma(E_cm) with a Gaussian of sigma_cm = SIGMA_LAB * F about E_cm = E_lab * F over the WHOLE grid
(+-16 sigma here; trapezoid, 0.05 keV lab steps) and writes expected/truth.dat: k  E_lab  sigma_conv (barn).

    python3 truth.py path/to/bare/output/AZUREOut_aa=1_R=3.extrap
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_models as mm

HERE = os.path.dirname(os.path.abspath(__file__))
bare = np.loadtxt(sys.argv[1])
e, s = bare[:, 0], bare[:, 3]
pts = np.loadtxt(os.path.join(HERE, "models", "points.txt"))
scm = mm.SIGMA_LAB * mm.F
with open(os.path.join(HERE, "expected", "truth.dat"), "w") as f:
    f.write("# k  E_lab(MeV)  sigma_conv(barn)  -- truth.py: full Gaussian (sigma_cm = %.6g MeV) on the bare grid\n" % scm)
    for k, elab in pts:
        e0 = elab * mm.F
        w = np.exp(-0.5 * ((e - e0) / scm) ** 2)
        v = np.trapezoid(w * s, e) / np.trapezoid(w, e)
        f.write(f"{k:5.1f}  {elab:.6f}  {v:.10e}\n")
        print(f"{k:5.1f}  {elab:.6f}  {v:.6e}")
