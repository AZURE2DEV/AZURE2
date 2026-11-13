#!/usr/bin/env python3
"""Writes models/data/yield_curve.dat, the synthetic data of the worked example models/yield_story.azr: the model of
that project (the levels of tests/upgrade/T2/T2.azr) per uC at 33 lab energies 0.8 .. 1.2 MeV (in keV), times
(1 + 0.04 sin(2.3 i)), four significant digits, 4 % errors.  Run once with the engine of the commit that added it:

    PYAZR_ENGINE_DIR=<build>/pyazr-out python3 tests/units_scales/make_story.py
"""
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
import numpy as np  # noqa: E402
from pyazr import Session  # noqa: E402

E_KEV = [800 + 12.5 * k for k in range(33)]
with tempfile.TemporaryDirectory() as work:
    os.makedirs(os.path.join(work, "data"))
    os.makedirs(os.path.join(work, "output"))
    os.makedirs(os.path.join(work, "checks"))
    shutil.copy(os.path.join(HERE, "models", "yield_story.azr"), os.path.join(work, "m.azr"))
    with open(os.path.join(work, "data", "yield_curve.dat"), "w") as f:
        for e in E_KEV:
            f.write(f"{e:g} 0 1 1\n")
    with Session(os.path.join(work, "m.azr"), cwd=work) as m:
        model = np.concatenate([np.asarray(v, float) for v in [d.model for d in m.evaluate(np.asarray(m.params_rwa, float)).data]])
with open(os.path.join(HERE, "models", "data", "yield_curve.dat"), "w") as f:
    f.write("# E_lab(keV) unused Y(counts/uC) dY: synthetic, the model of yield_story.azr x (1 + 0.04 sin(2.3 i)), "
            "4 % errors (make_story.py)\n")
    for i, (e, y) in enumerate(zip(E_KEV, model)):
        v = float(f"{y * (1 + 0.04 * math.sin(2.3 * i)):.4g}")
        f.write(f"{e:g} 0 {v:.4g} {0.04 * v:.3g}\n")
