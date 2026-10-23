#!/usr/bin/env python3
"""Compare the secondary-gamma test runs in WORK with expected/ (called by check.sh).

Single-level models: W = 4 pi (dsigma/dOmega_gamma) / sigma from the UPOS and angle-integrated segments of one run,
against expected/<model>.dat (truth.py).  twolevel: the UPOS cross sections against expected/twolevel_pre_fix.out.
Exit status = number of failures.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SINGLE = {"spin2_plus": 5, "spin2_minus": 5, "spin2_zero": 5, "spin2_multil": 5, "mixing": 6}
TOL = 1e-7


def rows(work, model, key):
    """Output rows of the run, in segment order: E_cm, E_x, angle, cross section, ..."""
    path = os.path.join(work, model, "output", f"AZUREOut_aa=2_R={key}.out")
    d = np.loadtxt(path, ndmin=2)
    return d[d[:, 0] > 0]


def main(work):
    failures = 0
    for model, key in SINGLE.items():
        try:
            d = rows(work, model, key)
        except OSError:
            print(f"  FAIL  {model}: no output")
            failures += 1
            continue
        n = 13
        upos, sigma = d[:n, 3], d[n, 3]
        w = 4 * np.pi * upos / sigma
        exp = np.loadtxt(os.path.join(HERE, "expected", model + ".dat"))
        dev = np.max(np.abs(w - exp[:, 1]) / np.abs(exp[:, 1]))
        if np.allclose(d[:n, 2], exp[:, 0]) and dev < TOL:
            print(f"  ok    {model}: W(theta_gamma) = substate calculation, max rel. deviation {dev:.1e}")
        else:
            print(f"  FAIL  {model}: max rel. deviation {dev:.3e} from the substate calculation "
                  f"(0/90 deg {w[0] / w[6]:.4f}, expected {exp[0, 1] / exp[6, 1]:.4f})")
            failures += 1
    try:
        d = rows(work, "twolevel", 6)
        ref = np.loadtxt(os.path.join(HERE, "expected", "twolevel_pre_fix.out"), ndmin=2)
        ref = ref[ref[:, 0] > 0]
        dev = np.max(np.abs(d[:, 3] - ref[:, 3]) / np.abs(ref[:, 3]))
        if d.shape == ref.shape and dev < 1e-9:
            print(f"  ok    twolevel: one exit channel spin, unchanged from the pre-fix engine (max rel. {dev:.1e})")
        else:
            print(f"  FAIL  twolevel: differs from the pre-fix engine, max rel. {dev:.3e}")
            failures += 1
    except OSError:
        print("  FAIL  twolevel: no output")
        failures += 1
    return failures


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
