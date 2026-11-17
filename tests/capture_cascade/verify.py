#!/usr/bin/env python3
"""Compare the capture-cascade test runs in WORK with expected/ (called by check.sh).

Every model holds four segments of one capture pair: the cascade (UPOS) distribution, the angle-integrated cross
section, the primary-photon distribution and an angle-integrated UPOS segment.  W = 4 pi (dsigma/dOmega) / sigma is
compared with expected/<model>.dat (the secondary gamma ray) and expected/<model>_prim.dat (the primary, where there is
one).  Exit status = number of failures.
"""
import glob
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = {"single_e1": 1, "multipoles": 1, "twoj_plus": 1, "twoj_minus": 1, "ec_s": 1, "ec_p": 1, "ecres_both": 3}
TOL = 1e-7
N = 13


def rows(work, model):
    d = np.loadtxt(glob.glob(os.path.join(work, model, "output", "AZUREOut_*.out"))[0], ndmin=2)
    return d[d[:, 0] > 0]


def main(work):
    failures = 0
    for model, ne in MODELS.items():
        try:
            d = rows(work, model)
        except (OSError, IndexError):
            print(f"  FAIL  {model}: no output")
            failures += 1
            continue
        upos, sigma = d[:ne * N, 3], d[ne * N:ne * N + ne, 3]
        prim, iupos = d[ne * N + ne:2 * ne * N + ne, 3], d[2 * ne * N + ne:, 3]
        w = 4 * np.pi * upos / np.repeat(sigma, N)
        exp = np.loadtxt(os.path.join(HERE, "expected", model + ".dat"))
        dev = np.max(np.abs(w - exp[:, 1]))
        span = f"{w.min():.4f}..{w.max():.4f}"
        if dev < TOL:
            print(f"  ok    {model}: secondary W(theta) = substate calculation, max deviation {dev:.1e} (W {span})")
        else:
            print(f"  FAIL  {model}: secondary W(theta) deviates by {dev:.3e} from the substate calculation (W {span})")
            failures += 1
        ppath = os.path.join(HERE, "expected", model + "_prim.dat")
        if os.path.exists(ppath):
            devp = np.max(np.abs(4 * np.pi * prim / np.repeat(sigma, N) - np.loadtxt(ppath)[:, 1]))
            if devp < TOL:
                print(f"  ok    {model}: primary W(theta) = substate calculation, max deviation {devp:.1e}")
            else:
                print(f"  FAIL  {model}: primary W(theta) deviates by {devp:.3e}")
                failures += 1
        devi = np.max(np.abs(iupos / sigma - 1.0))
        if len(iupos) != ne or devi > 1e-12:
            print(f"  FAIL  {model}: angle-integrated UPOS segment differs from the capture cross section ({devi:.1e})")
            failures += 1
    if failures == 0:
        print("  ok    angle-integrated UPOS segments = the capture cross section into B, every model")
    return failures


if __name__ == "__main__":
    sys.exit(int(main(sys.argv[1])))
