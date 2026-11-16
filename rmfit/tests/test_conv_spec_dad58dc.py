#!/usr/bin/env python3
"""build_conv_spec follows AZURE2 after dad58dc (2026-10-11): one normalised Gaussian per data point, sigma evaluated at
the point -- constant (lab, x cm_conv), `sigmalab` expression (at the point's LAB energy, x cm_conv), legacy expression
(at the point's c.m. energy).  General expressions (a0..aN, x, ^, sqrt).  No engine.
Run:  python3 rmfit/tests/test_conv_spec_dad58dc.py"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.model import build_conv_spec, parse_target_effect  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def fold(spec, f):
    return spec.W @ f(spec.grid_cm)


def main():
    rows = np.array([[7.0, 0, 1.0, 0.1], [8.0, 0, 1.0, 0.1], [9.0, 0, 1.0, 0.1]])
    cmc = 16.0 / 17.0
    const = parse_target_effect('1  "170"  20  1  0.01443  0  0  ""  0  0  0  0  ""  0  0  0.04  20  50')
    lab = parse_target_effect('1  "170"  40  1  0.01443  0  0  ""  0  0  0  1  "sqrt(a0+a1*x^3)"  2  2.083333e-04  5.140511e-06  0  0.04  20  50  sigmalab')
    legacy = parse_target_effect('1  "170"  40  1  0.01443  0  0  ""  0  0  0  1  "sqrt(a0+a1*x^3)"  2  2.083333e-04  5.140511e-06  0  0.04  20  50')
    check("sigmalab keyword parsed (and absent on the legacy line)", lab.sigma_lab and not legacy.sigma_lab and not const.sigma_lab)
    check("general expression: sqrt(a0+a1*x^3) at 8 MeV", np.isclose(lab.sigma_eq(8.0), np.sqrt(2.083333e-4 + 5.140511e-6 * 512)),
          lab.sigma_eq(8.0))
    for name, te in (("constant", const), ("sigmalab", lab), ("legacy", legacy)):
        sp = build_conv_spec(None, te, rows, cmc)
        check(f"{name}: kernel normalised (a flat cross section folds to itself)", np.allclose(fold(sp, np.ones_like), 1.0, atol=1e-12),
              fold(sp, np.ones_like))
        check(f"{name}: a linear cross section folds to its value at the point", np.allclose(fold(sp, lambda g: 3 * g), 3 * rows[:, 0] * cmc, rtol=0, atol=3 * 1e-4))   # discrete +-3 sigma cut on a sigma/10 grid: <= 0.06 keV here
    # widths: second moment of each row's kernel
    def rms(sp):
        e_cm = rows[:, 0] * cmc
        return np.sqrt(np.array([np.sum(sp.W[i] * (sp.grid_cm - e_cm[i]) ** 2) for i in range(3)]))
    trunc = np.sqrt(1 - 2 * 3 * np.exp(-4.5) / np.sqrt(2 * np.pi) / 0.9973002)    # rms of a +-3 sigma truncated Gaussian / sigma
    s_lab = np.sqrt(2.083333e-4 + 5.140511e-6 * rows[:, 0] ** 3)
    check("sigmalab: sigma = sigma_lab(E_lab of the point) x cm_conv", np.allclose(rms(build_conv_spec(None, lab, rows, cmc)), s_lab * cmc * trunc, rtol=2e-3),
          (rms(build_conv_spec(None, lab, rows, cmc)), s_lab * cmc * trunc))
    s_leg = np.sqrt(2.083333e-4 + 5.140511e-6 * (rows[:, 0] * cmc) ** 3)
    check("legacy: sigma = expression at the c.m. energy of the point", np.allclose(rms(build_conv_spec(None, legacy, rows, cmc)), s_leg * trunc, rtol=2e-3))
    check("constant: sigma x cm_conv", np.allclose(rms(build_conv_spec(None, const, rows, cmc)), 0.01443 * cmc * trunc, rtol=2e-3))
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all conv-spec (dad58dc) checks passed")


if __name__ == "__main__":
    main()
