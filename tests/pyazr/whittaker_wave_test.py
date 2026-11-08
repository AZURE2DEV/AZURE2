#!/usr/bin/env python3
"""pyazr.whittaker_wave: the engine's GSL-free Whittaker function through _azure2 (docs/upgrade/coulomb.md 1.4).

Checks, against mpmath's whitw at 30 digits (independent of the engine's integral representation):
  (a) W and rho W'/W at a few (l, eta_b, rho), relative 1e-11 (physics_whittaker measures 2e-12 on its grid);
  (b) the scaled form at l = 0, eta_b = 300, rho = 1e-3, where W ~ 10^-680 is outside the double range: log10_W
      against mpmath to 1e-11 absolute, W = 0;
  (c) broadcasting: l scalar, eta_b (3, 1), rho (4,) give arrays of shape (3, 4), element by element as scalars;
  (d) ValueError for rho <= 0, eta_b < 0, non-integer l, naming the element.

Run from anywhere:  python3 tests/pyazr/whittaker_wave_test.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))


def main():
    try:
        import numpy as np
        import mpmath as mp
        from pyazr import whittaker_wave
    except Exception as err:  # same skip as the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 0
    mp.mp.dps = 30
    failures = []

    def check(name, ok, detail=""):
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
        if not ok:
            failures.append(name)

    worstW = worstD = 0.0
    for l, eta, rho in [(0, 0.0, 0.5), (0, 0.7, 2.0), (1, 2.3, 0.3), (2, 10.0, 5.0), (3, 0.05, 20.0), (5, 30.0, 1.0)]:
        w = whittaker_wave(l, eta, rho)
        W = mp.whitw(-eta, l + 0.5, 2 * rho)
        D = rho + eta - mp.whitw(1 - eta, l + 0.5, 2 * rho) / W   # DLMF 13.15.23
        worstW = max(worstW, abs(float(w["W"]) / float(W) - 1.0))
        worstD = max(worstD, abs(float(w["rho_dlog"]) - float(D)) / max(abs(float(D)), 1.0))
    check(f"(a) W against mpmath: {worstW:.1e}", worstW <= 1e-11)
    check(f"(a) rho W'/W against mpmath: {worstD:.1e}", worstD <= 1e-11)

    w = whittaker_wave(0, 300.0, 1e-3)
    ref = float(mp.log10(mp.whitw(-300, 0.5, 2e-3)))
    check(f"(b) outside the double range: log10_W {float(w['log10_W']):.6f} against {ref:.6f}, W = {float(w['W'])}",
          abs(float(w["log10_W"]) - ref) <= 1e-11 * max(1.0, abs(ref)) and float(w["W"]) == 0.0)

    etas, rhos = np.array([[0.0], [1.5], [7.0]]), np.array([0.2, 1.0, 3.0, 9.0])
    grid = whittaker_wave(2, etas, rhos)
    same = grid["W"].shape == (3, 4)
    for i in range(3):
        for j in range(4):
            s = whittaker_wave(2, float(etas[i, 0]), float(rhos[j]))
            same = same and grid["W"][i, j] == s["W"] and grid["rho_dlog"][i, j] == s["rho_dlog"]
    check("(c) broadcasting (3, 1) x (4,): shape (3, 4), the same values as element by element", same)

    for label, args in [("rho = 0", (0, 1.0, [1.0, 0.0])), ("eta_b < 0", (0, -0.5, 1.0)), ("l = 1.5", (1.5, 1.0, 1.0))]:
        try:
            whittaker_wave(*args)
            check(f"(d) ValueError for {label}", False, "no exception")
        except ValueError as err:
            check(f"(d) ValueError for {label}: {err}", "whittaker_wave" in str(err))
    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
