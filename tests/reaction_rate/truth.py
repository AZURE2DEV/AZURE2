#!/usr/bin/env python3
"""N_A <sigma v> from a dense mode-3 extrapolation of sigma(E), independently of ReactionRate.cpp.

    python3 truth.py <m1> <m2> <T9 file> <AZUREOut_*.extrap ...>  >  rates.dat

Reads every .extrap (columns E_cm, E_x, angle, sigma [b], S), merges and sorts the points
(duplicates within 1e-9 dropped), keeps sigma >= 0 and finite, and integrates
sigma E exp(-E/kT) by the trapezoid rule on the merged grid, which is dense enough
(1 % of E per step in every decade, 4 eV across a narrow resonance) for 1e-5.
CODATA 2018 constants; masses in u as in the .azr.
"""
import math
import sys

NA = 6.02214076e23
C_CM_S = 2.99792458e10
KB = 8.617333262e-2      # MeV / GK
U_MEV = 931.49410242


def main():
    m1, m2 = float(sys.argv[1]), float(sys.argv[2])
    temps = [float(t) for t in open(sys.argv[3]).read().split()]
    pts = []
    for f in sys.argv[4:]:
        for line in open(f):
            t = line.split()
            if len(t) < 4:
                continue
            try:
                e, s = float(t[0]), float(t[3])
            except ValueError:
                continue
            if not (e > 0.0) or not math.isfinite(s):
                continue
            pts.append((e, max(s, 0.0)))
    pts.sort()
    grid = []
    for e, s in pts:
        if grid and abs(e - grid[-1][0]) <= 1e-9 * e:
            continue
        grid.append((e, s))
    mu = m1 * m2 / (m1 + m2) * U_MEV
    print("%12s %18s" % ("T9", "Rate"))
    for T in temps:
        kT = KB * T
        total = 0.0
        for (e0, s0), (e1, s1) in zip(grid[:-1], grid[1:]):
            f0 = s0 * e0 * math.exp(-e0 / kT)
            f1 = s1 * e1 * math.exp(-e1 / kT)
            total += 0.5 * (f0 + f1) * (e1 - e0)
        rate = 1e-24 * NA * C_CM_S * math.sqrt(8.0 / math.pi / mu) / kT ** 1.5 * total
        print("%12g %18.8e" % (T, rate))


if __name__ == "__main__":
    main()
