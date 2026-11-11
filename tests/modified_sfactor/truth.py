#!/usr/bin/env python3
"""The modified S factor S*(E) = E sigma(E) exp(2 pi eta + g E) from the engine's own cross section (mpmath, no
AZURE2 code).

eta = Z1 Z2 alpha sqrt(mu c^2 / (2 E)), non-relativistic, mu = m1 m2 / (m1 + m2) of the atomic masses the project's
<levels> print for the entrance pair (p + 15N: 1.00728 u, 14.99627 u, Z = 1, 7), alpha and the atomic mass unit as in
physics/include/azure2/physics/Constants.h (CODATA 2018), E the c.m. energy of the output row (MeV).

  truth.py check <AZUREOut .extrap of models/pa_extrap.azr>
      blocks: 1 sigma, 2 S, 3 S* (g = 0.46 /MeV), 4 S* (g = 0).  Prints the largest relative difference of S* (block
      3, column 5) from sigma (block 1, column 4) x E exp(2 pi eta + g E), of block 3's sigma column from block 1's,
      and whether block 4 equals block 2 character for character; exit 1 beyond 1e-10 or when the text differs.
  truth.py data <AZUREOut .extrap>
      writes models/data/sigma.dat and models/data/sstar.dat: block 1's sigma times 1 + 0.03 sin(3.7 i), 5 % errors,
      every second energy, and the same rows as S* with g = 0.46 /MeV (value and error times E exp(2 pi eta + g E)).
"""
import os
import sys

import mpmath as mp

mp.mp.dps = 40
U_MEV = mp.mpf("931.49410242")
ALPHA = 1 / mp.mpf("137.035999084")
Z1, Z2, M1, M2 = 1, 7, mp.mpf("1.00728"), mp.mpf("14.99627")
G = mp.mpf("0.46")
TOL = mp.mpf("1e-10")


def conversion(e_mev, g):
    e = mp.mpf(e_mev)
    mu = M1 * M2 / (M1 + M2)
    eta = Z1 * Z2 * ALPHA * mp.sqrt(mu * U_MEV / (2 * e))
    return e * mp.exp(2 * mp.pi * eta + g * e)


def blocks(path):
    out, cur = [], []
    for line in open(path):
        if line.strip():
            cur.append(line.rstrip("\n"))
        elif cur:
            out.append(cur)
            cur = []
    if cur:
        out.append(cur)
    return out


def check(path):
    b = blocks(path)
    if len(b) != 4 or len({len(x) for x in b}) != 1:
        print("expected 4 blocks of equal length, got", [len(x) for x in b])
        return 1
    worst_s, worst_sigma = mp.mpf(0), mp.mpf(0)
    for sig_row, s_row in zip(b[0], b[2]):
        e, sigma = sig_row.split()[0], mp.mpf(sig_row.split()[3])
        f = s_row.split()
        if f[0] != e:
            print("grids differ:", e, f[0])
            return 1
        want = sigma * conversion(e, G)
        worst_s = max(worst_s, abs(mp.mpf(f[4]) / want - 1))
        worst_sigma = max(worst_sigma, abs(mp.mpf(f[3]) / sigma - 1))
    same = b[3] == b[1]
    print("S* against sigma E exp(2 pi eta + g E): %.2e; its sigma column: %.2e; g = 0 equals S: %s"
          % (worst_s, worst_sigma, "yes" if same else "no"))
    return 0 if worst_s <= TOL and worst_sigma <= TOL and same else 1


def data(path):
    b = blocks(path)
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models", "data")
    with open(os.path.join(here, "sigma.dat"), "w") as fs, open(os.path.join(here, "sstar.dat"), "w") as fx:
        fs.write("# E_cm [MeV]   (unused)   sigma [b]   error [b]  (truth.py data: the model times 1 + 0.03 sin(3.7 i), 5 % errors)\n")
        fx.write("# E_cm [MeV]   (unused)   S* [MeV b]   error [MeV b]  (the rows of sigma.dat times E exp(2 pi eta + 0.46 E))\n")
        for i, row in enumerate(b[0][::2]):
            e, sigma = row.split()[0], mp.mpf(row.split()[3])
            v = sigma * (1 + mp.mpf("0.03") * mp.sin(mp.mpf("3.7") * i))
            dv = v * mp.mpf("0.05")
            c = conversion(e, G)
            fs.write("%s  0  %s  %s\n" % (e, mp.nstr(v, 17, min_fixed=0, max_fixed=0), mp.nstr(dv, 17, min_fixed=0, max_fixed=0)))
            fx.write("%s  0  %s  %s\n" % (e, mp.nstr(v * c, 17, min_fixed=0, max_fixed=0),
                                      mp.nstr(dv * c, 17, min_fixed=0, max_fixed=0)))
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("check", "data"):
        print(__doc__)
        sys.exit(2)
    sys.exit(check(sys.argv[2]) if sys.argv[1] == "check" else data(sys.argv[2]))
