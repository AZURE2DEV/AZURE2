#!/usr/bin/env python3
"""Reference values of the Coulomb functions of a neutral pair (eta = 0) from mpmath, for tests/coulomb_neutral_test.cpp.

At eta = 0 the regular and irregular functions are Riccati-Bessel functions, F_l(rho) = rho j_l(rho) and
G_l(rho) = -rho y_l(rho), here sqrt(pi rho / 2) J_{l+1/2}(rho) and -sqrt(pi rho / 2) Y_{l+1/2}(rho) by mpmath's
besselj and bessely (independent of the Coulomb routines of physics/coulomb and of mpmath's own coulombf).  The
derivatives follow from l and l + 1 by Abramowitz and Stegun 14.2.1 at eta = 0, u'_l = (l + 1)/rho u_l - u_{l+1},
evaluated in the working precision of mpmath.

Grid: l = 0 ... 8; rho = 1e-6 ... 5e4 (small rho for each l, every half integer from 0.5 to 10 around rho ~ l, and
large rho up to 5e4); 41 values of rho, 369 rows.  Columns: l rho F F' G G'.

Usage: make_neutral_coulomb.py OUTFILE [DIGITS]  (default 40 digits; run twice with different DIGITS to estimate the
reference's own error, which is printed as the largest relative change).
"""
import sys

import mpmath as mp

LS = range(0, 9)
RHOS = (["1e-6", "1e-5", "1e-4", "7.5e-4", "3e-3", "0.01", "0.03", "0.1", "0.25"]
        + [str(k / 2) for k in range(1, 21)]
        + ["12", "15", "20", "30", "50", "100", "397", "1000", "5178", "1e4", "2e4", "5e4"])


def fmt(x):
    return mp.nstr(x, 20, min_fixed=0, max_fixed=0)


def rows():
    for rho_s in RHOS:
        rho = mp.mpf(rho_s)
        c = mp.sqrt(mp.pi * rho / 2)
        F = [c * mp.besselj(l + mp.mpf(1) / 2, rho) for l in range(0, 10)]
        G = [-c * mp.bessely(l + mp.mpf(1) / 2, rho) for l in range(0, 10)]
        for l in LS:
            s = mp.mpf(l + 1) / rho
            yield l, rho_s, [F[l], s * F[l] - F[l + 1], G[l], s * G[l] - G[l + 1]]


def main():
    out = sys.argv[1]
    mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    table = list(rows())
    with open(out, "w") as f:
        f.write(f"# mpmath {mp.__version__}, {mp.mp.dps} digits; tests/data/make_neutral_coulomb.py\n")
        f.write("# l rho F Fp G Gp  (eta = 0: F = rho j_l, G = -rho y_l)\n")
        for l, r, v in table:
            f.write(f"{l} {r} " + " ".join(fmt(x) for x in v) + "\n")
    # The reference's own error: the same table at 20 more digits.
    mp.mp.dps += 20
    worst = 0
    for (l, r, v), (_, _, w) in zip(table, rows()):
        for a, b in zip(v, w):
            worst = max(worst, abs(a - b) / abs(b))
    # Independent check of the closed forms: mpmath's coulombf/coulombg at eta = 0 on a subset.
    cw = 0
    for l, r, v in table:
        rho = mp.mpf(r)
        if mp.mpf("1e-4") <= rho <= 100:
            cw = max(cw, abs(mp.coulombf(l, 0, rho) - v[0]) / abs(v[0]), abs(mp.coulombg(l, 0, rho) - v[2]) / abs(v[2]))
    print(f"{len(table)} rows; largest relative change at {mp.mp.dps} digits {mp.nstr(worst, 3)}; "
          f"coulombf/coulombg on 1e-4 <= rho <= 100 {mp.nstr(cw, 3)}")


if __name__ == "__main__":
    main()
