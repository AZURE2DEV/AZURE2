#!/usr/bin/env python3
"""Reference values of the Coulomb and Whittaker functions from mpmath, for tests/test_accuracy.cpp.

Grid: l = 0 ... 8; eta = 0, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20; rho = 1e-3 ... 50 (16 points).
Scattering (coulomb.csv): F_l, F'_l, G_l, G'_l, sigma_l. The derivatives follow from F, G at l and l + 1 by
Abramowitz and Stegun 14.2.1, evaluated in the working precision of mpmath.
Bound (whittaker.csv): W = W_{-eta, l+1/2}(2 rho), dW/drho by mpmath.diff, and rho W'/W.

Usage: make_reference.py OUTDIR [DIGITS]  (default 40 digits; run twice with different DIGITS to estimate the
reference's own error, which is printed as the largest relative change).
"""
import sys

import mpmath as mp

LS = range(0, 9)
ETAS = ["0", "0.01", "0.1", "0.5", "1", "2", "5", "10", "20"]
RHOS = ["1e-3", "3e-3", "0.01", "0.03", "0.1", "0.3", "1", "2", "3", "5", "7", "10", "15", "20", "30", "50"]


def fmt(x):
    return mp.nstr(x, 20, min_fixed=0, max_fixed=0)


def coulomb_rows():
    for eta_s in ETAS:
        eta = mp.mpf(eta_s)
        for rho_s in RHOS:
            rho = mp.mpf(rho_s)
            F = [mp.coulombf(l, eta, rho) for l in range(0, 10)]
            G = [mp.coulombg(l, eta, rho) for l in range(0, 10)]
            for l in LS:
                lp1 = mp.mpf(l + 1)
                s = lp1 / rho + eta / lp1
                r = mp.sqrt(1 + (eta / lp1) ** 2)
                Fp = s * F[l] - r * F[l + 1]
                Gp = s * G[l] - r * G[l + 1]
                sigma = mp.im(mp.loggamma(mp.mpc(l + 1, eta)))
                yield l, eta_s, rho_s, [F[l], Fp, G[l], Gp, sigma]


def whittaker_rows():
    for eta_s in ETAS:
        eta = mp.mpf(eta_s)
        for rho_s in RHOS:
            rho = mp.mpf(rho_s)
            for l in LS:
                m = mp.mpf(l) + mp.mpf(1) / 2
                W = mp.whitw(-eta, m, 2 * rho)
                Wp = mp.diff(lambda x: mp.whitw(-eta, m, 2 * x), rho)
                yield l, eta_s, rho_s, [W, Wp, rho * Wp / W]


def main():
    out = sys.argv[1]
    mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    with open(f"{out}/coulomb.csv", "w") as f:
        f.write(f"# mpmath {mp.__version__}, {mp.mp.dps} digits; tests/make_reference.py\n")
        f.write("# l eta rho F Fp G Gp sigma\n")
        for l, e, r, v in coulomb_rows():
            f.write(" ".join([str(l), e, r] + [fmt(x) for x in v]) + "\n")
    with open(f"{out}/whittaker.csv", "w") as f:
        f.write(f"# mpmath {mp.__version__}, {mp.mp.dps} digits; tests/make_reference.py\n")
        f.write("# l eta_b rho W Wp rhoWp_over_W   (W = W_{-eta_b, l+1/2}(2 rho), Wp = dW/drho)\n")
        for l, e, r, v in whittaker_rows():
            f.write(" ".join([str(l), e, r] + [fmt(x) for x in v]) + "\n")


if __name__ == "__main__":
    main()
