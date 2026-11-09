#!/usr/bin/env python3
"""Reference values of the Whittaker function W_l(rho) = W_{-eta_b, l+1/2}(2 rho) for physics_whittaker.

mpmath (40 digits) is the truth: log10 W and rho W'/W (primes d/drho) on l x eta_b x rho.  The derivative uses
DLMF 13.15.23, z dW_{k,m}/dz = (z/2 - k) W_{k,m}(z) - W_{k+1,m}(z), i.e. rho W'/W = rho + eta_b - W_{1-eta_b,
l+1/2}(2rho) / W_{-eta_b, l+1/2}(2rho); it is checked here against mpmath's numerical derivative before anything is
written.  log10 W rather than W, so that values far outside the double range (eta_b = 1000, rho = 1e-4) are kept.

The grid leaves out eta_b rho > 2000.  eta_b rho = Z1 Z2 alpha mu c^2 a / (hbar c) is fixed by the pair and the
channel radius a (about 100 for Z1 Z2 = 50, mu = 10 u, a = 6 fm), and mpmath's hyperu does not converge at the far
corner (l = 0, eta_b = 1000, rho = 200).

Run:  python3 make_whittaker_reference.py > whittaker.csv     (mpmath 1.4.1; about a minute)
"""
import sys

import mpmath as mp

mp.mp.dps = 40

LS = [0, 1, 2, 3, 5, 8, 10]
ETAS = ["0", "0.01", "0.1", "0.5", "1", "2", "2.5", "3", "5", "10", "30", "100", "300", "1000"]
RHOS = ["1e-4", "1e-3", "1e-2", "0.1", "0.5", "1", "2", "5", "10", "20", "50", "100", "200"]


def values(l, eta, rho):
    w = mp.whitw(-eta, l + mp.mpf(1) / 2, 2 * rho)
    w1 = mp.whitw(1 - eta, l + mp.mpf(1) / 2, 2 * rho)
    return mp.log10(w), rho + eta - w1 / w


def main():
    for l, eta, rho in ((0, mp.mpf("0.3"), mp.mpf("1.2")), (3, mp.mpf(5), mp.mpf("0.5")), (8, mp.mpf(30), mp.mpf(10))):
        _, d = values(l, eta, rho)
        num = rho * mp.diff(lambda r: mp.whitw(-eta, l + mp.mpf(1) / 2, 2 * r), rho) / \
            mp.whitw(-eta, l + mp.mpf(1) / 2, 2 * rho)
        if abs(d - num) > mp.mpf("1e-25") * max(1, abs(num)):
            sys.exit(f"DLMF 13.15.23 check failed at l={l} eta={eta} rho={rho}: {d} vs {num}")
    print("# l,eta_b,rho,log10_W,rho_dlog  (mpmath %s, %d digits; make_whittaker_reference.py)" % (mp.__version__, mp.mp.dps))
    for l in LS:
        for e in ETAS:
            for r in RHOS:
                if mp.mpf(e) * mp.mpf(r) > 2000:
                    continue
                lw, d = values(l, mp.mpf(e), mp.mpf(r))
                print(f"{l},{e},{r},{mp.nstr(lw, 22)},{mp.nstr(d, 22)}")


if __name__ == "__main__":
    main()
