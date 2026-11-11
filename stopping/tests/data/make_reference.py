#!/usr/bin/env python3
"""Reference values for the stopping module's tests (stopping/tests/*_test.cpp), written to reference.csv.

Independent of the C++ code: every value is computed here from the published formulas (and, for the nuclear
stopping, from the ZBL universal potential itself), with mpmath and scipy.  The coefficients are read from the
committed table stopping/data/electronic_coefficients.tsv.  Run from anywhere:

    python3 stopping/tests/data/make_reference.py      # rewrites stopping/tests/data/reference.csv

Rows: kind,<arguments...>,value.  See the tests for what each kind checks.
"""
import math
import os

import mpmath
from scipy import integrate, optimize

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = os.path.join(HERE, "..", "..", "data", "electronic_coefficients.tsv")
mpmath.mp.dps = 30

ALPHA = 1 / 137.035999084
HBARC = 197.3269804
ME = 0.51099895
U = 931.49410242


def rows():
    out = {}
    with open(TABLE) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Z\t") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            out[(int(c[0]), c[1])] = [float(x) for x in c[3:8]]
    return out


def fit_p(a, t_kev):
    scale = 1.0
    if t_kev < 10:
        scale, t_kev = math.sqrt(t_kev / 10), 10.0
    lo = a[1] * t_kev ** 0.45
    hi = a[2] / t_kev * math.log(1 + a[3] / t_kev + a[4] * t_kev)
    return scale * lo * hi / (lo + hi)


def fit_he(a, t_mev):
    scale = 1.0
    if t_mev < 1e-3:
        scale, t_mev = math.sqrt(t_mev / 1e-3), 1e-3
    lo = a[0] * (1000 * t_mev) ** a[1]
    hi = a[2] / t_mev * math.log(1 + a[3] / t_mev + a[4] * t_mev)
    return scale * lo * hi / (lo + hi)


def zbl_fit(eps):
    if eps > 30:
        return math.log(eps) / (2 * eps)
    return math.log(1 + 1.1383 * eps) / (2 * (eps + 0.01321 * eps ** 0.21226 + 0.19593 * eps ** 0.5))


def zbl_universal_integral(eps):
    """s_n(eps) = eps int_0^inf sin^2(theta/2) d(b^2) for the ZBL universal potential, classical scattering."""
    def phi(x):
        return (0.18175 * math.exp(-3.1998 * x) + 0.50986 * math.exp(-0.94229 * x)
                + 0.28022 * math.exp(-0.4029 * x) + 0.028171 * math.exp(-0.20162 * x))

    def theta(b):
        g = lambda x: 1 - phi(x) / (x * eps) - (b / x) ** 2
        x0 = optimize.brentq(g, 1e-12, 1e4 + b, xtol=1e-15, rtol=1e-14)
        # u = x0/x = 1 - s^2: theta = pi - 2 (b/x0) int_0^1 2 s ds / sqrt(g(x0/(1-s^2)))
        f = lambda s: 2 * s / math.sqrt(max(g(x0 / (1 - s * s)), 1e-300)) if s < 1 else 0.0
        val, _ = integrate.quad(f, 0, 1, epsabs=1e-13, epsrel=1e-11, limit=400)
        return math.pi - 2 * b / x0 * val

    f = lambda lnb: 2 * math.exp(2 * lnb) * math.sin(theta(math.exp(lnb)) / 2) ** 2
    val, _ = integrate.quad(f, math.log(1e-8), math.log(200.0), epsabs=1e-12, epsrel=1e-8, limit=400)
    return eps * val


BARKAS_F = [(0.02, 21.5), (0.03, 20.0), (0.04, 18.0), (0.05, 15.6), (0.06, 15.0), (0.07, 14.0), (0.08, 13.5),
            (0.09, 13.0), (0.1, 12.2), (0.2, 9.25), (0.3, 7.0), (0.4, 6.0), (0.5, 4.5), (0.6, 3.5), (0.7, 3.0),
            (0.8, 2.5), (0.9, 2.0), (1.0, 1.7), (1.2, 1.2), (1.3, 1.0), (1.4, 0.86), (1.5, 0.7), (1.6, 0.61),
            (1.7, 0.52), (1.8, 0.5), (1.9, 0.43), (2.0, 0.42), (2.1, 0.3), (2.4, 0.2), (3.0, 0.13), (3.08, 0.1),
            (3.1, 0.09), (3.3, 0.08), (3.5, 0.07), (3.8, 0.06), (4.0, 0.051), (4.1, 0.04), (4.8, 0.03),
            (5.0, 0.024), (5.1, 0.02), (6.0, 0.013), (6.5, 0.01), (7.0, 0.009), (7.1, 0.008), (8.0, 0.006),
            (9.0, 0.0032), (10.0, 0.0025)]


def barkas(z2, beta2):
    b = 1.8 if z2 == 1 else 0.6 if z2 == 2 else 1.8 if z2 <= 10 else 1.4 if z2 <= 17 else 1.8 if z2 == 18 \
        else 1.4 if z2 <= 25 else 1.35 if z2 <= 50 else 1.3
    x = beta2 / (ALPHA ** 2 * z2)
    w = b / math.sqrt(x)
    if w <= BARKAS_F[0][0]:
        f = BARKAS_F[0][1]
    elif w >= BARKAS_F[-1][0]:
        f = BARKAS_F[-1][1] * BARKAS_F[-1][0] / w
    else:
        k = next(i for i, (xx, _) in enumerate(BARKAS_F) if xx >= w)
        (x0, y0), (x1, y1) = BARKAS_F[k - 1], BARKAS_F[k]
        f = y0 + (y1 - y0) * (w - x0) / (x1 - x0)
    return f / (math.sqrt(z2) * x ** 1.5)


def shell(i, eta):
    e2 = eta ** -2
    return ((0.422377 * e2 + 0.0304043 * e2 ** 2 - 0.00038106 * e2 ** 3) * 1e-6 * i ** 2
            + (3.858019 * e2 - 0.1667989 * e2 ** 2 + 0.00157955 * e2 ** 3) * 1e-9 * i ** 3)


def bloch(y):
    return float(mpmath.digamma(1) - mpmath.re(mpmath.digamma(1 + 1j * y)))


def bethe(z1, z2, i, beta2):
    g2 = 1 / (1 - beta2)
    eta = math.sqrt(beta2 * g2)
    l0 = math.log(2 * ME * 1e6 * beta2 * g2 / i) - beta2
    ll = l0 - shell(i, max(eta, 0.13)) / z2 + z1 * barkas(z2, beta2) + z1 * z1 * bloch(z1 * ALPHA / math.sqrt(beta2))
    re = ALPHA * HBARC / ME * 1e-13
    k = 4 * math.pi * re * re * ME * 1e6 * 1e15
    return k * z1 * z1 * z2 / beta2 * ll


def main():
    tab = rows()
    out = ["# kind,arguments...,value  (written by make_reference.py; do not edit)"]
    for z in (1, 6, 13, 47, 79):
        for t in (1.0, 5.0, 10.0, 30.0, 100.0, 500.0, 1000.0):
            out.append(f"fit_p,{z},{t!r},{fit_p(tab[(z, 'p')], t)!r}")
    for z in (2, 8, 22, 79):
        for t in (5e-4, 1e-3, 0.01, 0.1, 1.0, 4.0):
            out.append(f"fit_he,{z},{t!r},{fit_he(tab[(z, 'He')], t)!r}")
    for y in (0.01, 0.1, 0.5, 1.0, 3.0):
        out.append(f"bloch,{y!r},{bloch(y)!r}")
    for i, eta in ((166.0, 0.1463), (790.0, 0.2), (19.2, 0.13), (718.0, 0.5)):
        out.append(f"shell,{i!r},{eta!r},{shell(i, eta)!r}")
    for z1, z2, i, b2 in ((1, 13, 166.0, 0.006376), (1, 79, 790.0, 0.02096), (2, 8, 95.0, 0.004), (1, 1, 19.2, 0.01)):
        out.append(f"bethe,{z1},{z2},{i!r},{b2!r},{bethe(z1, z2, i, b2)!r}")
    for eps in (0.003, 0.03, 0.3, 3.0, 30.0):
        out.append(f"zbl_fit,{eps!r},{zbl_fit(eps)!r}")
        out.append(f"zbl_universal,{eps!r},{zbl_universal_integral(eps)!r}")
    with open(os.path.join(HERE, "reference.csv"), "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
