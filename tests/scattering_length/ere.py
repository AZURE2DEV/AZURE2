#!/usr/bin/env python3
"""s-wave scattering length from AZURE2's own phase shifts, by the Coulomb-modified effective-range
expansion: K(E) = C0^2 k cot(delta) + 2 k eta h(eta) = -1/a0 + r0 k^2/2 + ...  (K = k cot(delta) for a
neutral pair), fitted as a polynomial in k^2 and read at k = 0.  Independent of the threshold code
(ScatteringLength.cpp): it uses only the phase shifts AZURE2 writes for a phase-shift test segment.

    ere.py FILE BLOCK Z1Z2 M1 M2 [DEGREE]

FILE: AZUREOut_*.extrap (columns E_cm, E_x, angle, phase in degrees), one block per test segment,
blocks separated by blank lines; BLOCK is 1-based.  Prints a0 in fm.  Standard library only."""
import math
import sys

HBARC, UCONV, ALPHA = 197.32696310, 931.4940880, 1.0 / 137.0359996790


def h_function(eta):
    """h(eta) = Re psi(i eta) - ln(eta), from Re psi(i eta) = -gamma + eta^2 sum 1/(n (n^2 + eta^2))."""
    n_max = 200000
    s = sum(1.0 / (n * (n * n + eta * eta)) for n in range(1, n_max + 1))
    s += 1.0 / (2.0 * n_max * n_max)  # tail of the sum, ~ int dn / n^3
    return -0.57721566490153286 + eta * eta * s - math.log(eta)


def fit_intercept(xs, ys, degree):
    """Least-squares polynomial of the given degree; returns the constant term."""
    n = degree + 1
    ata = [[sum(x ** (i + j) for x in xs) for j in range(n)] for i in range(n)]
    aty = [sum(y * x ** i for x, y in zip(xs, ys)) for i in range(n)]
    for c in range(n):  # Gauss-Jordan with partial pivoting
        p = max(range(c, n), key=lambda r: abs(ata[r][c]))
        ata[c], ata[p], aty[c], aty[p] = ata[p], ata[c], aty[p], aty[c]
        for r in range(n):
            if r != c:
                f = ata[r][c] / ata[c][c]
                ata[r] = [a - f * b for a, b in zip(ata[r], ata[c])]
                aty[r] -= f * aty[c]
    return aty[0] / ata[0][0]


def main():
    path, block, z1z2, m1, m2 = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), \
        float(sys.argv[5])
    degree = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    blocks, current = [], []
    for line in open(path):
        if line.strip():
            current.append(line.split())
        elif current:
            blocks.append(current)
            current = []
    if current:
        blocks.append(current)
    rows = blocks[block - 1]
    mu = m1 * m2 / (m1 + m2) * UCONV
    xs, ys = [], []
    for r in rows:
        e, delta = float(r[0]), math.radians(float(r[3]))
        k = math.sqrt(2.0 * mu * e) / HBARC
        if z1z2 == 0:
            K = k / math.tan(delta)
        else:
            eta = z1z2 * ALPHA * math.sqrt(mu / (2.0 * e))
            two_pi_eta = 2.0 * math.pi * eta
            c02 = two_pi_eta / math.expm1(two_pi_eta)
            K = c02 * k / math.tan(delta) + 2.0 * k * eta * h_function(eta)
        xs.append(k * k)
        ys.append(K)
    print("%.10e" % (-1.0 / fit_intercept(xs, ys, degree)))


if __name__ == "__main__":
    main()
