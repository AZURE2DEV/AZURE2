#!/usr/bin/env python3
"""Check (g) of tests/shift_sqrt_subpoints: the analytic plateau of the narrow resonance under narrow_sqrt's shift.

    python3 -I plateau.py NARROW_OUTPUT NORES_OUTPUT [TOLERANCE]

NARROW_OUTPUT is the menu-1 output of models/narrow_sqrt.azr, NORES_OUTPUT that of models/narrow_sqrt_nores.azr (the
same model with the 1/2+ level's proton width at zero: the non-resonant part).  Each of the four yields (a thin
target's quantity sigma, the yield over 1e-24 N) minus its non-resonant part is compared with

    sigma_r = (lambda_r^2 / 2) omega Gamma_gamma (m + M) / (M eps(E_r,lab) N) x F,
    F = [arctan(2 (E - E_r)/Gamma) - arctan(2 (E_back - E_r)/Gamma)] / pi,

E the point's true (shifted) c.m. energy, E_back = E - f N eps(E_r) its back edge (the arctan of the back edge is 1e-3
of the front one's, so its stopping is not refined), Gamma = Gamma_p.  The strength is omega Gamma_gamma Gamma_p /
(sum of the particle widths) = omega Gamma_gamma: the radiative widths do not enter the level matrix
(docs/source/theory/conventions.rst), so neither the total width of the Lorentzian nor the denominator of the
strength holds Gamma_gamma.  The level (models/narrow_sqrt.azr, <levels>): E_x 7.5567 MeV, S_p 7.297 MeV, J = 1/2,
p + 14N (j = 1/2, 1) with masses 1.00783 u and 14.0031 u, Gamma_p 9.74937 eV, five gamma widths summing to 38.5027
meV; the target: N = 3.2e19 atoms/cm^2 and the stopping formula of its <targets> record.  Exit 1 when a ratio is off
1 by more than TOLERANCE.
"""
import glob
import math
import sys

HBARC, AMU = 197.3269804, 931.49410242  # MeV fm, MeV/u (physics/Constants.h)
M_P, M_N14 = 1.00783, 14.0031
E_R = 7.5567 - 7.297                    # MeV c.m.
GAMMA_P = 9.74937e-6                    # MeV
GAMMA_G = (0.5787 + 8.87 + 6.5 + 22.5 + 0.054) * 1e-9  # MeV
OMEGA = 2.0 / (2.0 * 3.0)               # (2J+1)/((2j_p+1)(2j_N+1))
N = 3.2e19                              # atoms/cm^2
A = [0.146403, 0.00484558, 5.54014, 0.455903, 1763.46, 0.919014, 3951.93, 0.0550341,
     3.31007, 0.00495744, 0.540621, 0.77994, 1104.15, 0.967848, 2235.35, 0.0531612]


def eps(x):
    """The target's stopping formula, MeV cm^2/atom at lab energy x (MeV)."""
    def part(o):
        lo = A[o] * x ** A[o + 1] + A[o + 2] * x ** A[o + 3]
        hi = A[o + 4] / x ** A[o + 5] * math.log(A[o + 6] / x + A[o + 7] * x)
        return lo * hi / (lo + hi)
    return (0.5 * part(0) + 0.5 * part(8)) / 1e21


def values(out):
    files = glob.glob(out + "/AZUREOut_aa=1_*.out")
    return [(float(l.split()[0]), float(l.split()[3])) for l in open(files[0]) if len(l.split()) >= 5]


def main():
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 5e-4
    f = M_N14 / (M_P + M_N14)
    mu = M_P * M_N14 / (M_P + M_N14) * AMU
    lam2_half = 2.0 * math.pi ** 2 * HBARC ** 2 / (2.0 * mu * E_R) / 100.0  # b
    eps_r = eps(E_R / f)
    plateau = lam2_half * OMEGA * GAMMA_G / f / (N * eps_r)  # b, the thin target's sigma = Y / (1e-24 N)
    worst, bad = 0.0, 0
    for (e, y), (e2, y0) in zip(values(sys.argv[1]), values(sys.argv[2])):
        back = e - f * N * eps_r
        F = (math.atan(2.0 * (e - E_R) / GAMMA_P) - math.atan(2.0 * (back - E_R) / GAMMA_P)) / math.pi
        r = (y - y0) / (plateau * F)
        worst = max(worst, abs(r - 1.0))
        bad += abs(r - 1.0) > tol or abs(e - e2) > 1e-12
        print("E_cm %.6f MeV: engine %.10e, non-resonant %.4e, analytic %.10e (arctan %.6f), ratio %.6f" %
              (e, y, y0, plateau * F, F, r))
    print("largest |ratio - 1| %.2e, tolerance %.1e" % (worst, tol))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
