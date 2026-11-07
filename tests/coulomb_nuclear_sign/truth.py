"""Independent elastic cross section for tests/coulomb_nuclear_sign (run once; expected/dsdo.dat is committed).

p + 12C hard-sphere scattering (every partial wave l <= 4, no resonant width) in the Coulomb field, computed from
the textbook partial-wave sum with mpmath Coulomb functions, without any AZURE2 convention:

    f(theta) = f_C(theta) + (1 / 2ik) sum_l (2l+1) exp(2i w_l) [exp(-2i phi_l) - 1] P_l(cos theta),
    f_C(theta) = -eta / (2k sin^2(theta/2)) exp(-i eta ln sin^2(theta/2)),
    w_l = sigma_l - sigma_0,  phi_l = arg(G_l + i F_l) at rho = k a,

both amplitudes taken relative to the common factor exp(2i sigma_0).  For a spin-1/2 projectile on a spin-0
target with a J-independent collision matrix the spin-flip amplitude vanishes and dsigma/dOmega = |f|^2.  The
interference term 2 Re[f_C^* f_N] is what the test pins: with the opposite sign it changes dsigma/dOmega by
7-74 % on this grid.

Usage:  python3 truth.py AZUREOut_aa=1_R=1.extrap > expected/dsdo.dat
The c.m. energies and angles are read from columns 1 and 3 of an AZURE2 run of hard_sphere_pc.azr (kinematics
are not under test); masses, radius and constants are those of the model and of include/Constants.h.
"""
import sys

import mpmath as mp

mp.mp.dps = 30
HBARC = mp.mpf("197.3269804")  # MeV fm
UCONV = mp.mpf("931.49410242")  # MeV / u
ALPHA = 1 / mp.mpf("137.035999084")
M1, M2 = mp.mpf("1.00782"), mp.mpf("12.0")  # u, as in hard_sphere_pc.azr
Z1Z2 = 6
A = mp.mpf("5.0")  # fm
LMAX = 4


def dsdo(e_cm, theta_deg):
    mu = M1 * M2 / (M1 + M2) * UCONV  # MeV
    k = mp.sqrt(2 * mu * e_cm) / HBARC  # fm^-1
    eta = Z1Z2 * ALPHA * mp.sqrt(mu / (2 * e_cm))
    rho = k * A
    th = mp.radians(theta_deg)
    s2 = mp.sin(th / 2) ** 2
    fc = -eta / (2 * k * s2) * mp.expj(-eta * mp.log(s2))
    fn = mp.mpc(0)
    w = mp.mpf(0)
    for l in range(LMAX + 1):
        if l > 0:
            w += mp.atan(eta / l)
        F = mp.coulombf(l, eta, rho)
        G = mp.coulombg(l, eta, rho)
        phi = mp.atan2(F, G)
        fn += (2 * l + 1) * mp.expj(2 * w) * (mp.expj(-2 * phi) - 1) * mp.legendre(l, mp.cos(th))
    fn /= 2j * k
    return abs(fc + fn) ** 2 * 10  # fm^2 -> mb; /1000 below gives b


def main(path):
    print("# E_cm(MeV) theta_cm(deg) dsigma/dOmega_cm(b/sr)  truth.py, mpmath %s" % mp.__version__)
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 4 or p[0].startswith("#"):
                continue
            e, th = mp.mpf(p[0]), mp.mpf(p[2])
            if th <= 0:
                continue
            print("%.10e %.10e %.10e" % (float(e), float(th), float(dsdo(e, th) / 1000)))


if __name__ == "__main__":
    main(sys.argv[1])
