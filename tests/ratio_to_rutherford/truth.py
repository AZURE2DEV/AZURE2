"""Independent numbers for tests/ratio_to_rutherford (mpmath; no AZURE2 code or convention).

Hard-sphere elastic scattering in the Coulomb field (every partial wave l <= L, no resonant width), from the textbook
partial-wave sum (as tests/coulomb_nuclear_sign/truth.py):

    f_C(theta) = -eta / (2k sin^2(theta/2)) exp(-i eta ln sin^2(theta/2)),
    f_N(theta) = (1 / 2ik) sum_l (2l+1) exp(2i w_l) [exp(-2i phi_l) - 1] P_l(cos theta),
    w_l = sigma_l - sigma_0,  phi_l = arg(G_l + i F_l) at rho = k a,

both relative to the common factor exp(2i sigma_0).  Distinct spin-0 particles (alpha + 12C):
dsigma/dOmega = |f_C + f_N|^2, Rutherford |f_C|^2.  Identical spin-0 bosons (alpha + alpha, even l only):
dsigma/dOmega = |f_C(theta) + f_C(pi - theta) + 2 f_N(theta)|^2, Mott |f_C(theta) + f_C(pi - theta)|^2.
Non-relativistic kinematics, k = sqrt(2 mu E)/hbar c, eta = Z1 Z2 alpha sqrt(mu / 2E), E_cm = E_lab m2/(m1 + m2);
the masses, radii and L of the models, the constants of include/Constants.h.

Usage (check.sh calls these; the output of `data` is committed as models/data/aC_rth_cm.dat):
  truth.py extrap aC|aa FILE           the value of every row of an AZURE2 .extrap file (blocks in the models' order)
  truth.py data                        synthetic ratio data on the grid of models/data/aC_rth_cm.dat
  truth.py multiply IN OUT             ratio data (E_lab theta_cm R dR) times the Rutherford cross section (b/sr)
  truth.py chi2 IN SHIFT_MEV           chi^2 of the ratio data against the hard-sphere model, both forms, at a lab
                                       energy shift: "ratio <chi2> multiplied <chi2>"
  truth.py fold IN SIGMA_LAB_MEV       per row: ratio of the Gaussian folds and fold of the ratio (window +-5 sigma,
                                       normalised on it, as the engine's kernel), and the folded Rutherford value
"""
import sys

import mpmath as mp

mp.mp.dps = 30
HBARC = mp.mpf("197.3269804")  # MeV fm
UCONV = mp.mpf("931.49410242")  # MeV / u
ALPHA = 1 / mp.mpf("137.035999084")

PAIRS = {
    # m1, m2 (u), Z1 Z2, channel radius (fm), L, identical
    "aC": (mp.mpf("4.00260325"), mp.mpf("12.0"), 12, mp.mpf("5.0"), 4, False),
    "aa": (mp.mpf("4.00260325"), mp.mpf("4.00260325"), 4, mp.mpf("4.5"), 6, True),
}
# the blocks of each model's .extrap file, in its order: what each row holds
BLOCKS = {"aC": ["dsdo", "ratio", "ratio"], "aa": ["dsdo", "ratio"]}


def kinematics(pair, e_cm):
    m1, m2, z1z2, a, lmax, ident = PAIRS[pair]
    mu = m1 * m2 / (m1 + m2) * UCONV
    k = mp.sqrt(2 * mu * e_cm) / HBARC
    eta = z1z2 * ALPHA * mp.sqrt(mu / (2 * e_cm))
    return k, eta


def amplitudes(pair, e_cm, theta_deg):
    """(nuclear + Coulomb amplitude, Coulomb amplitude) in fm, symmetrised for identical particles."""
    m1, m2, z1z2, a, lmax, ident = PAIRS[pair]
    k, eta = kinematics(pair, e_cm)
    th = mp.radians(theta_deg)

    def fc(t):
        s2 = mp.sin(t / 2) ** 2
        return -eta / (2 * k * s2) * mp.expj(-eta * mp.log(s2))

    fn = mp.mpc(0)
    w = mp.mpf(0)
    for l in range(lmax + 1):
        if l > 0:
            w += mp.atan(eta / l)
        if ident and l % 2:
            continue
        F = mp.coulombf(l, eta, k * a)
        G = mp.coulombg(l, eta, k * a)
        phi = mp.atan2(F, G)
        fn += (2 * l + 1) * mp.expj(2 * w) * (mp.expj(-2 * phi) - 1) * mp.legendre(l, mp.cos(th))
    fn /= 2j * k
    if ident:
        coul = fc(th) + fc(mp.pi - th)
        return coul + 2 * fn, coul
    return fc(th) + fn, fc(th)


def dsdo(pair, e_cm, theta_deg):
    """b/sr (1 fm^2 = 0.01 b)."""
    return abs(amplitudes(pair, e_cm, theta_deg)[0]) ** 2 / 100


def coulomb(pair, e_cm, theta_deg):
    return abs(amplitudes(pair, e_cm, theta_deg)[1]) ** 2 / 100


def ratio(pair, e_cm, theta_deg):
    f, c = amplitudes(pair, e_cm, theta_deg)
    return abs(f) ** 2 / abs(c) ** 2


def e_cm_of_lab(pair, e_lab):
    m1, m2 = PAIRS[pair][0], PAIRS[pair][1]
    return e_lab * m2 / (m1 + m2)


def rows(path):
    with open(path) as f:
        for line in f:
            p = line.split()
            if p and not p[0].startswith("#"):
                yield [mp.mpf(x) for x in p]


def extrap(pair, path):
    """Every row of the .extrap file: E_cm theta_cm AZURE2 truth (blocks separated by blank lines)."""
    block, blank = 0, False
    with open(path) as f:
        for line in f:
            p = line.split()
            if not p:
                blank = True
                continue
            if blank:
                block += 1
                blank = False
            e, th, v = mp.mpf(p[0]), mp.mpf(p[2]), p[3]
            kind = BLOCKS[pair][block]
            t = dsdo(pair, e, th) if kind == "dsdo" else ratio(pair, e, th)
            print("%d %s %.12e %.12e %s %.15e" % (block, kind, float(e), float(th), v, float(t)))


def data():
    """Ratio data: the hard-sphere ratio times 1 + 0.02 sin(7.3 i) (a fixed pattern, not noise), 3 % errors."""
    print("# E_lab(MeV) theta_cm(deg) R dR: alpha + 12C hard-sphere ratio to Rutherford x (1 + 0.02 sin(7.3 i)), "
          "3 % errors (tests/ratio_to_rutherford/truth.py data)")
    i = 0
    for e_lab in ["2.0", "2.5", "3.0", "3.5", "4.0", "4.5", "5.0"]:
        for th in range(20, 170, 20):
            r = ratio("aC", e_cm_of_lab("aC", mp.mpf(e_lab)), mp.mpf(th)) * (1 + mp.mpf("0.02") * mp.sin(mp.mpf("7.3") * i))
            print("%s %d %.12e %.12e" % (e_lab, th, float(r), float(r * mp.mpf("0.03"))))
            i += 1


def multiply(src, dst):
    with open(dst, "w") as out:
        out.write("# the rows of %s times the Rutherford cross section (b/sr; truth.py multiply)\n" % src)
        for e_lab, th, r, dr in rows(src):
            c = coulomb("aC", e_cm_of_lab("aC", e_lab), th)
            out.write("%s %s %.17e %.17e\n" % (mp.nstr(e_lab, 17), mp.nstr(th, 17), float(r * c), float(dr * c)))


def chi2_values(src, shift):
    """Both forms at lab energy E + shift: the ratio form compares sigma(E')/sigma_R(E') with R, the multiplied form
    sigma(E') with R sigma_R(E), error dR sigma_R(E).  Returns (ratio, multiplied)."""
    shift = mp.mpf(shift)
    cr = cm = mp.mpf(0)
    for e_lab, th, r, dr in rows(src):
        e0, e1 = e_cm_of_lab("aC", e_lab), e_cm_of_lab("aC", e_lab + shift)
        f1, c1 = amplitudes("aC", e1, th)
        s1, c1 = abs(f1) ** 2 / 100, abs(c1) ** 2 / 100
        c0 = coulomb("aC", e0, th)
        cr += ((s1 / c1 - r) / dr) ** 2
        cm += ((s1 - r * c0) / (dr * c0)) ** 2
    return cr, cm


def chi2(src, shift):
    cr, cm = chi2_values(src, shift)
    print("ratio %.15e multiplied %.15e" % (float(cr), float(cm)))


@mp.workdps(20)
def fold(src, sigma_lab, rows_wanted=None):
    """Per row of `src` (or the rows numbered in rows_wanted): E_lab theta, the ratio of the Gaussian folds
    <sigma>/<sigma_R>, the fold of the ratio <sigma/sigma_R>, and <sigma_R> (b/sr).  The Gaussian has the c.m. width of
    the lab width sigma_lab and is cut at +-5 widths and normalised on that window, as the engine's kernel
    (TargetEffect::convolutionRange); Gauss-Legendre, 8 panels of 12 nodes, 20 digits."""
    m1, m2 = PAIRS["aC"][0], PAIRS["aC"][1]
    s = mp.mpf(sigma_lab) * m2 / (m1 + m2)  # c.m. width
    x, w = zip(*mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(3, mp.mp.prec))  # 3 * 2^(3-1) = 12 nodes
    out = []
    for n, (e_lab, th, r, dr) in enumerate(rows(src)):
        if rows_wanted is not None and n not in rows_wanted:
            continue
        e0 = e_cm_of_lab("aC", e_lab)
        sums = [mp.mpf(0)] * 4  # norm, sigma, sigma_R, sigma / sigma_R
        for p in range(8):
            lo = e0 - 5 * s + p * 10 * s / 8
            half = 10 * s / 16
            for xi, wi in zip(x, w):
                e = lo + half * (1 + xi)
                g = wi * half * mp.exp(-(e - e0) ** 2 / (2 * s * s))
                f, c = amplitudes("aC", e, th)
                sig, ref = abs(f) ** 2 / 100, abs(c) ** 2 / 100
                sums = [sums[0] + g, sums[1] + g * sig, sums[2] + g * ref, sums[3] + g * sig / ref]
        out.append((e_lab, th, sums[1] / sums[2], sums[3] / sums[0], sums[2] / sums[0]))
    return out


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "extrap":
        extrap(sys.argv[2], sys.argv[3])
    elif cmd == "data":
        data()
    elif cmd == "multiply":
        multiply(sys.argv[2], sys.argv[3])
    elif cmd == "chi2":
        chi2(sys.argv[2], sys.argv[3])
    elif cmd == "fold":
        for e_lab, th, rof, fr, ref in fold(sys.argv[2], sys.argv[3]):
            print("%s %s %.15e %.15e %.15e" % (mp.nstr(e_lab, 6), mp.nstr(th, 6), float(rof), float(fr), float(ref)))
    else:
        sys.exit(__doc__)
