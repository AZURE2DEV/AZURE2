#!/usr/bin/env python3
"""Independent reference for the global optical potentials (src/ThmOptical.cpp).

(1) Parameters.  The ten numbers V,R,a,W,RW,aW,WD,RD,aD,RC of each model at
a few (projectile, target, lab energy) points, evaluated from the RIPL-3
optical-model library (Capote et al., Nucl. Data Sheets 110 (2009) 3107;
www-nds.iaea.org/RIPL-3/optical/om-data/om-parameter-u.dat, read and
evaluated here by its documented coefficient formulas, standard and Koning
forms) for

    ancai06      RIPL 6200   (An & Cai 2006)
    kd03         RIPL 2405 (n), 5405 (p)   (Koning & Delaroche 2003, global)
    bg71         RIPL 7100 (t), 8100 (3He) (Becchetti & Greenlees)
    mcfadden66   RIPL 9100   (McFadden & Satchler 1966)
    avrigeanu94  RIPL 9600   (Avrigeanu, Hodgson & Avrigeanu 1994)

and from the FRONT21 front end of TWOFNR (J.A. Tostevin, Surrey; the
subroutines `deuteron`, option 3, and `liangchuntian`, transcribed here) for

    daehnick80   (Daehnick, Childs & Vrcelj 1980)   -- RIPL 6112-6116 only
                 carries a polynomial refit of its imaginary terms; its real
                 depth and the imaginary diffuseness agree with these formulas
                 for 40Ca, 90Zr, 120Sn and 208Pb (printed as a check);
    liang09      (Liang, Li & Cai 2009).

Radii R = r0 A^(1/3) (RIPL's convention), spin-orbit dropped, negative
imaginary depths set to 0 (as the engine does).

(2) Elastic scattering, d + 40Ca at 56 MeV with daehnick80: the nuclear
S matrix from scipy's DOP853 (rtol 1e-11, power-series start at 1e-3 fm,
matched to mpmath's Coulomb functions at 25 fm), the ratio to Rutherford
sigma/sigma_R at the angles of the Hatanaka et al. (1980) data (EXFOR
E0682-022), printed with the data's ratio.  The engine computes S_l by
Numerov + COUL (ThmDistortion::Wave).

Run:  python3 tests/reference/thm_optical_reference.py [om-parameter-u.dat]
It prints the tables hard-coded in thm_optical_test.cpp (about a minute).
"""
import math
import os
import re
import sys

# --------------------------------------------------------------------------
# RIPL-3 reader (the part of om-parameter-u.readme these entries use).


def _num(t):
    m = re.match(r'^([+-]?\d*\.?\d*)([+-]\d+)$', t)
    if m and m.group(1) not in ('', '+', '-'):
        return float(m.group(1)) * 10 ** int(m.group(2))
    return float(t)


def ripl_load(path, wanted):
    lines = open(path).read().split('\n')
    out = {}
    for i, line in enumerate(lines):
        s = line.strip()
        if not (re.fullmatch(r'\d+', s) and int(s) in wanted and int(s) not in out):
            continue
        if re.fullmatch(r'[\d.\s+-]*', lines[i + 1]):
            continue
        toks = []
        j = i + 7
        while len(toks) < 2000 and j < len(lines):
            toks += lines[j].split()
            j += 1
        k = [0]

        def nx():
            k[0] += 1
            return _num(toks[k[0] - 1])
        e = dict(emin=nx(), emax=nx(), zmin=nx(), zmax=nx(), amin=nx(), amax=nx())
        e['imodel'], e['zp'], e['ap'], e['irel'], e['idr'] = [nx() for _ in range(5)]
        e['terms'] = []
        for _ in range(6):
            rng = []
            for _ in range(abs(int(nx()))):
                ep = nx()
                rco = [nx() for _ in range(13)]
                aco = [nx() for _ in range(13)]
                pot = [nx() for _ in range(25)]
                rng.append(dict(epot=ep, rco=rco, aco=aco, pot=pot))
            e['terms'].append(rng)
        e['coul'] = [[nx() for _ in range(8)] for _ in range(int(nx()))]
        out[int(s)] = e
    return out


def ripl_eval(e, Z, A, E):
    """Spherical entry at lab energy E: [(depth, R, a) for terms 1..4], RC."""
    a13 = A ** (1.0 / 3.0)
    eta = (A - 2.0 * Z) / A
    rc = 0.0
    if e['coul']:
        c = e['coul'][0]
        rc = (c[1] / a13 + c[2] + c[3] * A ** (-2.0 / 3.0) + c[4] * A ** (-5.0 / 3.0) + c[7] * A) * a13
    out = []
    for t in range(4):
        rng = e['terms'][t]
        if not rng:
            out.append((0.0, 0.0, 0.0))
            continue
        r = rng[0]
        rco, aco, pot = r['rco'], r['aco'], r['pot']
        assert rco[12] == 0.0
        R = (abs(rco[0]) + rco[2] * eta + rco[3] / A + rco[4] / math.sqrt(A) + rco[5] * A ** (2 / 3)
             + rco[6] * A + rco[7] * A ** 2 + rco[8] * A ** 3 + rco[9] * a13 + rco[10] / a13
             + rco[1] * E + rco[11] * E * E) * a13
        a = (abs(aco[0]) + aco[1] * E + aco[2] * eta + aco[3] / A + aco[4] / math.sqrt(A)
             + aco[5] * A ** (2 / 3) + aco[6] * A + aco[7] * A ** 2 + aco[8] * A ** 3 + aco[9] * a13
             + aco[10] / a13 + aco[11] * E)
        if pot[23] == 0.0:  # standard form
            ec1 = 0.4 * Z / a13
            ec2 = 1.73 * Z / rc if rc else 0.0
            V = (pot[0] + pot[6] * eta + pot[7] * ec1 + pot[8] * A + pot[9] * a13 + pot[10] * A ** (-2 / 3)
                 + pot[11] * ec2 + (pot[1] + pot[12] * eta + pot[13] * A) * E + pot[2] * E * E
                 + pot[3] * E ** 3 + pot[5] * math.sqrt(E) + (pot[4] + pot[14] * eta + pot[15] * E) * math.log(E)
                 + pot[16] * ec1 / E ** 2)
        else:  # Koning form (pot(24) = 1), terms 1, 2, 4
            ef = pot[17] + pot[18] * A
            f = E - ef
            n = int(pot[12])
            if t == 0:
                b1 = pot[0] + pot[1] * A + pot[7] * eta
                b2 = pot[2] + pot[3] * A
                b3 = pot[4] + pot[5] * A
                b4 = pot[6]
                b5 = pot[8]
                b15 = pot[13] + pot[14] * A or 1.0
                V = b1 * (b15 - b2 * f + b3 * f * f - b4 * f ** 3)
                if b5:
                    V += b5 * b1 * (1.73 * Z / rc) * (b2 - 2 * b3 * f + 3 * b4 * f * f)
            elif t == 1:
                b6 = pot[0] + pot[1] * A + pot[7] * eta
                b7 = pot[2] + pot[3] * A + pot[8] * eta
                V = b6 * f ** n / (f ** n + b7 ** n)
            elif t == 3:
                b8 = pot[0] + pot[7] * eta + pot[6] * A + pot[8] * A ** (-1 / 3)
                b9 = pot[1] + pot[9] * A + pot[2] / (1 + math.exp((A - pot[3]) / pot[4]))
                b10 = pot[5] + pot[10] * A
                b12 = pot[11] or 1.0  # 0 in the file: the exponent of (E - EF), 1
                V = b8 * math.exp(-b9 * f ** b12) * f ** n / (f ** n + b10 ** n)
            else:
                raise ValueError("term %d of a Koning-form entry" % (t + 1))
        out.append((V, R, a))
    return out, rc


def ten_from_ripl(e, Z, A, E):
    (v, w, _vd, wd), rc = ripl_eval(e, Z, A, E)
    p = [v[0], v[1], v[2], max(0.0, w[0]), w[1], w[2], max(0.0, wd[0]), wd[1], wd[2], rc]
    for t in (1, 2):
        if p[3 * t] == 0.0:
            p[3 * t + 1] = p[3 * t + 2] = 0.0
    return p


# --------------------------------------------------------------------------
# FRONT21 (TWOFNR front end) transcriptions.


def front_daehnick(Z, A, E):
    a13 = A ** (1 / 3)
    n = A - Z
    bet = math.exp(-(E / 100.0) ** 2)
    v = 88.5 - 0.26 * E + 0.88 * Z / a13
    wv = (12.2 + 0.026 * E) * (1 - bet)
    ws = (12.2 + 0.026 * E) * bet
    ai = 0.53 + 0.07 * a13
    for m in (8, 20, 28, 50, 82, 126):
        ai -= 0.04 * math.exp(-((m - n) / 2.0) ** 2)
    return [v, 1.17 * a13, 0.709 + 0.0017 * E, wv, 1.325 * a13, ai, ws, 1.325 * a13, ai, 1.30 * a13]


def front_liang(Z, A, E):
    at13 = A ** (1 / 3)
    vd = 118.36 - 0.2071 * E + 6.3961e-5 * E * E + 26.001 * (A - 2 * Z) / A + 0.5668 * Z / at13
    wd = -6.8871 + 0.3115 * E - 6.8096e-4 * E * E
    wisd = 20.119 - 0.1626 * E - 5.4067 * (A - 2 * Z) / A + 1.2087 * at13
    p = [vd, (1.1657 + 0.0401 / at13) * at13, 0.6641 + 0.0305 * at13,
         max(0.0, wd), (1.4022 + 0.0418 / at13) * at13, 0.7732 + 0.0219 * at13,
         max(0.0, wisd), (1.1802 + 0.0587 / at13) * at13, 0.6292 + 0.0657 * at13, 1.289 * at13]
    for t in (1, 2):
        if p[3 * t] == 0.0:
            p[3 * t + 1] = p[3 * t + 2] = 0.0
    return p


# model, (Zp, Ap), (Zt, At), E_lab
CASES = [
    ("ancai06", (1, 2), (9, 19), 5.83), ("ancai06", (1, 2), (12, 24), 2.0),
    ("ancai06", (1, 2), (20, 40), 56.0), ("ancai06", (1, 2), (82, 208), 100.0),
    ("daehnick80", (1, 2), (20, 40), 56.0), ("daehnick80", (1, 2), (28, 58), 20.0),
    ("daehnick80", (1, 2), (82, 208), 80.0), ("daehnick80", (1, 2), (9, 19), 5.8),
    ("kd03", (0, 1), (10, 20), 2.6), ("kd03", (0, 1), (26, 56), 14.0),
    ("kd03", (0, 1), (82, 208), 0.5), ("kd03", (0, 1), (82, 208), 150.0),
    ("kd03", (1, 1), (20, 40), 30.0), ("kd03", (1, 1), (82, 208), 65.0), ("kd03", (1, 1), (12, 24), 3.0),
    ("bg71", (1, 3), (40, 90), 15.0), ("bg71", (2, 3), (28, 58), 30.0),
    ("liang09", (2, 3), (8, 16), 10.0), ("liang09", (2, 3), (50, 120), 100.0),
    ("mcfadden66", (2, 4), (8, 16), 20.0), ("mcfadden66", (2, 4), (82, 208), 24.7),
    ("avrigeanu94", (2, 4), (20, 40), 20.0), ("avrigeanu94", (2, 4), (82, 208), 50.0),
]
RIPL = {("ancai06", 2): 6200, ("kd03", 1): 5405, ("kd03", 0): 2405, ("bg71", 3): 7100,
        ("bg71", 3.5): 8100, ("mcfadden66", 4): 9100, ("avrigeanu94", 4): 9600}


def reference_ten(lib, model, zp, ap, zt, at, e):
    if model == "daehnick80":
        return front_daehnick(zt, at, e), "FRONT21"
    if model == "liang09":
        return front_liang(zt, at, e), "FRONT21"
    if model == "kd03":
        iref = 5405 if zp == 1 else 2405
    elif model == "bg71":
        iref = 8100 if zp == 2 else 7100
    else:
        iref = {"ancai06": 6200, "mcfadden66": 9100, "avrigeanu94": 9600}[model]
    return ten_from_ripl(lib[iref], zt, at, e), "RIPL %d" % iref


# --------------------------------------------------------------------------
# Elastic scattering.

HBARC = 197.32696310
UCONV = 931.4940880
FSTRUC = 1 / 137.0359996790
M_D, M_CA40 = 2.0135532134, 39.9516979  # nuclear masses (u); 40Ca: AME2020 atomic - 20 m_e + B_e
HATANAKA = [(14.72, 0.166), (18.92, 0.1175), (21.01, 0.1531), (23.11, 0.1052), (31.48, 0.02596),
            (34.1, 0.04067), (36.7, 0.03507), (44.51, 0.01292), (49.69, 0.01451), (57.43, 0.00614),
            (67.69, 0.00285), (77.87, 0.000875)]  # theta_cm (deg), dsigma/dOmega (b/sr), EXFOR E0682-022


def elastic(p, elab, angles, lmax=45):
    import mpmath as mp
    import numpy as np
    from scipy.integrate import solve_ivp
    mu = M_D * M_CA40 / (M_D + M_CA40) * UCONV
    ecm = elab * M_CA40 / (M_D + M_CA40)
    k = math.sqrt(2 * mu * ecm) / HBARC
    eta = 20 * FSTRUC * mu / (HBARC * k)
    fac = 2 * mu / HBARC ** 2
    zz = 20 * FSTRUC * HBARC
    V, R, a, W, RW, aW, WD, RD, aD, RC = p

    def U(r):
        vc = zz / r if r >= RC else zz * (3 - r * r / RC ** 2) / (2 * RC)
        f = lambda rr, aa: 1 / (1 + math.exp((r - rr) / aa))
        x = math.exp((r - RD) / aD)
        return vc - V * f(R, a) - 1j * W * f(RW, aW) - 4j * WD * x / (1 + x) ** 2

    rm = 25.0
    S = []
    for l in range(lmax + 1):
        def rhs(r, y):
            g = l * (l + 1) / r ** 2 + fac * U(r) - k * k
            return [y[1], g * y[0]]
        r0 = 1e-3
        sol = solve_ivp(rhs, (r0, rm), [complex(r0 ** (l + 1)), complex((l + 1) * r0 ** l)],
                        method="DOP853", rtol=1e-11, atol=1e-300)
        u, du = sol.y[0, -1], sol.y[1, -1]
        rho = k * rm
        F, G = float(mp.coulombf(l, eta, rho)), float(mp.coulombg(l, eta, rho))
        dF = float(mp.diff(lambda t: mp.coulombf(l, eta, t), rho)) * k
        dG = float(mp.diff(lambda t: mp.coulombg(l, eta, t), rho)) * k
        # u = A (H- - S H+), H+- = G +- iF
        Hp, Hm, dHp, dHm = G + 1j * F, G - 1j * F, dG + 1j * dF, dG - 1j * dF
        L = du / u
        S.append((dHm - L * Hm) / (dHp - L * Hp))
    sig0 = float(mp.arg(mp.gamma(1 + 1j * eta)))
    out = []
    for th in angles:
        t = math.radians(th)
        s2 = math.sin(t / 2) ** 2
        fc = -eta / (2 * k * s2) * complex(math.cos(-eta * math.log(s2) + 2 * sig0),
                                           math.sin(-eta * math.log(s2) + 2 * sig0))
        fn = 0j
        sig = sig0
        x = math.cos(t)
        pm, pl = 1.0, x
        for l in range(lmax + 1):
            if l > 0:
                sig += math.atan(eta / l)
            P = 1.0 if l == 0 else (x if l == 1 else None)
            if l >= 2:
                pm, pl = pl, ((2 * l - 1) * x * pl - (l - 1) * pm) / l
                P = pl
            fn += (2 * l + 1) * complex(math.cos(2 * sig), math.sin(2 * sig)) * (S[l] - 1) * P / (2j * k)
        ratio = abs(fc + fn) ** 2 / abs(fc) ** 2
        ruth = abs(fc) ** 2 * 10.0  # fm^2/sr -> mb/sr
        out.append((th, ratio, ruth))
    return S, out


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "om-parameter-u.dat"
    if not os.path.exists(path):
        sys.exit("needs RIPL-3 om-parameter-u.dat (www-nds.iaea.org/RIPL-3/optical/om-data/)")
    lib = ripl_load(path, {6200, 2405, 5405, 7100, 8100, 9100, 9600, 6112, 6113, 6114, 6115, 6116})
    print("// model, Zp, Ap, Zt, At, E_lab, V,R,a, W,RW,aW, WD,RD,aD, RC   (source)")
    for model, (zp, ap), (zt, at), e in CASES:
        p, src = reference_ten(lib, model, zp, ap, zt, at, e)
        print('    {"%s", %d, %d, %d, %d, %g, {%s}},  // %s' % (
            model, zp, ap, zt, at, e, ", ".join("%.10g" % v for v in p), src))
    print("\nDaehnick: FRONT21 formulas against RIPL 6112-6116 (V at 56 MeV, imaginary diffuseness)")
    for iref, (z, a) in zip((6112, 6113, 6114, 6115, 6116), ((20, 40), (28, 58), (40, 90), (50, 120), (82, 208))):
        (v, w, _x, wd), _rc = ripl_eval(lib[iref], z, a, 56.0)
        f = front_daehnick(z, a, 56.0)
        print("  %3d: V %.4f / %.4f   aI %.4f / %.4f   WV %.3f / %.3f   WD %.3f / %.3f" % (
            a, f[0], v[0], f[5], w[2], f[3], w[0], f[6], wd[0]))
    print("\nElastic d + 40Ca, 56 MeV, daehnick80: theta, sigma/sigma_R (DOP853), data/sigma_R (Hatanaka 1980)")
    p = front_daehnick(20, 40, 56.0)
    S, out = elastic(p, 56.0, [a for a, _ in HATANAKA])
    for l in (0, 5, 10, 20):
        print("  S_%d = %.10f %+.10fi" % (l, S[l].real, S[l].imag))
    for (th, ratio, ruth), (_a, d) in zip(out, HATANAKA):
        print("    {%.2f, %.8e, %.5e},  // data/sigma_R; model/data %.3f" % (th, ratio, d * 1000.0 / ruth,
                                                                         ratio * ruth / (d * 1000.0)))


if __name__ == "__main__":
    main()
