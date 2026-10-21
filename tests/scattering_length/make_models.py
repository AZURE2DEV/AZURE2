#!/usr/bin/env python3
"""Builds the scattering-length test models (committed; rerun only to change them).

he3he4.azr  3He+4He, Paneru et al., PRC 109, 015802 (2024) Table V: the 1/2+ background level
            (Ex 14.000 MeV, gamma~ 1.683 MeV^1/2, a = 4.2 fm), the only s-wave term.  a0 = 33.10 fm.
be7p.azr    7Be+p, Paneru et al., PRC 99, 045807 (2019) Table II, a = 4.3 fm: the 2- level (Ex 3.520,
            gamma~(S=2) 1.700) and the 1- background (Ex 9.000, gamma~_el(S=1) 1.433, gamma~ in
            7Be*(0.429)+p (S=1) -1.822, closed at the 7Be+p threshold).
n16o.azr    n+16O with one 1/2+ level (Ex 6.0 MeV in 17O, gamma~ 0.8 MeV^1/2, a = 4.0 fm): neutral, so
            a0 = a (1 - gamma~^2 / E~) in closed form (S = 0 for l = 0).
Widths are written as observed (Brune) partial widths in eV, Gamma~_c = 2 P_c gamma~_c^2 / (1 + sum gamma~^2 S'),
so that AZURE2 recovers the published gamma~ (parameters.out g_int); each model has
a phase-shift test segment at low energy for the effective-range check of check.sh.
Needs mpmath."""
import mpmath as mp
import os
mp.mp.dps = 30
hbarc, uconv, fstruc = mp.mpf("197.32696310"), mp.mpf("931.4940880"), 1 / mp.mpf("137.0359996790")
HERE = os.path.dirname(os.path.abspath(__file__))

CONFIG = """<config>
true #Perform A-Matrix Calculation
output/ #Full Path to Output Directory
checks/ #Full Path to Checks Directory
none #Compond Nucleus Check
none #Boundary Condition Check
none #Data Check
none #Lo-Matrix and Penetrability Check
none #Legendre Polynomial Check
none #Coulomb Amplitudes Check
none #Reaction Pathway Check
none #Angular Distributions Check
</config>
"""
TAIL = """<segmentsData>
</segmentsData>
<segmentsTest>
{tests}
</segmentsTest>
<targetInt>
</targetInt>
<lastRun>
9997
""
""
0 0
0 ""
-1. -1. -1.
</lastRun>
"""

def penetrability(E, m1, m2, z1z2, a, l=0):
    mu = m1 * m2 / (m1 + m2) * uconv
    k = mp.sqrt(2 * mu * E) / hbarc
    rho = k * a
    if z1z2 == 0:
        return rho if l == 0 else None
    eta = z1z2 * fstruc * mp.sqrt(mu / (2 * E))
    F, G = mp.coulombf(l, eta, rho), mp.coulombg(l, eta, rho)
    return rho / (F**2 + G**2)

def shift(E, m1, m2, z1z2, a, l=0):
    mu = m1 * m2 / (m1 + m2) * uconv
    k = mp.sqrt(2 * mu * E) / hbarc
    rho = k * a
    if z1z2 == 0:
        return mp.mpf(0)
    eta = z1z2 * fstruc * mp.sqrt(mu / (2 * E))
    F, G = mp.coulombf(l, eta, rho), mp.coulombg(l, eta, rho)
    Fp = mp.diff(lambda r: mp.coulombf(l, eta, r), rho)
    Gp = mp.diff(lambda r: mp.coulombg(l, eta, r), rho)
    return rho * (F * Fp + G * Gp) / (F**2 + G**2)

def widths_eV(channels):
    """Observed (Brune) partial widths of one level from its gamma~: Brune, PRC 66, 044611 (2002),
    Gamma_c = 2 P_c gamma~_c^2 / (1 + sum_c' gamma~_c'^2 dS_c'/dE), all at E~.
    channels: [(gamma, Ecm, m1, m2, z1z2, a), ...]"""
    den = 1 + sum(g**2 * mp.diff(lambda E: shift(E, m1, m2, z, a), Ec) for g, Ec, m1, m2, z, a in channels)
    return [mp.sign(g) * 2 * penetrability(Ec, m1, m2, z, a) * g**2 / den * 1e6 for g, Ec, m1, m2, z, a in channels]

def line(J, pi, Ex, ir, s, l, lid, gamma_eV, j1, pi1, j2, pi2, e2, m1, m2, z1, z2, entSep, sep, a):
    return (f"{J:g} {pi} {Ex:.6f} 1 1 {ir} {int(round(2*s))} {2*l} {lid} 1 1 {mp.nstr(gamma_eV, 12)} "
            f"{j1:g} {pi1} {j2:g} {pi2} {e2:g} {m1} {m2} {z1} {z2} {entSep} {sep} 0 0 0.0 0 {a} 0 0 0")

def write(name, levels, tests):
    with open(os.path.join(HERE, name), "w") as f:
        f.write(CONFIG + "<levels>\n" + "\n\n".join(levels) + "\n\n</levels>\n" + TAIL.format(tests="\n".join(tests)))

# 3He+4He
m1, m2, Sa, a = 3.0160293, 4.0026032, 1.5866, 4.2
g, = widths_eV([(1.683, 14.000 - Sa, m1, m2, 4, a)])
write("he3he4.azr",
      [line(0.5, 1, 14.000, 1, 0.5, 0, 1, g, 0.5, 1, 0, 1, 0, m1, m2, 2, 2, Sa, Sa, a)],
      ["1 1 1 0.05 1.05 0.02 0 0 0 2 0.5 0"])          # lab 3He energies (beam 3He on 4He)

# 7Be+p
m1, m2, Sp, a, Exs = 1.00783, 7.01693, 0.1375, 4.3, 0.42908
g2m, = widths_eV([(1.700, 3.520 - Sp, m1, m2, 4, a)])
g1m_el, g1m_in = widths_eV([(1.433, 9.000 - Sp, m1, m2, 4, a), (-1.822, 9.000 - Sp - Exs, m1, m2, 4, a)])
el = dict(j1=0.5, pi1=1, j2=1.5, pi2=-1, e2=0, m1=m1, m2=m2, z1=1, z2=4, entSep=Sp, sep=Sp, a=a)
inl = dict(el, j2=0.5, e2=Exs)
write("be7p.azr",
      [line(2, -1, 3.520, 1, 2, 0, 1, g2m, **el),
       line(1, -1, 9.000, 1, 1, 0, 2, g1m_el, **el) + "\n" + line(1, -1, 9.000, 2, 1, 0, 2, g1m_in, **inl)],
      ["1 1 1 0.02 0.30 0.01 0 0 0 2 2 0", "1 1 1 0.02 0.30 0.01 0 0 0 2 1 0"])   # lab proton energies

# n+16O (neutral)
m1, m2, Sn, a = 1.00866492, 15.99491462, 4.1431, 4.0
g, = widths_eV([(0.8, 6.0 - Sn, m1, m2, 0, a)])
write("n16o.azr",
      [line(0.5, 1, 6.0, 1, 0.5, 0, 1, g, 0.5, 1, 0, 1, 0, m1, m2, 0, 8, Sn, Sn, a)],
      ["1 1 1 0.01 0.41 0.01 0 0 0 2 0.5 0"])

# he3he4_fit.azr: the 3He+4He model with the level energy and width free, phase-shift data generated
# from he3he4.azr (data/he3he4_phase.dat: AZURE2 phase shifts at E_lab = 1.0 .. 5.0 MeV, +- 0.5 deg)
# and a measured a0 = 31.0 +- 0.3 fm, 2 fm below the model's 33.10, so the constraint pulls the fit.
src = open(os.path.join(HERE, "he3he4.azr")).read()
lev = [l for l in src.splitlines() if l.startswith("0.5 1 14.000000")][0].split()
lev[3] = "0"     # level energy free
lev[10] = "0"    # width free
fit = src.replace([l for l in src.splitlines() if l.startswith("0.5 1 14.000000")][0], " ".join(lev))
fit = fit.replace("<segmentsData>\n</segmentsData>\n",
                  "<segmentsData>\n1 1 1 0 10 0 180 2 0.5 0 1 0 0 0 0 0 data/he3he4_phase.dat 0 0\n</segmentsData>\n"
                  "<scatteringLengths>\n# active pairKey s value sigma [imValue imSigma]\n1 1 0.5 31.0 0.3\n</scatteringLengths>\n")
open(os.path.join(HERE, "he3he4_fit.azr"), "w").write(fit)

# be7p_fit.azr: the 7Be+p model with every level energy and width free, a second 1- level (Ex 5.0 MeV,
# gamma~ 0.5 in 7Be+p and 0.3 in 7Be*+p) so the 1- group has two levels (the Brune off-diagonal terms),
# phase-shift data generated from be7p.azr (data/be7p_phase_2m.dat, _1m.dat: E_lab 0.3 .. 1.9 MeV,
# +- 0.5 deg) and constraints on every kind of row: a01, a02, the coherent combination of pair 1, and
# both parts of the complex 1- length of 7Be*+p (its 7Be+p channel is open at that threshold).
m1, m2, Sp, a, Exs = 1.00783, 7.01693, 0.1375, 4.3, 0.42908
g5_el, g5_in = widths_eV([(0.5, 5.000 - Sp, m1, m2, 4, a), (0.3, 5.000 - Sp - Exs, m1, m2, 4, a)])
src = open(os.path.join(HERE, "be7p.azr")).read()
out = []
for l in src.splitlines():
    f = l.split()
    if len(f) == 31 and l[0].isdigit():
        f[3] = "0"
        f[10] = "0"
        l = " ".join(f)
    out.append(l)
fit = "\n".join(out) + "\n"
extra = (line(1, -1, 5.000, 1, 1, 0, 3, g5_el, **el) + "\n" + line(1, -1, 5.000, 2, 1, 0, 3, g5_in, **inl)).split("\n")
extra = [" ".join(x.split()[:3] + ["0"] + x.split()[4:10] + ["0"] + x.split()[11:]) for x in extra]
fit = fit.replace("\n</levels>", "\n" + "\n".join(extra) + "\n\n</levels>")
fit = fit.replace("<segmentsData>\n</segmentsData>\n",
                  "<segmentsData>\n1 1 1 0 10 0 180 2 2 0 1 0 0 0 0 0 data/be7p_phase_2m.dat 0 0\n"
                  "1 1 1 0 10 0 180 2 1 0 1 0 0 0 0 0 data/be7p_phase_1m.dat 0 0\n</segmentsData>\n"
                  "<scatteringLengths>\n"
                  "1 1 1 17.0 0.5\n"
                  "1 1 2 -4.0 0.5\n"
                  "1 1 -1 5.0 1.0\n"
                  "1 2 1 9.0 0.5 -0.3 0.05\n"
                  "0 2 0 25.0 1.0\n"
                  "</scatteringLengths>\n")
open(os.path.join(HERE, "be7p_fit.azr"), "w").write(fit)
