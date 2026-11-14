#!/usr/bin/env python3
"""pyazr.endf, reading (docs/upgrade/endf.md sections 2, 7.2, 8.3): numbers, records, interpolation, MF=3, MF=6,
curves at a data set's points, the module's own format check.

Hand values: every number below is computed here from the manual's definitions (ENDF-102, 2023: Table 16 and eqs.
1-6, pp. 25-26; eq. 6.1 p. 133; eq. 6.8 p. 143; eqs. 6.9-6.20 pp. 145-146), not by pyazr.endf.  The Coulomb terms
of LAW=5 and the lab energies are the engine's rows in pyazr.endf; here they are recomputed from the manual's
eqs. 6.9, 6.11, 6.12 and the non-relativistic E_lab = E_cm (m1 + m2)/m2.
The files are written by the test helper tests/pyazr/endf_testfiles.py; the library extract
tests/pyazr/data/endf/h-002_He_004_extract.endf (its README says what it is) is read as it is.
Exits 77 (skipped) without the compiled engine module.
Run from anywhere:  python3 tests/pyazr/endf_read_test.py
"""
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
try:
    import numpy as np
    import pyazr
    from pyazr import _azure2  # noqa: F401
    import pyazr.endf as E
    from pyazr.endf.format_check import check_format
except ImportError as e:
    print(f"skip: {e}")
    sys.exit(77)
from endf_testfiles import Writer, material

EXTRACT = os.path.join(HERE, "data", "endf", "h-002_He_004_extract.endf")
MN = 1.00866491595               # neutron mass (u), ENDF-102 p. 394
U_EV = 931.49410242e6            # u c^2 (eV), CODATA 2018
ALPHA = 1 / 137.035999084
HBARC_FM_EV = 197.3269804e6      # hbar c (eV fm)
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {str(detail)[:500]}"))
    if not ok:
        failures.append(name)


def close(a, b, rel=1e-9, absolute=0.0):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return bool(np.all(np.abs(a - b) <= rel * np.abs(b) + absolute))


def raises(call, cls, *parts):
    try:
        call()
    except cls as e:
        text = str(e)
        return all(p in text for p in parts) and all(getattr(e, k) for k in ("what", "allowed", "fix")), text
    except Exception as e:  # noqa: BLE001
        return False, f"{type(e).__name__}: {e}"
    return False, "nothing raised"


TMP = tempfile.mkdtemp(prefix="endf_read_")
AWR, AWI = 11.896910, 0.998623      # 12C (atomic) and p (nuclear, p. 24) in neutron masses
M1, M2 = AWI * MN, AWR * MN

# ------------------------------------------------------------------------------------------------- 1. numbers
print("numbers (p. 31)")
for text, value in (("-1.234567+6", -1.234567e6), (" 1.23456E+6", 1.23456e6), (" -9546287.6", -9546287.6),
                    (" 1159.42866", 1159.42866), ("   35070343", 35070343.0), ("  1.0e-5   ", 1e-5),
                    ("        1+3", 1000.0), ("           ", 0.0), (" 2.5-12    ", 2.5e-12), ("- 1.5 + 2  ", -150.0)):
    check(f"read {text!r}", E.read_float(text) == value, E.read_float(text))
ok, text = raises(lambda: E.read_float(" 1.0D+3"), pyazr.FormatError, "p. 31")
check("an exponent with D is refused (p. 31)", ok, text)
for name, digits in (("classic", 7), ("e", 6), ("widest", 9)):
    f = E.NUMBER_FORMS[name]
    for x in (1.234567891234, -9546287.61234, 3.5070343e7, 2.2e-12, -4.0198e6):
        s = f(x)
        check(f"{name} form of {x!r} fits 11 columns and keeps {digits} digits", len(s) == 11 and
              abs(E.read_float(s) - x) <= 0.5 * 10 ** (1 - digits) * abs(x) * 10, s)

# ------------------------------------------------------------------------------------------------- 2. laws
print("interpolation laws (Table 16, eqs. 1-6, pp. 25-26)")
xs, ys = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
nbt, laws = (2, 3, 4, 5, 6), (1, 2, 3, 4, 5)
hand = {1.5: 1.0, 2.5: 3.0, 3.5: 4 + 4 * math.log(3.5 / 3) / math.log(4 / 3),
        4.5: 8 * math.exp(math.log(2) * 0.5), 5.5: 16 * math.exp(math.log(2) * math.log(5.5 / 5) / math.log(6 / 5))}
got = E.interpolate(list(hand), xs, ys, nbt, laws)
check("INT 1-5 at mid-interval", close(got, list(hand.values()), 1e-14), (got, list(hand.values())))
check("a value at a point is the point's", close(E.interpolate([1.0, 4.0, 6.0], xs, ys, nbt, laws), [1, 8, 32]))
check("outside the table: NaN", np.all(np.isnan(E.interpolate([0.5, 6.5], xs, ys, nbt, laws))))
# INT=6: eq. 6 against eqs. 1-4 (the two forms of the manual)
e1, e2, s1, s2, t, e = 2.0, 4.0, 1.0, 3.0, 1.0, 3.0
b = math.log(s2 * e2 / (s1 * e1)) / (1 / math.sqrt(e1 - t) - 1 / math.sqrt(e2 - t))
a = math.exp(b / math.sqrt(e1 - t)) * s1 * e1
check("INT=6 against eqs. 1-4", close(E.interpolate(e, [e1, e2], [s1, s2], (2,), (6,), threshold=t),
                                         a / e * math.exp(-b / math.sqrt(e - t)), 1e-13))
check("INT=6 with a zero end is zero (p. 26)", E.interpolate(e, [e1, e2], [0.0, s2], (2,), (6,), threshold=t) == 0.0)
ok, text = raises(lambda: E.interpolate(1.5, [1, 2], [1, 2], (2,), (12,)), pyazr.FormatRefused, "INT=12")
check("INT=12 refused for a table of values", ok, text)

# ------------------------------------------------------------------------------------------------- 3. a material
print("a hand-written material: MF=1, MF=3, MF=6 LAW=2 (p. 143)")
QM = -7.5519e6
x3 = [1e6, 2e6, 3e6, 4e6, 5e6, 6e6]
y3 = [0.001, 0.002, 0.004, 0.008, 0.016, 0.032]
sections = {
    (3, 800): Writer().cont(6012.0, AWR, 0, 0, 0, 0).tab1(QM, QM, 0, 0, nbt, laws, x3, y3),
    (3, 801): Writer().cont(6012.0, AWR, 0, 0, 0, 0).tab1(QM, QM - 2.0e6, 0, 0, (6,), (2,), x3, y3),
    (6, 800): (Writer().cont(6012.0, AWR, 0, 2, 2, 0)
               .tab1(2004.0, 3.9737, 0, 2, (2,), (2,), [1e6, 6e6], [1.0, 0.5])
               .tab2(0.0, 0.0, 0, 0, (2,), (2,), 2)
               .list(0.0, 2e6, 0, 0, 2, [0.1, 0.2]).list(0.0, 4e6, 0, 0, 2, [0.3, -0.1])
               .tab1(5009.0, 8.9, 0, 4, (2,), (2,), [1e6, 6e6], [1.0, 1.0])),
    (6, 801): (Writer().cont(6012.0, AWR, 0, 2, 1, 0)
               .tab1(2004.0, 3.9737, 0, 2, (2,), (2,), [1e6, 6e6], [1.0, 1.0])
               .tab2(0.0, 0.0, 0, 0, (2,), (2,), 2)
               .list(0.0, 2e6, 12, 0, 2, [-1.0, 0.3, 1.0, 0.7]).list(0.0, 4e6, 14, 0, 2, [-1.0, 0.25, 1.0, 1.0])),
}
hand_file = material(os.path.join(TMP, "hand.endf"), sections)
m = E.read(hand_file)
check("read: the sections", sorted(m.sections) == [(1, 451), (3, 800), (3, 801), (6, 800), (6, 801)],
      sorted(m.sections))
h = m.header
check("header: ZA, AWR, AWI, NSUB, projectile", (h.za, h.awr, h.awi, h.nsub, h.projectile_za) ==
      (6012.0, AWR, AWI, 10010, 1001), h)
xs800 = m.cross_section(800)
check("MF=3/MT=800: QM, QI and the laws", (xs800.qm, xs800.qi, xs800.table.interpolation) == (QM, QM, laws))
check("MF=3/MT=800 at mid-interval (eV)", close(xs800([1e6 * k for k in hand]), [v * 1e-3 for v in hand.values()],
                                                1e-12))
mu = 0.5
p2 = (3 * mu * mu - 1) / 2
pa = 0.5 + 1.5 * 0.1 * mu + 2.5 * 0.2 * p2
pb = 0.5 + 1.5 * 0.3 * mu + 2.5 * (-0.1) * p2
d = m.distribution(800)
check("MF=6/MT=800: LCT, products", (d.lct, [(p.zap, p.law) for p in d.products]) == (2, [(2004, 2), (5009, 4)]))
from pyazr.endf.distributions import two_body_probability, elastic_cross_section
got = two_body_probability(d.product(2004), [2e6, 3e6, 4e6], [mu] * 3)
check("LAW=2 LANG=0 (eq. 6.8) at the panels and between them", close(got, [pa, (pa + pb) / 2, pb], 1e-14), got)
d1 = m.distribution(801).product(2004)
got = two_body_probability(d1, [2e6, 4e6], [mu, mu])
check("LAW=2 LANG=12 (linear in mu) and 14 (ln p linear in mu)",
      close(got, [0.3 + 0.4 * 0.75, 0.25 * math.exp(math.log(4) * 0.75)], 1e-14), got)
s = str(m.summary())
check("summary names the target, projectile and reactions", all(t in s for t in ("6-C - 12", "incident proton",
                                                                                   "(z,a0)", "(z,a1)", "LANG=0")), s)
check("the hand-written file passes the own check without a warning", not [
    f for f in check_format(hand_file).findings if f.severity != "note"], check_format(hand_file))


# ------------------------------------------------------------------------------------------------- 4. LAW=5
print("LAW=5 charged-particle elastic (eqs. 6.9-6.20, pp. 145-146)")


def eta_manual(e_lab_ev, z1=1, z2=6, m1=M1):
    return z1 * z2 * math.sqrt(ALPHA ** 2 * U_EV * m1 / (2 * e_lab_ev))      # eq. 6.12


def rutherford_manual(e_lab_ev, mu_, z1=1, z2=6, m1=M1, m2=M2):
    a = m2 / m1
    k = a / (1 + a) * math.sqrt(2 * U_EV * m1 * e_lab_ev) / HBARC_FM_EV * 10.0   # eq. 6.11: fm^-1 x 10 = b^-1/2
    return eta_manual(e_lab_ev, z1, z2, m1) ** 2 / (k * k * (1 - mu_) ** 2)       # eq. 6.9, b/sr


def legendre_p(l, x):
    return np.polynomial.legendre.legval(x, [0] * l + [1])


def law5_file(name, ltp, lidp, panels, nl, sigma_ni=None, za=6012.0, awr=AWR, nsub=10010, awi=AWI, spin=0.0):
    w = Writer().cont(za, awr, 0, 2, 1, 0).tab1(nsub // 10, awi, 0, 5, (2,), (2,), [1e6, 3e6], [1.0, 1.0])
    w.tab2(spin, 0.0, lidp, 0, (2,), (2,), len(panels))
    for energy, values in panels:
        w.list(0.0, energy, ltp, 0, nl, values)
    secs = {(6, 2): w, (3, 2): Writer().cont(za, awr, 0, 0, 0, 0).tab1(0.0, 0.0, 0, 0, (2,), (2,), [1e6, 3e6],
                                                                       sigma_ni or [1.0, 1.0])}
    return E.read(material(os.path.join(TMP, name), secs, za=za, awr=awr, awi=awi, nsub=nsub))


mu, e = 0.3, 2e6
ca, cb = [0.1, 0.05, 0.02], [0.2, 0.0, 0.04]
mm = law5_file("ltp2.endf", 2, 0, [(1e6, ca), (3e6, cb)], 2)
coul = E.CoulombTerms(mm)
sc = coul.sigma(np.array([e]), np.array([mu]))[0]
check("the Coulomb term (engine rows on the file's masses) against eqs. 6.9, 6.11, 6.12", close(
    sc, rutherford_manual(e, mu), 2e-6), (sc, rutherford_manual(e, mu)))
check("eta against eq. 6.12", close(coul.eta(np.array([e]))[0], eta_manual(e), 2e-6))
sr = lambda c: sum((2 * l + 1) / 2 * c[l] * legendre_p(l, mu) for l in range(3))   # noqa: E731, eq. 6.17
got = elastic_cross_section(mm.distribution(2).product(1001), [e], [mu], coul)[0] - sc
check("LTP=2 distinguishable: sigma - sigma_c = sigma_R/(1 - mu) (eqs. 6.15, 6.17)",
      close(got, (sr(ca) + sr(cb)) / 2 / (1 - mu), 1e-12), got)

ci_a, ci_b = [0.1, 0.05], [0.3, -0.02]
mi = law5_file("ltp2i.endf", 2, 1, [(1e6, ci_a), (3e6, ci_b)], 1, za=2004.0, awr=3.968219, nsub=20040, awi=3.968219)
ci = E.CoulombTerms(mi)
sri = lambda c: sum((4 * l + 1) / 2 * c[l] * legendre_p(2 * l, mu) for l in range(2))  # noqa: E731, eq. 6.18
got = elastic_cross_section(mi.distribution(2).product(2004), [e], [mu], ci)[0] - ci.sigma(
    np.array([e]), np.array([mu]), identical=True)[0]
check("LTP=2 identical: sigma - sigma_c = sigma_R/(1 - mu^2) (eqs. 6.16, 6.18)",
      close(got, (sri(ci_a) + sri(ci_b)) / 2 / (1 - mu * mu), 1e-12), got)

b_d, a_d = [0.3, 0.1, 0.05], [complex(0.2, -0.1), complex(-0.05, 0.3)]
vals = b_d + [v for z in a_d for v in (z.real, z.imag)]
ml = law5_file("ltp1.endf", 1, 0, [(1e6, vals), (3e6, vals)], 1)


def ltp1_part(energy, identical=False, b=b_d, a=a_d):
    eta = eta_manual(energy)
    terms = [(2 * l + 1) / 2 * a[l] * legendre_p(l, mu) for l in range(len(a))]
    if not identical:   # eq. 6.13
        inter = -2 * eta / (1 - mu) * (np.exp(1j * eta * math.log((1 - mu) / 2)) * sum(terms)).real
        return inter + sum((2 * l + 1) / 2 * b[l] * legendre_p(l, mu) for l in range(len(b)))
    s_ = sum(((1 + mu) * np.exp(1j * eta * math.log((1 - mu) / 2)) + (-1) ** l * (1 - mu) *
              np.exp(1j * eta * math.log((1 + mu) / 2))) * terms[l] for l in range(len(a)))   # eq. 6.14
    return -2 * eta / (1 - mu * mu) * s_.real + sum((4 * l + 1) / 2 * b[l] * legendre_p(2 * l, mu)
                                                     for l in range(len(b)))


cl = E.CoulombTerms(ml)
got = elastic_cross_section(ml.distribution(2).product(1001), [e], [mu], cl)[0] - cl.sigma(np.array([e]),
                                                                                            np.array([mu]))[0]
want = (ltp1_part(1e6) + ltp1_part(3e6)) / 2
check("LTP=1 distinguishable (eq. 6.13), NW = 4 NL + 3", close(got, want, 1e-5), (got, want))
b_i = [0.3, 0.1]
vals_i = b_i + [v for z in a_d for v in (z.real, z.imag)]
mli = law5_file("ltp1i.endf", 1, 1, [(1e6, vals_i), (3e6, vals_i)], 1, za=2004.0, awr=3.968219, nsub=20040,
                awi=3.968219)
cli = E.CoulombTerms(mli)
got = elastic_cross_section(mli.distribution(2).product(2004), [e], [mu], cli)[0] - cli.sigma(
    np.array([e]), np.array([mu]), identical=True)[0]
eta_i = lambda energy: 4 * math.sqrt(ALPHA ** 2 * U_EV * 3.968219 * MN / (2 * energy))  # noqa: E731
want = 0.0
for energy in (1e6, 3e6):
    eta = eta_i(energy)
    terms = [(2 * l + 1) / 2 * a_d[l] * legendre_p(l, mu) for l in range(2)]
    s_ = sum(((1 + mu) * np.exp(1j * eta * math.log((1 - mu) / 2)) + (-1) ** l * (1 - mu) *
              np.exp(1j * eta * math.log((1 + mu) / 2))) * terms[l] for l in range(2))
    want += 0.5 * (-2 * eta / (1 - mu * mu) * s_.real + sum((4 * l + 1) / 2 * b_i[l] * legendre_p(2 * l, mu)
                                                            for l in range(2)))
check("LTP=1 identical (eq. 6.14), NW = 3 NL + 3", close(got, want, 1e-5), (got, want))
for ltp, ppairs, want_p in ((12, [-1.0, 0.2, 1.0, 0.6], 0.2 + 0.4 * 0.65), (15, [-1.0, 0.2, 1.0, 0.6], 0.46),
                            (14, [-1.0, 0.25, 1.0, 1.0], 0.25 * math.exp(math.log(4) * 0.65))):
    mt_ = law5_file(f"ltp{ltp}.endf", ltp, 0, [(1e6, ppairs), (3e6, ppairs)], 2, sigma_ni=[0.5, 0.7])
    ct = E.CoulombTerms(mt_)
    got = elastic_cross_section(mt_.distribution(2).product(1001), [e], [mu], ct, mt_.cross_section(2))[0] - ct.sigma(
        np.array([e]), np.array([mu]))[0]
    check(f"LTP={ltp}: sigma - sigma_c = sigma_NI(E) P_NI(mu) (eqs. 6.19-6.20)", close(got, 0.6 * want_p, 1e-12),
          (got, 0.6 * want_p))

# ------------------------------------------------------------------------------------------------- 5. curves
print("curves at a data set's points (section 7.2 of the design)")
P_ENT = E.PairSpec(1, 6, 1.00782503, 12.0, 1, 12, j1=0.5)
P_A0 = E.PairSpec(2, 5, 4.00260325, 9.0133288, 4, 9)
P_A1 = E.PairSpec(2, 5, 4.00260325, 9.0133288, 4, 9, excitation=2.0)
to_cm = lambda elab: elab * P_ENT.m2_u / (P_ENT.m1_u + P_ENT.m2_u) * 1e-6     # noqa: E731, MeV
c = m.curves(E.CurveRequest("sigma", (to_cm(3e6), to_cm(3.5e6), to_cm(9e6)), (0, 0, 0), P_ENT, P_A0))
check("sigma of (p,a0): MT=800 at the lab energies of the points", close(c.values[:2], [0.004, hand[3.5] * 1e-3],
                                                                          1e-9) and c.source == "MF=3/MT=800", c)
check("a point outside the table: NaN, not covered, a note", np.isnan(c.values[2]) and not c.coverage[2] and c.notes)
c = m.curves(E.CurveRequest("sigma", (to_cm(4e6),), (0,), P_ENT, P_A1))
check("the exit pair with the residual at 2 MeV is MT=801 (QM - QI)", c.source == "MF=3/MT=801" and
      close(c.values, [0.008], 1e-9), c)
c = m.curves(E.CurveRequest("dsigma_domega", (to_cm(3e6),), (60.0,), P_ENT, P_A0))
check("dsigma/dOmega = sigma y p / 2 pi (eq. 6.1)", close(c.values, [0.004 * 0.8 * (pa + pb) / 2 / (2 * math.pi)],
                                                          1e-9) and c.unit == "b/sr", c)
ecm = to_cm(3e6)
eta_cm = 6 * math.sqrt(ALPHA ** 2 * U_EV * P_ENT.m1_u * P_ENT.m2_u / (P_ENT.m1_u + P_ENT.m2_u) / (2 * ecm * 1e6))
c = m.curves(E.CurveRequest("sfactor", (ecm,), (0,), P_ENT, P_A0))
check("sfactor = sigma E exp(2 pi eta)", close(c.values, [0.004 * ecm * math.exp(2 * math.pi * eta_cm)], 2e-6), c)
for request, part in ((E.CurveRequest("analyzing_power", (1.0,), (90,), P_ENT, P_A0), "analyzing_power"),
                      (E.CurveRequest("sigma", (1.0,), (0,), P_ENT, P_A1.__class__(2, 5, 4.0026, 9.01, 4, 9,
                                                                                    excitation=1.0)), "1 MeV"),
                      (E.CurveRequest("sigma", (1.0,), (0,), E.PairSpec(2, 6, 4.0026, 12, 4, 12), P_A0), "2004"),
                      (E.CurveRequest("sigma", (1.0,), (0,), P_ENT, E.PairSpec(0, 7, 0, 13, 0, 13, photon=True)),
                       "one final state")):
    ok, text = raises(lambda: m.curves(request), pyazr.FormatRefused, part)
    check(f"refused with its three parts: {part}", ok, text)
cr = mm.curves(E.CurveRequest("ratio_to_rutherford", (to_cm(e),), (math.degrees(math.acos(mu)),), P_ENT, P_ENT))
r_pair = rutherford_manual(e, mu, m1=P_ENT.m1_u, m2=P_ENT.m2_u)
check("ratio_to_rutherford: (sigma_c + sigma_R/(1 - mu)) / Rutherford of the data set's pair", close(
    cr.values, [(rutherford_manual(e, mu) + (sr(ca) + sr(cb)) / 2 / (1 - mu)) / r_pair], 5e-6), cr)

# a session's data set: T1's elastic data at 90 deg lab against the same request built by hand
shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T1"), os.path.join(TMP, "T1"),
                ignore=shutil.ignore_patterns("output*"))
os.makedirs(os.path.join(TMP, "T1", "output"), exist_ok=True)
s1 = pyazr.Session(os.path.join(TMP, "T1", "T1.azr"))
ev = s1.evaluate(s1.parameters.x0)
data = ev.data[0]
c_session = mm.curves(data, session=s1)
pair = [p for p in s1.pairs if p.key == 1][0]
hand_request = E.CurveRequest("dsigma_domega", tuple(data.energies), tuple(data.angles), E.PairSpec.of(pair),
                              E.PairSpec.of(pair), tuple(data.extra["energy_lab"]))
c_hand = mm.curves(hand_request)
check("Material.curves(data, session=s) is the curve of the request made from the data set by hand",
      np.array_equal(c_session.values, c_hand.values, equal_nan=True) and c_session.name == data.name, c_session)
check("the curve covers the points inside 1-3 MeV lab only", np.array_equal(
    c_session.coverage, (np.asarray(data.extra["energy_lab"]) >= 1.0) & (np.asarray(data.extra["energy_lab"]) <= 3.0)))
s1.close()

# ------------------------------------------------------------------------------------------------- 6. columns
print("the fields' column rules (pp. 30-35): refused by the reader or reported by the check")
with open(hand_file) as f:
    lines = f.read().splitlines()
k800 = next(i for i, ln in enumerate(lines) if ln[70:75] == " 3800" and ln[75:80].strip() == "2")
k801 = next(i for i, ln in enumerate(lines) if ln[70:75] == " 3801" and ln[75:80].strip() == "2")


def variant(name, change):
    v = list(lines)
    change(v)
    path = os.path.join(TMP, name)
    with open(path, "w") as f:
        f.write("\n".join(v) + "\n")
    return path


def codes(path):
    return {(f.severity, f.code) for f in check_format(path).findings}


p = variant("long.endf", lambda v: v.__setitem__(k800 + 2, v[k800 + 2] + "  junk"))
check("a line beyond column 80: read as if not there", close(E.read(p).cross_section(800)(3.5e6), hand[3.5] * 1e-3))
check("... and reported as a warning (p. 30)", ("warning", "BEYOND_COL_80") in codes(p), codes(p))
p = variant("seq.endf", lambda v: v.__setitem__(k800 + 2, v[k800 + 2][:75] + "   17"))
check("a wrong sequence number: read; reported as a warning (p. 33)", E.read(p).cross_section(800) is not None and
      ("warning", "SEQUENCE") in codes(p), codes(p))
p = variant("np_more.endf", lambda v: v.__setitem__(k800, v[k800][:55] + f"{7:11d}" + v[k800][66:]))
ok, text = raises(lambda: E.read(p).cross_section(800), pyazr.FormatError, "p. ")
check("NP larger than the lines: refused by the reader (file:line, p. 34-35)", ok and ":" in text, text)
check("... and an error of the check", ("error", "RECORD") in codes(p), codes(p))
p = variant("np_less.endf", lambda v: (v.__setitem__(k801, v[k801][:55] + f"{2:11d}" + v[k801][66:]),
                                        v.__setitem__(k801 + 1, f"{2:11d}{2:11d}".ljust(66) + v[k801 + 1][66:])))
ok, text = raises(lambda: E.read(p).cross_section(801), pyazr.FormatError, "lines after the last record")
check("NP smaller than the lines: refused by the reader", ok, text)
check("... and an error of the check", ("error", "COUNT") in codes(p), codes(p))
p = variant("np_5.endf", lambda v: (v.__setitem__(k801, v[k801][:55] + f"{5:11d}" + v[k801][66:]),
                                    v.__setitem__(k801 + 1, f"{5:11d}{2:11d}".ljust(66) + v[k801 + 1][66:])))
check("NP one smaller, the same lines: read (5 points); the check warns of the values after the last",
      len(E.read(p).cross_section(801).table.x) == 5 and ("warning", "COUNT") in codes(p), codes(p))
p = variant("d.endf", lambda v: v.__setitem__(k800 + 2, v[k800 + 2][:11] + " 2.000000D6" + v[k800 + 2][22:]))
ok, text = raises(lambda: E.read(p).cross_section(800), pyazr.FormatError, "p. 31", "d.endf")
check("a number with D: refused, with file, line, columns and p. 31", ok and f":{k800 + 3}:12-22" in text, text)
p = variant("short.endf", lambda v: v.__setitem__(k800 + 2, v[k800 + 2][:60]))
ok, text = raises(lambda: E.read(p), pyazr.FormatError, "MAT, MF, MT")
check("a line without its control numbers: refused by read()", ok, text)
check("... and an error of the check", ("error", "SHORT_LINE") in codes(p) or ("error", "CONTROL") in codes(p),
      codes(p))
k_dir = next(i for i, ln in enumerate(lines) if ln[22:44] == f"{3:11d}{800:11d}")
p = variant("dir.endf", lambda v: v.__setitem__(k_dir, v[k_dir][:44] + f"{99:11d}" + v[k_dir][55:]))
check("a directory NC that is not the section's lines: an error of the check (p. 41)",
      ("error", "DIRECTORY") in codes(p), codes(p))
p = variant("nosend.endf", lambda v: v.__delitem__(next(i for i, ln in enumerate(v) if ln[70:75] == " 3  0"
                                                        and ln[75:80].strip() == "99999")))
check("a section without SEND: an error of the check (p. 34)", ("error", "NO_SEND") in codes(p), codes(p))
unseq = material(os.path.join(TMP, "unseq.endf"), sections, sequence=False)
check("a file without sequence numbers: a note only", codes(unseq) == {("note", "UNSEQUENCED")}, codes(unseq))
for form in ("e", "widest"):
    other = E.read(material(os.path.join(TMP, f"{form}.endf"), sections, form=form))
    check(f"the {form} form read back gives the same table", close(other.cross_section(800).table.y, y3, 1e-6))

# ------------------------------------------------------------------------------------------------- 7. extract
print("the ENDF/B-VIII.1 extract (3He + 4He, MAT 228)")
x = E.read(EXTRACT)
s = str(x.summary())
check("summary: MAT 228, helion sublibrary, LRF=7 KRM=4, MT=600", all(t in s for t in (
    "MAT 228", "incident helion", "R-Matrix Limited (KRM=4 full R-matrix, IFG=1, KRL=0, 2 pairs, 10 spin groups)",
    "(z,p0)", "97 energies")), s)
t = x.cross_section(600)
check("MF=3/MT=600 at its own energies gives its values", close(t(list(t.table.x)), list(t.table.y), 0))
check("MF=3/MT=600: QM = QI = -4.0198 MeV, 97 points from 0.05 MeV", (t.qm, t.qi, len(t.table.x), t.table.x[0]) ==
      (-4019800.0, -4019800.0, 97, 50000.0))
found = {(f.severity, f.code, f.mf, f.mt) for f in check_format(EXTRACT).findings}
check("the own check of the extract: the directory names the sections left out, LRP=2 is a warning", found == {
    ("note", "UNSEQUENCED", 0, 0), ("warning", "LRP", 1, 451), ("error", "DIRECTORY", 1, 451)}, found)

if os.environ.get("ENDF_TEST_KEEP"):  # keep the files written, for a run of the real codes on them
    print(f"files kept in {TMP}")
else:
    shutil.rmtree(TMP, ignore_errors=True)
print(f"{len(failures)} failure(s)" + (": " + ", ".join(failures) if failures else ""))
sys.exit(1 if failures else 0)
