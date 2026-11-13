#!/usr/bin/env python3
"""Units of format 3 and the quantities that are scales of the engine's observables (project_format3.md sections 3.3,
7.2; data_kinds.md sections 4.1, 4.6), on 12C(p,g)13N with the levels of tests/13N:

1. units of a data file: sigma in nb against the same data in b (values exact in nb), energies in keV against MeV
   (exact in keV): every model value and chi^2 term identical to the bit; an angle column in rad against deg and an
   E_cm column against E_lab (both converted by this script) to 1e-12; an areal density in atoms/b against atoms/cm2
   to 1e-12;
2. quantity sfactor: the chi^2 of S data (keV*b) equals the chi^2 of the same data converted to sigma by this script
   (mpmath, sigma = S exp(-2 pi eta)/E with the <levels> masses and charges) to 1e-12; with a fixed 5 keV shift the
   two differ by what the script computes from the engine's sigma at the shifted energies; the analytic gradient
   against central differences; pyazr.physics.sfactor_conversion against the script;
3. quantity yield on an infinite target: a file per uC with charge state 2 equals the file per particle times
   1e-6 / (2 e) (chi^2 to 1e-12); counts per uC with a fixed efficiency likewise;
4. pyazr.sfactor_to_cross_section (atomic masses through the engine's function) against the script, and against the
   mass-number form it replaced (which is 6.7 % off for p + 12C at 100 keV).

Needs the compiled module and mpmath; skips cleanly without them.

Run from anywhere:  python3 tests/pyazr/units_scales_test.py
"""
import math
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


try:
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    import mpmath
    import pyazr
    from pyazr import Session, _azure2
except Exception as err:
    print(f"skip: the engine or mpmath is not importable ({type(err).__name__}: {err})")
    sys.exit(77)

mpmath.mp.dps = 40
# The pair p + 12C of tests/13N's <levels> (masses in u as written there) and the constants of physics/Constants.h.
M1, M2, Z1, Z2 = mpmath.mpf("1.00728"), mpmath.mpf("12"), 1, 6
UCONV, ALPHA = mpmath.mpf("931.49410242"), 1 / mpmath.mpf("137.035999084")
E_CHARGE = mpmath.mpf("1.602176634e-19")


def conversion(e_cm):
    """E exp(2 pi eta) (MeV), eta = Z1 Z2 alpha sqrt(mu c^2 / 2E), mu = m1 m2 / (m1 + m2) (u)."""
    e = mpmath.mpf(e_cm)
    mu = M1 * M2 / (M1 + M2)
    eta = Z1 * Z2 * ALPHA * mpmath.sqrt(mu * UCONV / (2 * e))
    return e * mpmath.exp(2 * mpmath.pi * eta)


def cm_of_lab(e_lab):
    return mpmath.mpf(e_lab) * M2 / (M1 + M2)


LEVELS = re.search(r"<levels>.*?</levels>\n", open(os.path.join(ROOT, "tests", "13N", "13N.azr")).read(), re.S).group(0)
CONFIG = ("<config>\ntrue #Perform A-Matrix Calculation\noutput/ #Output\nchecks/ #Checks\n" + "none #check\n" * 8 +
          "format=3\n</config>\n")


def project_text(dataset, targets=""):
    return CONFIG + LEVELS + "<data>\n" + dataset + "</data>\n<extrapolations>\n</extrapolations>\n<targets>\n" + targets + \
        "</targets>\n"


def dataset(file, quantity="sigma", unit="b", columns="E_lab[MeV] _ value error", extra="", reaction="12C(p,g)13N"):
    return (f"dataset D\n    reaction       {reaction}\n    quantity       {quantity}\n    unit           {unit}\n"
            f"    file           data/{file}\n    columns        {columns}\n    norm           1 fixed\n" + extra)


def write_project(work, name, text, files):
    d = os.path.join(work, name)
    os.makedirs(os.path.join(d, "output"))
    os.makedirs(os.path.join(d, "checks"))
    os.makedirs(os.path.join(d, "data"))
    for fname, rows in files.items():
        with open(os.path.join(d, "data", fname), "w") as f:
            for row in rows:
                f.write(" ".join(row) + "\n")
    with open(os.path.join(d, "m.azr"), "w") as f:
        f.write(text)
    return d


def terms_of(d):
    """(model values, chi^2 terms) at the project's start parameters."""
    with Session(os.path.join(d, "m.azr"), cwd=d) as m:
        x = np.asarray(m.params_rwa, float)
        model = np.concatenate([np.asarray(v, float) for v in [d.model for d in m.evaluate(x).data]])
        return model, np.asarray([m.evaluate(x).chi2], float)


def same_bits(a, b):
    return a.shape == b.shape and a.tobytes() == b.tobytes()


# beam energies (lab): 300 .. 1000 keV in steps of 12.5 keV (exact in keV and, as decimals, in MeV)
E_KEV = [300 + 12.5 * k for k in range(57)]
E_MEV = ["%.4f" % (e / 1000) for e in E_KEV]
TARGET_INF = ("target T\n    kind           infinite\n    active         12C\n    composition    C 1\n"
              "    stopping       formula \"a0*x^a1\" a0=4e-20 a1=-0.5 MeV*cm2/atom\n    straggling     none\n")

with tempfile.TemporaryDirectory() as work:
    # the model at the start parameters, as data of a dummy file: the cross sections the data below are built from
    probe = write_project(work, "probe", project_text(dataset("p.dat")), {"p.dat": [(e, "0", "1", "1") for e in E_MEV]})
    sigma, _ = terms_of(probe)
    # values exact in nb: sigma in nb rounded to a multiple of 2^-10 nb (a dyadic number), 5 % errors likewise
    wobble = [1 + 0.05 * math.sin(3.1 * k) for k in range(len(E_KEV))]
    nb = [round(s * 1e9 * w * 1024) / 1024 for s, w in zip(sigma, wobble)]
    dnb = [round(0.05 * v * 1024) / 1024 or 1.0 / 1024 for v in nb]

    def decimal(x):
        """The exact decimal of a double (a dyadic rational) as text."""
        return str(mpmath.mpf(x)) if x == int(x) else format(mpmath.mpf(x), ".40g").rstrip("0")

    print("1. units of the columns of a data file")
    rows_b = [(e, "0", decimal(mpmath.mpf(v) / 10**9), decimal(mpmath.mpf(dv) / 10**9)) for e, v, dv in zip(E_MEV, nb, dnb)]
    rows_nb = [(e, "0", repr(v), repr(dv)) for e, v, dv in zip(E_MEV, nb, dnb)]
    rows_kev = [(repr(e), "0", repr(v), repr(dv)) for e, v, dv in zip(E_KEV, nb, dnb)]
    in_b = write_project(work, "b", project_text(dataset("d.dat")), {"d.dat": rows_b})
    in_nb = write_project(work, "nb", project_text(dataset("d.dat", unit="nb")), {"d.dat": rows_nb})
    in_kev = write_project(work, "kev", project_text(dataset("d.dat", unit="nb", columns="E_lab[keV] _ value error")),
                           {"d.dat": rows_kev})
    mb_, cb = terms_of(in_b)
    mn, cn = terms_of(in_nb)
    mk, ck = terms_of(in_kev)
    check("nb file = b file: model values and chi^2 terms identical to the bit", same_bits(mb_, mn) and same_bits(cb, cn),
          f"chi2 {cb.sum()!r} {cn.sum()!r}")
    check("keV energies = MeV energies: identical to the bit", same_bits(mn, mk) and same_bits(cn, ck),
          f"chi2 {cn.sum()!r} {ck.sum()!r}")
    check("the data are not trivial (chi^2 > 1)", cb.sum() > 1.0, f"{cb.sum()!r}")
    # E_cm[keV]: the c.m. energies by this script (non-relativistic, the <levels> masses)
    rows_cm = [(format(cm_of_lab(mpmath.mpf(e) / 1000) * 1000, ".20g"), "0", repr(v), repr(dv))
               for e, v, dv in zip(E_KEV, nb, dnb)]
    in_cm = write_project(work, "cm", project_text(dataset("d.dat", unit="nb", columns="E_cm[keV] _ value error")),
                          {"d.dat": rows_cm})
    mc, cc = terms_of(in_cm)
    check("E_cm[keV] column = E_lab column to 1e-12 (chi^2)", rel(cc.sum(), cn.sum()) <= 1e-12,
          f"{cc.sum()!r} {cn.sum()!r} ({rel(cc.sum(), cn.sum()):.1e})")
    # an angle column in rad: dsigma/dOmega of the capture gamma at 55 deg lab
    rows_deg = [(e, "55", repr(v / (4 * math.pi)), repr(dv / (4 * math.pi))) for e, v, dv in zip(E_MEV, nb, dnb)]
    rows_rad = [(e, repr(math.radians(55.0)), repr(v / (4 * math.pi)), repr(dv / (4 * math.pi))) for e, v, dv in zip(E_MEV, nb, dnb)]
    in_deg = write_project(work, "deg", project_text(dataset("d.dat", "dsigma_domega", "nb/sr", "E_lab[MeV] theta_lab[deg] value error")),
                           {"d.dat": rows_deg})
    in_rad = write_project(work, "rad", project_text(dataset("d.dat", "dsigma_domega", "nb/sr", "E_lab[MeV] theta_lab[rad] value error")),
                           {"d.dat": rows_rad})
    cd_, cr_ = terms_of(in_deg)[1].sum(), terms_of(in_rad)[1].sum()
    check("theta_lab[rad] = theta_lab[deg] to 1e-12 (chi^2)", rel(cr_, cd_) <= 1e-12, f"{cr_!r} {cd_!r} ({rel(cr_, cd_):.1e})")
    # a target layer's areal density in atoms/b and in atoms/cm2 (the target-averaged cross section, eq. 2)
    layer = ("target L\n    kind           layer\n    active         12C\n    composition    C 1\n"
             "    areal_density  {}\n    stopping       formula \"a0*x^a1\" a0=4e-20 a1=-0.5 MeV*cm2/atom\n"
             "    straggling     none\n")
    folded = "    target         L\n    quadrature     points 40\n"
    in_cm2 = write_project(work, "cm2", project_text(dataset("d.dat", unit="nb", extra=folded), layer.format("3e18 atoms/cm2 fixed")),
                           {"d.dat": rows_nb})
    in_atb = write_project(work, "atb", project_text(dataset("d.dat", unit="nb", extra=folded), layer.format("3e-6 atoms/b fixed")),
                           {"d.dat": rows_nb})
    c2, c3 = terms_of(in_cm2)[1].sum(), terms_of(in_atb)[1].sum()
    check("areal density 3e-6 atoms/b = 3e18 atoms/cm2 to 1e-12 (chi^2)", rel(c3, c2) <= 1e-12 and rel(c2, cn.sum()) > 1e-6,
          f"{c3!r} {c2!r} (unfolded {cn.sum()!r})")

    print("2. quantity sfactor")
    s_kevb = [mpmath.mpf(v) * 1e-9 * conversion(cm_of_lab(e)) * 1000 for v, e in zip(nb, E_MEV)]
    ds_kevb = [mpmath.mpf(dv) * 1e-9 * conversion(cm_of_lab(e)) * 1000 for dv, e in zip(dnb, E_MEV)]
    rows_s = [(e, "0", format(s, ".17g"), format(ds, ".17g")) for e, s, ds in zip(E_MEV, s_kevb, ds_kevb)]
    # the same S data as sigma by this script: sigma = S / (E exp(2 pi eta)), b from keV b
    rows_sig = [(e, "0", format(mpmath.mpf(s[2]) / 1000 / conversion(cm_of_lab(e)), ".17g"),
                 format(mpmath.mpf(s[3]) / 1000 / conversion(cm_of_lab(e)), ".17g")) for e, s in zip(E_MEV, rows_s)]
    in_s = write_project(work, "s", project_text(dataset("s.dat", "sfactor", "keV*b")), {"s.dat": rows_s})
    in_sig = write_project(work, "sig", project_text(dataset("s.dat")), {"s.dat": rows_sig})
    ms, cs = terms_of(in_s)
    msig, csig = terms_of(in_sig)
    check("chi^2 of S data = chi^2 of the script's sigma to 1e-12", rel(cs.sum(), csig.sum()) <= 1e-12,
          f"{cs.sum()!r} {csig.sum()!r} ({rel(cs.sum(), csig.sum()):.1e})")
    check("the model is S (MeV b) = conversion x sigma", max(rel(ms[i], msig[i] * float(conversion(cm_of_lab(E_MEV[i]))))
                                                             for i in range(len(ms))) <= 1e-13)
    shift = "    shift          offset 5 keV fixed\n"
    in_s5 = write_project(work, "s5", project_text(dataset("s.dat", "sfactor", "keV*b", extra=shift)), {"s.dat": rows_s})
    in_sig5 = write_project(work, "sig5", project_text(dataset("s.dat", extra=shift)), {"s.dat": rows_sig})
    _, cs5 = terms_of(in_s5)
    msig5, csig5 = terms_of(in_sig5)
    # the script: S model at the shifted energy E' = E + 5 keV (lab), sigma(E') from the engine's sigma form
    expected = mpmath.mpf(0)
    for i, (e, s, ds) in enumerate(zip(E_MEV, s_kevb, ds_kevb)):
        model = conversion(cm_of_lab(mpmath.mpf(e) + mpmath.mpf("0.005"))) * msig5[i] * 1000
        expected += ((model - s) / ds) ** 2
    expected = float(expected)
    check("with a 5 keV shift: chi^2 of S data = the script's (the factor at the shifted energy) to 1e-10",
          rel(cs5.sum(), expected) <= 1e-10, f"{cs5.sum()!r} {expected!r}")
    check("  and the sigma form differs from it (the data converted at the unshifted energy) by %.4g" %
          (csig5.sum() - cs5.sum()), rel(csig5.sum(), cs5.sum()) > 1e-4)
    e_test = np.array([0.1, 0.4, 1.2])
    import pyazr.physics
    conv = pyazr.physics.sfactor_conversion(1, 6, 1.00728, 12.0, e_test)
    check("pyazr.physics.sfactor_conversion = the script to 1e-13", all(rel(conv[i], float(conversion(e_test[i]))) <= 1e-13 for i in range(3)),
          str(conv))
    # the gradient of the S-factor chi^2 (norm and the level parameters) against central differences
    in_sg = write_project(work, "sg", project_text(dataset("s.dat", "sfactor", "keV*b").replace(
        "    norm           1 fixed\n", "    norm           1 fit prior 1 5%\n    shift          offset 2 keV fit\n")),
        {"s.dat": rows_s})
    with Session(os.path.join(in_sg, "m.azr"), cwd=in_sg) as m:
        x = np.asarray(m.params_rwa, float)
        c, g = m.objective().value_and_grad(x)
        worst = 0.0
        for i in range(len(x)):
            h = 1e-6 * max(abs(x[i]), 1e-3)
            xp, xm = x.copy(), x.copy()
            xp[i] += h
            xm[i] -= h
            fd = (m.objective()(xp) - m.objective()(xm)) / (2 * h)
            worst = max(worst, abs(g[i] - fd) / max(np.max(np.abs(g)), 1e-300))
        check(f"gradient: {len(x)} free parameters (with norm and shift), analytic = central differences to 1e-6 of the "
              f"largest", len(x) >= 3 and worst <= 1e-6, f"worst {worst:.2e}, chi2 {c!r}")

    print("3. quantity yield on an infinite target, per particle and per uC")
    yprobe = write_project(work, "yprobe", project_text(dataset("y.dat", "yield", "1/particle", extra="    target T\n    quadrature     points 40\n"),
                                                        TARGET_INF), {"y.dat": [(e, "0", "1", "1") for e in E_MEV]})
    yields, _ = terms_of(yprobe)
    ypart = [float(y) * w for y, w in zip(yields, wobble)]
    rows_yp = [(e, "0", repr(y), repr(0.05 * y)) for e, y in zip(E_MEV, ypart)]
    per_uc = 1e-6 / (2 * float(E_CHARGE))  # particles per uC at charge state 2
    rows_yq = [(e, "0", repr(y * per_uc), repr(0.05 * y * per_uc)) for e, y in zip(E_MEV, ypart)]
    rows_yc = [(e, "0", repr(y * per_uc * 0.25), repr(0.05 * y * per_uc * 0.25)) for e, y in zip(E_MEV, ypart)]
    folded_inf = "    target T\n    quadrature     points 40\n"
    in_yp = write_project(work, "yp", project_text(dataset("y.dat", "yield", "1/particle", extra=folded_inf), TARGET_INF), {"y.dat": rows_yp})
    in_yq = write_project(work, "yq", project_text(dataset("y.dat", "yield", "1/uC", extra=folded_inf + "    charge_state 2\n"),
                                                   TARGET_INF), {"y.dat": rows_yq})
    in_yc = write_project(work, "yc", project_text(dataset("y.dat", "yield", "counts/uC",
                                                           extra=folded_inf + "    charge_state 2\n    efficiency 0.25 fixed\n"),
                                                   TARGET_INF), {"y.dat": rows_yc})
    cyp, cyq, cyc = (terms_of(d)[1].sum() for d in (in_yp, in_yq, in_yc))
    check("per uC (q = 2) = per particle x 1e-6/(2e): chi^2 to 1e-12", rel(cyq, cyp) <= 1e-12 and cyp > 1.0,
          f"{cyq!r} {cyp!r} ({rel(cyq, cyp):.1e})")
    check("counts per uC with efficiency 0.25: chi^2 to 1e-12", rel(cyc, cyp) <= 1e-12, f"{cyc!r} {cyp!r}")

    print("4. pyazr.sfactor_to_cross_section")
    s_kev_b, e_kev = 20.0, 100.0  # 20 keV b at 100 keV (c.m.)
    got = pyazr.sfactor_to_cross_section(s_kev_b, e_kev, 1, 6, 1.00728, 12.0)
    want = float(mpmath.mpf(s_kev_b) / 1000 / conversion(mpmath.mpf(e_kev) / 1000))
    old = s_kev_b * math.exp(-2 * math.pi * 0.157488 * 6 * math.sqrt(12 / 13 / 0.1)) / e_kev  # mass numbers, 0.157488
    check("sigma from S with the <levels> masses = the script to 1e-13", rel(got, want) <= 1e-13, f"{got!r} {want!r}")
    check("the mass-number form is %.1f %% off at 100 keV" % (100 * (old / want - 1)), rel(old, want) > 0.05)

if failures:
    print(f"\n{len(failures)} FAILED")
    sys.exit(1)
print("\nall passed")
