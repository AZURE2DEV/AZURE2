#!/usr/bin/env python3
"""The limits of the yield of a layer (README.md), each closed form computed here from models/truth.json and never from
the engine; the engine's model values are read from the output of a calculation (menu 1) of models/layer_yield.azr.

    python3 -I tests/layer_yield/check.py <truth.json> <output directory>      (prints a table; exit 1 on a failure)

Frames and units (README.md, "Formulas"): E_lab and E_cm in MeV, sigma in b (c.m.), epsilon in lab MeV cm^2 per active
nucleus, N in atoms/cm^2, Y in reactions per incident particle; per uC: Y 1e-6 C / (q e).  Only the standard library.
"""
import json
import math
import os
import sys

HBARC = 197.3269804      # MeV fm (physics/Constants.h)
UMEV = 931.49410242      # MeV per u
E_CHARGE = 1.602176634e-19  # C

failures = []
rows_out = []


def finite(r):
    """A relative difference, infinite when it is not a number (an engine value of NaN must fail, never be skipped)."""
    return r if math.isfinite(r) else math.inf


def check(name, engine, closed, rel, tol, limit):
    ok = rel <= tol
    rows_out.append((name, closed, engine, rel, tol, limit, ok))
    if not ok:
        failures.append(name)


def blocks(path):
    """The rows of an AZUREOut file, one list per data set (blank lines separate them), each row a list of floats."""
    out, cur = [], []
    with open(path) as f:
        for line in f:
            if line.strip():
                cur.append([float(x) for x in line.split()])
            elif cur:
                out.append(cur)
                cur = []
    if cur:
        out.append(cur)
    return out


def f_cm(p):
    return p["m2"] / (p["m1"] + p["m2"])


def mu(p):
    return p["m1"] * p["m2"] / (p["m1"] + p["m2"]) * UMEV


def sigma_hard_sphere(p, e_cm):
    """s-wave hard sphere of radius a, spin-0 target, spin-1/2 projectile (g = 1): 4 pi sin^2(ka)/k^2, in b."""
    k = math.sqrt(2.0 * mu(p) * e_cm) / HBARC
    return 4.0 * math.pi * math.sin(k * p["R"]) ** 2 / k ** 2 * 0.01


def sigma_bw(p, L, e_cm):
    """Single-level Breit-Wigner of observed constant widths, c.m.: (pi/k^2) omega Gp Ga / ((E-Er)^2 + G^2/4), in b."""
    k2 = 2.0 * mu(p) * e_cm / HBARC ** 2
    omega = (2 * L["J"] + 1) / ((2 * p["j1"] + 1) * (2 * p["j2"] + 1))
    gp, ga = L["Gp_eV"] * 1e-6, L["Ga_eV"] * 1e-6
    g = gp + ga
    return math.pi / k2 * 0.01 * omega * gp * ga / ((e_cm - L["Er_cm"]) ** 2 + 0.25 * g * g)


def exit_pow(E0, N, a0, a1):
    """E_out of a layer with epsilon = a0 x^a1 (lab): int_{E_out}^{E0} dE/eps = N exactly."""
    b = 1.0 - a1
    return (E0 ** b - a0 * b * N) ** (1.0 / b)


def mean_energy_pow(E0, Eo, a0, a1):
    """E_bar = int E dE/eps / int dE/eps over [E_out, E0] for epsilon = a0 x^a1 (exact)."""
    b, c = 1.0 - a1, 2.0 - a1
    return (b / c) * (E0 ** c - Eo ** c) / (E0 ** b - Eo ** b)


def resonance_yield(p, L, eps, E0, Delta):
    """Eq. (R3) of README.md: the Breit-Wigner integral over the layer, lambda_r and widths constant across it,
    epsilon constant (lab), E_cm = f E_lab: Y = 1e-24 (2 pi/k_r^2) omega gamma / (f eps)
                                               [atan(2 (f E0 - Er)/G) - atan(2 (f (E0 - Delta) - Er)/G)]."""
    f = f_cm(p)
    k2 = 2.0 * mu(p) * L["Er_cm"] / HBARC ** 2
    omega = (2 * L["J"] + 1) / ((2 * p["j1"] + 1) * (2 * p["j2"] + 1))
    gp, ga = L["Gp_eV"] * 1e-6, L["Ga_eV"] * 1e-6
    G = gp + ga
    wg = omega * gp * ga / G
    pref = 1e-24 * (2.0 * math.pi / k2 * 0.01) * wg / (f * eps)
    return pref * (math.atan(2 * (f * E0 - L["Er_cm"]) / G) - math.atan(2 * (f * (E0 - Delta) - L["Er_cm"]) / G)), \
        pref * math.pi


def main():
    truth = json.load(open(sys.argv[1]))
    out = sys.argv[2]
    lv = truth["levels"]
    pN, n12 = lv["p15N"], lv["n12C"]
    L = lv["level"]
    S = truth["sets"]
    nb = blocks(os.path.join(out, "AZUREOut_aa=3_R=3.out"))   # const, const_spread, const_iso
    pb = blocks(os.path.join(out, "AZUREOut_aa=1_R=2.out"))   # thin .. deriv_bare in file order
    names = ["thin", "thin_sigma", "thin_bare", "thick", "thick_inf", "res", "res_q1", "res_q2", "deriv", "deriv_bare",
             "res_bare"]
    if len(nb) != 3 or len(pb) != len(names):
        print(f"FAIL: {len(nb)} neutron and {len(pb)} proton blocks in the output (3 and {len(names)} expected)")
        return 1
    P = dict(zip(names, pb))
    model = lambda rows: [r[3] for r in rows]

    # 1. constant sigma: Y = 1e-24 N sigma(E_bar).  The hard sphere of a = 0.1 fm is constant to its slope
    #    dln(sigma)/dln(E) = -(ka)^2/3 (< 1e-6), which turns the quadrature's offset of the weighted mean energy (the
    #    weights 1/eps_lin are off by about (h/E)^2/4 with h the node spacing, 1.02 Delta/10) into a difference of
    #    s (Delta/E_bar) (h/E_out)^2/4: the tolerance is that plus the output's 11 digits (2e-10).
    a0, a1 = truth["eps_inv"]["a0"], truth["eps_inv"]["a1"]
    N = S["const"]["N"]
    for b, name, w, points, limit in ((nb[0], "const", 1.0, 10, "11 digits + sigma slope x E_bar offset"),
                                      (nb[2], "const_iso", 0.5, 10, "as const; the abundance 0.5 (Weighted)"),
                                      (nb[1], "const_spread", 1.0, 100, "as const + the 0.3 keV Gaussian on the grid")):
        worst, wy, wc, wt = -1.0, 0, 0, 0
        for E0, y in zip(S["const"]["E0"], model(b)):
            Eo = exit_pow(E0, N, a0, a1)
            Eb = mean_energy_pow(E0, Eo, a0, a1)
            c = w * 1e-24 * N * sigma_hard_sphere(n12, f_cm(n12) * Eb)
            ka = math.sqrt(2.0 * mu(n12) * f_cm(n12) * Eb) / HBARC * n12["R"]
            h = 1.02 * (E0 - Eo) / points
            tol = ka * ka / 3 * (E0 - Eo) / Eb * (h / Eo) ** 2 / 4 + 2e-10
            if name == "const_spread":
                tol += 1e-7  # the Gaussian average of sigma/eps on a grid of step 0.2 keV (measured 2.8e-8)
            r = finite(abs(y / c - 1))
            if r - tol >= worst - wt:
                worst, wy, wc, wt = r, y, c, tol
        check(name, wy, wc, worst, wt, limit)

    # 2. thin layer: Y -> 1e-24 N sigma(E0 - Delta/2); the second-order remainder from the bare sigma at +-5 keV
    ap, bp = truth["eps_pow"]["a0"], truth["eps_pow"]["a1"]
    N = S["thin"]["N"]
    bare = model(P["thin_bare"])
    worst = (0.0, 0, 0, 0.0)
    worst_avg = (0.0, 0, 0)
    for i, (E0, y, avg) in enumerate(zip(S["thin"]["E0"], model(P["thin"]), model(P["thin_sigma"]))):
        Eo = exit_pow(E0, N, ap, bp)
        D = E0 - Eo
        Em, Eb = E0 - 0.5 * D, mean_energy_pow(E0, Eo, ap, bp)
        sm, s0, sp = bare[3 * i:3 * i + 3]
        h = 0.005
        d1, d2 = (sp - sm) / (2 * h), (sp - 2 * s0 + sm) / h ** 2
        remainder = (d2 * D * D / 24 + d1 * (Eb - Em)) / s0     # sigma''Delta^2/24 + sigma' (E_bar - E_mid)
        c = 1e-24 * N * s0
        r = finite(abs(y / c - 1 - remainder))
        if r >= worst[0]:
            worst = (r, y, c, remainder)
        ra = finite(abs(y / (1e-24 * N) / avg - 1))
        if ra >= worst_avg[0]:
            worst_avg = (ra, y / (1e-24 * N), avg)
    check("thin: Y vs 1e-24 N sigma(E0-D/2) (+ remainder %.1e)" % worst[3], worst[1], worst[2], worst[0], 1e-7,
          "bare sigma's 11 digits; the third-order remainder")
    check("thin: Y/(1e-24 N) vs the thin-target sigma", worst_avg[1], worst_avg[2], worst_avg[0], 1e-6,
          "two grids (yield +2 %), two E_out (quadrature / Newton-Simpson)")

    # 3. thick layer: Y(E0) = Y_inf(E0) - Y_inf(E_out); and -> Y_inf(E0)
    N = S["thick"]["N"]
    inf = model(P["thick_inf"])
    n = len(S["thick"]["E0"])
    worst, worst_lim = (0.0, 0, 0), (0.0, 0, 0)
    for i, (E0, y) in enumerate(zip(S["thick"]["E0"], model(P["thick"]))):
        c = inf[i] - inf[n + i]
        r = finite(abs(y / c - 1))
        if r >= worst[0]:
            worst = (r, y, c)
        r = finite(abs(y / inf[i] - 1))
        if r >= worst_lim[0]:
            worst_lim = (r, y, inf[i])
    check("thick: Y vs Y_inf(E0) - Y_inf(E_out)", worst[1], worst[2], worst[0], 1e-5,
          "two adaptive grids at the resonance (1.5e-6 at 4x)")
    rows_out.append(("thick: Y vs Y_inf(E0) (the cut below E_out)", worst_lim[2], worst_lim[1], worst_lim[0], None,
                     "Y_inf(E_out)/Y_inf(E0), the yield below E_out", True))

    # 4. narrow resonance, constant epsilon: the arctan form (eq. R3), relative to the plateau
    eps = truth["eps_const"]
    N = S["res"]["N"]
    worst = (0.0, 0, 0)
    plateau = None
    for E0, y in zip(S["res"]["E0"], model(P["res"])):
        c, plateau = resonance_yield(pN, L, eps, E0, N * eps)
        r = finite(abs(y - c) / plateau)
        if r >= worst[0]:
            worst = (r, y, c)
    check("res: Y vs the arctan form (|dY| / plateau %.4e)" % plateau, worst[1], worst[2], worst[0], 3e-4,
          "R-matrix vs Breit-Wigner line shape (1.2e-4); grid (1e-5)")

    # the same layer with its areal density a fit parameter (grid laid out for N0 + 3 sigma = 1.3 N0, cut at E_out),
    # E_0 such that E_out lies on the resonance
    worst = (0.0, 0, 0)
    for E0, y in zip(S["deriv"]["E0"], model(P["deriv"])):
        c, plateau = resonance_yield(pN, L, eps, E0, N * eps)
        r = finite(abs(y - c) / plateau)
        if r >= worst[0]:
            worst = (r, y, c)
    check("deriv: N fitted (grid to 1.3 N), E_out on the resonance", worst[1], worst[2], worst[0], 3e-4,
          "as res")

    # the Breit-Wigner of eq. (R1) against the engine's bare sigma across the resonance (what limits limit 4)
    dev = [finite(abs(y / sigma_bw(pN, L, f_cm(pN) * e) - 1)) for e, y in zip(S["res_bare"]["E0"], model(P["res_bare"]))]
    rows_out.append(("res_bare: engine sigma vs Breit-Wigner, |E-Er| <= 1 G", 0.0, 0.0, max(dev[3:6]), None,
                     "R-matrix vs Breit-Wigner line shape (no tolerance: a measurement)", True))
    rows_out.append(("res_bare: engine sigma vs Breit-Wigner, |E-Er| = 30 G", 0.0, 0.0, max(dev[0], dev[-1]), None,
                     "R-matrix vs Breit-Wigner line shape (no tolerance: a measurement)", True))

    # 5. per unit charge: q = 1 and 2 against the per-particle yield x 1e-6 / (q e)
    idx = [S["res"]["E0"].index(e) for e in S["res_q"]["E0"]]
    yres = model(P["res"])
    for q, b in ((1, P["res_q1"]), (2, P["res_q2"])):
        worst = (0.0, 0, 0)
        for j, y in zip(idx, model(b)):
            c = yres[j] * 1e-6 / (q * E_CHARGE)
            r = finite(abs(y / c - 1))
            if r >= worst[0]:
                worst = (r, y, c)
        check(f"res_q{q}: Y per uC vs Y x 1e-6/({q} e)", worst[1], worst[2], worst[0], 2e-10, "output digits (11)")
    r = max(finite(abs(a / b / 2 - 1)) for a, b in zip(model(P["res_q1"]), model(P["res_q2"])))
    check("res_q1 / res_q2 = 2", 2.0, 2.0 * (1 + r), r, 2e-10, "output digits (11)")

    print(f"{'limit':58s} {'closed form':>17s} {'engine':>17s} {'rel. diff':>9s} {'tol':>8s}  limited by")
    for name, c, e, r, tol, limit, ok in rows_out:
        t = "-" if tol is None else f"{tol:.1e}"
        print(f"{name:58s} {c:17.10e} {e:17.10e} {r:9.2e} {t:>8s}  {limit}{'' if ok else '   FAIL'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
