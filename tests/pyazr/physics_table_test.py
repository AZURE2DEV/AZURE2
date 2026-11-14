#!/usr/bin/env python3
"""pyazr.physics: every row of the binding's physics table (api/src/pybind_physics.cpp) against a reference.

Each row is tested by a ``row_test("<name>", ...)`` (tests/consistency check 41 requires one per row):

- against mpmath (30 digits) where the function is a special function or a closed form of the engine's constants
  (Coulomb functions, phase, Whittaker function, channel functions, kinematics, widths, Coulomb cross sections,
  S-factor conversions);
- against the engine's own use on tests/13N/13N.azr: the session's channel functions (exact Coulomb routine) and its
  Wigner limits (``Parameter.wigner_limit``);
- against the Python formula the row replaces (nds.py), over the range given; the
  largest relative difference per formula is printed in a table (``python-formula`` lines).

The calling convention of the generated functions (broadcasting, scalars, errors="mask", DomainError with the
argument, its value and the domain) is tested at the end.  Run from anywhere:  python3 tests/pyazr/physics_table_test.py
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []
tested = set()
formula_table = []


def row_test(name, what, ok, detail=""):
    tested.add(name)
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}: {what}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(f"{name}: {what}")


def check(what, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {what}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(what)


def rel(a, b, floor=0.0):
    """|a - b| / max(|b|, floor): relative, or absolute below `floor` (a quantity that is 0, as S_0 of a neutral
    channel)."""
    a, b = float(a), float(b)
    scale = max(abs(b), floor)
    return abs(a - b) / scale if scale > 0 else abs(a)


def main():
    try:
        import numpy as np
        import mpmath as mp
        import pyazr
        from pyazr import physics as ph
        from pyazr.sources import nds
    except Exception as err:  # the skip of the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 77
    mp.mp.dps = 30
    U, HC, ALPHA = mp.mpf("931.49410242"), mp.mpf("197.3269804"), 1 / mp.mpf("137.035999084")

    def formula(name, replaced_by, values, reference, where, floor=0.0):
        """Largest relative difference of a Python formula from the row that replaces it (absolute below `floor`),
        recorded, printed and checked: the two are the same formula, so they agree to rounding (1e-12)."""
        same_length = len(values) == len(reference) and len(values) > 0
        d = max(rel(v, r, floor) for v, r in zip(values, reference)) if same_length else float("inf")
        formula_table.append((name, replaced_by, d, where))
        check(f"python-formula  {name} = {replaced_by} to 1e-12 (max rel diff {d:.2e} over {where})", d <= 1e-12,
              f"{len(values)} values against {len(reference)}")
        return d

    # pairs (z1, z2, m1, m2) and their energy windows (MeV) for the channel functions
    pairs = [("p+12C", 1, 6, 1.00727646688, 12.0, 0.05, 20.0), ("a+12C", 2, 6, 4.00150617913, 12.0, 0.3, 20.0),
             ("n+16O", 0, 8, 1.00866491595, 15.99491461957, 0.01, 20.0), ("12C+12C", 6, 6, 12.0, 12.0, 2.0, 20.0)]
    mu_of = lambda m1, m2: mp.mpf(m1) * m2 / (mp.mpf(m1) + m2)

    def mp_eta(z1, z2, m1, m2, e):
        return mp.sqrt(U / 2) * ALPHA * z1 * z2 * mp.sqrt(mu_of(m1, m2) / mp.mpf(e))

    def mp_rho(m1, m2, a, e):
        return mp.sqrt(2 * U) / HC * mp.mpf(a) * mp.sqrt(mu_of(m1, m2) * mp.mpf(e))

    def mp_waves(l, eta, rho):
        F, G = mp.coulombf(l, eta, rho), mp.coulombg(l, eta, rho)
        F1, G1 = mp.coulombf(l + 1, eta, rho), mp.coulombg(l + 1, eta, rho)
        lp = l + 1
        R, S = mp.sqrt(1 + (eta / lp) ** 2), lp / rho + eta / lp  # DLMF 33.4.1
        return F, S * F - R * F1, G, S * G - R * G1

    def mp_channel(l, z1, z2, m1, m2, a, e):
        eta, rho = mp_eta(z1, z2, m1, m2, e), mp_rho(m1, m2, a, e)
        F, Fp, G, Gp = mp_waves(l, eta, rho)
        return rho / (F * F + G * G), rho * (F * Fp + G * Gp) / (F * F + G * G), mp.atan2(F, G)

    # ------------------------------------------------------------------ kinematics
    print("kinematics")
    e = np.geomspace(0.01, 20.0, 40)
    for label, z1, z2, m1, m2, _, _ in pairs:
        eta = ph.sommerfeld_eta(z1, z2, m1, m2, e)
        ref = [mp_eta(z1, z2, m1, m2, x) for x in e]
        row_test("sommerfeld_eta", f"{label} = Z1 Z2 alpha (mu c^2 / 2E)^1/2 (mpmath) to 1e-15",
                 all(rel(eta[i], ref[i]) <= 1e-15 or ref[i] == 0 and eta[i] == 0 for i in range(e.size)))
        k = ph.wave_number(m1, m2, e)
        row_test("wave_number", f"{label} = (2 mu c^2 E)^1/2 / hbar c (mpmath) to 1e-15",
                 all(rel(k[i], mp_rho(m1, m2, 1, e[i])) <= 1e-15 for i in range(e.size)))
    # relativistic: sqrt(s) - m1 - m2 with s = (m1 + m2)^2 + 2 m2 E_lab, masses in MeV
    elab = np.geomspace(0.05, 50.0, 30)
    ecm = ph.lab_to_cm_energy(4.00150617913, 12.0, elab)
    row_test("lab_to_cm_energy", "a on 12C = E_lab m2 / (m1 + m2) to 1e-15",
             all(rel(ecm[i], elab[i] * 12.0 / 16.00150617913) <= 1e-15 for i in range(elab.size)))
    M1, M2 = 4.00150617913 * U, 12 * U
    ecm_r = ph.lab_to_cm_energy(4.00150617913, 12.0, elab, relativistic=True)
    row_test("lab_to_cm_energy", "relativistic = sqrt(s) - m1 - m2 (mpmath) to 1e-9",
             all(rel(ecm_r[i], mp.sqrt((M1 + M2) ** 2 + 2 * M2 * elab[i]) - M1 - M2) <= 1e-9 for i in range(elab.size)))
    back = ph.cm_to_lab_energy(4.00150617913, 12.0, ecm)
    row_test("cm_to_lab_energy", "inverts lab_to_cm_energy to 1e-15",
             all(rel(back[i], elab[i]) <= 1e-15 for i in range(elab.size)))
    back_r = ph.cm_to_lab_energy(4.00150617913, 12.0, ecm_r, relativistic=True)
    row_test("cm_to_lab_energy", "relativistic: inverts lab_to_cm_energy to 1e-9",
             all(rel(back_r[i], elab[i]) <= 1e-9 for i in range(elab.size)))
    formula("nds.cm_to_lab", "cm_to_lab_energy", ph.cm_to_lab_energy(1.00727646688, 12.0, e),
            [nds.cm_to_lab(x, 12.0, 1.00727646688) for x in e], "p+12C, E_cm = 0.01-20 MeV")

    # angles: forward map theta_cm -> theta_lab written here (non-relativistic, gamma = V / v3'), then the rows invert it
    def mp_forward(m1, m2, m3, m4, q, elab, tcm):
        m1, m2, m3, m4, q, elab = map(mp.mpf, (m1, m2, m3, m4, q, elab))
        V = mp.sqrt(2 * m1 * elab) / (m1 + m2)                       # c.m. velocity (units with masses in u)
        e_after = elab * m2 / (m1 + m2) + q                          # c.m. kinetic energy after
        v3 = mp.sqrt(2 * (e_after * m4 / (m3 + m4)) / m3)            # speed of 3 in the c.m.
        t = mp.radians(tcm)
        tl = mp.atan2(v3 * mp.sin(t), v3 * mp.cos(t) + V)
        # dOmega_lab / dOmega_cm = d cos(theta_lab) / d cos(theta_cm)
        c = mp.diff(lambda x: mp.cos(mp.atan2(v3 * mp.sin(mp.acos(x)), v3 * x + V)), mp.cos(t))
        return mp.degrees(tl), c

    tcm = np.linspace(5.0, 175.0, 18)
    for (m1, m2) in ((1.00727646688, 12.0), (4.00150617913, 15.99491461957), (1.0, 1.0 + 1e-9)):
        lab = [float(mp_forward(m1, m2, m1, m2, 0, 1.0, t)[0]) for t in tcm]
        jac = [mp_forward(m1, m2, m1, m2, 0, 1.0, t)[1] for t in tcm]
        r = ph.elastic_lab_to_cm_angle(m1, m2, lab)
        row_test("elastic_lab_to_cm_angle", f"m1/m2 = {m1 / m2:.3f}: theta_cm and J from the forward map (mpmath) "
                 "to 1e-9", all(abs(r["theta_cm"][i] - tcm[i]) <= 1e-9 * 180 and rel(r["jacobian"][i], jac[i]) <= 1e-8
                                 for i in range(tcm.size)), (r, jac[:3]))
        ret = nds.elastic_lab_to_cm(np.asarray(lab), m2, m1)
        formula("nds.elastic_lab_to_cm (theta_cm)", "elastic_lab_to_cm_angle", r["theta_cm"], ret[0],
                f"m1/m2 = {m1 / m2:.3f}, theta_lab {min(lab):.1f}-{max(lab):.1f} deg")
        formula("nds.elastic_lab_to_cm (J)", "elastic_lab_to_cm_angle", r["jacobian"], ret[1],
                f"m1/m2 = {m1 / m2:.3f}, theta_lab {min(lab):.1f}-{max(lab):.1f} deg")
    edge = ph.elastic_lab_to_cm_angle(4.0, 12.0, [0.0, 180.0])
    row_test("elastic_lab_to_cm_angle", "theta_lab = 0, 180 deg (m1/m2 = 1/3): theta_cm 0, 180, J = 1/(1 +- 1/3)^2",
             abs(edge["theta_cm"][0]) < 1e-12 and abs(edge["theta_cm"][1] - 180) < 1e-12
             and rel(edge["jacobian"][0], 0.5625) < 1e-12 and rel(edge["jacobian"][1], 2.25) < 1e-12, edge)
    _, st = ph.elastic_lab_to_cm_angle(12.0, 4.0, [10.0, 30.0], errors="mask")
    row_test("elastic_lab_to_cm_angle", "a beam heavier than the target: status -2 beyond asin(m2/m1) = 19.47 deg",
             list(st) == [0, -2], st)
    # 12C(a,p)15N: Q = -4.966 MeV; 7Li(p,a)4He: Q = 17.347 MeV
    # (the Marion-Young map assumes m1 + m2 = m3 + m4; the masses here keep it)
    for (m1, m2, m3, m4, q, el) in ((4.00150617913, 12.0, 1.00727646688, 16.00150617913 - 1.00727646688, -4.966, 10.0),
                                     (1.00727646688, 7.01435, 4.00150617913, 8.02162646688 - 4.00150617913, 17.347,
                                      0.5)):
        fw = [mp_forward(m1, m2, m3, m4, q, el, t) for t in tcm]
        r = ph.reaction_lab_to_cm_angle(m1, m2, m3, m4, q, el, [float(f[0]) for f in fw])
        row_test("reaction_lab_to_cm_angle", f"Q = {q} MeV, E_lab = {el} MeV: theta_cm and J from the forward map "
                 "(mpmath) to 1e-8", all(abs(r["theta_cm"][i] - tcm[i]) <= 1e-8 * 180 and rel(r["jacobian"][i],
                                                                                              fw[i][1]) <= 1e-7
                                         for i in range(tcm.size)), r)
    _, st = ph.reaction_lab_to_cm_angle(4.0015, 12.0, 1.00728, 15.0, -4.966, [1.0, 10.0], 30.0, errors="mask")
    row_test("reaction_lab_to_cm_angle", "below the exit threshold: status -2", list(st) == [-2, 0], st)

    # ------------------------------------------------------------------ Coulomb functions
    print("Coulomb functions")
    grid = [(l, eta, rho) for l in (0, 1, 3, 6) for eta in (0.0, 0.5, 2.0, 10.0) for rho in (0.2, 1.0, 5.0, 20.0)]
    L, ETA, RHO = (np.array(c, float) for c in zip(*grid))
    ref = [mp_waves(l, mp.mpf(eta), mp.mpf(rho)) for l, eta, rho in grid]
    # where each algorithm is accurate (physics/coulomb/README.md, "Accuracy"): coulsm and coul give a wrong F for
    # eta = 0, l >= 2 at small rho (F enters P, S, phi only through F/G), coulfg is JWKB inside the turning point
    turning = ETA + np.sqrt(ETA ** 2 + L * (L + 1))
    regions = {"coulsm": (ETA > 0, "eta > 0 (eta = 0: G, G' only)"), "coul": (ETA > 0, "eta > 0 (eta = 0: G, G' only)"),
               "coulfg": (RHO >= turning, "rho >= the turning point"), "coulcc": (ETA >= 0, "everywhere")}
    for method, tol in (("coulsm", 1e-9), ("coul", 1e-9), ("coulfg", 1e-9), ("coulcc", 1e-9)):
        w, st = ph.coulomb_waves(L, ETA, RHO, method=method, errors="mask")
        region, text = regions[method]
        worst = max(max(rel(w[k][i], ref[i][j]) for j, k in enumerate(("F", "Fp", "G", "Gp"))
                        if region[i] or (method in ("coulsm", "coul") and k in ("G", "Gp")))
                    for i in range(len(grid)) if region[i] or method in ("coulsm", "coul"))
        row_test("coulomb_waves", f"method {method}: F, F', G, G' against mpmath on l = 0-6, eta = 0-10, rho = 0.2-20, "
                 f"{text}: largest rel diff {worst:.1e} (<= {tol:g}), {int(np.sum(st != 0))} not computed",
                 worst <= tol and np.all(st == 0), worst)
    # the recorded COULCC defect (docs/upgrade/HANDOFF.md): l = 1, eta = 0.493, rho = 8.55
    cc = ph.coulomb_waves(1, 0.493, 8.55, method="coulcc")
    sm = ph.coulomb_waves(1, 0.493, 8.55)
    rf = mp_waves(1, mp.mpf("0.493"), mp.mpf("8.55"))
    print(f"  note  COULCC at l = 1, eta = 0.493, rho = 8.55: F = {cc['F']:.15g}, G = {cc['G']:.15g}; coulsm F = "
          f"{sm['F']:.15g}, G = {sm['G']:.15g}; mpmath F = {float(rf[0]):.15g}, G = {float(rf[2]):.15g}")
    sig = ph.coulomb_phase(L, ETA)
    row_test("coulomb_phase", "sigma_l = Im ln Gamma(l + 1 + i eta) (mpmath) to 1e-12 absolute",
             all(abs(sig[i] - float(mp.im(mp.loggamma(l + 1 + 1j * mp.mpf(eta))))) <= 1e-12 for i, (l, eta, _) in
                 enumerate(grid)))
    row_test("coulomb_phase", "equals the sigma of coulomb_waves to 1e-14 absolute",
             np.all(np.abs(sig - w["sigma"]) <= 1e-14))

    wg = [(l, eb, rho) for l in (0, 1, 2, 4) for eb in (0.0, 0.3, 2.0, 20.0) for rho in (0.05, 0.5, 2.0, 10.0)]
    Lw, Ew, Rw = (np.array(c, float) for c in zip(*wg))
    ww = ph.whittaker_wave(Lw, Ew, Rw)
    refw = [mp.whitw(-mp.mpf(eb), l + mp.mpf(1) / 2, 2 * mp.mpf(rho)) for l, eb, rho in wg]
    row_test("whittaker_wave", "W against mpmath whitw on l = 0-4, eta_b = 0-20, rho = 0.05-10 to 1e-11",
             all(rel(ww["W"][i], refw[i]) <= 1e-11 for i in range(len(wg))))
    row_test("whittaker_wave", "rho W'/W against mpmath to 1e-9",
             all(rel(ww["rho_dlog"][i], Rw[i] * mp.diff(lambda r: mp.whitw(-mp.mpf(Ew[i]), Lw[i] + 0.5, 2 * r),
                                                           mp.mpf(Rw[i])) / refw[i]) <= 1e-9 for i in range(len(wg))))

    # ------------------------------------------------------------------ channel functions
    print("channel functions")
    worst = {"P": 0.0, "S": 0.0, "phi": 0.0}
    for label, z1, z2, m1, m2, e_lo, e_hi in pairs:
        energies = np.geomspace(e_lo, e_hi, 9)
        for l in range(4):
            for a in (3.5, 5.5):
                P = ph.penetrability(l, z1, z2, m1, m2, a, energies)
                S = ph.shift(l, z1, z2, m1, m2, a, energies)
                phi = ph.hard_sphere_phase(l, z1, z2, m1, m2, a, energies)
                for i, x in enumerate(energies):
                    rP, rS, rphi = mp_channel(l, z1, z2, m1, m2, a, x)
                    worst["P"] = max(worst["P"], rel(P[i], rP))
                    worst["S"] = max(worst["S"], abs(S[i] - float(rS)))
                    worst["phi"] = max(worst["phi"], abs(phi[i] - float(rphi)))
    rng = "p+12C 0.05-20, a+12C 0.3-20, n+16O 0.01-20, 12C+12C 2-20 MeV; l = 0-3; a = 3.5, 5.5 fm"
    row_test("penetrability", f"against mpmath, {rng}: largest rel diff {worst['P']:.1e} (<= 1e-9)", worst["P"] <= 1e-9)
    # S = rho (F F' + G G') / (F^2 + G^2) cancels: with F, G good to 5e-10 (coulsm), S is good to a few 1e-9 absolute
    # (2.6e-8 relative at S = -0.06, a + 12C, l = 1, 20 MeV, a = 5.5 fm)
    row_test("shift", f"against mpmath, same range: largest abs diff {worst['S']:.1e} (<= 5e-9)", worst["S"] <= 5e-9)
    row_test("hard_sphere_phase", f"atan2(F, G) against mpmath, same range: largest abs diff {worst['phi']:.1e} rad "
             "(<= 1e-9, the precision of coulsm)", worst["phi"] <= 1e-9)

    dS_worst = 0.0
    for label, z1, z2, m1, m2, e_lo, e_hi in pairs:
        for l in (0, 2):
            for x in np.geomspace(max(e_lo, 0.1), e_hi, 4):
                d = ph.shift_derivative(l, z1, z2, m1, m2, 5.0, x)
                ref = mp.diff(lambda E: mp_channel(l, z1, z2, m1, m2, 5.0, E)[1], mp.mpf(x))
                dS_worst = max(dS_worst, rel(d, ref, 1e-6))
    # the engine's Richardson difference with h = 1e-3 E carries the few-1e-9 noise of S: 9.5e-6 relative at p + 12C,
    # l = 0, 3.42 MeV (the Python formula, mpmath with h = 1e-3 MeV, is right there to 1e-11; docs/upgrade/HANDOFF.md)
    row_test("shift_derivative", f"against mpmath's derivative of S, four pairs, l = 0, 2, a = 5 fm, E = 0.1-20 MeV: "
             f"largest rel diff {dS_worst:.1e} (<= 2e-5, the engine's difference)", dS_worst <= 2e-5)

    # M_l(p) = (B - 1) j_l(pa) - pa j_l'(pa), j_l by mpmath (spherical Bessel through besselj of order l + 1/2)
    jl = lambda l, x: mp.sqrt(mp.pi / (2 * x)) * mp.besselj(l + mp.mpf(1) / 2, x)
    worst_M = 0.0
    for l in range(4):
        for p in (0.05, 0.3, 1.0, 2.5):
            for B in (-2.0, 0.0, 1.5):
                x = mp.mpf(p) * 5
                ref = (B - 1) * jl(l, x) - x * mp.diff(lambda y: jl(l, y), x)
                worst_M = max(worst_M, abs(ph.thm_vertex_function(l, p, 5.0, B) - float(ref)))
    row_test("thm_vertex_function", f"(B - 1) j_l(pa) - pa j_l'(pa) (mpmath), l = 0-3, p = 0.05-2.5 fm^-1, a = 5 fm, "
             f"B = -2-1.5: largest abs diff {worst_M:.1e} (<= 1e-9)", worst_M <= 1e-9)
    # ------------------------------------------------------------------ widths
    print("widths")
    for label, z1, z2, m1, m2, _, _ in pairs:
        gW = ph.wigner_limit(m1, m2, np.array([3.0, 5.0]))
        ref = [3 * HC ** 2 / (2 * mu_of(m1, m2) * U * a * a) for a in (3, 5)]
        row_test("wigner_limit", f"{label} = 3 (hbar c)^2 / (2 mu c^2 a^2) (mpmath) to 1e-15",
                 all(rel(gW[i], ref[i]) <= 1e-15 for i in range(2)))
        Gsp = ph.single_particle_width(1, z1, z2, m1, m2, 5.0, np.array([0.5, 3.0]) + (2.0 if z1 * z2 > 30 else 0))
        refs = [2 * mp_channel(1, z1, z2, m1, m2, 5.0, x)[0] * 3 * HC ** 2 / (2 * mu_of(m1, m2) * U * 25)
                for x in np.array([0.5, 3.0]) + (2.0 if z1 * z2 > 30 else 0)]
        row_test("single_particle_width", f"{label}, l = 1, a = 5 fm = 2 P gamma_W^2 (mpmath) to 1e-9",
                 all(rel(Gsp[i], refs[i]) <= 1e-9 for i in range(2)), (Gsp, refs))
    eg = np.geomspace(0.01, 20.0, 15)
    for L_ in (1, 2, 3):
        P = ph.photon_penetrability(eg, L_)
        row_test("photon_penetrability", f"L = {L_}: (E / hbar c)^(2L+1) (mpmath) to 1e-14",
                 all(rel(P[i], (mp.mpf(eg[i]) / HC) ** (2 * L_ + 1)) <= 1e-14 for i in range(eg.size)))

    # ------------------------------------------------------------------ the engine's own use: 13N
    print("the engine's own use (tests/13N/13N.azr, exact Coulomb routine)")
    azr = os.path.join(ROOT, "tests", "13N", "13N.azr")
    s = pyazr.Session(azr, use_coulomb_table=False)
    try:
        for pair in s.pairs:
            if pair.ptype != 0:
                continue
            energies = np.array([0.1, 0.4, 1.0, 2.5])
            for l in range(3):
                cf = s.channel_functions(pair.key, energies, l)
                args = (l, pair.Z1, pair.Z2, pair.M1, pair.M2, pair.channel_radius, energies)
                row_test("penetrability", f"pair {pair.key}, l = {l}: = the session's P bit for bit",
                         np.array_equal(ph.penetrability(*args), cf["P"]), (ph.penetrability(*args), cf["P"]))
                row_test("shift", f"pair {pair.key}, l = {l}: = the session's S bit for bit",
                         np.array_equal(ph.shift(*args), cf["S"]))
                row_test("hard_sphere_phase", f"pair {pair.key}, l = {l}: = -(the session's delta_hs) to 1e-15 abs",
                         np.all(np.abs(ph.hard_sphere_phase(*args) + cf["delta_hs"]) <= 1e-15))
            for p in s.parameters.of_kind("width"):
                if p.wigner_limit and p.pair == pair.number:
                    row_test("wigner_limit", f"pair {pair.key}: = the engine's Parameter.wigner_limit bit for bit",
                             ph.wigner_limit(pair.M1, pair.M2, pair.channel_radius) == p.wigner_limit,
                             (ph.wigner_limit(pair.M1, pair.M2, pair.channel_radius), p.wigner_limit))
                    break
        cf = s.channel_functions(1, np.array([[0.5, -0.2], [1.0, 2.0]]), 0)
        check("channel_functions keeps the shape of e_cm and marks the closed channel",
              cf["P"].shape == (2, 2) and cf["open"].tolist() == [[True, False], [True, True]], cf["open"])
    finally:
        s.close()

    # ------------------------------------------------------------------ resonance quantities, Coulomb scattering
    print("resonance quantities, Coulomb scattering, S factor")
    row_test("statistical_factor", "7/4 (p + 15N, 3-), 10 (12C + 12C, 2+), 2 (a + a, 0+) by hand",
             list(ph.statistical_factor([6, 4, 0], [1, 0, 0], [1, 0, 0], [False, True, True])) == [1.75, 10.0, 2.0])
    r = ph.resonance_strength(1.75, 1.702, 12.92, 14.622)
    row_test("resonance_strength", "T2 level 5: 1.75 x 1.702 x 12.92 / 14.622 = 2.6318028 keV, d/dGamma = -value/Gamma",
             rel(r["value"], 1.75 * 1.702 * 12.92 / 14.622) <= 1e-15 and rel(r["d_gamma"], -r["value"] / 14.622) <= 1e-14)
    b = ph.branching_ratio(12.92, 14.622)
    row_test("branching_ratio", "12.92 / 14.622 and its derivatives",
             rel(b["value"], 12.92 / 14.622) <= 1e-15 and rel(b["d_gamma_b"], 1 / 14.622) <= 1e-15
             and rel(b["d_gamma"], -12.92 / 14.622 ** 2) <= 1e-15)
    _, st = ph.branching_ratio([1.0, 1.0], [2.0, 0.0], errors="mask")
    row_test("branching_ratio", "Gamma = 0: status 2 (the argument gamma)", list(st) == [0, 2], st)
    _, st = ph.resonance_strength(1.0, 1.0, 1.0, 0.0, errors="mask")
    row_test("resonance_strength", "Gamma = 0: status 4 (the argument gamma)", st == 4, st)

    def mp_rutherford(z1, z2, m1, m2, e, th):
        eta, k = mp_eta(z1, z2, m1, m2, e), mp_rho(m1, m2, 1, e)
        return (eta / (2 * k * mp.sin(mp.radians(th) / 2) ** 2)) ** 2 / 100  # fm^2 -> b
    e3, th3 = np.array([1.5, 2.25, 3.75]), np.array([20.0, 90.0, 160.0])
    ru = ph.rutherford_cross_section(2, 6, 4.00260325, 12.0, e3, th3)
    row_test("rutherford_cross_section", "(eta / 2k sin^2(theta/2))^2 (mpmath), a + 12C, to 1e-13; d/dE = -2 sigma/E",
             all(rel(ru["value"][i], mp_rutherford(2, 6, 4.00260325, 12.0, e3[i], th3[i])) <= 1e-13
                 and rel(ru["d_dE"][i], -2 * ru["value"][i] / e3[i]) <= 1e-13 for i in range(3)))

    def mp_mott(z, m, e, th):
        eta, k = mp_eta(z, z, m, m, e), mp_rho(m, m, 1, e)
        t = mp.radians(th)
        f = lambda x: -eta / (2 * k * mp.sin(x / 2) ** 2) * mp.exp(-1j * eta * mp.log(mp.sin(x / 2) ** 2))
        return abs(f(t) + f(mp.pi - t)) ** 2 / 100  # spin 0: symmetric
    mo = ph.mott_cross_section(2, 4.00260325, 0, e3, th3)
    row_test("mott_cross_section", "|f_C(theta) + f_C(pi - theta)|^2 (mpmath), a + a, to 1e-12",
             all(rel(mo["value"][i], mp_mott(2, 4.00260325, e3[i], th3[i])) <= 1e-12 for i in range(3)), mo)
    eS = np.array([0.05, 0.3, 1.2, 4.0])
    conv = ph.sfactor_conversion(1, 6, 1.00728, 12.0, eS)
    row_test("sfactor_conversion", "E exp(2 pi eta) (mpmath), p + 12C, to 1e-13",
             all(rel(conv[i], eS[i] * mp.exp(2 * mp.pi * mp_eta(1, 6, 1.00728, 12.0, eS[i]))) <= 1e-13
                 for i in range(4)))
    _, st = ph.sfactor_conversion(0, 6, 1.00866, 12.0, 1.0, errors="mask")
    row_test("sfactor_conversion", "an uncharged pair: status -2", st in (0, -2), st)
    mconv = ph.modified_sfactor_conversion(6, 6, 12.0, 12.0, 0.46, eS + 2)
    row_test("modified_sfactor_conversion", "E exp(2 pi eta + g E) (mpmath), 12C + 12C, g = 0.46, to 1e-13",
             all(rel(mconv[i], (eS[i] + 2) * mp.exp(2 * mp.pi * mp_eta(6, 6, 12.0, 12.0, eS[i] + 2) + mp.mpf("0.46") *
                                                    (eS[i] + 2))) <= 1e-13 for i in range(4)))

    # ------------------------------------------------------------------ the calling convention
    print("the calling convention")
    P = ph.penetrability(0, 1, 6, 1.00728, 12.0, 5.0, 1.0)
    check("scalars in: a float out", isinstance(P, float), type(P))
    grid2 = ph.penetrability(np.arange(3)[:, None], 1, 6, 1.00728, 12.0, 5.0, np.array([0.5, 1.0, 2.0, 4.0]))
    check("arrays broadcast: l (3, 1) x e_cm (4,) gives (3, 4)", grid2.shape == (3, 4), grid2.shape)
    check("a broadcast element = the scalar call", grid2[1, 2] == ph.penetrability(1, 1, 6, 1.00728, 12.0, 5.0, 2.0))
    v, st = ph.penetrability(0, 1, 6, 1.00728, 12.0, 5.0, [1.0, -1.0, np.nan], errors="mask")
    check("errors='mask': NaN and the status of the argument (e_cm is argument 7)",
          v[0] > 0 and np.isnan(v[1]) and np.isnan(v[2]) and st.tolist() == [0, 7, 7], (v, st))
    check("status_text names the argument and its domain",
          ph.status_text("penetrability", 7) == "e_cm outside its domain (finite > 0)",
          ph.status_text("penetrability", 7))
    try:
        ph.penetrability(0, 1, 6, 1.00728, 12.0, 5.0, [1.0, -1.0])
        check("errors='raise': DomainError", False)
    except pyazr.DomainError as err:
        msg = str(err)
        check("errors='raise': DomainError (a ValueError) naming the argument, value, element and domain",
              isinstance(err, ValueError) and all(t in msg for t in ("e_cm = -1.0 MeV", "element 1", "e_cm finite > 0",
                                                                     "help(pyazr.physics.penetrability)")), msg)
    try:
        ph.penetrability(0.5, 1, 6, 1.00728, 12.0, 5.0, 1.0)
        check("a half-integer l: DomainError", False)
    except pyazr.DomainError as err:
        check("a half-integer l: DomainError naming l and 'integer in [0, 10000]'",
              "l = 0.5" in str(err) and "integer in [0, 10000]" in str(err), str(err))
    try:
        ph.coulomb_waves(0, 1.0, 1.0, method="gsl")
        check("an unknown method: DomainError", False)
    except pyazr.DomainError as err:
        check("an unknown method: DomainError listing the methods", "coulsm, coul, coulfg, coulcc" in str(err), str(err))
    try:
        ph.penetrability(0, 1, 6, 1.00728, 12.0, 5.0, "one")
        check("a string energy: TypeError", False)
    except TypeError as err:
        check("a string energy: TypeError naming e_cm", "e_cm" in str(err), str(err))
    rows = pyazr.kinds()["physics"]
    check("pyazr.kinds()['physics'] = the table, one function per row", set(rows) == set(ph.functions) == set(ph.table()))
    check("every generated function has the row's signature and a docstring with units",
          all(list(__import__("inspect").signature(ph.functions[n]).parameters)[:-1] ==
              [a["name"] for a in rows[n]["arguments"]] and "Returns" in ph.functions[n].__doc__ for n in rows))
    import importlib.util
    spec = importlib.util.spec_from_file_location("gen_physics_docs", os.path.join(ROOT, "tools", "gen_physics_docs.py"))
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    with open(gen.PAGE) as f:
        check("docs/source/user_guide/physics_functions.rst is what tools/gen_physics_docs.py makes of the rows",
              f.read() == gen.generate(rows), "run tools/gen_physics_docs.py")
    import json
    check("the rows are strict JSON (no NaN or infinity: the web UI's /api/physics)",
          json.loads(json.dumps(rows, allow_nan=False)) == rows)
    # what the web UI's Physics tab needs of a row to draw a first plot without typing (web_ui.md section 6.6)
    def inside(a, x):
        lo, hi = a["lower"], a["upper"]
        return ((lo is None or (x > lo if a["lower_open"] else x >= lo))
                and (hi is None or (x < hi if a["upper_open"] else x <= hi)))
    bad = []
    for n, r in rows.items():
        names = {a["name"]: a for a in r["arguments"]}
        if not r["group"] or r["axis"] not in names or names[r["axis"]]["sweep"] is None:
            bad.append(f"{n}: group {r['group']!r}, axis {r['axis']!r} without a sweep")
        for a in r["arguments"]:
            if a["default"] is None and a["example"] is None and a["name"] != r["axis"]:
                bad.append(f"{n}.{a['name']}: neither a default nor an example")
            if a["example"] is not None and not inside(a, a["example"]):
                bad.append(f"{n}.{a['name']}: example {a['example']} outside {a['domain']}")
            w = a["sweep"]
            if w is not None and not (w["lower"] < w["upper"] and w["spacing"] in ("linear", "log")
                                      and inside(a, w["lower"]) and inside(a, w["upper"])
                                      and (w["spacing"] == "linear" or w["lower"] > 0)):
                bad.append(f"{n}.{a['name']}: sweep {w} not inside {a['domain']}")
    check("every row has a group, its axis a sweep, every argument a default or an example, all inside the domain",
          not bad, bad)
    missing = sorted(set(rows) - tested)
    check("every row has a row_test here (check 41 asks for one in some test)", not missing, missing)

    print("\nlargest relative difference of each Python formula from the row that replaces it:")
    for name, by, d, where in formula_table:
        print(f"  {name:42s} {by:28s} {d:9.2e}  {where}")
    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
