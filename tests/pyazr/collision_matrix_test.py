#!/usr/bin/env python3
"""pyazr.physics.collision_matrix: the level matrix and the collision matrix of one J^pi block at one energy from
arrays, without a session or a project (docs/upgrade/physics_functions.md; physics/CollisionMatrix.h).

Checked against values that do not come from the function:

- A is the inverse of the level matrix written here in numpy from its definition, and U follows from A by the
  formula of the docstring (two levels, two channels; the proton channel's P, S, phi and Coulomb phase are the rows
  of pyazr.physics for p + 12C);
- U is unitary and symmetric for real parameters, T = exp(2i omega) - U;
- one level in one channel: U = exp(2i (omega - phi + atan(P gamma^2 / (E_l - E - gamma^2 (S - B)))));
- one level in two channels with B = S: |U_12|^2 is the Breit-Wigner form with Gamma_c = 2 P_c gamma_c^2;
- no level: hard-sphere scattering; a closed channel couples to nothing;
- Brune's parameters (PRC 66, 044611, Eq. (33)) and Park's (PRC 104, 064612, Eqs. (21), (22)): A is the inverse of
  the alternative level matrix written here in numpy; and the three parametrizations describe one physics: the
  formal levels of a block are taken to Brune's and Park's parameters by numpy (Brune's eigenvalue problem with
  shift functions written here), and U is the same in all three at every energy;
- two levels closer than 1e-6 MeV: the Taylor form against the divided difference;
- one Park level: A in closed form with J = 1 - sum gamma^2 dS/dE;
- a photon channel: its width is not in the level matrix, |U_p,gamma|^2 is the one-level capture form; the
  Reich-Moore term adds 2 gamma^2 to the total width;
- the errors name the argument.

Run from anywhere:  python3 tests/pyazr/collision_matrix_test.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(what, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {what}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(what)


def raises(kind, text, call):
    try:
        call()
    except kind as err:
        return text in str(err), str(err)
    except Exception as err:  # another class: reported
        return False, f"{type(err).__name__}: {err}"
    return False, "nothing raised"


def main():
    try:
        import numpy as np
        import pyazr
        from pyazr import physics as ph
    except Exception as err:  # the skip of the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 77

    # The channels: p + 12C, l = 0, a = 3.4 fm at 0.4 MeV from the rows; a second channel with made-up functions.
    pair = dict(z1=1, z2=6, m1_u=1.00782503, m2_u=12.0)
    e_cm, radius = 0.4, 3.4
    eta = ph.sommerfeld_eta(e_cm=e_cm, **pair)
    P = np.array([ph.penetrability(0, radius=radius, e_cm=e_cm, **pair), 0.3])
    S = np.array([ph.shift(0, radius=radius, e_cm=e_cm, **pair), -0.1])
    phi = np.array([ph.hard_sphere_phase(0, radius=radius, e_cm=e_cm, **pair), 0.2])
    omega = np.array([ph.coulomb_phase(0, eta) - ph.coulomb_phase(0, eta), 0.35])
    B = np.array([-0.4, 0.05])
    E = 2.35  # on the scale of the level energies
    levels = np.array([2.37, 3.5])
    gamma = np.array([[0.61, 0.12], [-0.83, 0.47]])

    print("two levels, two channels against numpy")
    out = ph.collision_matrix(levels, gamma, E, penetrability=P, shift=S, boundary=B, hard_sphere_phase=phi,
                              coulomb_phase=omega)
    A, U, T = out["A"], out["U"], out["T"]
    L = S - B + 1j * P
    M = np.diag(levels - E).astype(complex) - (gamma * L) @ gamma.T
    d = np.abs(A @ M - np.eye(2)).max()
    check("A is the inverse of (E_l - E) delta - sum_c gamma gamma L_c", d < 1e-12, f"largest difference {d:.2e}")
    Omega = np.exp(1j * (omega - phi))
    want = np.outer(Omega, Omega) * (np.eye(2) + 2j * np.sqrt(np.outer(P, P)) * (gamma.T @ np.linalg.inv(M) @ gamma))
    d = np.abs(U - want).max()
    check("U = Omega Omega' (1 + 2i sqrt(P P') gamma^T A gamma')", d < 1e-12, f"largest difference {d:.2e}")
    d = np.abs(U @ U.conj().T - np.eye(2)).max()
    check("U is unitary for real widths and boundary values", d < 1e-12, f"largest difference {d:.2e}")
    d = np.abs(U - U.T).max()
    check("U is symmetric", d < 1e-13, f"largest difference {d:.2e}")
    d = np.abs(T - (np.diag(np.exp(2j * omega)) - U)).max()
    check("T = exp(2i omega) delta - U", d < 1e-14, f"largest difference {d:.2e}")
    check("the shapes and types: A (2, 2), U (2, 2), T (2, 2), complex; no other key",
          sorted(out) == ["A", "T", "U"] and all(out[k].shape == (2, 2) and out[k].dtype == complex for k in out),
          f"{ {k: (v.shape, v.dtype) for k, v in out.items()} }")

    print("one level")
    g, e1 = 0.61, 2.37
    one = ph.collision_matrix([e1], [[g]], E, penetrability=P[:1], shift=S[:1], boundary=B[:1],
                              hard_sphere_phase=phi[:1], coulomb_phase=[0.2])
    delta = -phi[0] + np.arctan2(P[0] * g * g, e1 - E - g * g * (S[0] - B[0]))
    d = abs(one["U"][0, 0] - np.exp(2j * (0.2 + delta)))
    check("one channel: U = exp(2i (omega - phi + atan(P g^2 / (E_l - E - g^2 (S - B)))))", d < 1e-13, f"{d:.2e}")
    neutral = ph.collision_matrix([e1], [[g]], E, penetrability=P[:1], shift=S[:1], boundary=B[:1],
                                  hard_sphere_phase=phi[:1])
    d = abs(neutral["U"][0, 0] - np.exp(2j * delta))
    check("coulomb_phase not given is zero", d < 1e-13, f"{d:.2e}")
    g2 = np.array([[0.61, 0.12]])
    bw = ph.collision_matrix([e1], g2, E, penetrability=P, shift=S, boundary=S, hard_sphere_phase=phi)
    width = 2.0 * P * g2[0] ** 2
    want = width[0] * width[1] / ((e1 - E) ** 2 + 0.25 * width.sum() ** 2)
    d = abs(abs(bw["U"][0, 1]) ** 2 - want) / want
    check("two channels, B = S: |U_12|^2 = Gamma_1 Gamma_2 / ((E_l - E)^2 + Gamma^2 / 4)", d < 1e-12, f"{d:.2e}")

    print("limits")
    none = ph.collision_matrix(np.zeros(0), np.zeros((0, 2)), E, penetrability=P, shift=S, boundary=B,
                               hard_sphere_phase=phi, coulomb_phase=omega)
    d = np.abs(none["U"] - np.diag(Omega ** 2)).max()
    check("no level: U = exp(2i (omega - phi)) on the diagonal (hard sphere), A of shape (0, 0)",
          d < 1e-15 and none["A"].shape == (0, 0), f"{d:.2e}, A {none['A'].shape}")
    closed = ph.collision_matrix(levels, gamma, E, penetrability=[P[0], 0.0], shift=S, boundary=B,
                                 hard_sphere_phase=phi, coulomb_phase=omega)
    Uc = closed["U"]
    check("a closed channel (P = 0) couples to nothing and the open one is unitary by itself",
          abs(Uc[0, 1]) == 0.0 and abs(Uc[1, 0]) == 0.0 and abs(abs(Uc[0, 0]) - 1.0) < 1e-13,
          f"U = {Uc}")

    print("Brune's and Park's parameters: three levels, two channels")
    # Shift functions of the two channels on the scale of the level energies, written here (any smooth real S_c(E)
    # serves: the equivalence of the parametrizations is algebra on S_c).  P and the phases are those of above.
    def S_of(e):
        return np.array([-0.9 + 0.31 * e - 0.021 * e * e, -1.4 + 0.17 * e + 0.013 * e * e])

    def dS_of(e):
        return np.array([0.31 - 0.042 * e, 0.17 + 0.026 * e])

    e_formal = np.array([2.37, 3.5, 6.2])
    g_formal = np.array([[0.61, 0.12], [-0.83, 0.47], [1.3, -0.9]])
    B3 = S_of(2.0)

    def formal_matrix(e_ref):  # the real symmetric matrix whose eigenvalue E = e_ref is a Brune energy
        return np.diag(e_formal) - (g_formal * (S_of(e_ref) - B3)) @ g_formal.T

    e_brune, g_brune = np.zeros(3), np.zeros((3, 2))
    for i in range(3):  # Brune (2002): E~_i is the i-th eigenvalue of the formal matrix with B_c = S_c(E~_i)
        e_i = e_formal[i]
        for _ in range(200):
            w, v = np.linalg.eigh(formal_matrix(e_i))
            step, e_i = abs(w[i] - e_i), w[i]
            if step < 1e-15:
                break
        e_brune[i], g_brune[i] = e_i, v[:, i] @ g_formal  # gamma~_ic = a_i^T gamma_c, a_i^T a_i = 1
    S_lev = np.array([S_of(e) for e in e_brune])
    dS_lev = np.array([dS_of(e) for e in e_brune])
    norm = 1.0 / (1.0 + (g_brune ** 2 * dS_lev).sum(axis=1))  # Park's J in Brune's amplitudes
    g_park = g_brune * np.sqrt(norm)[:, None]
    J_park = 1.0 - (g_park ** 2 * dS_lev).sum(axis=1)
    check("the numpy transformation moved the levels (the comparison below is not of equal inputs)",
          np.abs(e_brune - e_formal).max() > 1e-2 and np.abs(g_brune - g_formal).max() > 1e-2,
          f"E~ = {e_brune}")
    d = np.abs(J_park - norm).max()
    check("Park's J = 1 - sum gamma_P^2 S' equals 1 / (1 + sum gamma~^2 S')", d < 1e-14, f"{d:.2e}")

    def brune_matrix(e):
        s_e = S_of(e)
        m = np.zeros((3, 3), complex)
        for l in range(3):
            for k in range(3):
                if l == k:
                    m[l, l] = e_brune[l] - e - (g_brune[l] ** 2 * (s_e - S_lev[l] + 1j * P)).sum()
                else:
                    chord = (S_lev[l] * (e - e_brune[k]) - S_lev[k] * (e - e_brune[l])) / (e_brune[l] - e_brune[k])
                    m[l, k] = -(g_brune[l] * g_brune[k] * (s_e + 1j * P - chord)).sum()
        return m

    def park_matrix(e):
        s_e = S_of(e)
        m = np.zeros((3, 3), complex)
        for l in range(3):
            for k in range(3):
                kernel = (g_park[l] * g_park[k] * (s_e - S_lev[l] + 1j * P)).sum()
                if l == k:
                    overlap = J_park[l]
                else:
                    overlap = -(g_park[l] * g_park[k] * (S_lev[l] - S_lev[k])).sum() / (e_brune[l] - e_brune[k])
                m[l, k] = (e_brune[l] - e) * overlap - kernel
        return m

    worst = dict(brune_A=0.0, park_A=0.0, brune_U=0.0, park_U=0.0, unitary=0.0)
    energies = [1.9, 2.3, 2.36, 3.1, 3.52, 4.4, 6.0, 7.5]
    for e in energies:
        common = dict(penetrability=P, shift=S_of(e), hard_sphere_phase=phi, coulomb_phase=omega)
        formal = ph.collision_matrix(e_formal, g_formal, e, boundary=B3, **common)
        brune = ph.collision_matrix(e_brune, g_brune, e, parametrization="brune", shift_at_levels=S_lev, **common)
        park = ph.collision_matrix(e_brune, g_park, e, parametrization="park", shift_at_levels=S_lev,
                                   park_norms=J_park, **common)
        worst["brune_A"] = max(worst["brune_A"], np.abs(brune["A"] @ brune_matrix(e) - np.eye(3)).max())
        worst["park_A"] = max(worst["park_A"], np.abs(park["A"] @ park_matrix(e) - np.eye(3)).max())
        worst["brune_U"] = max(worst["brune_U"], np.abs(brune["U"] - formal["U"]).max())
        worst["park_U"] = max(worst["park_U"], np.abs(park["U"] - formal["U"]).max())
        worst["unitary"] = max(worst["unitary"], np.abs(brune["U"] @ brune["U"].conj().T - np.eye(2)).max(),
                               np.abs(park["U"] @ park["U"].conj().T - np.eye(2)).max())
    check("brune: A is the inverse of Brune's Eq. (33) written in numpy", worst["brune_A"] < 1e-12,
          f"largest difference {worst['brune_A']:.2e}")
    check("park: A is the inverse of Park's Eqs. (21), (22) written in numpy", worst["park_A"] < 1e-12,
          f"largest difference {worst['park_A']:.2e}")
    check(f"brune: U equals U of the formal levels at {len(energies)} energies", worst["brune_U"] < 1e-12,
          f"largest difference {worst['brune_U']:.2e}")
    check(f"park: U equals U of the formal levels at {len(energies)} energies", worst["park_U"] < 1e-12,
          f"largest difference {worst['park_U']:.2e}")
    check("brune, park: U is unitary", worst["unitary"] < 1e-12, f"largest difference {worst['unitary']:.2e}")
    print(f"    largest differences: A brune {worst['brune_A']:.1e}, A park {worst['park_A']:.1e}; "
          f"U brune - formal {worst['brune_U']:.1e}, U park - formal {worst['park_U']:.1e}")

    print("coincident levels and one Park level")
    # Two levels 0.5e-6 MeV apart (below the 1e-6 MeV of physics/LevelPairShift.h): the chord of S between them is
    # its Taylor form with the mean derivative; the divided difference of numpy differs from it by the third
    # derivative of S (0 here: S is quadratic) and by its own rounding, 2 eps |S| / gap = 4e-10.
    e_pair = np.array([2.5, 2.5 + 0.5e-6])
    g_pair = np.array([[0.4, 0.2], [0.3, -0.5]])
    S_pair = np.array([S_of(x) for x in e_pair])
    dS_pair = np.array([dS_of(x) for x in e_pair])
    e = 2.9
    common = dict(penetrability=P, shift=S_of(e), hard_sphere_phase=phi, coulomb_phase=omega)
    near = ph.collision_matrix(e_pair, g_pair, e, parametrization="brune", shift_at_levels=S_pair,
                               shift_derivative_at_levels_per_MeV=dS_pair, **common)
    m = np.zeros((2, 2), complex)
    for l in range(2):
        k = 1 - l
        m[l, l] = e_pair[l] - e - (g_pair[l] ** 2 * (S_of(e) - S_pair[l] + 1j * P)).sum()
        chord = (S_pair[l] * (e - e_pair[k]) - S_pair[k] * (e - e_pair[l])) / (e_pair[l] - e_pair[k])
        m[l, k] = -(g_pair[l] * g_pair[k] * (S_of(e) + 1j * P - chord)).sum()
    d = np.abs(near["U"] - np.outer(Omega, Omega) * (
        np.eye(2) + 2j * np.sqrt(np.outer(P, P)) * (g_pair.T @ np.linalg.inv(m) @ g_pair))).max()
    check("brune, levels 0.5 eV apart: U of the Taylor form equals U of the divided difference to 1e-8", d < 1e-8,
          f"largest difference {d:.2e}")
    ok, got = raises(ValueError, "reads shift_derivative_at_levels_per_MeV",
                     lambda: ph.collision_matrix(e_pair, g_pair, e, parametrization="brune",
                                                 shift_at_levels=S_pair, **common))
    check("coincident levels without the derivative: ValueError naming it", ok, got)
    gp, el = np.array([[0.61, 0.12]]), 2.37
    s_l, ds_l = S_of(el), dS_of(el)
    j = 1.0 - (gp[0] ** 2 * ds_l).sum()
    one = ph.collision_matrix([el], gp, e, parametrization="park", shift_at_levels=[s_l], park_norms=[j], **common)
    want = 1.0 / ((el - e) * j - (gp[0] ** 2 * (S_of(e) - s_l + 1j * P)).sum())
    d = abs(one["A"][0, 0] - want) / abs(want)
    check("one Park level: A = 1 / ((E_l - E) J - sum gamma^2 (S - S_l + i P))", d < 1e-14, f"{d:.2e}")

    print("photon channels")
    # One level, a proton channel and a photon channel.  The photon's width is not in the level matrix:
    # |U_p,gamma|^2 = Gamma_p Gamma_gamma / ((E_l - E)^2 + Gamma_p^2 / 4) with B = S, Gamma_gamma = 2 P_gamma g^2.
    e, el = 2.35, 2.37
    k_gamma = 2.3 / 197.3269804  # E_gamma / (hbar c), fm^-1, for an E1 photon of 2.3 MeV
    Pg = np.array([P[0], k_gamma ** 3])
    gg = np.array([[0.61, 0.9]])
    photon = dict(penetrability=Pg, shift=[S[0], 0.0], hard_sphere_phase=[phi[0], 0.0], coulomb_phase=[0.2, 0.0],
                  photon_channels=[False, True])
    cap = ph.collision_matrix([el], gg, e, boundary=[S[0], 0.0], **photon)
    width = 2.0 * Pg * gg[0] ** 2
    want = width[0] * width[1] / ((el - e) ** 2 + 0.25 * width[0] ** 2)
    d = abs(abs(cap["U"][0, 1]) ** 2 - want) / want
    check("formal, B = S: |U_p,gamma|^2 = Gamma_p Gamma_gamma / ((E_l - E)^2 + Gamma_p^2 / 4)", d < 1e-12, f"{d:.2e}")
    d = abs(cap["A"][0, 0] - 1.0 / (el - e - 1j * Pg[0] * gg[0, 0] ** 2))
    check("the photon width is not in A: A = 1 / (E_l - E - i P_p gamma_p^2)", d < 1e-12 * abs(cap["A"][0, 0]),
          f"{d:.2e}")
    d = abs(cap["T"][0, 1] + cap["U"][0, 1]) + abs(cap["U"][0, 1] - cap["U"][1, 0])
    check("T_p,gamma = -U_p,gamma, U symmetric", d < 1e-15, f"{d:.2e}")
    capb = ph.collision_matrix([el], gg, e, parametrization="brune", shift_at_levels=[[S[0], 0.0]], **photon)
    d = np.abs(capb["U"] - cap["U"]).max()
    check("brune with S_l = S: the same U as formal with B = S", d < 1e-15, f"{d:.2e}")
    # Two levels with a photon channel against numpy: the photon column is left out of the sum over channels.
    g2p = np.array([[0.61, 0.9], [-0.83, 0.4]])
    two = ph.collision_matrix(levels, g2p, E, boundary=[B[0], 0.0], **photon)
    m = np.diag(levels - E).astype(complex) - np.outer(g2p[:, 0], g2p[:, 0]) * (S[0] - B[0] + 1j * Pg[0])
    Om = np.exp(1j * np.array([0.2 - phi[0], 0.0]))
    want = np.outer(Om, Om) * (np.eye(2) + 2j * np.sqrt(np.outer(Pg, Pg)) * (g2p.T @ np.linalg.inv(m) @ g2p))
    d = max(np.abs(two["A"] @ m - np.eye(2)).max(), np.abs(two["U"] - want).max() / np.abs(want).max())
    check("two levels: A and U against numpy with the photon column out of the level matrix", d < 1e-12, f"{d:.2e}")
    # Reich-Moore: i gamma^2 of the photon channel on the diagonal, the photon's factor 1 (a session's use_rmc).
    rm = dict(photon, penetrability=[P[0], 1.0])
    rmc = ph.collision_matrix([el], gg, e, boundary=[S[0], 0.0], use_rmc=True, **rm)
    wp, wg = 2.0 * P[0] * gg[0, 0] ** 2, 2.0 * gg[0, 1] ** 2
    want = wp * wg / ((el - e) ** 2 + 0.25 * (wp + wg) ** 2)
    d = abs(abs(rmc["U"][0, 1]) ** 2 - want) / want
    check("use_rmc: |U_p,gamma|^2 = Gamma_p Gamma_gamma / ((E_l - E)^2 + (Gamma_p + Gamma_gamma)^2 / 4), "
          "Gamma_gamma = 2 gamma^2", d < 1e-12, f"{d:.2e}")

    print("errors")
    kw = dict(penetrability=P, shift=S, boundary=B, hard_sphere_phase=phi, coulomb_phase=omega)
    nb = dict(penetrability=P, shift=S, hard_sphere_phase=phi, coulomb_phase=omega)
    S2 = np.zeros((2, 2))
    for what, kind, text, call in [
            ("reduced_widths of the wrong shape", ValueError, "reduced_widths must have the shape (levels, channels)",
             lambda: ph.collision_matrix(levels, gamma[:1], E, **kw)),
            ("a per-channel array of the wrong length", ValueError, "shift must hold one value per channel (2",
             lambda: ph.collision_matrix(levels, gamma, E, **dict(kw, shift=S[:1]))),
            ("a value that is not finite", ValueError, "boundary of channel 1 is not finite",
             lambda: ph.collision_matrix(levels, gamma, E, **dict(kw, boundary=[0.0, float("nan")]))),
            ("a negative penetrability", ValueError, "penetrability of channel 0 is negative",
             lambda: ph.collision_matrix(levels, gamma, E, **dict(kw, penetrability=[-0.1, 0.3]))),
            ("an energy that is not finite", ValueError, "energy is not finite",
             lambda: ph.collision_matrix(levels, gamma, float("inf"), **kw)),
            ("a singular level matrix", pyazr.EngineError, "singular at energy = 1.0 MeV",
             lambda: ph.collision_matrix([1.0, 1.0], [[1.0], [1.0]], 1.0, penetrability=[0.0], shift=[1.0],
                                         boundary=[0.0], hard_sphere_phase=[0.0])),
            ("an unknown parametrization", ValueError, "allowed: 'formal', 'brune', 'park'",
             lambda: ph.collision_matrix(levels, gamma, E, parametrization="observed", **kw)),
            ("formal without boundary", ValueError, "boundary (B_c, one value per channel) is required",
             lambda: ph.collision_matrix(levels, gamma, E, **nb)),
            ("brune with boundary", ValueError, "boundary is not an argument of parametrization 'brune'",
             lambda: ph.collision_matrix(levels, gamma, E, parametrization="brune", shift_at_levels=S2, **kw)),
            ("brune without shift_at_levels", ValueError, "shift_at_levels (S_c at each level energy",
             lambda: ph.collision_matrix(levels, gamma, E, parametrization="brune", **nb)),
            ("shift_at_levels of the wrong shape", ValueError, "shift_at_levels must have the shape (levels, channels)",
             lambda: ph.collision_matrix(levels, gamma, E, parametrization="brune", shift_at_levels=S2[:1], **nb)),
            ("park without park_norms", ValueError, "park_norms (J_l = 1 - sum_c gamma_lc^2 dS_c/dE",
             lambda: ph.collision_matrix(levels, gamma, E, parametrization="park", shift_at_levels=S2, **nb)),
            ("park_norms with brune", ValueError, "park_norms is not an argument of parametrization 'brune'",
             lambda: ph.collision_matrix(levels, gamma, E, parametrization="brune", shift_at_levels=S2,
                                         park_norms=[1.0, 1.0], **nb)),
            ("shift_at_levels with formal", ValueError, "shift_at_levels is not an argument of parametrization 'formal'",
             lambda: ph.collision_matrix(levels, gamma, E, shift_at_levels=S2, **kw)),
            ("a photon channel with a shift", ValueError, "shift of channel 1 must be 0: it is a photon channel",
             lambda: ph.collision_matrix(levels, gamma, E, photon_channels=[False, True], **kw)),
            ("photon_channels of the wrong length", ValueError, "photon_channels must hold one True or False per",
             lambda: ph.collision_matrix(levels, gamma, E, photon_channels=[True], **kw)),
            ("use_rmc with brune", ValueError, "use_rmc needs parametrization 'formal'",
             lambda: ph.collision_matrix(levels, gamma, E, parametrization="brune", shift_at_levels=S2,
                                         use_rmc=True, **nb)),
            ("use_rmc without a photon channel", ValueError, "photon_channels names none",
             lambda: ph.collision_matrix(levels, gamma, E, use_rmc=True, **kw))]:
        ok, got = raises(kind, text, call)
        check(f"{what}: {kind.__name__} with '{text}'", ok, got)
    check("collision_matrix is exported by pyazr.physics and is not a row of the table",
          "collision_matrix" in ph.__all__ and "collision_matrix" not in ph.functions)

    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
