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

    print("errors")
    kw = dict(penetrability=P, shift=S, boundary=B, hard_sphere_phase=phi, coulomb_phase=omega)
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
                                         boundary=[0.0], hard_sphere_phase=[0.0]))]:
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
