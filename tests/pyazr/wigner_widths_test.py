#!/usr/bin/env python3
"""Session.wigner_widths: Gamma_W = 2 P_l(E_r) gamma_W^2 in eV, against an independent calculation.

gamma_W^2 = 3 hbar^2 / (2 mu a^2) (the engine's Wigner limit), P_l the penetrability at the pair's channel radius,
E_r the pair's c.m. energy at the level's physical (observed) energy (docs/upgrade/pyazr3.md section 11).

On T1, T2, T3 (tests/upgrade) and tests/13N, in Brune's, the formal and Park's parametrization:
(a) every entry against 2 P gamma_W^2 computed here: P = rho / (F_l^2 + G_l^2) with mpmath's Coulomb functions,
    eta = Z1 Z2 alpha (mu c^2 / 2E)^(1/2), rho = (2 mu c^2 E)^(1/2) a / hbar c, mu from the pair's masses, the
    constants written below; to 1e-8 (a session with a Coulomb table interpolates F and G in it; measured: 2.5e-10
    with the table, 3.6e-10 with AZURE2_COULOMB_TABLE unset).  A pair with a hybrid nuclear potential (alpha
    + 4He of T3) has the hybrid model's penetrability, not the Coulomb one: its entries are checked against 2 P
    gamma_W^2 with P of Session.channel_functions, and differ from the Coulomb value (by 0.9 % there);
(b) a neutral channel (7Be + n of T3) against the closed forms P_0 = rho, P_1 = rho^3 / (1 + rho^2), to 1e-12;
(c) the same numbers in the three parametrizations (6e-7 measured: the formal basis goes through the engine's
    transformation of the level energy);
(d) which channels have an entry: the free particle channels open at the level's physical energy, no photon, fixed
    or closed channel; eps is ignored;
(e) dimensionless_widths: wigner_width is the entry, theta2 = |Gamma| / Gamma_W, theta2 / theta2_formal of a Brune
    level is 1 / (1 + sum_c gamma_c^2 dS_c/dE) over its particle channels (the level's norm; dS/dE of
    pyazr.physics.shift_derivative), the relation the earlier probe had folded into Gamma_W;
(f) the numbers of the design note's table (section 11.3), at 1e-6: T2 width_1_3 is 4.020247e6 eV (the probe gave
    3.208103e6), T1 width_3_1 is 1.210273e6 eV in the formal basis too (the probe gave 1.827438e6);
(g) a pair with a spin-orbit potential (tests/spin_orbit_potential): the entry is 2 P gamma_W^2 with P of
    channel_functions at the level's J;
(h) dimensionless_widths of a FIXED particle channel (pyazr3.md section 13): open at the level's physical energy, its
    wigner_width is 2 P gamma_W^2 of mpmath (or of channel_functions on a hybrid pair) as a free one's, to 1e-8, and
    theta2 = |value| / wigner_width; closed, both None;
(i) level_energy, is_open and a photon's e_gamma are taken at the physical energy (the row of transform): in Brune's
    and Park's parametrization equal bit for bit to the fit vector's energy (what they were before), in the formal
    one every row equal to Brune's (level_energy and e_gamma to 1e-6, is_open the same);
(j) the entries equal 2 P gamma_W^2 with P of channel_functions at E_r = E_x - (S + e), the route of pyazr
    3.0.0.dev5, bit for bit (the record's P is the same function at the same channel energy).

Run from anywhere:  python3 tests/pyazr/wigner_widths_test.py
"""
import os
import shutil
import sys
import tempfile
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


try:
    import mpmath as mp
    import numpy as np
    import pyazr
    from pyazr import physics
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(77)

mp.mp.dps = 30
U, HC, ALPHA = mp.mpf("931.49410242"), mp.mpf("197.3269804"), 1 / mp.mpf("137.035999084")  # Constants.h
TABLE = bool(os.environ.get("AZURE2_COULOMB_TABLE"))
TOL = 1e-8


def mp_wigner_width(pair, l, e):
    """2 P_l gamma_W^2 in eV for the pair at c.m. energy e (MeV), from the pair's masses, charges and radius."""
    mu = mp.mpf(pair.M1) * pair.M2 / (mp.mpf(pair.M1) + pair.M2)
    a = mp.mpf(pair.channel_radius)
    eta = mp.sqrt(U / 2) * ALPHA * pair.Z1 * pair.Z2 * mp.sqrt(mu / mp.mpf(e))
    rho = mp.sqrt(2 * U * mu * mp.mpf(e)) / HC * a
    P = rho / (mp.coulombf(l, eta, rho) ** 2 + mp.coulombg(l, eta, rho) ** 2)
    return float(2 * P * 3 * HC ** 2 / (2 * mu * U * a * a) * 1e6), float(rho)


PROJECTS = {"T1": ("tests/upgrade/T1", "T1.azr"), "T2": ("tests/upgrade/T2", "T2.azr"),
            "T3": ("tests/upgrade/T3", "T3.azr"), "13N": ("tests/13N", "13N.azr")}
MODES = {"Brune": dict(use_brune=True), "formal": dict(use_brune=False), "Park": dict(use_park=True)}
tmp = tempfile.mkdtemp(prefix="wigner_widths_")
try:
    by_mode, rows_by_name = {}, {}
    worst_all = 0.0
    for name, (directory, azr) in PROJECTS.items():
        for mode, options in MODES.items():
            work = os.path.join(tmp, f"{name}_{mode}")
            shutil.copytree(os.path.join(ROOT, directory), work, ignore=shutil.ignore_patterns("output", "expected"))
            os.makedirs(os.path.join(work, "output"), exist_ok=True)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s = pyazr.Session(os.path.join(work, azr), **options)
            with s:
                x = s.parameters.x0
                got = s.wigner_widths(x)
                physical = s.transform(x)
                ex = {(p.jgroup, p.level): float(physical[p.index]) for p in s.parameters.of_kind("energy")}
                expected, neutral, hybrid = {}, [], {}
                for p in s.parameters.of_kind("width"):
                    pair = s.pairs.by_number(p.pair)
                    if p.fixed or p.free_index is None or pair.is_photon:
                        continue
                    e = ex[(p.jgroup, p.level)] - pair.sep_energy - pair.excitation
                    if e > 0:
                        expected[p.free_index], rho = mp_wigner_width(pair, p.L, e)
                        if s.nuclear_potential(pair.number).enabled:   # the hybrid model's P, not the Coulomb one
                            P = float(s.channel_functions(pair.key, e, p.L)["P"])
                            hybrid[p.free_index] = (expected[p.free_index], 2 * P * p.wigner_limit * 1e6)
                            expected[p.free_index] = hybrid[p.free_index][1]
                        if pair.Z1 * pair.Z2 == 0 and p.L in (0, 1):
                            closed = rho if p.L == 0 else rho ** 3 / (1 + rho ** 2)
                            mu = pair.M1 * pair.M2 / (pair.M1 + pair.M2)
                            gw2 = 3 * float(HC) ** 2 / (2 * mu * float(U) * pair.channel_radius ** 2)
                            neutral.append((p.name, got.get(p.free_index), 2 * closed * gw2 * 1e6))
                worst = max((abs(got[k] / v - 1) for k, v in expected.items() if k in got), default=float("inf"))
                worst_all = max(worst_all, worst)
                check(f"(a) {name}, {mode}: {len(expected) - len(hybrid)} entries = 2 P gamma_W^2 of mpmath"
                      + (f" ({len(hybrid)} more on a pair with a nuclear potential = 2 P gamma_W^2 of its channel "
                         f"functions, {max(abs(c / h - 1) for c, h in hybrid.values()):.1e} from the Coulomb value)"
                         if hybrid else "") + f", largest relative difference {worst:.1e} (<= {TOL:g})",
                      set(got) == set(expected) and worst <= TOL
                      and all(abs(c / h - 1) > 1e-6 for c, h in hybrid.values()),
                      f"entries {sorted(got)} against {sorted(expected)}")
                for label, value, want in neutral:
                    check(f"(b) {name}, {mode}: neutral channel {label} = the closed form of P_l "
                          f"({want:.9g} eV) to 1e-12", value is not None and abs(value / want - 1) <= 1e-12,
                          f"{value!r}")
                by_mode[(name, mode)] = {s.parameters.free.names[k]: v for k, v in got.items()}
                # (d) what has an entry
                rows = {p.free_index: p for p in s.parameters.of_kind("width") if p.free_index is not None}
                photons = [k for k in got if rows[k].radiation_type in ("E", "M")]
                check(f"(d) {name}, {mode}: no photon and no fixed channel has an entry; eps is ignored",
                      not photons and all(not rows[k].fixed for k in got) and s.wigner_widths(x, eps=1e-2) == got
                      and s.wigner_widths() == got)
                # (e) the table of dimensionless widths
                table = s.dimensionless_widths(x)
                particle = [c for c in table.particles if c.free_index in got and not c.fixed]
                ok = all(c.wigner_width == got[c.free_index] for c in particle) and all(
                    c.theta2 is None or abs(c.theta2 - abs(c.value) / got[c.free_index]) <= 1e-15 * c.theta2
                    for c in particle)
                check(f"(e) {name}, {mode}: ChannelWidth.wigner_width is the entry, theta2 = |Gamma| / Gamma_W", ok)
                # (h) fixed particle channels
                n_open, n_closed, worst_fixed, fixed_ok = 0, 0, 0.0, True
                for c in table.particles:
                    if not c.fixed:
                        continue
                    pair = s.pairs.by_number(c.pair)
                    q = next(q for q in s.parameters.of_kind("width") if q.index == c.index)
                    e = ex[(c.jgroup, c.level)] - pair.sep_energy - pair.excitation
                    if e > 0 and q.wigner_limit:
                        n_open += 1
                        if s.nuclear_potential(pair.number).enabled:
                            want = 2 * float(s.channel_functions(pair.key, e, c.L)["P"]) * q.wigner_limit * 1e6
                        else:
                            want = mp_wigner_width(pair, c.L, e)[0]
                        if c.wigner_width is None:
                            fixed_ok = False
                            continue
                        worst_fixed = max(worst_fixed, abs(c.wigner_width / want - 1))
                        fixed_ok = fixed_ok and c.is_open and (c.value is None or c.theta2 == abs(c.value) / c.wigner_width)
                    else:
                        n_closed += 1
                        fixed_ok = fixed_ok and c.wigner_width is None and c.theta2 is None
                check(f"(h) {name}, {mode}: {n_open} fixed open channels have wigner_width = 2 P gamma_W^2 (largest "
                      f"relative difference {worst_fixed:.1e}, <= {TOL:g}) and theta2 = |value| / wigner_width; "
                      f"{n_closed} fixed closed ones neither", fixed_ok and worst_fixed <= TOL)
                # (i) the physical energy; Brune and Park: the fit vector's, as before
                rows_by_name[(name, mode)] = {c.name: c for c in table}
                energy_ok = all(c.level_energy == ex[(c.jgroup, c.level)] for c in table)
                if mode != "formal":
                    before = {}
                    for q in s.parameters.of_kind("energy"):
                        free_e = not q.fixed and q.free_index is not None and q.free_index < x.size
                        before[(q.jgroup, q.level)] = (float(x[q.free_index]) if free_e
                                                       else (0.0 if q.level_energy is None else q.level_energy))
                    same = all(c.level_energy == before[(c.jgroup, c.level)] for c in table) and all(
                        c.is_open == (before[(c.jgroup, c.level)] > c.threshold) for c in table.particles) and all(
                        c.e_gamma == before[(c.jgroup, c.level)] - s.pairs.by_number(c.pair).excitation
                        for c in table.photons)
                    check(f"(i) {name}, {mode}: level_energy is the physical energy, equal bit for bit to the fit "
                          f"vector's, and is_open, e_gamma as before", energy_ok and same)
                else:
                    check(f"(i) {name}, formal: level_energy is the physical energy (the row of transform)", energy_ok)
                # (j) the route of 3.0.0.dev5: P of channel_functions at E_x - (S + e)
                same_route = True
                for k, v in got.items():
                    q = rows[k]
                    pair = s.pairs.by_number(q.pair)
                    potential = s.nuclear_potential(pair.number)
                    spin_orbit = potential.enabled and potential.spin_orbit is not None
                    P = float(s.channel_functions(pair.key, ex[(q.jgroup, q.level)] - (pair.sep_energy + pair.excitation),
                                                  q.L, J=q.J if spin_orbit else None)["P"])
                    same_route = same_route and v == 2.0 * P * q.wigner_limit * 1e6
                check(f"(j) {name}, {mode}: {len(got)} entries equal the channel_functions route bit for bit", same_route)
                if mode == "Brune":
                    with pyazr.locking.session_lock(s):
                        amplitudes = np.asarray(s._engine.params_all_rwa(), float)   # the binding: fixed amplitudes
                    amplitudes[s.parameters.free_mask] = x
                    worst_norm, n_norm = 0.0, 0
                    for c in particle:
                        if not c.theta2 or not c.theta2_formal or c.theta2_formal < 1e-12:
                            continue
                        total, usable = 0.0, True
                        for q in s.parameters.of_kind("width"):
                            if (q.jgroup, q.level) != (c.jgroup, c.level):
                                continue
                            qp = s.pairs.by_number(q.pair)
                            if qp.is_photon:
                                continue
                            if s.nuclear_potential(qp.number).enabled:   # dS/dE of the Coulomb functions is not its
                                usable = False
                                continue
                            e = ex[(q.jgroup, q.level)] - qp.sep_energy - qp.excitation
                            if e <= 0:
                                usable = usable and amplitudes[q.index] == 0.0   # a closed channel: not checked here
                                continue
                            total += amplitudes[q.index] ** 2 * float(physics.shift_derivative(
                                q.L, qp.Z1, qp.Z2, qp.M1, qp.M2, qp.channel_radius, e))
                        if usable:
                            n_norm += 1
                            worst_norm = max(worst_norm, abs(c.theta2 / c.theta2_formal * (1 + total) - 1))
                    check(f"(e) {name}, Brune: theta2 / theta2_formal = 1 / (1 + sum gamma^2 dS/dE) for {n_norm} "
                          f"channels, largest relative difference {worst_norm:.1e} (<= 2e-5, the engine's dS/dE)",
                          n_norm > 0 and worst_norm <= 2e-5)
    for name in PROJECTS:
        ref = by_mode[(name, "Brune")]
        worst = max(abs(by_mode[(name, mode)][k] / v - 1) for mode in ("formal", "Park") for k, v in ref.items()
                    if k in by_mode[(name, mode)])
        same = all(set(by_mode[(name, mode)]) == set(ref) for mode in ("formal", "Park"))
        check(f"(c) {name}: the same channels and numbers in the three parametrizations, largest relative difference "
              f"{worst:.1e} (<= 5e-6)", same and worst <= 5e-6)
    for name in PROJECTS:
        brune, formal = rows_by_name[(name, "Brune")], rows_by_name[(name, "formal")]
        de = max(abs(formal[k].level_energy - c.level_energy) for k, c in brune.items())
        dg = max((abs(formal[k].e_gamma - c.e_gamma) for k, c in brune.items() if c.is_photon), default=0.0)
        dw = max((abs(formal[k].wigner_width / c.wigner_width - 1) for k, c in brune.items()
                  if c.wigner_width and formal[k].wigner_width), default=0.0)
        same_open = all(formal[k].is_open == c.is_open and (formal[k].wigner_width is None) == (c.wigner_width is None)
                        for k, c in brune.items())
        check(f"(i) {name}: the formal basis's rows are Brune's: level_energy {de:.1e} MeV, e_gamma {dg:.1e} MeV "
              f"(<= 1e-6), wigner_width {dw:.1e} (free and fixed, <= 5e-6), is_open the same",
              set(formal) == set(brune) and de <= 1e-6 and dg <= 1e-6 and dw <= 5e-6 and same_open)
    t2, t1f = by_mode[("T2", "Brune")]["width_1_3"], by_mode[("T1", "formal")]["width_3_1"]
    check(f"(f) T2 width_1_3: {t2:.7g} eV = 4.020247e6 (the probe gave 3.208103e6); T1 formal width_3_1: {t1f:.7g} eV "
          f"= 1.210273e6 (the probe gave 1.827438e6)",
          abs(t2 / 4.020247e6 - 1) <= 1e-6 and abs(t1f / 1.210273e6 - 1) <= 1e-6)
    print(f"largest relative difference from mpmath over all projects: {worst_all:.2e} "
          f"({'Coulomb table' if TABLE else 'exact Coulomb routine'})")

    # (g) a pair with a spin-orbit potential
    work = os.path.join(tmp, "spin_orbit")
    shutil.copytree(os.path.join(ROOT, "tests", "spin_orbit_potential"), work,
                    ignore=shutil.ignore_patterns("output", "expected"))
    os.makedirs(os.path.join(work, "output"), exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = pyazr.Session(os.path.join(work, "spin_orbit_potential.azr"))
    with s:
        x = s.parameters.x0
        got = s.wigner_widths(x)
        n, worst = 0, 0.0
        for p in s.parameters.of_kind("width"):
            if p.free_index not in got or p.fixed:
                continue
            pair = s.pairs.by_number(p.pair)
            potential = s.nuclear_potential(pair.number)
            if not (potential.enabled and potential.spin_orbit):
                continue
            ex = float(s.transform(x)[next(q.index for q in s.parameters.of_kind("energy")
                                             if (q.jgroup, q.level) == (p.jgroup, p.level))])
            P = float(s.channel_functions(pair.key, ex - pair.sep_energy - pair.excitation, p.L, J=p.J)["P"])
            n += 1
            worst = max(worst, abs(got[p.free_index] / (2 * P * p.wigner_limit * 1e6) - 1))
        check(f"(g) spin-orbit pair: {n} entries = 2 P gamma_W^2 with P of channel_functions at the level's J "
              f"({worst:.1e})", n > 0 and worst <= 1e-15, f"{len(got)} entries, {n} on a spin-orbit pair")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok")
