#!/usr/bin/env python3
"""pyazr.physics.optical_scattering: elastic scattering of a pair in an optical potential (physics/transfer through
api/src/pybind_optical.cpp), without a session or a project.

Checked against values that do not come from the function:

- FRESCO 3.5 (the reference files of physics/tests/data/fresco/, written by make_fresco_benchmarks.py; FRESCO is not
  run here): p + 15N at 5 MeV (spin 1/2, spin-orbit), d + 15N at 10 MeV (spin 1, spin-orbit) and 6Li + 12C at 20 MeV
  (radii by both masses): every S-matrix element, sigma / sigma_R at every angle and the reaction cross section, to
  the tolerances of the C++ test of the library (physics/tests/transfer_fresco_test.cpp: 2e-5, 1e-4, 1e-5);
- closed forms: no potential gives S = 1 and no reaction cross section; a real potential gives |S| = 1; a neutral
  pair in a Woods-Saxon well with l_max = 0 has dsigma/dOmega = sin^2(delta_0) / k^2 at every angle, delta_0 from
  the S the function returns;
- k and eta against the rows wave_number and sommerfeld_eta;
- the errors name the argument.

Run from anywhere:  python3 tests/pyazr/optical_scattering_test.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
FRESCO = os.path.join(ROOT, "physics", "tests", "data", "fresco")
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


def reference(case):
    """The lines of FRESCO's reference run of a case: S {(l, two_j): S}, ratio [(theta, value)], reaction (mb)."""
    S, ratio, reaction, run = {}, [], None, None
    with open(os.path.join(FRESCO, case, "reference.txt")) as f:
        for line in f:
            parts = line.split()
            if not parts or parts[0].startswith("#"):
                continue
            if parts[0] == "run":
                run = parts[1]
            elif run != "reference":
                continue
            elif parts[0] == "S":
                S[(int(parts[1]), int(parts[2]))] = complex(float(parts[4]), float(parts[5]))
            elif parts[0] == "ratio":
                ratio.append((float(parts[1]), float(parts[2])))
            elif parts[0] == "reaction":
                reaction = float(parts[1])
    return S, ratio, reaction


def main():
    try:
        import numpy as np
        import pyazr  # noqa: F401
        from pyazr import physics as ph
    except Exception as err:  # the skip of the other pyazr tests
        print(f"skip: pyazr engine not importable ({type(err).__name__}: {err})")
        return 77

    # The inputs of the three FRESCO jobs (job.nin; masses are mass numbers there).
    p15n = dict(radius_convention="target_mass", target_mass_number=15, coulomb_r0_fm=1.25,
                real_volume=(50.0, 1.2, 0.65), imaginary_volume=(1.0, 1.3, 0.5), imaginary_surface=(6.0, 1.3, 0.5),
                real_spin_orbit=(6.0, 1.1, 0.6))
    d15n = dict(radius_convention="target_mass", target_mass_number=15, coulomb_r0_fm=1.3,
                real_volume=(90.0, 1.15, 0.75), imaginary_volume=(2.0, 1.35, 0.6),
                imaginary_surface=(10.0, 1.35, 0.6), real_spin_orbit=(6.0, 1.07, 0.66))
    li6c12 = dict(radius_convention="both_masses", target_mass_number=12, projectile_mass_number=6,
                  coulomb_r0_fm=0.73, real_volume=(120.0, 0.80, 0.85), imaginary_volume=(35.0, 0.95, 0.85))
    cases = [("p15N_elastic", p15n, dict(z1=1, z2=7, m1_u=1.0, m2_u=15.0, e_cm=5.0 * 15.0 / 16.0),
              dict(projectile_spin=0.5, l_max=22, step_fm=0.05, matching_radius_fm=30.0)),
             ("d15N_elastic", d15n, dict(z1=1, z2=7, m1_u=2.0, m2_u=15.0, e_cm=10.0 * 15.0 / 17.0),
              dict(projectile_spin=1, l_max=32, step_fm=0.05, matching_radius_fm=30.0)),
             ("li6c12_elastic", li6c12, dict(z1=3, z2=6, m1_u=6.0, m2_u=12.0, e_cm=20.0 * 12.0 / 18.0),
              dict(projectile_spin=1, l_max=39, step_fm=0.025, matching_radius_fm=30.0))]
    print("against FRESCO 3.5")
    for case, potential, pair, grid in cases:
        S_ref, ratio_ref, reaction_ref = reference(case)
        angles = [min(t, 180.0) for t, _ in ratio_ref]
        out = ph.optical_scattering(potential=potential, angles_deg=angles, **pair, **grid)
        ours = {(int(l), int(j)): s for l, j, s in zip(out["l"], out["two_j"], out["S"])}
        if not out["spin_orbit"]:  # one wave per l (two_j = 2l): FRESCO lists every j
            ours = {(l, j): ours[(l, 2 * l)] for (l, j) in S_ref if (l, 2 * l) in ours}
        common = [key for key in S_ref if key in ours]
        d = max(abs(ours[key] - S_ref[key]) for key in common)
        check(f"{case}: {len(common)} S-matrix elements, largest |dS| {d:.1e}", len(common) > 40 and d < 2e-5,
              f"{len(common)} elements, {d:.2e}")
        d = np.abs(out["ratio_to_rutherford"] / np.array([v for _, v in ratio_ref]) - 1.0).max()
        check(f"{case}: sigma / sigma_R at {len(angles)} angles, largest relative difference {d:.1e}",
              len(angles) > 20 and d < 1e-4, f"{d:.2e}")
        d = abs(out["reaction_mb"] / reaction_ref - 1.0)
        check(f"{case}: reaction cross section {out['reaction_mb']:.5f} mb against {reaction_ref:.5f} "
              f"(relative {d:.1e})", d < 1e-5, f"{d:.2e}")

    print("closed forms")
    pair = dict(z1=1, z2=7, m1_u=1.0, m2_u=15.0, e_cm=4.6875)
    out = ph.optical_scattering(potential={}, angles_deg=[30.0, 90.0, 150.0], l_max=8, **pair)
    d = max(np.abs(out["S"] - 1.0).max(), np.abs(out["ratio_to_rutherford"] - 1.0).max(), abs(out["reaction_mb"]))
    check("no potential (a point charge): S = 1, sigma / sigma_R = 1, no reaction cross section", d < 1e-9,
          f"{d:.2e}")
    d = max(abs(out["k"] - ph.wave_number(1.0, 15.0, 4.6875)), abs(out["eta"] - ph.sommerfeld_eta(1, 7, 1.0, 15.0,
                                                                                                4.6875)))
    check("k and eta are those of the rows wave_number and sommerfeld_eta", d < 1e-13, f"{d:.2e}")
    d = np.abs(out["coulomb_phase"] - [ph.coulomb_phase(l, out["eta"]) for l in range(9)]).max()
    check("coulomb_phase is the row coulomb_phase for l = 0 ... l_max", d < 1e-12 and out["coulomb_phase"].size == 9,
          f"{d:.2e}")
    real = ph.optical_scattering(potential=dict(p15n, imaginary_volume=(0, 0, 0), imaginary_surface=(0, 0, 0)),
                                 projectile_spin=0.5, l_max=12, **pair)
    d = np.abs(np.abs(real["S"]) - 1.0).max()
    check("a real potential: |S| = 1 for every partial wave, no reaction cross section",
          d < 1e-6 and abs(real["reaction_mb"]) < 1e-3 and np.abs(real["S"] - 1.0).max() > 0.1,
          f"largest ||S| - 1| {d:.2e}, reaction {real['reaction_mb']:.2e} mb")
    well = dict(real_volume=(30.0, 2.5, 0.3))
    n = ph.optical_scattering(0, 6, 1.008665, 12.0, 0.05, well, angles_deg=[20.0, 90.0, 160.0], l_max=0)
    delta = np.angle(n["S"][0]) / 2.0
    want = 10.0 * np.sin(delta) ** 2 / n["k"] ** 2  # dsigma/dOmega = sin^2(delta) / k^2, fm^2 -> mb
    d = np.abs(n["dsigma_domega_mb_sr"] / want - 1.0).max()
    check("a neutral pair, l_max = 0: dsigma/dOmega = sin^2(delta_0) / k^2 at every angle (mb/sr)",
          d < 1e-12 and "ratio_to_rutherford" not in n and abs(np.sin(delta)) > 1e-2, f"{d:.2e}")

    print("errors")
    for what, kind, text, call in [
            ("an unknown member of the potential", ValueError, "the potential has no member 'depth'",
             lambda: ph.optical_scattering(potential=dict(depth=50.0), **pair)),
            ("a term that is not three numbers", ValueError, "must be (depth_MeV, r0_fm, diffuseness_fm)",
             lambda: ph.optical_scattering(potential=dict(real_volume=(50.0, 1.2)), **pair)),
            ("an unknown radius convention", ValueError, "allowed: 'absolute'",
             lambda: ph.optical_scattering(potential=dict(radius_convention="A13"), **pair)),
            ("an angle outside (0, 180]", ValueError, "is outside 0 < theta <= 180",
             lambda: ph.optical_scattering(potential={}, angles_deg=[0.0], **pair)),
            ("a spin that is not 0, 1/2 or 1", ValueError, "projectile_spin = 0.3 is not 0, 0.5 or 1",
             lambda: ph.optical_scattering(potential={}, projectile_spin=0.3, **pair)),
            ("an energy that is not positive", ValueError, "E must be positive",
             lambda: ph.optical_scattering(potential={}, **dict(pair, e_cm=0.0)))]:
        ok, got = raises(kind, text, call)
        check(f"{what}: {kind.__name__} with '{text}'", ok, got)
    check("optical_scattering is exported by pyazr.physics and is not a row of the table",
          "optical_scattering" in ph.__all__ and "optical_scattering" not in ph.functions)

    if failures:
        print(f"FAIL: {len(failures)} check(s)")
        return 1
    print("ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())
