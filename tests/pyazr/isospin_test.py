#!/usr/bin/env python3
"""Isospin-labelled levels across compound nuclei (the <isospin> block), through pyazr.

The model is tests/isospin_levels/models/iso_free.azr (make_models.py): a toy
A=4 system in Hale's charge-independent scheme, 4Li the reference nucleus, 4H
and 4He shifted, one 2- T=1 level in all three and one 2- T=0 level in 4He,
every reference parameter free, priors on the two shifts (0.3 MeV) and on the
d+d mix amplitude (0.1).  The data files are synthetic (smooth curves).

Checks:
  1. the isospin parameters: kind, names, prior centre and sigma, their place
     right after the level block; the derived level slots are fixed;
  2. penalties()["isospin"] is the prior sum, and the engine's chi-squared
     carries exactly it (against the same model with the priors removed);
  3. chi2_and_grad against central differences of calculate_chi2_rwa, every
     free parameter (owner widths, reference energies, shifts, mix);
  4. residual_jacobian columns against central differences of the residuals;
  5. the physical-parameter chi-squared is refused (NaN);
  6. the file side: AzrModel.isospin, engine_level_keys per nucleus,
     renumbering keeps the block's IDs in step, join_levels and set_isospin_*
     rebuild the model from separate levels;
  7. save_fit of moved shifts and mix: the block carries them, the snapshot
     re-evaluates equal (save_fit's own verification) and reloads equal.

Run from anywhere:  python3 tests/pyazr/isospin_test.py
"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
MODELS = os.path.join(ROOT, "tests", "isospin_levels", "models")
sys.path.insert(0, ROOT)

try:
    from pyazr import azure2                                  # noqa: E402
    from pyazr.azrfile import AzrModel                        # noqa: E402
except ImportError as e:  # as the other pyazr tests: no package dependencies -> nothing to test
    print(f"pyazr not importable ({e}); nothing to test")
    sys.exit(0)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def write_data(work):
    # lab-energy grids of the data segments, smooth made-up cross sections, 1 % errors
    for name, e0, e1, de in (("nt.dat", 2.0, 12.0, 0.25), ("ptn.dat", 1.2, 10.0, 0.25),
                             ("ddp.dat", 0.1, 4.0, 0.1), ("ddn.dat", 0.1, 4.0, 0.1)):
        e = np.arange(e0, e1 + 1e-9, de)
        y = 0.05 + 0.02 * np.exp(-((e - 0.5 * (e0 + e1)) / (0.3 * (e1 - e0))) ** 2)
        with open(os.path.join(work, name), "w") as f:
            for a, b in zip(e, y):
                f.write(f"{a:.6f} 0 {b:.12e} {0.01 * b:.12e}\n")


def main():
    work = tempfile.mkdtemp(prefix="isospin_")
    try:
        return run(work)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def run(work):
    os.makedirs(os.path.join(work, "output"), exist_ok=True)
    os.makedirs(os.path.join(work, "checks"), exist_ok=True)
    write_data(work)
    model = os.path.join(work, "iso_free.azr")
    shutil.copy(os.path.join(MODELS, "iso_free.azr"), model)
    # the same model without priors: their sum must be the engine's difference
    text = open(model).read()
    noprior = os.path.join(work, "iso_noprior.azr")
    with open(noprior, "w") as f:
        f.write(text.replace("-0.86 0.3", "-0.86 0").replace("-0.64 0.3", "-0.64 0")
                    .replace("0.02 1 0.1", "0.02 1 0"))

    with azure2(model, cwd=work, use_brune=False) as s:
        # ---- 1. parameters ---------------------------------------------------------
        print("1. parameters")
        iso = s.parameters.isospin
        names = [p.name for p in iso]
        check("three isospin parameters", names == ["isospin_shift_1_4", "isospin_shift_2_4", "isospin_mix_1"], str(names))
        check("kind 'isospin', free", all(p.kind == "isospin" and not p.fixed for p in iso), str(iso))
        check("prior centres and widths",
              [(p.prior_centre, p.prior_sigma) for p in iso] == [(-0.86, 0.3), (-0.64, 0.3), (0.0, 0.1)],
              str([(p.prior_centre, p.prior_sigma) for p in iso]))
        nlev = sum(1 for p in s.parameters if p.kind in ("energy", "width"))
        check("right after the level block", [p.index for p in iso] == [nlev, nlev + 1, nlev + 2], str([p.index for p in iso]))
        free_lev = [p.name for p in s.parameters if p.kind in ("energy", "width") and not p.fixed]
        check("only the reference slots are free",
              free_lev == ["energy_1", "width_1_1", "energy_4", "width_4_1", "width_4_3"], str(free_lev))

        x = np.asarray(s.params_rwa, float)

        # ---- 2. penalty ---------------------------------------------------------------
        print("\n2. penalty")
        pen = s.penalties(x)["isospin"]
        want = ((-0.8 + 0.86) / 0.3) ** 2 + ((-0.6 + 0.64) / 0.3) ** 2 + (0.02 / 0.1) ** 2
        check("penalties()['isospin'] is the prior sum", abs(pen - want) < 1e-12, f"{pen} vs {want}")
        chi = float(np.sum(s.calculate_chi2_rwa(x)))
        with azure2(noprior, cwd=work, use_brune=False) as s0:
            chi0 = float(np.sum(s0.calculate_chi2_rwa(np.asarray(s0.params_rwa, float))))
        check("the engine's chi-squared carries it", abs((chi - chi0) - want) < 1e-8 * max(1.0, chi), f"{chi - chi0} vs {want}")

        # ---- 3. gradient --------------------------------------------------------------
        print("\n3. gradient")
        c, g = s.chi2_and_grad(x)
        check("chi2_and_grad value", abs(c - chi) < 1e-8 * max(1.0, chi), f"{c} vs {chi}")
        worst = 0.0
        for i in range(x.size):
            h = 1e-6 * (abs(x[i]) + 1.0)
            xp, xm = x.copy(), x.copy()
            xp[i] += h
            xm[i] -= h
            fd = (float(np.sum(s.calculate_chi2_rwa(xp))) - float(np.sum(s.calculate_chi2_rwa(xm)))) / (2 * h)
            worst = max(worst, abs(g[i] - fd) / (abs(g[i]) + abs(fd) + 1.0))
        check(f"analytic gradient == central differences over {x.size} free parameters (worst {worst:.1e})", worst < 1e-5)

        # ---- 4. residual Jacobian -----------------------------------------------------
        print("\n4. residual Jacobian")
        r, J = s.residual_jacobian(x)
        r, J = np.asarray(r, float), np.asarray(J, float)
        worst = 0.0
        for i in range(x.size):
            h = 1e-6 * (abs(x[i]) + 1.0)
            xp, xm = x.copy(), x.copy()
            xp[i] += h
            xm[i] -= h
            fd = (np.asarray(s.residual_jacobian(xp)[0], float) - np.asarray(s.residual_jacobian(xm)[0], float)) / (2 * h)
            col = J[:, i]
            worst = max(worst, float(np.max(np.abs(col - fd)) / (np.max(np.abs(fd)) + 1e-12)))
        check(f"residual_jacobian == central differences, every column (worst {worst:.1e})", worst < 1e-5)

        # ---- 5. physical parameters refused -------------------------------------------
        print("\n5. physical parameters")
        try:
            v = float(np.sum(s.sess.calculate_chi2_physical(np.asarray(s.params, float))))
            check("calculate_chi2_physical refused (NaN)", np.isnan(v), str(v))
        except Exception as e:  # an exception is a refusal too
            check("calculate_chi2_physical refused", True, str(e))

        # ---- 7. save_fit -------------------------------------------------------------
        print("\n7. save_fit")
        xs = x.copy()
        ish = {p.name: p.free_index for p in s.parameters.isospin}
        xs[ish["isospin_shift_1_4"]] = -0.83
        xs[ish["isospin_shift_2_4"]] = -0.61
        xs[ish["isospin_mix_1"]] = 0.035
        chi_s = float(np.sum(s.calculate_chi2_rwa(xs)))
        snap = os.path.join(work, "iso_saved.azr")
        try:
            s.save_fit(snap, xs)          # verify=True re-evaluates the snapshot
            saved_ok = True
        except Exception as e:
            saved_ok = False
            check("save_fit verified the snapshot", False, str(e)[:300])
        if saved_ok:
            check("save_fit verified the snapshot", True)
            blk = AzrModel.from_file(snap).isospin
            check("the block carries the fitted shifts and mix",
                  abs(blk["shifts"][(1, 4)][1] + 0.83) < 1e-9 and abs(blk["shifts"][(2, 4)][1] + 0.61) < 1e-9
                  and abs(blk["mixes"][0][4] - 0.035) < 1e-9, str(blk))
            with azure2(snap, cwd=work, use_brune=False) as s2:
                chi2 = float(np.sum(s2.calculate_chi2_rwa(np.asarray(s2.params_rwa, float))))
            check("the snapshot reloads with the same chi-squared", abs(chi2 - chi_s) < 1e-6 * max(1.0, chi_s), f"{chi2} vs {chi_s}")

    # ---- 6. file side -----------------------------------------------------------------
    print("\n6. file side")
    m = AzrModel.from_file(os.path.join(MODELS, "iso_mix.azr"))
    iso = m.isospin
    check("block parsed", iso["reference"] == (3, 4, 1) and iso["levels"] == {1: 1.0, 2: 0.0}
          and iso["shifts"][(1, 4)] == (2, -0.86, False, -0.86, 0.0) and iso["mixes"] == [(1, 5, 1, 1.0, 0.05, False, 0.0)], str(iso))
    keys = m.engine_level_keys()
    check("one engine level per nucleus of the T=1 level", sorted(keys) == [(1, 1), (2, 1), (3, 1), (3, 2)], str(sorted(keys)))
    check("each key holds that nucleus's lines only",
          [c.pair for c in keys[(1, 1)].channels] == [1] and [c.pair for c in keys[(2, 1)].channels] == [2]
          and [c.pair for c in keys[(3, 1)].channels] == [3, 4, 5], str({k: [c.pair for c in v.channels] for k, v in keys.items()}))
    keys[(2, 1)].set_energy(4.2)
    check("set_energy on a member touches only its lines",
          [c.levelE for c in m.levels[0].channels] == [5.0, 4.2, 24.1736, 24.1736, 24.1736], str([c.levelE for c in m.levels[0].channels]))
    # renumbering: IDs 7 and 9 in the file come back as 1 and 2, the block with them
    text = open(os.path.join(MODELS, "iso_mix.azr")).read()
    m7 = AzrModel.from_file(os.path.join(MODELS, "iso_mix.azr"))
    for c in m7.levels[0].channels:
        c._set("levelID", 7)
    for c in m7.levels[1].channels:
        c._set("levelID", 9)
    m7._suffix = m7._suffix.replace("level 1 1", "level 7 1").replace("level 2 0", "level 9 0").replace("mix 1 5", "mix 7 5")
    m7.write(os.path.join(work, "renum.azr"))
    r = AzrModel.from_file(os.path.join(work, "renum.azr")).isospin
    check("renumbering keeps the block in step", r["levels"] == {1: 1.0, 2: 0.0} and r["mixes"][0][0] == 1, str(r))
    # rebuild iso_mix.azr from separate levels: one per nucleus, no block
    b = AzrModel.from_file(os.path.join(MODELS, "iso_mix.azr"))
    lv1 = b.levels[0]
    parts = {}
    for c in lv1.channels:
        parts.setdefault(b.nucleus_of_pair(c.pair), []).append(c)
    from pyazr.azrfile import AzrLevel
    b.levels = [AzrLevel(ch) for ch in parts.values()] + b.levels[1:]
    lines = b._suffix.splitlines()
    i, j = lines.index("<isospin>"), lines.index("</isospin>")
    b._suffix = "\n".join(lines[:i] + lines[j + 1:]) + "\n"
    b._renumber()
    check("no block after removing it", not b.has_isospin_block)
    joined = b.join_levels(*b.levels[:3], T=1)
    b.set_isospin_level(b.levels[1], 0)
    b.set_isospin_reference("4Li", pair=1)
    b.set_isospin_shift("4H", pair=2, value=-0.86, vary=False, centre=-0.86, sigma=0)
    b.set_isospin_shift("4He", pair=3, value=-0.64, vary=False, centre=-0.64, sigma=0)
    b.add_isospin_mix(joined, pair=5, L=1, S=1, value=0.05, vary=False, sigma=0)
    b.write(os.path.join(work, "rebuilt.azr"))
    rb = AzrModel.from_file(os.path.join(work, "rebuilt.azr"))
    check("rebuilt block equals the original's", rb.isospin == iso, f"{rb.isospin} vs {iso}")
    same_lines = sorted(c.to_line().split() for lv in rb.levels for c in lv.channels) == \
        sorted(c.to_line().split() for lv in AzrModel.from_file(os.path.join(MODELS, "iso_mix.azr")).levels for c in lv.channels)
    check("rebuilt levels equal the original's (same lines, one ID per isospin level)", same_lines and len(rb.levels) == 2)
    # soft tie: the line, and a session in which the T=1 level's derived slots are free
    sm = AzrModel.from_file(os.path.join(MODELS, "iso.azr"))
    sm.set_isospin_soft(1, 0.05, 0.02)
    check("set_isospin_soft writes the line", sm.isospin["soft"] == {1: (0.05, 0.02)}, str(sm.isospin["soft"]))
    sm.set_output_dir(os.path.join(work, "output"))
    sm.write(os.path.join(work, "soft.azr"))
    with azure2(os.path.join(work, "soft.azr"), cwd=work, use_brune=False, data_mode=False) as ss:
        free = [p.name for p in ss.parameters if p.kind in ("energy", "width") and not p.fixed]
    check("soft slots are free parameters", {"energy_2", "width_2_1", "energy_3", "width_3_1", "width_3_2", "width_3_3"} <= set(free),
          str(free))

    print()
    if failures:
        print(f"isospin: {len(failures)} FAILED: {failures}")
        return 1
    print("isospin: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
