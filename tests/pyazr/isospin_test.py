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
  5. the physical-parameter chi-squared is refused (NaN).

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

    print()
    if failures:
        print(f"isospin: {len(failures)} FAILED: {failures}")
        return 1
    print("isospin: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
