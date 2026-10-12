#!/usr/bin/env python3
"""A parameter file's FIXED level parameters reach the engine (azure2.load_sav).

AZURE2's external parameter file sets every parameter it names, fixed level
energies and widths included (AZUREParams::ReadUserParameters).  pyazr's
free-vector calls take the fixed values from the session, i.e. from the .azr,
so `x = full_rwa_from_sav(p)[~fixed]; calculate_chi2_rwa(x)` evaluated a
different model from the CLI whenever the file's fixed values differed from
the .azr's -- silently.  11B+alpha, 10-12-26: 10 fixed widths off by up to
1.6 %, chi2 9925.03 from pyazr against 9868.04 from the CLI at the same file.

Pinned on tests/identical_pp_res with one width fixed (AzrModel), scaled in
the file:

  * the free-vector-only route ignores the change (the old behaviour);
  * load_sav applies it, returns the free vector, and the engine's
    params_all_rwa carries the file's value;
  * the chi2 after load_sav is AZURE2's own mode-1 chi2 with the same file
    (when a built AZURE2 binary is found next to the module);
  * a file whose #parametrization tag is another basis is refused.

Needs the compiled engine; skips cleanly without it.
Run from anywhere:  python3 tests/pyazr/sav_fixed_params_test.py
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2, AzrModel
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(0)


def write_sav(path, names, values, basis):
    with open(path, "w") as fh:
        fh.write(f"{'#parametrization':>28s} {float(basis): .7e} {0.0: .7e}\n")
        for n, v in zip(names, values):
            fh.write(f"{n:>28s} {float(v): .7e} {0.0: .7e}\n")


with tempfile.TemporaryDirectory() as tmp:
    work = os.path.join(tmp, "identical_pp_res")
    shutil.copytree(os.path.join(ROOT, "tests", "identical_pp_res"), work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
        os.makedirs(os.path.join(work, junk))
    # The project ships with every width free: fix one, so the file has a
    # fixed level parameter to carry.
    model = AzrModel.from_file(os.path.join(work, "identical_pp_res.azr"))
    model.levels[1].channels[0].channel_fixed = True
    azr = os.path.join(work, "fixed_width.azr")
    model.write(azr)
    sav = os.path.join(work, "fixed.sav")

    with azure2(azr, cwd=work) as m:
        names = [p.name for p in m.parameters]
        fixed = np.asarray(m.fixed_params, bool)
        full0 = np.asarray(m.sess.params_all_rwa(), float)
        x0 = np.asarray(m.params_rwa, float)
        chi2_azr = m.calculate_chi2_rwa(x0)[0]

        cand = [i for i, p in enumerate(m.parameters)
                if p.kind == "width" and fixed[i] and full0[i] != 0.0]
        check("the model has a fixed nonzero width to move", bool(cand))
        if not cand:
            sys.exit(1)
        k = cand[0]
        full = full0.copy()
        full[k] *= 0.7
        write_sav(sav, names, full, m.basis)
        print(f"  (moved fixed {names[k]}: {full0[k]:.6g} -> {full[k]:.6g})")

        print("1. the free vector alone cannot carry a fixed value")
        x_old = m.full_rwa_from_sav(sav)[~fixed]
        chi2_old = m.calculate_chi2_rwa(x_old)[0]
        # (the file holds 7 digits, so the free values move by rounding only)
        check("free-vector-only route = the .azr's model",
              np.isclose(chi2_old, chi2_azr, rtol=1e-6), f"{chi2_old} vs {chi2_azr}")

        print("\n2. load_sav takes the file's fixed values into the engine")
        x = m.load_sav(sav)
        check("returns the free vector", np.array_equal(x, x_old))
        check("params_rwa set to it", np.array_equal(np.asarray(m.params_rwa), x))
        live = np.asarray(m.sess.params_all_rwa(), float)
        want = m.read_sav(sav)[names[k]]
        check("engine holds the file's fixed value", live[k] == want,
              f"{live[k]} vs {want}")
        chi2_new = m.calculate_chi2_rwa(x)[0]
        check("the chi2 moves with it", not np.isclose(chi2_new, chi2_old, rtol=1e-6),
              f"{chi2_new} vs {chi2_old}")
        print(f"        chi2: .azr fixed values {chi2_old:.6f}, file's {chi2_new:.6f}")

        print("\n3. a file in another basis is refused")
        other = os.path.join(work, "other_basis.sav")
        write_sav(other, names, full, 0 if m.basis != 0 else 1)
        try:
            m.load_sav(other)
            check("refused", False, "no error")
        except ValueError:
            check("refused", True)

    print("\n4. the same file through AZURE2 itself")
    binary = None
    for cand_bin in sorted(glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))):
        if os.path.isfile(cand_bin) and os.access(cand_bin, os.X_OK):
            binary = cand_bin
            break
    if binary is None:
        print("  skip  no AZURE2 binary to cross-check against")
    else:
        # The pyazr defaults: Brune, ignore externals, long-wavelength on.
        subprocess.run([binary, "--no-gui", "--no-readline", "--use-brune",
                        "--ignore-externals", "fixed_width.azr"],
                       cwd=work, input="1\nfixed.sav\n\n7\n", text=True,
                       capture_output=True, timeout=900)
        total = None
        out = os.path.join(work, "output", "chiSquared.out")
        if os.path.exists(out):
            for line in open(out):
                if line.startswith("Total-Chi-Squared:"):
                    total = float(line.split()[1])
        check("AZURE2 produced a total", total is not None, "no chiSquared.out line")
        if total is not None:
            print(f"        AZURE2 {total}, pyazr {chi2_new:.6f}")
            # chiSquared.out prints 6 significant digits.
            check("pyazr after load_sav = AZURE2 mode 1",
                  np.isclose(total, chi2_new, rtol=2e-6), f"{total} vs {chi2_new}")

print()
if failures:
    print(f"{len(failures)} check(s) failed: {', '.join(failures)}")
    sys.exit(1)
print("all checks passed")
