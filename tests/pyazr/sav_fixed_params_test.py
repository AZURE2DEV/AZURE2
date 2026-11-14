#!/usr/bin/env python3
"""A parameter file's FIXED level parameters reach the engine (Session(file, parameter_file=path)).

AZURE2's external parameter file sets every parameter it names, fixed level
energies and widths included (AZUREParams::ReadUserParameters).  pyazr's
free-vector calls take the fixed values from the session, i.e. from the .azr,
so `x = parameters.read(p); evaluate(x).chi2` evaluates a
different model from the CLI whenever the file's fixed values differed from
the .azr's -- silently.  11B+alpha, 10-12-26: 10 fixed widths off by up to
1.6 %, chi2 9925.03 from pyazr against 9868.04 from the CLI at the same file.

Pinned on tests/identical_pp_res with one width fixed (Project.set_channel), scaled
in the file:

  * the free-vector-only route ignores the change (the old behaviour);
  * a session made with parameter_file= applies it, its x0 is the free vector,
    and the engine's params_all_rwa carries the file's value;
  * the chi2 of that session is AZURE2's own mode-1 chi2 with the same file
    (when a built AZURE2 binary is found next to the module);
  * a file whose #parametrization tag is another basis is refused;
  * the file's names are matched against the run's as AZURE2 matches them
    (experiment::MatchParameterFileNames, bound): a name the project does not
    have warns with W_PARAMETER_FILE_NAMES and its list, a file with no name of
    the run is refused with E_PARAMETER_FILE_NO_MATCH.

Needs the compiled engine; skips cleanly without it.
Run from anywhere:  python3 tests/pyazr/sav_fixed_params_test.py

This test reads the binding itself (`_azure2.Session`), so it uses the private handle `Session._engine`.
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
    from pyazr import Session, Project
    from pyazr.parameters import read_parameter_file
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


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
    model = Project.open(os.path.join(work, "identical_pp_res.azr"))
    level = model.levels()[1]
    c = level.channels[0]
    model.set_channel(level, c.pair, c.l, c.s, fixed=True)
    azr = model.save(os.path.join(work, "fixed_width.azr"))
    sav = os.path.join(work, "fixed.sav")

    with Session(azr, cwd=work) as m:
        P = m.parameters
        names = P.names
        fixed = ~P.free_mask
        full0 = np.asarray(m._engine.params_all_rwa(), float)
        x0 = P.x0
        chi2_azr = m.evaluate(x0).chi2

        cand = [i for i, p in enumerate(P)
                if p.kind == "width" and fixed[i] and full0[i] != 0.0]
        check("the model has a fixed nonzero width to move", bool(cand))
        if not cand:
            sys.exit(1)
        k = cand[0]
        full = full0.copy()
        full[k] *= 0.7
        write_sav(sav, names, full, m.basis)
        print(f"  (moved fixed {names[k]}: {full0[k]:.6g} -> {full[k]:.6g})")

        print("1. the fit vector alone cannot carry a fixed value")
        x_old = P.read(sav)
        chi2_old = m.evaluate(x_old).chi2
        # (the file holds 7 digits, so the free values move by rounding only)
        check("free-vector-only route = the .azr's model",
              np.isclose(chi2_old, chi2_azr, rtol=1e-6), f"{chi2_old} vs {chi2_azr}")

    print("\n2. a session made with parameter_file= takes the file's fixed values into the engine")
    with Session(azr, cwd=work, parameter_file=sav) as m:
        x = m.parameters.x0
        check("its x0 is the file's free vector", np.array_equal(x, x_old))
        live = np.asarray(m._engine.params_all_rwa(), float)
        want = read_parameter_file(sav)[1][names[k]]
        check("engine holds the file's fixed value", live[k] == want,
              f"{live[k]} vs {want}")
        chi2_new = m.evaluate(x).chi2
        check("the chi2 moves with it", not np.isclose(chi2_new, chi2_old, rtol=1e-6),
              f"{chi2_new} vs {chi2_old}")
        print(f"        chi2: .azr fixed values {chi2_old:.6f}, file's {chi2_new:.6f}")
        basis = m.basis

    print("\n3. a file in another basis is refused")
    other = os.path.join(work, "other_basis.sav")
    write_sav(other, names, full, 0 if basis != 0 else 1)
    try:
        Session(azr, cwd=work, parameter_file=other).close()
        check("refused", False, "no error")
    except ValueError:
        check("refused", True)

    print("\n3b. the file's names against the run's (the engine's rule)")
    import warnings
    extra = os.path.join(work, "extra_name.sav")
    write_sav(extra, names + ["no_such_parameter"], list(full) + [1.0], basis)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        Session(azr, cwd=work, parameter_file=extra).close()
    said = " ".join(str(w.message) for w in caught)
    check("a name the project lacks warns, naming it",
          "W_PARAMETER_FILE_NAMES" in said and "no_such_parameter" in said, said[:300])
    foreign = os.path.join(work, "foreign.sav")
    write_sav(foreign, ["segment_1_norm", "other_width"], [1.0, 2.0], basis)
    try:
        Session(azr, cwd=work, parameter_file=foreign).close()
        check("a file with no name of the run is refused", False, "no error")
    except ValueError as err:
        check("a file with no name of the run is refused", "E_PARAMETER_FILE_NO_MATCH" in str(err), str(err)[:300])

    print("\n4. the same file through AZURE2 itself")
    binary = os.environ.get("AZURE2_BIN")  # ctest's: the binary this build made
    for cand_bin in [] if binary else sorted(glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))):
        if os.path.isfile(cand_bin) and os.access(cand_bin, os.X_OK):
            binary = cand_bin
            break
    if binary is None:
        print("  skip  no AZURE2 binary to cross-check against")
    else:
        # The pyazr defaults: Brune, ignore externals, long-wavelength on.
        subprocess.run([binary, "calculate", "--parameter-file=fixed.sav", "--no-gui", "--use-brune",
                        "--ignore-externals", "fixed_width.azr"],
                       cwd=work, stdin=subprocess.DEVNULL, text=True,
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
            check("pyazr with parameter_file= = AZURE2 mode 1",
                  np.isclose(total, chi2_new, rtol=2e-6), f"{total} vs {chi2_new}")

print()
if failures:
    print(f"{len(failures)} check(s) failed: {', '.join(failures)}")
    sys.exit(1)
print("all checks passed")
