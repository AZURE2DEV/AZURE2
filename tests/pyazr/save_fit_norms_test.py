#!/usr/bin/env python3
"""save_fit writes the fitted direct-segment norms where a run on its own
reads them, and says what that costs.

A data set's ``norm`` field (format 3: ``<value> fit [prior <centre> <width>]``)
is both the start value of a run and, through its prior, the centre of the
data set's normalization prior.  save_fit used to put the fitted norms in the
companion ``.sav`` only, so the ``_fit.azr`` run alone (CLI, GUI, a fresh
session) used the old norms and gave another chi-squared (18O 2010
cross-check: 7483 against 1872).  Now, on tests/15N_p_a (three data sets with
free norms):

  2. save_fit (norms="fitted", the default) at a vector with moved norms:
     the norm fields carry them and the prior of the one norm with an
     error (schardt, 15 %) keeps its loaded centre; the
     snapshot reopened alone gives the fit's chi-squared (rel 1e-9) and the
     fit's prior terms, and so does the AZURE2 binary run on it
     (chiSquared.out, rel 1e-6).  Between October 3 and this change the
     prior moved to the fitted norm (Total-Norm-Chi-Squared 0).
  (A Python list handed to calculate_chi2_rwa, as here, used to crash the
  process: the GIL-free binding reference-counted the converted array.)
  3. norms="nominal": the norm fields keep the loaded values (the prior
     centres), the .sav has the fitted norms, and the snapshot alone does not
     reproduce the fit (the documented trade-off); without a .sav it is
     refused.
  4. close_session=True (one engine in memory while verifying): the session
     is closed before the snapshot is reopened, the snapshot is verified and
     reproduces the fit's chi-squared, and the closed session refuses work.

Needs the compiled engine (and the binary for the CLI check of 2); skips
those parts cleanly.

Run from anywhere:  python3 tests/pyazr/save_fit_norms_test.py
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "15N_p_a")
AZR = "15N_p_a.azr"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")


def fresh_copy(dst):
    shutil.copytree(SOURCE, dst)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(dst, junk), ignore_errors=True)
        os.makedirs(os.path.join(dst, junk))


def norm_fields(path):
    """{segment key: the data set's norm field} of the project at PATH (keys in data-set order)."""
    from pyazr import Project
    p = Project.open(path)
    sets = [r for r in p.records().values() if r.record_kind == "dataset"]
    return {k: r["norm"] for k, r in enumerate(sets, 1)}


def norm_values(path):
    return {k: float(v.split()[0]) for k, v in norm_fields(path).items()}


try:
    from pyazr import azure2
except Exception as err:                                   # engine not built
    print(f"skip 2-3: engine not available ({type(err).__name__}: {err})")
    azure2 = None

binary = os.environ.get("AZURE2_BIN")
if not binary:
    cands = [c for c in glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))
             if os.path.isfile(c) and os.access(c, os.X_OK)]
    binary = max(cands, key=os.path.getmtime) if cands else None


def cli_chi2(work, azr):
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    subprocess.run([binary, "--no-gui", "--no-readline", azr], cwd=work,
                   input="1\n\n\n7\n", text=True, capture_output=True, timeout=600)
    text = open(os.path.join(work, "output", "chiSquared.out")).read()
    tot = [l.split() for l in text.splitlines() if l.startswith("Total-Chi-Squared:")][0]
    return float(tot[1]), float(tot[3])


if azure2 is not None:
    with tempfile.TemporaryDirectory() as tmp:
        work = os.path.join(tmp, "w")
        fresh_copy(work)
        with azure2(os.path.join(work, AZR), cwd=work) as s:
            norms = s.parameters.norms
            idx = [p.free_index for p in norms if not p.fixed]
            keys = [p.segment_key for p in norms if not p.fixed]
            check("three free norms", len(idx) == 3, keys)
            x = list(s.params_rwa)
            want = {}
            for k, i, f in zip(keys, idx, (0.93, 1.07, 1.11)):
                x[i] = x[i] * f
                want[k] = x[i]
            chi_fit = s.calculate_chi2_rwa(x)[0]     # x is a list: this used to segfault
            chi_old = s.calculate_chi2_rwa(s.params_rwa)[0]
            check("the moved norms change chi2", rel(chi_fit, chi_old) > 1e-3, (chi_fit, chi_old))

            # 2. ---------------------------------------------------------------
            print("2. norms='fitted' (default)")
            fit_azr, fit_sav = s.save_fit(os.path.join(work, "fit.azr"), x)
            got = norm_values(fit_azr)
            check("the norm fields carry the fitted norms",
                  all(abs(got[k] - want[k]) <= 1e-15 * abs(want[k]) for k in want), got)
            prior = norm_fields(fit_azr)[3].split("prior")[1].split()
            check("the prior of schardt (key 3) keeps its loaded centre",
                  prior == ["0.9950877877494786", "15%"], norm_fields(fit_azr)[3])
            pen_fit = s.penalties(x)
            pen_fit_norm = float(pen_fit["norm"].sum())
            check("the fit pays a norm prior", pen_fit_norm > 0.1, pen_fit)
            with azure2(fit_azr, cwd=work) as t:
                chi_alone = t.calculate_chi2_rwa(t.params_rwa)[0]
                pen = t.penalties(t.params_rwa)
            check("the snapshot alone gives the fit's chi2", rel(chi_alone, chi_fit) < 1e-9,
                  (chi_alone, chi_fit))
            print(f"        chi2 fit {chi_fit:.6f}, snapshot alone {chi_alone:.6f}, "
                  f"at the loaded norms {chi_old:.6f}; norm prior {pen_fit_norm:.6f}")
            check("and the fit's priors (centred where the fit had them)",
                  rel(float(pen["norm"].sum()), pen_fit_norm) < 1e-9
                  and float(pen["shift"].sum()) == 0.0, (pen, pen_fit))
            if binary:
                c, n = cli_chi2(work, "fit.azr")
                check("the binary run on the snapshot gives it too", rel(c, chi_fit) < 1e-6,
                      (c, chi_fit))
                check("with the fit's Total-Norm-Chi-Squared", rel(n, pen_fit_norm) < 1e-5,
                      (n, pen_fit_norm))
            else:
                print("        skip: no AZURE2 binary for the CLI check")

            # 3. ---------------------------------------------------------------
            print("3. norms='nominal'")
            nom_azr, nom_sav = s.save_fit(os.path.join(work, "nom.azr"), x, norms="nominal")
            check("the norm fields keep the loaded values",
                  norm_fields(nom_azr) == norm_fields(os.path.join(work, AZR)))
            sav = {l.split()[0]: float(l.split()[1]) for l in open(nom_sav)}
            names = {p.segment_key: p.name for p in norms}
            check("the .sav has the fitted norms",
                  all(rel(sav[names[k]], want[k]) < 1e-7 for k in want), sav)
            with azure2(nom_azr, cwd=work) as t:
                chi_nom = t.calculate_chi2_rwa(t.params_rwa)[0]
            check("alone it is not the fit (the trade-off)", rel(chi_nom, chi_fit) > 1e-3,
                  (chi_nom, chi_fit))
            try:
                s.save_fit(os.path.join(work, "bad.azr"), x, norms="nominal", param_sav=False)
                check("nominal without a .sav is refused", False, "no error")
            except ValueError:
                check("nominal without a .sav is refused",
                      not os.path.exists(os.path.join(work, "bad.azr")))

        # 4. -------------------------------------------------------------------
        print("4. close_session=True")
        s = azure2(os.path.join(work, AZR), cwd=work)
        last_azr, _ = s.save_fit(os.path.join(work, "last.azr"), x, close_session=True)
        check("the session is closed", not s.is_alive())
        check("the snapshot is written and verified (same norms as 2)",
              norm_fields(last_azr) == norm_fields(fit_azr))
        with azure2(last_azr, cwd=work) as t:
            chi_last = t.calculate_chi2_rwa(t.params_rwa)[0]
        check("alone it gives the fit's chi2", rel(chi_last, chi_fit) < 1e-9, (chi_last, chi_fit))
        try:
            s.calculate_chi2_rwa(x)
            check("a closed session refuses further work", False, "no error")
        except RuntimeError:
            check("a closed session refuses further work", True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
if azure2 is None:
    # 77 is SKIP_RETURN_CODE in tests/pyazr/CMakeLists.txt: ctest reports a
    # test whose engine part could not run as Skipped, not as Passed.
    print("save_fit norm checks skipped (no engine)")
    sys.exit(77)
print("all save_fit norm checks passed")
