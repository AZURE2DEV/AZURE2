#!/usr/bin/env python3
"""save_fit writes the fitted direct-segment norms where a run on its own
reads them, and says what that costs.

A ``<segmentsData>`` norm field is both the start value of a run and the
centre of the segment's normalization prior.  save_fit used to put the
fitted norms in the companion ``.sav`` only, so the ``_fit.azr`` run alone
(CLI, GUI, a fresh session) used the old norms and gave another chi-squared
(18O 2010 cross-check: 7483 against 1872).  Now, on tests/15N_p_a (three
segments with free norms):

  1. AzrModel.segment_values / set_segment_values: the engine's segment keys,
     only the norm/shift fields of the named lines change (the rest of the
     line byte for byte), an unknown key is refused.  Pure Python.
  1b. AzrModel.prior_centres / set_prior_centres (explicit prior centres,
     "segment_N_norm prior_centre c" rows of <parameterSettings>): the block
     is created after </targetInt> when the file has none, rows are kept in
     key order, CLASSIC removes a row, a non-positive norm centre and an
     unknown key are refused.  Pure Python.
  2. save_fit (norms="fitted", the default) at a vector with moved norms:
     the segment lines carry them and the prior of the one norm with an
     error (segment 3, 15 %) keeps its loaded centre as an explicit row; the
     snapshot reopened alone gives the fit's chi-squared (rel 1e-9) and the
     fit's prior terms, and so does the AZURE2 binary run on it
     (chiSquared.out, rel 1e-6).  Between October 3 and this change the
     prior moved to the fitted norm (Total-Norm-Chi-Squared 0).
  (A Python list handed to calculate_chi2_rwa, as here, used to crash the
  process: the GIL-free binding reference-counted the converted array.)
  3. norms="nominal": the segment lines keep the loaded values (the prior
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
import importlib.util                                          # noqa: E402
_spec = importlib.util.spec_from_file_location(
    "azrfile", os.path.join(ROOT, "pyazr", "azrfile.py"))
_azrfile = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_azrfile)
AzrModel = _azrfile.AzrModel


def fresh_copy(dst):
    shutil.copytree(SOURCE, dst)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(dst, junk), ignore_errors=True)
        os.makedirs(os.path.join(dst, junk))


def segment_lines(path):
    lines, inside = [], False
    for line in open(path).read().splitlines():
        if line.strip() == "<segmentsData>":
            inside = True
        elif line.strip() == "</segmentsData>":
            inside = False
        elif inside and line.strip():
            lines.append(line)
    return lines


# 1. -------------------------------------------------------------------------
print("1. AzrModel segment values")
with tempfile.TemporaryDirectory() as tmp:
    src = os.path.join(SOURCE, AZR)
    m = AzrModel.from_file(src)
    v = m.segment_values()
    check("three segments, keys 1-3", sorted(v) == [1, 2, 3], v)
    check("norms and shifts read", abs(v[3][0] - 0.9950877877494786) < 1e-15 and v[3][1] == 0.0, v)
    m.set_segment_values({2: (0.75, None)})
    out = os.path.join(tmp, "x.azr")
    m.write(out)
    a, b = segment_lines(src), segment_lines(out)
    check("only line 2 changed", a[0] == b[0] and a[2] == b[2] and a[1] != b[1])
    check("only its norm field", a[1].split()[:8] == b[1].split()[:8]
          and a[1].split()[9:] == b[1].split()[9:] and float(b[1].split()[8]) == 0.75, b[1])
    others = open(src).read().replace(a[1], b[1])
    check("the rest of the file byte for byte", open(out).read() == others)
    try:
        m.set_segment_values({9: (1.0, None)})
        check("an unknown key is refused", False, "no error")
    except KeyError:
        check("an unknown key is refused", True)

print("1b. AzrModel prior centres")
with tempfile.TemporaryDirectory() as tmp:
    src = os.path.join(SOURCE, AZR)
    plain = os.path.join(ROOT, "tests", "7Li_p_a", "7Li_p_a.azr")
    text = open(plain).read()
    check("tests/7Li_p_a has no <parameterSettings>", "<parameterSettings>" not in text)
    m = AzrModel.from_file(plain)
    check("no centres in a classic file", m.prior_centres() == {})
    m.set_prior_centres({1: (0.0006, 0.002)})
    out = os.path.join(tmp, "p.azr")
    m.write(out)
    got = open(out).read()
    check("block created after </targetInt>",
          "</targetInt>\n<parameterSettings>\n# Prior centres" in got
          and "segment_1_norm prior_centre 0.0006\nsegment_1_energy_shift prior_centre 0.002\n"
              "</parameterSettings>\n" in got)
    check("read back", AzrModel.from_file(out).prior_centres() == {1: (0.0006, 0.002)})
    m.set_prior_centres({1: (None, _azrfile.CLASSIC)})
    check("CLASSIC removes one row", m.prior_centres() == {1: (0.0006, None)})
    m.set_prior_centres({1: (_azrfile.CLASSIC, None)})
    m.write(out)
    check("and the last one: no row and no comment left",
          AzrModel.from_file(out).prior_centres() == {}
          and "prior_centre" not in open(out).read() and "# Prior centres" not in open(out).read())
    for bad, exc in (({1: (0.0, None)}, ValueError), ({1: (-1.0, None)}, ValueError),
                     ({2: (1.0, None)}, KeyError)):
        try:
            m.set_prior_centres(bad)
            check(f"refused {bad}", False, "no error")
        except exc:
            check(f"refused {bad}", True)
    m15 = AzrModel.from_file(src)
    m15.set_prior_centres({3: (1.0, None), 1: (0.5, None)})
    rows = [l for l in m15.to_text().splitlines() if l.startswith("segment_") and "prior_centre" in l]
    check("rows in key order, inside the existing block",
          rows == ["segment_1_norm prior_centre 0.5", "segment_3_norm prior_centre 1"]
          and m15.to_text().count("<parameterSettings>") == 1, rows)

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
            snap = AzrModel.from_file(fit_azr)
            got = snap.segment_values()
            check("the segment lines carry the fitted norms",
                  all(abs(got[k][0] - want[k]) <= 1e-15 * abs(want[k]) for k in want), got)
            check("the prior of segment 3 keeps its loaded centre (explicit row)",
                  snap.prior_centres() == {3: (0.9950877877494786, None)}, snap.prior_centres())
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
            check("the segment lines keep the prior centres",
                  segment_lines(nom_azr) == segment_lines(os.path.join(work, AZR)))
            sav = {l.split()[0]: float(l.split()[1]) for l in open(nom_sav)}
            check("the .sav has the fitted norms",
                  all(rel(sav[f"segment_{k}_norm"], want[k]) < 1e-7 for k in want), sav)
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
              AzrModel.from_file(last_azr).segment_values() == AzrModel.from_file(fit_azr).segment_values())
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
    print("save_fit norm check 1 passed; 2-3 skipped (no engine)")
    sys.exit(77)
print("all save_fit norm checks passed")
