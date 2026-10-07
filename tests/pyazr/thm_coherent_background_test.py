#!/usr/bin/env python3
"""THM coherent background in pyazr: cbackground= on an experiment line.

A THM-only complex amplitude c(E) = c0 + c1 E times the entrance vertex M_l
is added to the resonant HOES amplitude of a (J^pi, entrance (s,l), exit
(s',l')) combination before squaring (ThmExperiment.h ThmCoherentBackground;
docs/source/theory/thm_implementation.rst, "Coherent background").  Re and Im
of c are ordinary fit parameters, the last block of the parameter vector.
Checked on tests/18O_p_a_thm with its six R-matrix parameters freed:

  1. AzrModel: set_thm_experiment(cbackground=...) / thm_experiments /
     set_thm_cbackground round trip; what AZURE2 refuses raises and leaves the
     model unchanged (syntax, an exit pair no segment has, a J^pi the model
     lacks, entranceL=coherent in either order).
  2. The engine: the parameters (kind "cbkg", names, values, fixed flags),
     calculate_chi2_rwa == the CLI, thm_experiments reports the values; a
     background fixed at zero gives residuals and Jacobian bit for bit those
     of none.
  3. Derivatives: residual_jacobian against central differences of the
     residuals for all ten free parameters (six R-matrix, four cbkg of a
     linear term), chi2_and_grad against 2 J^T r and differences of chi2.
  4. save_fit writes the fitted values into cbackground= and reads back.

Needs the compiled engine and an AZURE2 binary for 2-4; skips cleanly without.

Run from anywhere:  python3 tests/pyazr/thm_coherent_background_test.py
"""
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "18O_p_a_thm")
AZR = "18O_p_a_thm.azr"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
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
    for f in glob.glob(os.path.join(dst, "run.log")):
        os.remove(f)


def refuses(what, fn, frag):
    try:
        fn()
    except ValueError as err:
        check(f"refused: {what} ({err})", frag in str(err), str(err))
        return
    check(f"refused: {what}", False, "no ValueError")


with tempfile.TemporaryDirectory() as tmp:
    # -- 1. AzrModel ------------------------------------------------------------
    print("1. AzrModel: cbackground=")
    proj = os.path.join(tmp, "model")
    fresh_copy(proj)
    m = AzrModel.from_file(os.path.join(proj, AZR))
    m.set_thm_experiment("E", [1, 2], cbackground="1/2+:2:0.5,0,0.5,1:linear=0.25,-0.5f,0,1e-3")
    rec = m.thm_experiments()["E"]
    check("record: canonical value", rec.get("cbackground") ==
          "1/2+:2:1/2,0,1/2,1:linear=0.25,-0.5f,0,0.001", rec)
    one = os.path.join(proj, "one.azr")
    m.write(one)
    text = open(one).read()
    check("one canonical line", "experiment[E] segments=1,2 cbackground=1/2+:2:1/2,0,1/2,1:linear="
          "0.25,-0.5f,0,0.001\n" in text, text[-200:])
    m2 = AzrModel.from_file(one)
    check("reloaded: same record", m2.thm_experiments() == m.thm_experiments())
    m2.write(os.path.join(proj, "again.azr"))
    check("untouched: byte-identical write", open(os.path.join(proj, "again.azr")).read() == text)
    m.set_thm_experiment("E", [1, 2], cbackground=["1/2+:2", "1/2+:2:1/2,0,1/2,1"])
    check("a list of terms is joined with ';'",
          m.thm_experiments()["E"]["cbackground"] == "1/2+:2;1/2+:2:1/2,0,1/2,1")
    m2.set_thm_cbackground("E", "1/2+:2=1.5,2f")
    check("set_thm_cbackground replaces the value in place",
          m2.thm_experiments()["E"]["cbackground"] == "1/2+:2=1.5,2f" and
          m2.thm_experiments()["E"]["segments"] == [1, 2])
    before = str(m2._suffix) + str(m2._prefix)
    for what, fn, frag in [
        ("no parity", lambda: m2.set_thm_experiment("E", [1], cbackground="1/2:2"), "J^pi '1/2'"),
        ("bad exit key", lambda: m2.set_thm_experiment("E", [1], cbackground="1/2+:x"), "exit pair key"),
        ("two values for linear", lambda: m2.set_thm_experiment("E", [1], cbackground="1/2+:2:linear=1,2"),
         "expected 4 values"),
        ("bad channels", lambda: m2.set_thm_experiment("E", [1], cbackground="1/2+:2:1/2,0"), "channels"),
        ("empty term", lambda: m2.set_thm_experiment("E", [1], cbackground="1/2+:2;"), "an empty term"),
        ("exit pair of no segment", lambda: m2.set_thm_experiment("E", [1], cbackground="1/2+:1"),
         "no segment of the experiment has exit pair 1"),
        ("J^pi not in the model", lambda: m2.set_thm_experiment("E", [1], cbackground="3/2-:2"),
         "no J^pi = 3/2- group"),
        ("set_thm_cbackground: bad value", lambda: m2.set_thm_cbackground("E", "1/2+:2=a,b"), "value 'a'"),
        ("entranceL=coherent after", lambda: m2.set_thm_option("entranceL", "coherent"), "entranceL=coherent"),
        ("with lineshape=on", lambda: m2.set_thm_experiment("E", [1], beam="18O", target="d", spectator="n",
                                                             Ebeam=54, lineshape=True, cbackground="1/2+:2"),
         "lineshape=on and cbackground= cannot be combined"),
    ]:
        refuses(what, fn, frag)
    check("model unchanged after the refusals", str(m2._suffix) + str(m2._prefix) == before)
    m3 = AzrModel.from_file(os.path.join(proj, AZR))
    m3.set_thm_option("entranceL", "coherent")
    refuses("entranceL=coherent before", lambda: m3.set_thm_experiment("E", [1], cbackground="1/2+:2"),
            "entranceL=coherent")
    try:
        AzrModel.from_file(os.path.join(proj, AZR)).set_thm_cbackground("E", "1/2+:2")
        check("set_thm_cbackground without the key raises KeyError", False)
    except KeyError:
        check("set_thm_cbackground without the key raises KeyError", True)

    # -- the engine -----------------------------------------------------------
    try:
        import numpy as np
        from pyazr import azure2
    except Exception as err:                                   # engine not built
        print(f"skip the engine part: engine not available ({type(err).__name__}: {err})")
        sys.exit(1 if failures else 77)
    binary = os.environ.get("AZURE2_BIN")
    if not binary:
        cands = [c for c in glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))
                 if os.path.isfile(c) and os.access(c, os.X_OK)]
        binary = max(cands, key=os.path.getmtime) if cands else None
    if binary is None:
        print("skip the engine part: no AZURE2 binary to compare against")
        sys.exit(1 if failures else 77)

    def project(name, block):
        d = os.path.join(tmp, name)
        fresh_copy(d)
        mm = AzrModel.from_file(os.path.join(d, AZR))
        for lv in mm.levels:                                   # free all six parameters
            lv.set_fixed(False)
            for c in lv.channels:
                c.channel_fixed = False
        mm.write(os.path.join(d, "run.azr"))
        if block:
            with open(os.path.join(d, "run.azr"), "a") as f:
                f.write("\n<thm>\n" + block + "\n</thm>\n")
        return d

    block = ("experiment[E] segments=1,2 background=linear "
             "cbackground=1/2+:2:linear=0.3,-0.2,0.05,0.1")
    cli_dir, py_dir = project("cli", block), project("py", block)
    proc = subprocess.run([binary, "--no-gui", "--no-readline", "run.azr"], cwd=cli_dir,
                          input="1\n\n\n7\n", text=True, capture_output=True, timeout=600)
    hit = re.findall(r"Total Chi-Squared:\s*([0-9.eE+-]+)", proc.stdout.replace("\r", "\n"))
    if not hit:
        print("FAIL: the CLI run produced no chi-squared")
        print(proc.stdout[-2000:])
        sys.exit(1)
    cli_chi2 = float(hit[-1])

    # An experiment whose name holds "norm" and "shift": its cbkg parameters
    # are neither (they were taken for both by a substring test).
    named = project("named", "experiment[normshift] segments=1,2 cbackground=1/2+:2=0.3,-0.2")
    with azure2(os.path.join(named, "run.azr"), cwd=named) as s:
        cbi = [p.free_index for p in s.parameters.cbkg]
        norms = {int(round(i)) for i in s.norm_indices()}
        shifts = {int(round(i)) for i in s.shift_indices()}
        check("cbkg of experiment[normshift]: not norm_indices / shift_indices",
              len(cbi) == 2 and not set(cbi) & norms and not set(cbi) & shifts,
              f"cbkg {cbi} norms {sorted(norms)} shifts {sorted(shifts)}")

    print("2. the engine: parameters, chi2, report, zero background")
    names = ["cbkg_E_1/2+_2_1/2,0,1/2,1_" + p for p in ("re0", "im0", "re1", "im1")]
    with azure2(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
        cb = s.parameters.cbkg
        check("four cbkg parameters, last, in order", [p.name for p in cb] == names and
              [p.index for p in cb] == list(range(len(s.parameters) - 4, len(s.parameters))),
              [p.name for p in cb])
        check("cbkg records: value, J^pi, entrance channel, exit pair",
              [p.value for p in cb] == [0.3, -0.2, 0.05, 0.1] and
              all(p.jpi == "0.5+" and p.L == 0 and p.S == 0.5 and p.pair == 2 for p in cb),
              [str(p) for p in cb])
        x = np.asarray(s.params_rwa, float)
        check("ten free parameters", x.size == 10, x.size)
        chi2 = float(np.sum(s.calculate_chi2_rwa(x)))
        check(f"calculate_chi2_rwa == CLI (rel 1e-9): {chi2!r}", rel(chi2, cli_chi2) < 1e-9,
              f"{chi2} vs {cli_chi2}")
        rep = s.thm_experiments(x)["E"]
        check("thm_experiments reports the values", rep.get("cbkg") == dict(zip(names, [0.3, -0.2, 0.05, 0.1])),
              rep.get("cbkg"))
        x2 = x.copy()
        x2[-2] = 0.5
        rep = s.thm_experiments(x2)["E"]
        check("... at the parameters asked for", rep["cbkg"][names[2]] == 0.5, rep["cbkg"])

        print("3. derivatives (R-matrix and cbkg columns)")
        r, J = s.residual_jacobian(x)
        check("residual_jacobian: residuals == residuals()", np.max(np.abs(r - s.residuals(x))) == 0.0)
        for c in range(x.size):
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (s.residuals(xp) - s.residuals(xm)) / (2 * h)
            err = np.linalg.norm(J[:, c] - fd) / max(np.linalg.norm(fd), 1e-300)
            check(f"J[:, {c}] ({s.parameters.free[c].name}) vs central differences (rel {err:.1e})",
                  err < 1e-4)
        c2, g = s.chi2_and_grad(x)
        check("chi2_and_grad value == chi2", rel(c2, chi2) < 1e-9, f"{c2} vs {chi2}")
        gj = 2.0 * J.T @ r
        check("gradient == 2 J^T r", np.allclose(g, gj, rtol=1e-8, atol=1e-10 * np.max(np.abs(gj))),
              f"max diff {np.max(np.abs(g - gj))}")
        for c in range(x.size - 4, x.size):
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (s.calculate_chi2_rwa(xp)[0] - s.calculate_chi2_rwa(xm)[0]) / (2 * h)
            check(f"d chi2/d x[{c}] vs central differences (rel {rel(g[c], fd):.1e})",
                  rel(g[c], fd) < 1e-4, f"{g[c]} vs {fd}")

        print("4. save_fit carries the values in cbackground=")
        xs = x.copy()
        xs[-4:] = [0.31, -0.19, 0.04, 0.12]
        chi_s = float(s.calculate_chi2_rwa(xs)[0])
        path, sav = s.save_fit(os.path.join(py_dir, "fit.azr"), x=xs, param_sav=False)
    fit = AzrModel.from_file(path)
    check("the snapshot's cbackground carries the values",
          fit.thm_experiments()["E"]["cbackground"] == "1/2+:2:1/2,0,1/2,1:linear=0.31,-0.19,0.04,0.12",
          fit.thm_experiments()["E"])
    with azure2(path, cwd=py_dir) as s:
        y = np.asarray(s.params_rwa, float)
        check("reloaded: same parameter vector", np.allclose(y, xs, rtol=1e-12, atol=0), f"{y} vs {xs}")
        check("reloaded: same chi2", rel(float(s.calculate_chi2_rwa(y)[0]), chi_s) < 1e-9)

    zero_dir = project("zero", "experiment[E] segments=1,2 background=linear "
                               "cbackground=1/2+:2:linear=0f,0f,0f,0f")
    none_dir = project("none", "experiment[E] segments=1,2 background=linear")
    with azure2(os.path.join(none_dir, "run.azr"), cwd=none_dir) as s:
        x = np.asarray(s.params_rwa, float)
        r0, J0 = s.residual_jacobian(x)
    with azure2(os.path.join(zero_dir, "run.azr"), cwd=zero_dir) as s:
        check("fixed zero background: no free cbkg parameter",
              len(s.parameters.cbkg) == 4 and all(p.fixed for p in s.parameters.cbkg))
        r1, J1 = s.residual_jacobian(x)
    check("fixed zero background: residuals bit for bit", np.array_equal(r0, r1))
    check("fixed zero background: Jacobian bit for bit", np.array_equal(J0, J1))

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM coherent-background checks passed")
