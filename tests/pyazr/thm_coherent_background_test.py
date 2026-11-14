#!/usr/bin/env python3
"""THM coherent background in pyazr: cbackground= on an experiment line.

A THM-only complex amplitude c(E) = c0 + c1 E times the entrance vertex M_l
is added to the resonant HOES amplitude of a (J^pi, entrance (s,l), exit
(s',l')) combination before squaring (ThmExperiment.h ThmCoherentBackground;
docs/source/theory/thm_implementation.rst, "Coherent background").  Re and Im
of c are ordinary fit parameters, the last block of the parameter vector.
Checked on tests/18O_p_a_thm with its six R-matrix parameters freed:

  1. Project: set_thm_experiment(cbackground=...) / thm() round trip; what
     AZURE2 refuses is refused, by the edit or when a session is made, and
     the project is unchanged (syntax, an exit pair no data set has, a J^pi
     the model lacks, entranceL=coherent in either order).
  2. The engine: the parameters (kind "cbkg", names, values, fixed flags),
     evaluate().chi2 == the CLI, thm_experiments reports the values; a
     background fixed at zero gives residuals and Jacobian bit for bit those
     of none.
  3. Derivatives: objective().residuals_and_jacobian against central differences of the
     residuals for all ten free parameters (six R-matrix, four cbkg of a
     linear term), objective().value_and_grad against 2 J^T r and differences of chi2.
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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thm_records import datasets, refusal, set_experiment, thm_datasets, thm_set  # noqa: E402
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
try:
    from pyazr import Project
except Exception as err:                                   # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


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
    # -- 1. Project -------------------------------------------------------------
    print("1. Project: cbackground=")
    proj = os.path.join(tmp, "model")
    fresh_copy(proj)
    src = os.path.join(proj, AZR)
    m = Project.open(src)
    set_experiment(m, "E", [1, 2], cbackground="1/2+:2:0.5,0,0.5,1:linear=0.25,-0.5f,0,1e-3")
    rec = m.thm()["experiments"]["E"]
    check("record: the value as written", rec == {"cbackground": "1/2+:2:0.5,0,0.5,1:linear=0.25,-0.5f,0,1e-3"}, rec)
    one = os.path.join(proj, "one.azr")
    m.save(one)
    text = open(one).read()
    check("one line", "experiment[E] cbackground=1/2+:2:0.5,0,0.5,1:linear=0.25,-0.5f,0,1e-3\n" in text,
          text[text.find("<thm>"):text.find("</thm>")])
    m2 = Project.open(one)
    check("reloaded: same record", m2.thm() == m.thm())
    m2.save(os.path.join(proj, "again.azr"))
    check("untouched: byte-identical write", open(os.path.join(proj, "again.azr")).read() == text)
    set_experiment(m, "E", [1, 2], cbackground=["1/2+:2", "1/2+:2:1/2,0,1/2,1"])
    check("a list of terms is joined with ';'",
          m.thm()["experiments"]["E"]["cbackground"] == "1/2+:2;1/2+:2:1/2,0,1/2,1")
    m2.set_thm_experiment("E", cbackground="1/2+:2=1.5,2f")
    check("set_thm_experiment(cbackground=) replaces the value in place, the data sets stay",
          m2.thm()["experiments"]["E"]["cbackground"] == "1/2+:2=1.5,2f"
          and [m2[n]["thm"] for n in thm_datasets(m2)] == ["E", "E"])
    before = m2.to_bytes()
    bad = os.path.join(proj, "refused.azr")

    def refuses_both(what, make, frag):  # refused by the edit, or (rules that need the model) when a session is made
        msg = refusal(make, bad)
        lines = [ln for ln in (msg or "").splitlines() if frag in ln]
        check(f"refused: {what} ({lines[0].strip()[:110] if lines else ''})", bool(lines), msg and msg[-400:])

    def cb(value, **keys):
        return lambda: set_experiment(Project.open(one), "E", [1], cbackground=value, **keys)

    def coherent_after():
        p = Project.open(one)
        p.set_thm(entranceL="coherent")
        return p

    for what, make, frag in [
        ("no parity", cb("1/2:2"), "J^pi '1/2'"),
        ("bad exit key", cb("1/2+:x"), "exit pair key"),
        ("two values for linear", cb("1/2+:2:linear=1,2"), "expected 4 values"),
        ("bad channels", cb("1/2+:2:1/2,0"), "channels"),
        ("empty term", cb("1/2+:2;"), "an empty term"),
        ("exit pair of no data set", cb("1/2+:1"), "no segment of the experiment has exit pair 1"),
        ("J^pi not in the model", cb("3/2-:2"), "no J^pi = 3/2- group"),
        ("bad value", lambda: Project.open(one).set_thm_experiment("E", cbackground="1/2+:2=a,b"), "value 'a'"),
        ("entranceL=coherent after", coherent_after, "entranceL=coherent"),
        ("with lineshape=on", cb("1/2+:2", beam="18O", target="d", spectator="n", Ebeam=54, lineshape=True),
         "lineshape=on and cbackground= cannot be combined"),
    ]:
        refuses_both(what, make, frag)
    check("project unchanged after the refusals", m2.to_bytes() == before)

    def coherent_before():
        p = Project.open(src)
        p.set_thm(entranceL="coherent")
        return set_experiment(p, "E", [1], cbackground="1/2+:2")
    refuses_both("entranceL=coherent before", coherent_before, "entranceL=coherent")

    # -- the engine -----------------------------------------------------------
    import numpy as np
    from pyazr import Session
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
        mm = Project.open(os.path.join(d, AZR))
        for lv in mm.levels():                                 # free all six parameters
            mm.set_level(lv, fixed=False)
            for c in lv.channels:
                mm.set_channel(lv, c.pair, c.l, c.s, fixed=False)
        text = mm.text
        if block:  # the block's one experiment joined by both data sets
            text = thm_set(text, block, datasets(text))
        with open(os.path.join(d, "run.azr"), "w") as f:
            f.write(text)
        return d

    block = ("experiment[E] background=linear "
             "cbackground=1/2+:2:linear=0.3,-0.2,0.05,0.1")
    cli_dir, py_dir = project("cli", block), project("py", block)
    proc = subprocess.run([binary, "calculate", "--no-gui", "run.azr"], cwd=cli_dir,
                          stdin=subprocess.DEVNULL, text=True, capture_output=True, timeout=600)
    hit = re.findall(r"Total Chi-Squared:\s*([0-9.eE+-]+)", proc.stdout.replace("\r", "\n"))
    if not hit:
        print("FAIL: the CLI run produced no chi-squared")
        print(proc.stdout[-2000:])
        sys.exit(1)
    cli_chi2 = float(hit[-1])

    # An experiment whose name holds "norm" and "shift": its cbkg parameters
    # are neither (they were taken for both by a substring test).
    named = project("named", "experiment[normshift] cbackground=1/2+:2=0.3,-0.2")
    with Session(os.path.join(named, "run.azr"), cwd=named) as s:
        cbi = [p.free_index for p in s.parameters.of_kind("cbkg")]
        norms = {p.free_index for p in s.parameters.of_kind("norm").free}
        shifts = {p.free_index for p in s.parameters.of_kind("shift").free}
        check("cbkg of experiment[normshift]: not the free norms / shifts",
              len(cbi) == 2 and not set(cbi) & norms and not set(cbi) & shifts,
              f"cbkg {cbi} norms {sorted(norms)} shifts {sorted(shifts)}")

    print("2. the engine: parameters, chi2, report, zero background")
    names = ["cbkg_E_1/2+_2_1/2,0,1/2,1_" + p for p in ("re0", "im0", "re1", "im1")]
    with Session(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
        cb = s.parameters.of_kind("cbkg")
        check("four cbkg parameters, last, in order", [p.name for p in cb] == names and
              [p.index for p in cb] == list(range(len(s.parameters) - 4, len(s.parameters))),
              [p.name for p in cb])
        check("cbkg records: value, J^pi, entrance channel, exit pair",
              [p.value for p in cb] == [0.3, -0.2, 0.05, 0.1] and
              all(p.jpi == "1/2+" and p.L == 0 and p.S == 0.5 and p.pair == 2 for p in cb),
              [str(p) for p in cb])
        x = np.asarray(s.parameters.x0, float)
        check("ten free parameters", x.size == 10, x.size)
        chi2 = s.evaluate(x).chi2
        check(f"evaluate().chi2 == CLI (rel 1e-9): {chi2!r}", rel(chi2, cli_chi2) < 1e-9,
              f"{chi2} vs {cli_chi2}")
        rep = s.thm_experiments(x)["E"]
        check("thm_experiments reports the values", rep.get("cbkg") == dict(zip(names, [0.3, -0.2, 0.05, 0.1])),
              rep.get("cbkg"))
        x2 = x.copy()
        x2[-2] = 0.5
        rep = s.thm_experiments(x2)["E"]
        check("... at the parameters asked for", rep["cbkg"][names[2]] == 0.5, rep["cbkg"])

        print("3. derivatives (R-matrix and cbkg columns)")
        r, J = s.objective().residuals_and_jacobian(x)
        check("objective().residuals_and_jacobian: residuals == residuals()", np.max(np.abs(r - s.objective().residuals(x))) == 0.0)
        for c in range(x.size):
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (s.objective().residuals(xp) - s.objective().residuals(xm)) / (2 * h)
            err = np.linalg.norm(J[:, c] - fd) / max(np.linalg.norm(fd), 1e-300)
            check(f"J[:, {c}] ({s.parameters.free[c].name}) vs central differences (rel {err:.1e})",
                  err < 1e-4)
        c2, g = s.objective().value_and_grad(x)
        check("value_and_grad value == the objective", rel(c2, s.objective()(x)) < 1e-9, f"{c2} vs {chi2}")
        gj = 2.0 * J.T @ r
        check("gradient == 2 J^T r", np.allclose(g, gj, rtol=1e-8, atol=1e-10 * np.max(np.abs(gj))),
              f"max diff {np.max(np.abs(g - gj))}")
        for c in range(x.size - 4, x.size):
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (s.objective()(xp) - s.objective()(xm)) / (2 * h)
            check(f"d f/d x[{c}] vs central differences (rel {rel(g[c], fd):.1e})",
                  rel(g[c], fd) < 1e-4, f"{g[c]} vs {fd}")

        print("4. save_fit carries the values in cbackground=")
        xs = x.copy()
        xs[-4:] = [0.31, -0.19, 0.04, 0.12]
        chi_s = s.evaluate(xs).chi2
        path, sav = s.save_fit(os.path.join(py_dir, "fit.azr"), x=xs, param_sav=False)
    fit = Project.open(path)
    check("the snapshot's cbackground carries the values",
          fit.thm()["experiments"]["E"]["cbackground"] == "1/2+:2:1/2,0,1/2,1:linear=0.31,-0.19,0.04,0.12",
          fit.thm()["experiments"]["E"])
    with Session(path, cwd=py_dir) as s:
        y = np.asarray(s.parameters.x0, float)
        check("reloaded: same parameter vector", np.allclose(y, xs, rtol=1e-12, atol=0), f"{y} vs {xs}")
        check("reloaded: same chi2", rel(s.evaluate(y).chi2, chi_s) < 1e-9)

    zero_dir = project("zero", "experiment[E] background=linear "
                               "cbackground=1/2+:2:linear=0f,0f,0f,0f")
    none_dir = project("none", "experiment[E] background=linear")
    with Session(os.path.join(none_dir, "run.azr"), cwd=none_dir) as s:
        x = np.asarray(s.parameters.x0, float)
        r0, J0 = s.objective().residuals_and_jacobian(x)
    with Session(os.path.join(zero_dir, "run.azr"), cwd=zero_dir) as s:
        check("fixed zero background: no free cbkg parameter",
              len(s.parameters.of_kind("cbkg")) == 4 and all(p.fixed for p in s.parameters.of_kind("cbkg")))
        r1, J1 = s.objective().residuals_and_jacobian(x)
    check("fixed zero background: residuals bit for bit", np.array_equal(r0, r1))
    check("fixed zero background: Jacobian bit for bit", np.array_equal(J0, J1))

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM coherent-background checks passed")
