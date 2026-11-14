#!/usr/bin/env python3
"""The fixed-angle THM observable in pyazr (``theta=`` on an experiment line).

With theta=thmin-thmax the model of every segment of the experiment is the
HOES dsigma/dOmega of the exit pair averaged over the window (ThmAngular.h;
docs/source/theory/thm_implementation.rst, "Fixed-angle observable").
Checked on tests/7Li_p_a (7Li+p -> a+a, two entrance l per channel spin,
identical exit, 30 keV folding):

  1. Project: an experiment with theta=... round trip and canonical line;
     the refusals AZURE2 makes (malformed or reversed windows, beyond 180
     deg; together with entranceL=coherent when a session is made).
  2. The engine against the CLI with theta=50-70: evaluate().chi2
     (rel 1e-9) and the output files written by write_outputs.
  3. 0-180: 4 pi times the model is the angle-integrated model (rel 1e-12),
     the residuals (profiled norm) agree (1e-9 of their scale), and so does
     objective().residuals_and_jacobian; theta=all residuals are bit for bit those without it.
  4. Identical exit: 30-60 and 120-150 give the same model (1e-12).

Needs the compiled engine and an AZURE2 binary for 2-4; skips them cleanly.

Run from anywhere:  python3 tests/pyazr/thm_fixed_angle_test.py
"""
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "7Li_p_a")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thm_records import datasets, session_refusal, set_experiment, thm_set  # noqa: E402
AZR = "7Li_p_a.azr"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
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


with tempfile.TemporaryDirectory() as tmp:
    # -- 1. Project -----------------------------------------------------------
    print("1. Project: theta")
    proj = os.path.join(tmp, "model")
    fresh_copy(proj)
    src = os.path.join(proj, AZR)
    m = Project.open(src)
    set_experiment(m, "A", [1], theta=(50, 70))
    rec = {"theta": "50-70"}
    check("record carries theta", m.thm()["experiments"] == {"A": rec}, m.thm())
    one = os.path.join(proj, "one.azr")
    m.save(one)
    check("canonical line ends in theta=50-70", "experiment[A] theta=50-70\n" in open(one).read())
    check("reloaded: same record", Project.open(one).thm()["experiments"] == {"A": rec})
    set_experiment(m, "A", [1], theta="0-0")
    check("a single angle (0-0)", m.thm()["experiments"]["A"]["theta"] == "0-0")
    set_experiment(m, "A", [1], theta="all")
    check("theta=all", m.thm()["experiments"]["A"]["theta"] == "all")
    set_experiment(m, "A", [1], background="none")
    check("no theta: no key", "theta" not in m.thm()["experiments"]["A"])
    for value in ("70-50", "0-190", "forward", "60", "-10-20"):
        try:
            set_experiment(Project.open(src), "A", [1], theta=value)
            check(f"refused: theta={value}", False, "no ValueError")
        except ValueError as err:
            check(f"refused like AZURE2: theta={value}", "expected all or thmin-thmax" in str(err),
                  str(err))
    # theta with entranceL=coherent, in either order: the edits are taken, a session refuses the project
    m = Project.open(src)
    m.set_thm(entranceL="coherent")
    set_experiment(m, "A", [1], theta="50-70")
    bad, out = session_refusal(m, os.path.join(proj, "coherent_theta.azr"))
    check("refused by the engine: theta with entranceL=coherent", bad and "entranceL=coherent" in out, out[-400:])
    set_experiment(m, "A", [1], theta="all")
    bad, out = session_refusal(m, os.path.join(proj, "coherent_all.azr"))
    check("theta=all with entranceL=coherent is fine", not bad, out[-400:])

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
        text = open(os.path.join(d, AZR)).read()
        if block:  # the block's one experiment joined by the data set
            text = thm_set(text, block, datasets(text))
        with open(os.path.join(d, "run.azr"), "w") as f:
            f.write(text)
        return d

    def output_files(work):
        out = os.path.join(work, "output")
        return {os.path.basename(f): open(f).read()
                for f in sorted(glob.glob(os.path.join(out, "*")))
                if os.path.isfile(f) and not os.path.basename(f).startswith("intEC")}

    def model(d, x):
        with Session(os.path.join(d, "run.azr"), cwd=d) as s:
            return np.concatenate([np.asarray(v) for v in [d.model for d in s.evaluate(x).data]])

    block = "experiment[A] theta=50-70"
    cli_dir, py_dir = project("cli", block), project("py", block)
    proc = subprocess.run([binary, "calculate", "--no-gui", "run.azr"], cwd=cli_dir,
                          stdin=subprocess.DEVNULL, text=True, capture_output=True, timeout=600)
    hit = re.findall(r"Total Chi-Squared:\s*([0-9.eE+-]+)", proc.stdout.replace("\r", "\n"))
    from_cli = output_files(cli_dir)
    if not hit:
        print("FAIL: the CLI run produced no chi-squared")
        print(proc.stdout[-2000:])
        sys.exit(1)
    cli_chi2 = float(hit[-1])

    print("2. the engine against the CLI (theta=50-70)")
    with Session(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
        x = np.asarray(s.parameters.x0, float)
        chi2 = s.evaluate(x).chi2
        check(f"evaluate().chi2 == CLI (rel 1e-9): {chi2!r}",
              abs(chi2 / cli_chi2 - 1) < 1e-9, f"{chi2} vs {cli_chi2}")
        s.write_outputs(x)
    from_py = output_files(py_dir)
    for name in sorted(n for n in from_cli if n.startswith("AZUREOut_") or n in
                       ("chiSquared.out", "normalizations.out", "thm_experiments.out")):
        check(f"write_outputs: {name} identical to the CLI's", from_py.get(name) == from_cli[name])

    print("3. 0-180 is the angle-integrated model / 4 pi")
    nokey_dir = project("nokey", "experiment[A] background=none")
    all_dir = project("all", "experiment[A] theta=all")
    full_dir = project("full", "experiment[A] theta=0-180")
    m_nokey, m_full = model(nokey_dir, x), model(full_dir, x)
    worst = float(np.max(np.abs(4 * math.pi * m_full / m_nokey - 1)))
    check(f"4 pi model(0-180) == model (worst rel {worst:.1e}, {m_full.size} points)", worst < 1e-12)
    with Session(os.path.join(nokey_dir, "run.azr"), cwd=nokey_dir) as s:
        r_nokey = s.objective().residuals(x)
        _, J_nokey = s.objective().residuals_and_jacobian(x)
    with Session(os.path.join(all_dir, "run.azr"), cwd=all_dir) as s:
        check("theta=all: residuals bit for bit those without the key",
              np.array_equal(s.objective().residuals(x), r_nokey))
    with Session(os.path.join(full_dir, "run.azr"), cwd=full_dir) as s:
        r_full = s.objective().residuals(x)
        _, J_full = s.objective().residuals_and_jacobian(x)
    d_r = float(np.max(np.abs(r_full - r_nokey)) / np.max(np.abs(r_nokey)))
    check(f"residuals with the profiled norm agree ({d_r:.1e})", d_r < 1e-9)
    d_J = float(np.max(np.abs(J_full - J_nokey)) / np.max(np.abs(J_nokey)))
    check(f"objective().residuals_and_jacobian agrees ({d_J:.1e}, {J_full.shape})", d_J < 1e-6)

    print("4. identical exit: symmetric about 90 deg")
    f_dir = project("fwd", "experiment[A] theta=30-60")
    b_dir = project("bwd", "experiment[A] theta=120-150")
    m_f, m_b = model(f_dir, x), model(b_dir, x)
    worst = float(np.max(np.abs(m_f / m_b - 1)))
    check(f"model(30-60) == model(120-150) (worst rel {worst:.1e})", worst < 1e-12)
    shape = float(np.max(np.abs(4 * math.pi * m_f / m_nokey - 1)))
    check(f"and it is not the angle-integrated one (max rel {shape:.2f})", shape > 0.01)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM fixed-angle checks passed")
