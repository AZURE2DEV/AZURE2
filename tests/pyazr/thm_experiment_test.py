#!/usr/bin/env python3
"""THM experiments in pyazr: the <thm> experiment[<name>] lines and the data sets that name them.

The data sets of an experiment share one profiled norm and an optional
background b0 + b1 E + b2 E^2 (ThmExperiment.h; docs/source/theory/
thm_implementation.rst, "THM experiments").  In format 3 a data set joins an
experiment by its field `thm <name>`; the experiment's keys are its <thm> line.
Checked on tests/18O_p_a_thm (two THM data sets with profiled norms) with its
six R-matrix parameters freed:

  1. Project: thm() / set_thm_experiment(datasets=) / remove_thm_experiment
     round trip through a written file; everything AZURE2 refuses raises and
     leaves the project unchanged; other options live beside an experiment; an
     untouched project with an experiment is written back byte for byte.
  2. The engine, linear background: evaluate().chi2 == the CLI (rel 1e-9),
     write_output_files reproduces the CLI's files (thm_experiments.out
     included), data[i].norm is the shared n*, thm_background agrees with
     thm_experiments.out, residuals add up to chi2.
  3. Derivatives through the profile of norm and background:
     objective().residuals_and_jacobian against central differences of the residuals,
     objective().value_and_grad against 2 J^T r and against differences of chi2.
  4. A one-data-set experiment without background gives residuals identical
     to no experiment (bit for bit); a session refuses a bad theta window.

Needs the compiled engine and an AZURE2 binary; skips cleanly without them.

Run from anywhere:  python3 tests/pyazr/thm_experiment_test.py
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
from thm_records import set_experiment, thm_datasets, thm_set  # noqa: E402
AZR = "18O_p_a_thm.azr"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
try:
    import numpy as np
    from pyazr import Session, Project
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


def refuses(what, fn):
    try:
        fn()
    except (ValueError, KeyError) as err:
        return str(err)
    check(f"refused: {what}", False, "no error")
    return None


with tempfile.TemporaryDirectory() as tmp:
    # -- 1. Project -------------------------------------------------------------
    print("1. Project: experiment lines and the data sets naming them")
    proj = os.path.join(tmp, "model")
    fresh_copy(proj)
    src = os.path.join(proj, AZR)
    m = Project.open(src)
    members = thm_datasets(m)
    check("two THM data sets, no <thm> section in the plain project", len(members) == 2 and m.thm() is None,
          (members, m.thm()))
    m.set_thm_experiment("A", datasets=members, background="linear", beam="18O", target="d",
                         spectator="n", Ebeam=54)
    rec = {"background": "linear", "beam": "18O", "target": "d", "spectator": "n", "Ebeam": "54"}
    check("set_thm_experiment -> thm()", m.thm() == {"options": {}, "experiments": {"A": rec}}, m.thm())
    check("datasets= : both data sets name it", [m[n]["thm"] for n in members] == ["A", "A"],
          [m[n]["thm"] for n in members])
    one = os.path.join(proj, "one.azr")
    m.save(one)
    text = open(one).read()
    check("one canonical line in <thm>",
          "<thm>\nexperiment[A] background=linear beam=18O target=d spectator=n Ebeam=54\n</thm>" in text,
          text[text.find("<thm>") - 50:text.find("</thm>") + 10])
    m2 = Project.open(one)
    check("reloaded: same experiment", m2.thm() == {"options": {}, "experiments": {"A": rec}})
    m2.save(os.path.join(proj, "again.azr"))
    check("untouched project with an experiment: byte-identical write",
          open(os.path.join(proj, "again.azr")).read() == text)

    def coulomb_twice(p):  # the engine's global rule: the a + A Coulomb wave counted twice
        p.set_thm_experiment("A", distortion="coulomb")
        p.set_thm(coulombIntegral=1)

    before = m2.to_bytes()
    for what, fn, frag in [
        ("a data set that does not exist", lambda: m2.set_thm_experiment("B", datasets=["nope"]), "no record 'nope'"),
        ("a data set naming no experiment line", lambda: m2.set(members[1], thm="B"), "E_THM_EXPERIMENT_UNKNOWN"),
        ("unknown nuclide", lambda: m2.set_thm_experiment("A", beam="8Be", target="d",
                                                          spectator="n", Ebeam=5), "unknown nuclide"),
        ("partial kinematics", lambda: m2.set_thm_experiment("A", target=None), "all four"),
        ("bad background", lambda: m2.set_thm_experiment("A", background="cubic"), "cubic"),
        ("bad name", lambda: m2.set_thm_experiment("A B", background="none"), "a name is"),
        ("Ebeam <= 0", lambda: m2.set_thm_experiment("A", Ebeam=0), "Ebeam"),
        ("distortion=on", lambda: m2.set_thm_experiment("A", distortion="on"), "expected none, coulomb, optical"),
        ("theta=20", lambda: m2.set_thm_experiment("A", theta="20"), "expected all or thmin-thmax"),
        ("unknown key", lambda: m2.set_thm_experiment("A", foo="1"), "unknown key"),
        ("coulombIntegral=1 with R(E)", lambda: coulomb_twice(Project.open(one)), "already contains it"),
        ("an experiment its data sets name removed", lambda: m2.remove_thm_experiment("A"),
         "E_THM_EXPERIMENT_UNKNOWN"),
    ]:
        msg = refuses(what, fn)
        if msg is not None:
            check(f"refused: {what} ({msg.splitlines()[-1][:120]})", frag in msg, msg)
    check("project unchanged after the refusals", m2.to_bytes() == before)

    # Other options live beside experiment lines.
    m3 = Project.open(one)
    m3.set_thm(vertex="onshell")
    check("an option set beside an experiment", m3.thm() == {"options": {"vertex": "onshell"}, "experiments": {"A": rec}},
          m3.thm())
    for n in members:
        m3.set(n, thm="global")
    m3.remove_thm_experiment("A")
    check("remove_thm_experiment (its data sets back to `thm global`) keeps the other options",
          m3.thm() == {"options": {"vertex": "onshell"}, "experiments": {}}, m3.thm())
    m3.set_thm(vertex=None)
    check("all default and no experiment: an empty <thm>", m3.thm() == {"options": {}, "experiments": {}}, m3.thm())
    m5 = Project.open(one)
    set_experiment(m5, "A", members, background="none")
    check("the experiment redefined with one key (the others removed)",
          m5.thm()["experiments"] == {"A": {"background": "none"}}, m5.thm())

    # -- the engine -----------------------------------------------------------
    binary = os.environ.get("AZURE2_BIN")
    if not binary:
        cands = [c for c in glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))
                 if os.path.isfile(c) and os.access(c, os.X_OK)]
        binary = max(cands, key=os.path.getmtime) if cands else None
    if binary is None:
        print("skip the engine part: no AZURE2 binary to compare against")
        sys.exit(1 if failures else 77)

    def project(name, experiments=(), bad_block=None):
        """A copy with every level parameter free and EXPERIMENTS [(name, members, keys)]; BAD_BLOCK: a <thm>
        section written as text (with both data sets in its experiment), for what Project refuses."""
        d = os.path.join(tmp, name)
        fresh_copy(d)
        mm = Project.open(os.path.join(d, AZR))
        for lv in mm.levels():                                 # free all six parameters
            mm.set_level(lv, fixed=False)
            for c in lv.channels:
                mm.set_channel(lv, c.pair, c.l, c.s, fixed=False)
        for exp, ms, keys in experiments:
            set_experiment(mm, exp, ms, **keys)
        mm.save(os.path.join(d, "run.azr"))
        if bad_block:
            text = open(os.path.join(d, "run.azr")).read()
            with open(os.path.join(d, "run.azr"), "w") as f:
                f.write(thm_set(text, bad_block, thm_datasets(mm)))
        return d

    def output_files(work):
        out = os.path.join(work, "output")
        return {os.path.basename(f): open(f).read()
                for f in sorted(glob.glob(os.path.join(out, "*")))
                if os.path.isfile(f) and not os.path.basename(f).startswith("intEC")}

    linear = [("A", [1, 2], dict(background="linear"))]
    cli_dir, py_dir = project("cli", linear), project("py", linear)
    proc = subprocess.run([binary, "--no-gui", "--no-readline", "run.azr"], cwd=cli_dir,
                          input="1\n\n\n7\n", text=True, capture_output=True, timeout=600)
    hit = re.findall(r"Total Chi-Squared:\s*([0-9.eE+-]+)", proc.stdout.replace("\r", "\n"))
    from_cli = output_files(cli_dir)
    if not hit or "thm_experiments.out" not in from_cli:
        print("FAIL: the CLI run produced no chi-squared / thm_experiments.out")
        print(proc.stdout[-2000:])
        sys.exit(1)
    cli_chi2 = float(hit[-1])
    fx = from_cli["thm_experiments.out"]
    file_val = {}
    for k, v in re.findall(r"^(norm|b0|b1)\s+(\S+)", fx, re.M):
        file_val.setdefault(k, float(v))              # the value rows, not the covariance
    print(f"CLI: chi2 {cli_chi2!r}, {file_val}")

    print("2. the engine: shared norm and linear background")
    with Session(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
        x = np.asarray(s.parameters.x0, float)
        check("six free parameters", x.size == 6, x.size)
        chi2 = s.evaluate(x).chi2
        check(f"evaluate().chi2 == CLI (rel 1e-9): {chi2!r}", rel(chi2, cli_chi2) < 1e-9,
              f"{chi2} vs {cli_chi2}")
        r = s.objective().residuals(x)
        check("sum r^2 == chi2", rel(float(np.sum(r ** 2)), chi2) < 1e-12)
        check("evaluate().data[i].chi2 sums to chi2", rel(float(np.array([d.chi2 for d in s.evaluate(x).data]).sum()), chi2) < 1e-12)
        bg = s.thm_background("A", x)
        norms = np.array([d.norm for d in s.evaluate(x).data])
        check("data[i].norm: both segments carry the shared n*",
              norms[0] == norms[1] == bg["norm"], f"{norms} vs {bg['norm']}")
        check("thm_background == thm_experiments.out (10 digits)",
              rel(bg["norm"], file_val["norm"]) < 1e-9 and rel(bg["b"][0], file_val["b0"]) < 1e-9
              and rel(bg["b"][1], file_val["b1"]) < 1e-9, f"{bg} vs {file_val}")
        check("thm_background: status, chi2, covariance shape",
              bg["status"] == "profiled" and rel(bg["chi2"], chi2) < 1e-12
              and bg["cov"].shape == (3, 3) and np.all(bg["sigma"] > 0), bg)
        try:
            s.thm_background("nope", x)
            check("unknown experiment raises KeyError", False)
        except KeyError:
            check("unknown experiment raises KeyError", True)
        # The model as the output shows it: evaluate().data[i].model + b(E).
        mod = np.concatenate([np.asarray(v) for v in [d.model for d in s.evaluate(x).data]])
        ecm = np.concatenate([np.asarray(v) for v in s.energies])
        d = np.concatenate([np.asarray(v) for v in s.cross])
        e = np.concatenate([np.asarray(v) for v in s.cross_err])
        n = bg["norm"]
        f = mod + bg["b"][0] + bg["b"][1] * ecm
        mine = np.where(e != 0, (f - d * n) / np.where(e != 0, e * n, 1), 0.0)
        check("residuals = (m + b(E) - d n)/(e n)", np.max(np.abs(mine - r)) < 1e-9 * max(1, np.max(np.abs(r))),
              f"max diff {np.max(np.abs(mine - r))}")
        s.write_output_files(x)
    from_py = output_files(py_dir)
    for name in sorted(n for n in from_cli if n.startswith("AZUREOut_") or n in
                       ("chiSquared.out", "normalizations.out", "thm_experiments.out")):
        check(f"write_output_files: {name} identical to the CLI's", from_py.get(name) == from_cli[name])

    print("3. derivatives through the profile (norm + linear background)")
    with Session(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
        x = np.asarray(s.parameters.x0, float)
        chi2 = s.evaluate(x).chi2
        r, J = s.objective().residuals_and_jacobian(x)
        check("objective().residuals_and_jacobian: residuals == residuals()", np.max(np.abs(r - s.objective().residuals(x))) == 0.0)
        for c in range(x.size):
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (s.objective().residuals(xp) - s.objective().residuals(xm)) / (2 * h)
            err = np.linalg.norm(J[:, c] - fd) / max(np.linalg.norm(fd), 1e-300)
            check(f"J[:, {c}] ({s.parameters.free[c].name}) vs central differences "
                  f"(rel {err:.1e})", err < 1e-4)
        c2, g = s.objective().value_and_grad(x)
        check("value_and_grad value == the objective", rel(c2, s.objective()(x)) < 1e-9, f"{c2} vs {chi2}")
        gj = 2.0 * J.T @ r
        check("gradient == 2 J^T r", np.allclose(g, gj, rtol=1e-8, atol=1e-10 * np.max(np.abs(gj))),
              f"max diff {np.max(np.abs(g - gj))}")
        for c in range(x.size):
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (s.objective()(xp) - s.objective()(xm)) / (2 * h)
            check(f"d f/d x[{c}] vs central differences (rel {rel(g[c], fd):.1e})",
                  rel(g[c], fd) < 1e-4 or abs(g[c] - fd) < 1e-6 * np.max(np.abs(g)), f"{g[c]} vs {fd}")

    print("4. one-data-set experiments == no experiment; a session refuses a bad theta window")
    plain_dir = project("plain")
    one_dir = project("one", [("A", [1], dict(background="none")), ("B", [2], dict(background="none"))])
    with Session(os.path.join(plain_dir, "run.azr"), cwd=plain_dir) as s:
        x = np.asarray(s.parameters.x0, float)
        r0, J0 = s.objective().residuals_and_jacobian(x)
    with Session(os.path.join(one_dir, "run.azr"), cwd=one_dir) as s:
        r1, J1 = s.objective().residuals_and_jacobian(x)
        check("the experiments are reported",
              sorted(s.thm_experiments(x)) == ["A", "B"])
    check("residuals bit for bit", np.array_equal(r0, r1))
    check("Jacobian bit for bit", np.array_equal(J0, J1))
    bad_dir = project("bad", bad_block="experiment[A] background=linear theta=x")
    try:
        with Session(os.path.join(bad_dir, "run.azr"), cwd=bad_dir):
            pass
        check("a session refuses a bad theta", False, "it loaded")
    except Exception as err:
        check("a session refuses a bad theta", True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM-experiment checks passed")
