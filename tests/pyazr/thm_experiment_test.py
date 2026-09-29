#!/usr/bin/env python3
"""THM experiments in pyazr: the <thm> experiment[<name>] lines.

The segments of an experiment share one profiled norm and an optional
background b0 + b1 E + b2 E^2 (ThmExperiment.h; docs/source/theory/
thm_implementation.rst, "THM experiments").  Checked on tests/18O_p_a_thm (two
THM segments with free norms) with its six R-matrix parameters freed:

  1. AzrModel: thm_experiments / set_thm_experiment / clear_thm_experiment
     round trip through a written file; everything AZURE2 refuses raises and
     leaves the model unchanged; other options and comments survive; an
     untouched model with experiment lines is written back byte for byte.
  2. The engine, linear background: calculate_chi2_rwa == the CLI (rel 1e-9),
     write_output_files reproduces the CLI's files (thm_experiments.out
     included), segment_norms is the shared n*, thm_background agrees with
     thm_experiments.out, residuals add up to chi2.
  3. Derivatives through the profile of norm and background:
     residual_jacobian against central differences of the residuals,
     chi2_and_grad against 2 J^T r and against differences of chi2.
  4. A one-segment experiment without background gives residuals identical
     to no experiment (bit for bit); a session refuses a reserved key.

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
# AzrModel is pure Python: load azrfile.py directly, as roundtrip_test does,
# so these checks run where numpy (needed by the pyazr package) is missing.
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


def refuses(what, fn):
    try:
        fn()
    except ValueError as err:
        return str(err)
    check(f"refused: {what}", False, "no ValueError")
    return None


with tempfile.TemporaryDirectory() as tmp:
    # -- 1. AzrModel ------------------------------------------------------------
    print("1. AzrModel: experiment lines")
    proj = os.path.join(tmp, "model")
    fresh_copy(proj)
    src = os.path.join(proj, AZR)
    m = AzrModel.from_file(src)
    check("no experiments in the plain project", m.thm_experiments() == {})
    m.set_thm_experiment("A", [1, 2], background="linear", beam="18O", target="d",
                         spectator="n", Ebeam=54)
    rec = {"segments": [1, 2], "background": "linear", "beam": "18O", "target": "d",
           "spectator": "n", "Ebeam": 54.0}
    check("set_thm_experiment -> thm_experiments", m.thm_experiments() == {"A": rec},
          m.thm_experiments())
    one = os.path.join(proj, "one.azr")
    m.write(one)
    text = open(one).read()
    check("one canonical line after </targetInt>",
          "</targetInt>\n<thm>\nexperiment[A] segments=1,2 background=linear beam=18O "
          "target=d spectator=n Ebeam=54\n</thm>" in text, text[-300:])
    m2 = AzrModel.from_file(one)
    check("reloaded: same record", m2.thm_experiments() == {"A": rec})
    m2.write(os.path.join(proj, "again.azr"))
    check("untouched model with an experiment: byte-identical write",
          open(os.path.join(proj, "again.azr")).read() == text)

    before = m2.thm_experiments()
    snapshot = str(m2._suffix) + str(m2._prefix)
    for what, fn, frag in [
        ("segment 3 does not exist", lambda: m2.set_thm_experiment("B", [3]), "has only 2"),
        ("segment in two experiments", lambda: m2.set_thm_experiment("B", [2]), "already in"),
        ("unknown nuclide", lambda: m2.set_thm_experiment("A", [1], beam="8Be", target="d",
                                                          spectator="n", Ebeam=5), "unknown nuclide"),
        ("partial kinematics", lambda: m2.set_thm_experiment("A", [1], beam="18O"), "all four"),
        ("bad background", lambda: m2.set_thm_experiment("A", [1], background="cubic"), "cubic"),
        ("bad name", lambda: m2.set_thm_experiment("A B", [1]), "a name is"),
        ("Ebeam <= 0", lambda: m2.set_thm_experiment("A", [1], beam="18O", target="d",
                                                     spectator="n", Ebeam=0), "Ebeam"),
    ]:
        msg = refuses(what, fn)
        if msg is not None:
            check(f"refused: {what} ({msg})", frag in msg, msg)
    check("model unchanged after the refusals",
          m2.thm_experiments() == before and str(m2._suffix) + str(m2._prefix) == snapshot)

    # Engine-side refusals seen by the parser too.
    for line, frag in [("experiment[A] segments=1 ps=file.dat", "not implemented yet"),
                       ("experiment[A] segments=1 theta=20", "not implemented yet"),
                       ("experiment[A] segments=1 foo=1", "unknown key"),
                       ("experiment[A] segments=1\nexperiment[B] segments=1", "already in"),
                       ("experiment[A] segments=1 beam=p", "all four")]:
        path = os.path.join(proj, "bad.azr")
        with open(path, "w") as f:
            f.write(open(src).read().rstrip("\n") + "\n<thm>\n" + line + "\n</thm>\n")
        msg = refuses(line, lambda: AzrModel.from_file(path).thm_experiments())
        if msg is not None:
            check(f"file refused like AZURE2: {line!r}", frag in msg, msg)

    # Other options and comments live beside experiment lines.
    m3 = AzrModel.from_file(one)
    m3.set_thm_option("vertex", "onshell")
    check("an option set beside an experiment", m3.thm_options() == {"vertex": "onshell"}
          and m3.thm_experiments() == {"A": rec})
    m3.clear_thm_experiment("A")
    check("clear_thm_experiment keeps the other options",
          m3.thm_experiments() == {} and m3.thm_options() == {"vertex": "onshell"})
    m3.clear_thm_option("vertex")
    m3.write(os.path.join(proj, "cleared.azr"))
    check("all default and no experiment: no block",
          "<thm>" not in open(os.path.join(proj, "cleared.azr")).read())
    m4 = AzrModel.from_file(one)
    m4.clear_thm_experiment("A")
    m4.write(os.path.join(proj, "cleared2.azr"))
    check("clearing the only experiment removes the block",
          "<thm>" not in open(os.path.join(proj, "cleared2.azr")).read())
    m5 = AzrModel.from_file(one)
    m5.set_thm_experiment("A", "1-2")
    check("set replaces the record (text segment list accepted)",
          m5.thm_experiments() == {"A": {"segments": [1, 2], "background": "none"}},
          m5.thm_experiments())

    # -- the engine -----------------------------------------------------------
    try:
        import numpy as np
        from pyazr import azure2
    except Exception as err:                                   # engine not built
        print(f"skip the engine part: engine not available ({type(err).__name__}: {err})")
        sys.exit(1 if failures else 0)
    binary = os.environ.get("AZURE2_BIN")
    if not binary:
        cands = [c for c in glob.glob(os.path.join(ROOT, "build*", "src", "AZURE2*"))
                 if os.path.isfile(c) and os.access(c, os.X_OK)]
        binary = max(cands, key=os.path.getmtime) if cands else None
    if binary is None:
        print("skip the engine part: no AZURE2 binary to compare against")
        sys.exit(1 if failures else 0)

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

    def output_files(work):
        out = os.path.join(work, "output")
        return {os.path.basename(f): open(f).read()
                for f in sorted(glob.glob(os.path.join(out, "*")))
                if os.path.isfile(f) and not os.path.basename(f).startswith("intEC")}

    block = "experiment[A] segments=1,2 background=linear"
    cli_dir, py_dir = project("cli", block), project("py", block)
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
    with azure2(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
        x = np.asarray(s.params_rwa, float)
        check("six free parameters", x.size == 6, x.size)
        chi2 = float(np.sum(s.calculate_chi2_rwa(x)))
        check(f"calculate_chi2_rwa == CLI (rel 1e-9): {chi2!r}", rel(chi2, cli_chi2) < 1e-9,
              f"{chi2} vs {cli_chi2}")
        r = s.residuals(x)
        check("sum r^2 == chi2", rel(float(np.sum(r ** 2)), chi2) < 1e-12)
        check("segment_chi2 sums to chi2", rel(float(s.segment_chi2(x).sum()), chi2) < 1e-12)
        bg = s.thm_background("A", x)
        norms = s.segment_norms(x)
        check("segment_norms: both segments carry the shared n*",
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
        # The model as the output shows it: calculate_rwa + b(E).
        mod = np.concatenate([np.asarray(v) for v in s.calculate_rwa(x)])
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
    with azure2(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
        x = np.asarray(s.params_rwa, float)
        chi2 = float(np.sum(s.calculate_chi2_rwa(x)))
        r, J = s.residual_jacobian(x)
        check("residual_jacobian: residuals == residuals()", np.max(np.abs(r - s.residuals(x))) == 0.0)
        for c in range(x.size):
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (s.residuals(xp) - s.residuals(xm)) / (2 * h)
            err = np.linalg.norm(J[:, c] - fd) / max(np.linalg.norm(fd), 1e-300)
            check(f"J[:, {c}] ({s.parameters.free[c].name}) vs central differences "
                  f"(rel {err:.1e})", err < 1e-4)
        c2, g = s.chi2_and_grad(x)
        check("chi2_and_grad value == chi2", rel(c2, chi2) < 1e-9, f"{c2} vs {chi2}")
        gj = 2.0 * J.T @ r
        check("gradient == 2 J^T r", np.allclose(g, gj, rtol=1e-8, atol=1e-10 * np.max(np.abs(gj))),
              f"max diff {np.max(np.abs(g - gj))}")
        for c in range(x.size):
            h = 1e-5 * (abs(x[c]) + 1.0)
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            fd = (s.calculate_chi2_rwa(xp)[0] - s.calculate_chi2_rwa(xm)[0]) / (2 * h)
            check(f"d chi2/d x[{c}] vs central differences (rel {rel(g[c], fd):.1e})",
                  rel(g[c], fd) < 1e-4 or abs(g[c] - fd) < 1e-6 * np.max(np.abs(g)), f"{g[c]} vs {fd}")

    print("4. one-segment experiments == no experiment; a session refuses a reserved key")
    plain_dir = project("plain", "")
    one_dir = project("one", "experiment[A] segments=1\nexperiment[B] segments=2")
    with azure2(os.path.join(plain_dir, "run.azr"), cwd=plain_dir) as s:
        x = np.asarray(s.params_rwa, float)
        r0, J0 = s.residual_jacobian(x)
    with azure2(os.path.join(one_dir, "run.azr"), cwd=one_dir) as s:
        r1, J1 = s.residual_jacobian(x)
        check("the experiments are reported",
              sorted(s.thm_experiments(x)) == ["A", "B"])
    check("residuals bit for bit", np.array_equal(r0, r1))
    check("Jacobian bit for bit", np.array_equal(J0, J1))
    bad_dir = project("bad", "experiment[A] segments=1,2 ps=x")
    try:
        with azure2(os.path.join(bad_dir, "run.azr"), cwd=bad_dir):
            pass
        check("a session refuses a reserved key", False, "it loaded")
    except Exception as err:
        check("a session refuses a reserved key", True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM-experiment checks passed")
