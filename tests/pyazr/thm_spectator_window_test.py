#!/usr/bin/env python3
"""The spectator-momentum window of a THM experiment in pyazr (``ps=``).

The HOES model at E is the average of the cross section over the spectator
momentum p_s in the window, weight |phi(p_s)|^2 p_s^2, each node adding
T_s = p_s^2/2mu_sx to E + B in the entrance vertex (ThmLineshape.h
ThmSpectatorWindow; docs/source/theory/thm_implementation.rst,
"Spectator-momentum window").  Checked on tests/18O_p_a_thm, whose reaction
2H(18O,a15N)n has a deuteron Trojan horse (x = p, s = n):

  1. AzrModel: set_thm_experiment(..., ps=..., psNodes=...) round trip; the
     refusals AZURE2 makes (no kinematics, bad values, psNodes alone, a window
     together with spectatorEnergy for the same pair).
  2. The engine against the CLI with a Hulthen window: calculate_chi2_rwa
     (rel 1e-9) and the output files written by write_output_files.
  3. thm_vertex: the Gauss-Legendre nodes and Hulthen weights, T_s and rho
     against an independent evaluation here (1e-12), the window-averaged
     |M_0|^2 against M_0 = B sin(rho)/rho - cos(rho) to 1e-5 of its largest
     value (the engine's forward difference for j_l' costs ~1e-6); a delta
     experiment has one node.
  4. Linearity: a two-node window equals w_1 m(T_1) + w_2 m(T_2) from two
     sessions with spectatorEnergy = T_k (rel 1e-12); ps=delta residuals are
     bit for bit those without the key.

Needs the compiled engine and an AZURE2 binary for 2-4; skips them cleanly.

Run from anywhere:  python3 tests/pyazr/thm_spectator_window_test.py
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
SOURCE = os.path.join(ROOT, "tests", "18O_p_a_thm")
AZR = "18O_p_a_thm.azr"
KIN = "beam=18O target=d spectator=n Ebeam=54"

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
# AzrModel is pure Python: load azrfile.py directly, so these checks run
# where numpy (needed by the pyazr package) is missing.
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


# Independent kinematics: mu_sx from the built-in nuclear masses of p and n
# (AMU as the engine's table), mu_xA from the .azr pair (the engine's uconv,
# hbar c: include/Constants.h).
AMU, UCONV, HBARC = 931.49410242, 931.4940880, 197.32696310
MU_SX = 1.0072764675 * 1.0086649159 / (1.0072764675 + 1.0086649159) * AMU
MU_XA = 1.00727647 * 17.99477097 / (1.00727647 + 17.99477097) * UCONV
B, RADIUS = 2.224566, 5.1
HA, HB = 0.2317, 1.202                                         # Hulthen, fm^-1


with tempfile.TemporaryDirectory() as tmp:
    # -- 1. AzrModel ----------------------------------------------------------
    print("1. AzrModel: ps and psNodes")
    proj = os.path.join(tmp, "model")
    fresh_copy(proj)
    src = os.path.join(proj, AZR)
    m = AzrModel.from_file(src)
    m.set_thm_experiment("A", [1, 2], beam="18O", target="d", spectator="n", Ebeam=54,
                         ps="hulthen:0-40", psNodes=8)
    rec = {"segments": [1, 2], "background": "none", "beam": "18O", "target": "d",
           "spectator": "n", "Ebeam": 54.0, "ps": "hulthen:0-40", "psNodes": 8}
    check("record carries ps and psNodes", m.thm_experiments() == {"A": rec}, m.thm_experiments())
    one = os.path.join(proj, "one.azr")
    m.write(one)
    check("canonical line ends in ps=... psNodes=...",
          f"experiment[A] segments=1,2 {KIN} ps=hulthen:0-40 psNodes=8\n" in open(one).read())
    check("reloaded: same record", AzrModel.from_file(one).thm_experiments() == {"A": rec})
    m = AzrModel.from_file(src)
    m.set_thm_experiment("A", [1, 2], ps="delta")
    check("ps=delta needs no kinematics", m.thm_experiments()["A"].get("ps") == "delta")
    for kw, frag in [(dict(ps="hulthen:0-40"), "needs the kinematics"),
                     (dict(ps="hulthen:40-20", beam="18O", target="d", spectator="n", Ebeam=54),
                      "expected delta, hulthen"),
                     (dict(ps="hulthen:1.2,0.2:0-40", beam="18O", target="d", spectator="n",
                           Ebeam=54), "Hulthen a,b"),
                     (dict(ps="gauss:0:0-40", beam="18O", target="d", spectator="n", Ebeam=54),
                      "FWHM"),
                     (dict(psNodes=8), "psNodes= needs a ps window"),
                     (dict(ps="hulthen:0-40", psNodes=0, beam="18O", target="d", spectator="n",
                           Ebeam=54), "1 to 64")]:
        try:
            AzrModel.from_file(src).set_thm_experiment("A", [1, 2], **kw)
            check(f"refused: {kw}", False, "no ValueError")
        except ValueError as err:
            check(f"refused like AZURE2: {frag}", frag in str(err), str(err))
    for opt in ("spectatorEnergy", "spectatorEnergy[1]"):
        m = AzrModel.from_file(src)
        m.set_thm_option(opt, 0.4)
        try:
            m.set_thm_experiment("A", [1, 2], beam="18O", target="d", spectator="n", Ebeam=54,
                                 ps="gauss:60:0-40")
            check(f"refused: window with {opt}", False, "no ValueError")
        except ValueError as err:
            check(f"refused: window with {opt}", "spectatorEnergy both set" in str(err), str(err))
    path = os.path.join(proj, "bad.azr")
    with open(path, "w") as f:
        f.write(open(src).read().rstrip("\n") + "\n<thm>\nexperiment[A] segments=1 theta=1\n</thm>\n")
    try:
        AzrModel.from_file(path).thm_experiments()
        check("theta still reserved", False, "no ValueError")
    except ValueError as err:
        check("theta still reserved", "reserved for a later version" in str(err), str(err))

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
        text = open(os.path.join(d, AZR)).read()
        if block:
            text = text.rstrip("\n") + "\n<thm>\n" + block + "\n</thm>\n"
        with open(os.path.join(d, "run.azr"), "w") as f:
            f.write(text)
        return d

    def output_files(work):
        out = os.path.join(work, "output")
        return {os.path.basename(f): open(f).read()
                for f in sorted(glob.glob(os.path.join(out, "*")))
                if os.path.isfile(f) and not os.path.basename(f).startswith("intEC")}

    block = f"experiment[A] segments=1,2 {KIN} ps=hulthen:0-40"
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

    print("2. the engine against the CLI")
    with azure2(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
        x = np.asarray(s.params_rwa, float)
        chi2 = float(np.sum(s.calculate_chi2_rwa(x)))
        check(f"calculate_chi2_rwa == CLI (rel 1e-9): {chi2!r}", rel(chi2, cli_chi2) < 1e-9,
              f"{chi2} vs {cli_chi2}")
        s.write_output_files(x)

        print("3. thm_vertex against an independent evaluation")
        grid = np.linspace(0.45, 0.95, 11)
        r = s.thm_vertex("A", grid, x)
        t, w = np.polynomial.legendre.leggauss(16)
        p = 20.0 * (t + 1.0)
        q2 = (p / HBARC) ** 2
        ww = w * (1.0 / (HA ** 2 + q2) - 1.0 / (HB ** 2 + q2)) ** 2 * p * p
        ww /= ww.sum()
        check("16 nodes on [0, 40] MeV/c", r["p_s"].size == 16
              and np.max(np.abs(r["p_s"] - p)) < 1e-12, f"{r['p_s']} vs {p}")
        check("Hulthen weights |phi|^2 p^2, normalized (1e-12)",
              np.max(np.abs(r["weights"] - ww)) < 1e-12 and abs(r["weights"].sum() - 1) < 1e-14,
              f"{r['weights']} vs {ww}")
        check(f"mu_sx = {r['mu_sx']:.6f} MeV", rel(r["mu_sx"], MU_SX) < 1e-12, r["mu_sx"])
        check("T_s = p_s^2/2mu_sx", np.max(np.abs(r["T_s"] - p * p / (2 * MU_SX))) < 1e-12)
        check("B and radius", rel(r["B"], B) < 1e-12 and rel(r["radius"], RADIUS) < 1e-12,
              f"{r['B']} {r['radius']}")
        rho = np.sqrt(2 * MU_XA * (grid[:, None] + B + p[None, :] * p[None, :] / (2 * MU_SX))) \
            * RADIUS / HBARC
        d_rho = max(np.max(np.abs(np.asarray(row) / want - 1)) for row, want in zip(r["rho"], rho))
        check(f"rho = p_xA a/hbar c at every node ({d_rho:.1e})", d_rho < 1e-12)
        chans = r["channels"]
        check("one entrance channel, l = 0, two levels",
              len(chans) == 1 and chans[0]["l"] == 0 and len(chans[0]["levels"]) == 2, chans)
        levels = chans[0]["levels"]
        check("vertex=constant: one boundary for both levels",
              levels[0]["boundary"] == levels[1]["boundary"], levels)
        Bc = levels[0]["boundary"]

        def m0(rh):                        # (B - 1) j0 - rho j0' = B sin(rho)/rho - cos(rho)
            return Bc * np.sin(rh) / rh - np.cos(rh)

        want = (ww[None, :] * m0(rho) ** 2).sum(axis=1)
        rho0 = np.sqrt(2 * MU_XA * (grid + B)) * RADIUS / HBARC
        scale = np.max(want)
        d_avg = np.max(np.abs(levels[0]["M2"] - want)) / scale
        d_qf = np.max(np.abs(levels[0]["M2_qf"] - m0(rho0) ** 2)) / scale
        check(f"<|M_0|^2> over the window ({d_avg:.1e}) and at p_s = 0 ({d_qf:.1e})",
              d_avg < 1e-5 and d_qf < 1e-5, f"{levels[0]['M2']} vs {want}")
        check("the window changes |M_0|^2 (ratio "
              f"{np.min(levels[0]['M2'] / levels[0]['M2_qf']):.3f}-"
              f"{np.max(levels[0]['M2'] / levels[0]['M2_qf']):.3f})",
              np.max(np.abs(levels[0]["M2"] / levels[0]["M2_qf"] - 1)) > 0.05)
        try:
            s.thm_vertex("nope", grid, x)
            check("unknown experiment raises", False)
        except Exception as err:
            check(f"unknown experiment raises ({err})", "no THM experiment" in str(err))
    from_py = output_files(py_dir)
    for name in sorted(n for n in from_cli if n.startswith("AZUREOut_") or n in
                       ("chiSquared.out", "normalizations.out", "thm_experiments.out")):
        check(f"write_output_files: {name} identical to the CLI's", from_py.get(name) == from_cli[name])
    check("thm_experiments.out lists the nodes", from_cli["thm_experiments.out"].count("ps_node") == 16)

    print("4. linearity: two nodes == the weighted sum of two spectatorEnergy runs")
    two_dir = project("two", f"experiment[A] segments=1,2 {KIN} ps=gauss:50:10-40 psNodes=2")
    with azure2(os.path.join(two_dir, "run.azr"), cwd=two_dir) as s:
        x = np.asarray(s.params_rwa, float)
        m_two = np.concatenate([np.asarray(v) for v in s.calculate_rwa(x)])
        r = s.thm_vertex("A", [0.6], x)
    parts = []
    for k in range(2):
        d = project(f"se{k}", f"spectatorEnergy={float(r['T_s'][k])!r}\n"
                              f"experiment[A] segments=1,2 {KIN}")
        with azure2(os.path.join(d, "run.azr"), cwd=d) as s:
            parts.append(np.concatenate([np.asarray(v) for v in s.calculate_rwa(x)]))
    combo = r["weights"][0] * parts[0] + r["weights"][1] * parts[1]
    worst = float(np.max(np.abs(m_two / combo - 1)))
    check(f"model = w1 m(T1) + w2 m(T2) at {m_two.size} points (worst {worst:.1e}; "
          f"T = {r['T_s'][0]:.4f}, {r['T_s'][1]:.4f} MeV, w = {r['weights'][0]:.4f}, "
          f"{r['weights'][1]:.4f})", worst < 1e-12)
    nokey_dir = project("nokey", f"experiment[A] segments=1,2 {KIN}")
    delta_dir = project("delta", f"experiment[A] segments=1,2 {KIN} ps=delta")
    with azure2(os.path.join(nokey_dir, "run.azr"), cwd=nokey_dir) as s:
        r_nokey = s.residuals(x)
    with azure2(os.path.join(delta_dir, "run.azr"), cwd=delta_dir) as s:
        check("ps=delta: residuals bit for bit those without the key",
              np.array_equal(s.residuals(x), r_nokey))
        r = s.thm_vertex("A", [0.6, 0.8], x)
        lv = r["channels"][0]["levels"][0]
        check("delta: one node at p_s = 0, M2 == M2_qf",
              r["window"] == "delta" and list(r["p_s"]) == [0.0] and list(r["weights"]) == [1.0]
              and np.array_equal(lv["M2"], lv["M2_qf"]), r)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM spectator-window checks passed")
