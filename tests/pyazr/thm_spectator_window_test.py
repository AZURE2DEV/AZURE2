#!/usr/bin/env python3
"""The spectator-momentum window of a THM experiment in pyazr (``ps=``).

At fixed E the spectator direction fixes q = |p_s| = |k_sF - beta k_aA|, and
the HOES model at E is the average of the cross section over the accepted
directions with the event weight |phi(q)|^2 d cos(theta_cm) (the three-body
phase space at fixed E, = |phi|^2 q dq), each node adding T_s = q^2/2mu_sx to
E + B in the entrance vertex (ThmLineshape.h ThmSpectatorWindow;
docs/source/theory/thm_implementation.rst, "Spectator-momentum window").
Until October 2026 the weight was |phi|^2 p^2 dp; the checks of 3 and 4
changed with it.  Checked on tests/18O_p_a_thm, whose reaction 2H(18O,a15N)n
has a deuteron Trojan horse (x = p, s = n):

  1. AzrModel: set_thm_experiment(..., ps=..., psNodes=...) round trip; the
     refusals AZURE2 makes (no kinematics, bad values, psNodes alone, a window
     together with spectatorEnergy for the same pair).
  2. The engine against the CLI with a Hulthen window: calculate_chi2_rwa
     (rel 1e-9) and the output files written by write_output_files.
  3. thm_vertex: at every energy the Gauss-Legendre nodes in cos(theta_cm) on
     the part of the sphere with q in [0, 40] MeV/c, the Hulthen weights
     |phi(q)|^2, T_s and rho against an independent evaluation of the
     three-body kinematics here (1e-10), the window-averaged |M_0|^2 against
     M_0 = B sin(rho)/rho - cos(rho) to 1e-5 of its largest value (the
     engine's forward difference for j_l' costs ~1e-6); a delta experiment
     has one node.
  4. One direction: ps=hulthen:30-30 equals spectatorEnergy = 30^2/2mu_sx at
     every point (rel 1e-9); ps=delta residuals are bit for bit those without
     the key.

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
# Three-body kinematics (ThmDistortion::Setup): 18O on a deuteron target at
# 54 MeV, the Trojan horse d = p + n is the target, so x = k^_sF.k^_aA is
# -cos(theta_cm); masses of the engine's nuclide table, reduced masses in
# uconv as the engine's channels.
M_N, M_P, M_D, M_O18 = 1.0086649159, 1.0072764675, 2.0135532134, 17.9947732059
E_AA = 54.0 * M_D / (M_O18 + M_D)
K_AA = math.sqrt(2 * M_D * M_O18 / (M_D + M_O18) * UCONV * E_AA) / HBARC
KB = M_N / M_D * K_AA
MU_SF = M_N * (M_P + M_O18) / (M_N + M_P + M_O18) * UCONV
B_MASS = (M_P + M_N - M_D) * AMU


def nodes(e, pmin, pmax, n):
    """The accepted directions at E: q (MeV/c) and normalized weights."""
    import numpy as np
    ks = math.sqrt(2 * MU_SF * (E_AA - B_MASS - e)) / HBARC
    x_lo = max(-1.0, (ks * ks + KB * KB - (pmax / HBARC) ** 2) / (2 * ks * KB))
    x_hi = min(1.0, (ks * ks + KB * KB - (pmin / HBARC) ** 2) / (2 * ks * KB))
    t, w = np.polynomial.legendre.leggauss(n)
    x = 0.5 * (x_lo + x_hi) + 0.5 * (x_hi - x_lo) * t
    q = np.sqrt(ks * ks + KB * KB - 2 * ks * KB * x) * HBARC
    q2 = (q / HBARC) ** 2
    ww = w * (1.0 / (HA ** 2 + q2) - 1.0 / (HB ** 2 + q2)) ** 2
    order = np.argsort(q)
    return q[order], (ww / ww.sum())[order]


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
        f.write(open(src).read().rstrip("\n") + "\n<thm>\nexperiment[A] segments=1 theta=10-0\n</thm>\n")
    try:
        AzrModel.from_file(path).thm_experiments()
        check("a reversed theta window is refused", False, "no ValueError")
    except ValueError as err:
        check("a reversed theta window is refused", "expected all or thmin-thmax" in str(err), str(err))

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
        d_q = d_w = d_t = d_rho = 0.0
        ww_rows, rho_rows = [], []
        n_ok = True
        for i, e in enumerate(grid):
            q, ww = nodes(e, 0.0, 40.0, 16)
            got_q = np.asarray(r["p_s"][i])
            order = np.argsort(got_q)
            n_ok = n_ok and got_q.size == 16
            if got_q.size != 16:
                continue
            d_q = max(d_q, np.max(np.abs(got_q[order] - q)))
            d_w = max(d_w, np.max(np.abs(np.asarray(r["weights"][i])[order] - ww)))
            d_t = max(d_t, np.max(np.abs(np.asarray(r["T_s"][i])[order] - q * q / (2 * MU_SX))))
            rho = np.sqrt(2 * MU_XA * (e + B + q * q / (2 * MU_SX))) * RADIUS / HBARC
            d_rho = max(d_rho, np.max(np.abs(np.asarray(r["rho"][i])[order] / rho - 1)))
            ww_rows.append(ww)
            rho_rows.append(rho)
        check(f"16 nodes at every energy, q in [|k_sF - beta k_aA|, 40] MeV/c ({d_q:.1e} MeV/c)",
              n_ok and d_q < 1e-9, r["p_s"])
        check(f"weights |phi(q)|^2 d cos(theta_cm), normalized ({d_w:.1e})",
              n_ok and d_w < 1e-10 and all(abs(np.sum(w) - 1) < 1e-14 for w in r["weights"]))
        check("theta_cm reported for every node", len(r["theta_cm"]) == grid.size
              and all(np.all((np.asarray(t) >= 0) & (np.asarray(t) <= 180)) for t in r["theta_cm"]))
        check(f"mu_sx = {r['mu_sx']:.6f} MeV", rel(r["mu_sx"], MU_SX) < 1e-12, r["mu_sx"])
        check(f"T_s = p_s^2/2mu_sx ({d_t:.1e} MeV)", d_t < 1e-10)
        check("B and radius", rel(r["B"], B) < 1e-12 and rel(r["radius"], RADIUS) < 1e-12,
              f"{r['B']} {r['radius']}")
        check(f"rho = p_xA a/hbar c at every node ({d_rho:.1e})", d_rho < 1e-10)
        chans = r["channels"]
        check("one entrance channel, l = 0, two levels",
              len(chans) == 1 and chans[0]["l"] == 0 and len(chans[0]["levels"]) == 2, chans)
        levels = chans[0]["levels"]
        check("vertex=constant: one boundary for both levels",
              levels[0]["boundary"] == levels[1]["boundary"], levels)
        Bc = levels[0]["boundary"]

        def m0(rh):                        # (B - 1) j0 - rho j0' = B sin(rho)/rho - cos(rho)
            return Bc * np.sin(rh) / rh - np.cos(rh)

        want = np.array([(w * m0(rh) ** 2).sum() for w, rh in zip(ww_rows, rho_rows)])
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
    check("thm_experiments.out lists the nodes at the lowest and highest point",
          from_cli["thm_experiments.out"].count("ps_node") == 32)

    print("4. one direction: ps=hulthen:30-30 == spectatorEnergy = 30^2/2mu_sx")
    one_dir = project("one", f"experiment[A] segments=1,2 {KIN} ps=hulthen:30-30")
    with azure2(os.path.join(one_dir, "run.azr"), cwd=one_dir) as s:
        x = np.asarray(s.params_rwa, float)
        m_one = np.concatenate([np.asarray(v) for v in s.calculate_rwa(x)])
        r = s.thm_vertex("A", [0.6], x)
    d = project("se30", f"spectatorEnergy={30.0 * 30.0 / (2 * MU_SX)!r}\n"
                        f"experiment[A] segments=1,2 {KIN}")
    with azure2(os.path.join(d, "run.azr"), cwd=d) as s:
        m_se = np.concatenate([np.asarray(v) for v in s.calculate_rwa(x)])
    worst = float(np.max(np.abs(m_one / m_se - 1)))
    check(f"one node at q = 30 MeV/c, weight 1 ({r['p_s'][0]})",
          len(r["p_s"][0]) == 1 and abs(r["p_s"][0][0] - 30.0) < 1e-9 and r["weights"][0][0] == 1.0)
    check(f"model == spectatorEnergy session at {m_one.size} points (worst {worst:.1e})", worst < 1e-9)
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
              r["window"] == "delta" and [list(v) for v in r["p_s"]] == [[0.0], [0.0]]
              and [list(v) for v in r["weights"]] == [[1.0], [1.0]]
              and np.array_equal(lv["M2"], lv["M2_qf"]), r)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM spectator-window checks passed")
