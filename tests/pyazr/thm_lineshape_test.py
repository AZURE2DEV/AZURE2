#!/usr/bin/env python3
"""The Coulomb line shape of the THM spectator in pyazr (``lineshape=on``).

Each level's exit amplitude carries N_C = exp(pi zeta/2) (E_l - E - i G_l/2)^(-i
zeta), |N_C|^2 = exp[2 zeta arctan(2 (E_l - E)/G_l)] (Mukhamedzhanov, Kadyrov
& Pang, EPJA 56 (2020) 233, eq. 62), zeta = eta_sB - eta_0 (ThmLineshape.h;
docs/source/theory/thm_implementation.rst, "Coulomb line shape").  Checked on
tests/18O_p_a_thm with a made-up charged spectator, 18O(3He,a15N)d at 115 MeV:

  1. AzrModel: set_thm_experiment(..., lineshape=True) round trip; lineshape
     without kinematics and a bad value are refused like AZURE2 refuses them.
  2. The engine against the CLI: calculate_chi2_rwa (rel 1e-9) and the output
     files written by write_output_files (thm_experiments.out included).
  3. thm_lineshape: E_sF, eta_0, zeta and the eta_sb estimate against the
     formulas evaluated here (1e-12), the level widths against parameters.out,
     NC2 against exp[2 zeta arctan(2 (E_l - E)/G_l)].
  4. One level without folding: calculate_rwa(on) / calculate_rwa(off) equals
     the NC2 thm_lineshape reports at the data energies (1e-9); lineshape=off
     residuals are bit for bit those without the key; a formal (non-Brune)
     session and an unknown or line-shape-less experiment are refused.

Needs the compiled engine and an AZURE2 binary for 2-4; skips them cleanly.

Run from anywhere:  python3 tests/pyazr/thm_lineshape_test.py
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
KIN = "beam=18O target=3He spectator=d Ebeam=115"

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


# Independent kinematics (nuclear masses in u: the built-in table for 3He, d,
# p, 18O; the .azr pairs for p + 18O and a + 15N).
U, HBARC, ALPHA, AMU = 931.4940954, 197.3269804, 1 / 137.035999084, 931.49410242
M3HE, MD, MP, M18 = 3.0149322434, 2.0135532134, 1.0072764675, 17.9947732059
MF, MA, M15 = 1.00727647 + 17.99477097, 4.00150618, 14.99626884
E_AA = 115 * M3HE / (M18 + M3HE)
BIND = (MP + MD - M3HE) * AMU
MU_SF = MD * MF / (MD + MF) * U
MU_SB = MD * M15 / (MD + M15) * U
Q = 7.993600 - 4.013799           # entrance minus exit threshold (the .azr pairs)


def kin(E):
    esf = E_AA - BIND - E
    k = math.sqrt(2 * MU_SF * esf) / HBARC
    eta0 = 9 * ALPHA * MU_SF / (HBARC * k)
    zeta = 7 * ALPHA * MU_SB / (HBARC * k) - eta0
    vs = math.sqrt(2 * esf / MU_SF)
    mubb = MA * M15 / (MA + M15) * U
    vb = math.sqrt(2 * (E + Q) / mubb) * M15 / (MA + M15)
    eta_sb = 2 * ALPHA / max(vs, vb)
    return esf, eta0, zeta, eta_sb


with tempfile.TemporaryDirectory() as tmp:
    # -- 1. AzrModel ----------------------------------------------------------
    print("1. AzrModel: lineshape key")
    proj = os.path.join(tmp, "model")
    fresh_copy(proj)
    src = os.path.join(proj, AZR)
    m = AzrModel.from_file(src)
    m.set_thm_experiment("A", [1, 2], beam="18O", target="3He", spectator="d", Ebeam=115,
                         lineshape=True)
    rec = {"segments": [1, 2], "background": "none", "beam": "18O", "target": "3He",
           "spectator": "d", "Ebeam": 115.0, "lineshape": True}
    check("record carries lineshape", m.thm_experiments() == {"A": rec}, m.thm_experiments())
    one = os.path.join(proj, "one.azr")
    m.write(one)
    check("canonical line ends in lineshape=on",
          f"experiment[A] segments=1,2 {KIN} lineshape=on\n" in open(one).read())
    check("reloaded: same record", AzrModel.from_file(one).thm_experiments() == {"A": rec})
    try:
        AzrModel.from_file(src).set_thm_experiment("B", [1], lineshape=True)
        check("refused: lineshape without kinematics", False, "no ValueError")
    except ValueError as err:
        check(f"refused: lineshape without kinematics ({err})", "needs the kinematics" in str(err))
    for line, frag in [("experiment[A] segments=1 lineshape=yes", "expected on or off"),
                       ("experiment[A] segments=1 lineshape=on", "needs the kinematics")]:
        path = os.path.join(proj, "bad.azr")
        with open(path, "w") as f:
            f.write(open(src).read().rstrip("\n") + "\n<thm>\n" + line + "\n</thm>\n")
        try:
            AzrModel.from_file(path).thm_experiments()
            check(f"file refused: {line!r}", False, "no ValueError")
        except ValueError as err:
            check(f"file refused like AZURE2: {line!r}", frag in str(err), str(err))
    path = os.path.join(proj, "off.azr")
    with open(path, "w") as f:
        f.write(open(src).read().rstrip("\n") + f"\n<thm>\nexperiment[A] segments=1 {KIN} lineshape=off\n</thm>\n")
    check("lineshape=off: no lineshape in the record",
          "lineshape" not in AzrModel.from_file(path).thm_experiments()["A"])

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

    def project(name, block, one_level=False):
        d = os.path.join(tmp, name)
        fresh_copy(d)
        text = open(os.path.join(d, AZR)).read()
        if one_level:                                          # level 1 only, no folding
            text = re.sub(r"<targetInt>.*?</targetInt>", "<targetInt>\n</targetInt>", text, flags=re.S)
            text = "\n".join(l for l in text.split("\n")
                             if not (len(l.split()) > 30 and l.split()[2] == "8.805800"))
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

    block = f"experiment[A] segments=1,2 {KIN} lineshape=on"
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
    widths = {}
    level = None
    for line in from_cli["parameters.out"].split("\n"):
        mm = re.search(r"E_level =\s*(\S+)", line)
        if mm:
            level = float(mm.group(1))
        mm = re.search(r"G  =\s*(\S+) keV", line)
        if mm and level is not None:
            widths[level] = widths.get(level, 0.0) + float(mm.group(1)) / 1000

    print("2. the engine against the CLI")
    with azure2(os.path.join(py_dir, "run.azr"), cwd=py_dir) as s:
        x = np.asarray(s.params_rwa, float)
        chi2 = float(np.sum(s.calculate_chi2_rwa(x)))
        check(f"calculate_chi2_rwa == CLI (rel 1e-9): {chi2!r}", rel(chi2, cli_chi2) < 1e-9,
              f"{chi2} vs {cli_chi2}")
        s.write_output_files(x)

        print("3. thm_lineshape against the formulas")
        grid = np.linspace(0.45, 0.95, 11)
        r = s.thm_lineshape("A", grid, x)
        want = np.array([kin(E) for E in grid])
        # E_sF to round-off; the Sommerfeld parameters to 1e-7, the engine's
        # hbar c, u and alpha (Constants.h) being of an older CODATA than these.
        d_esf = np.max(np.abs(r["E_sF"] / want[:, 0] - 1))
        d_eta = np.max(np.abs(r["eta_0"] / want[:, 1] - 1))
        check(f"E_sF ({d_esf:.1e}), eta_0 ({d_eta:.1e})", d_esf < 1e-12 and d_eta < 1e-7,
              f"{r['E_sF']} {want[:, 0]} {r['eta_0']} {want[:, 1]}")
        check("E_aA and B", rel(r["E_aA"], E_AA) < 1e-12 and rel(r["B"], BIND) < 1e-12,
              f"{r['E_aA']} {r['B']}")
        check("one exit pair, a + 15N (key 2)", len(r["exits"]) == 1 and r["exits"][0]["pair"] == 2
              and (r["exits"][0]["Zb"], r["exits"][0]["ZB"]) == (2, 7), r["exits"])
        ex = r["exits"][0]
        d_z = np.max(np.abs(ex["zeta"] / want[:, 2] - 1))
        check(f"zeta = eta_sB - eta_0 < 0 ({d_z:.1e})", d_z < 1e-7 and np.all(ex["zeta"] < 0),
              f"{ex['zeta']} vs {want[:, 2]}")
        d_b = np.max(np.abs(ex["eta_sb"] / want[:, 3] - 1))
        check(f"eta_sb estimate ({d_b:.1e})", d_b < 1e-7,
              f"{ex['eta_sb']} vs {want[:, 3]}")
        check("two levels", [(l["jgroup"], l["level"]) for l in ex["levels"]] == [(1, 1), (1, 2)],
              ex["levels"])
        for l in ex["levels"]:
            ex_level = l["E_level"] + 7.9936
            G = widths.get(round(ex_level, 4))
            check(f"level {l['level']}: Gamma = parameters.out ({l['Gamma']:.6e} MeV)",
                  G is not None and rel(l["Gamma"], G) < 1e-7, f"{l['Gamma']} vs {G} ({widths})")
            nc2 = np.exp(2 * ex["zeta"] * np.arctan(2 * (l["E_level"] - grid) / l["Gamma"]))
            check(f"level {l['level']}: NC2 = exp[2 zeta arctan(2(E_l-E)/G_l)]",
                  np.max(np.abs(l["NC2"] / nc2 - 1)) < 1e-12)
        try:
            s.thm_lineshape("nope", grid, x)
            check("unknown experiment raises", False)
        except Exception as err:
            check(f"unknown experiment raises ({err})", "no THM experiment" in str(err))
    from_py = output_files(py_dir)
    for name in sorted(n for n in from_cli if n.startswith("AZUREOut_") or n in
                       ("chiSquared.out", "normalizations.out", "thm_experiments.out")):
        check(f"write_output_files: {name} identical to the CLI's", from_py.get(name) == from_cli[name])
    check("thm_experiments.out has the zeta row", "zeta[2]" in from_cli["thm_experiments.out"])

    print("4. one level: calculate_rwa(on)/calculate_rwa(off) = |N_C|^2")
    on_dir = project("on1", f"experiment[A] segments=1,2 {KIN} lineshape=on", one_level=True)
    off_dir = project("off1", f"experiment[A] segments=1,2 {KIN}", one_level=True)
    offk_dir = project("offk1", f"experiment[A] segments=1,2 {KIN} lineshape=off", one_level=True)
    with azure2(os.path.join(off_dir, "run.azr"), cwd=off_dir) as s:
        x = np.asarray(s.params_rwa, float)
        m_off = np.concatenate([np.asarray(v) for v in s.calculate_rwa(x)])
        r_off = s.residuals(x)
        try:
            s.thm_lineshape("A", [0.6], x)
            check("an experiment without the line shape raises", False)
        except Exception as err:
            check(f"an experiment without the line shape raises ({err})", "no line shape" in str(err))
    with azure2(os.path.join(offk_dir, "run.azr"), cwd=offk_dir) as s:
        check("lineshape=off: residuals bit for bit those without the key",
              np.array_equal(s.residuals(x), r_off))
    with azure2(os.path.join(on_dir, "run.azr"), cwd=on_dir) as s:
        m_on = np.concatenate([np.asarray(v) for v in s.calculate_rwa(x)])
        ecm = np.concatenate([np.asarray(v) for v in s.energies])
        r = s.thm_lineshape("A", ecm, x)
        nc2 = r["exits"][0]["levels"][0]["NC2"]
        worst = float(np.max(np.abs((m_on / m_off) / nc2 - 1)))
        check(f"ratio == NC2 at {ecm.size} points, {nc2.min():.3f}-{nc2.max():.3f} (worst {worst:.1e})",
              worst < 1e-9)
    try:
        with azure2(os.path.join(on_dir, "run.azr"), cwd=on_dir, use_brune=False):
            pass
        check("a formal (non-Brune) session refuses lineshape=on", False, "it loaded")
    except Exception as err:
        check("a formal (non-Brune) session refuses lineshape=on", True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM line-shape checks passed")
