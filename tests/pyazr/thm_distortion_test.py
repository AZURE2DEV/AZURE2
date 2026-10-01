#!/usr/bin/env python3
"""The THM distortion factor R(E) in pyazr (``distortion=`` on an experiment).

The model of every segment of the experiment is multiplied, before the
folding, by R(E) = rho(E)/rho(E_ref), rho = |M|^2/|M_PW|^2 (dwpw) or |M|^2
(dw), M the zero-range prior-form DWBA transfer amplitude (ThmDistortion.h;
docs/source/theory/thm_implementation.rst, "Distortion factor R(E)").
Checked on tests/18O_p_a_thm with a made-up charged spectator,
18O(3He,a15N)d at 115 MeV:

  1. AzrModel: set_thm_experiment(..., distortion=..., opticalAA=...) round
     trip and canonical line; values AZURE2 refuses raise ValueError the same
     way; the Coulomb consistency rules (coulombIntegral=1 with R(E) on a
     distorted a + A wave, a ps window with distortionRatio=dw) refused as
     AZURE2 refuses them.  Pure Python: runs without numpy.
  2. The engine against the CLI: calculate_chi2_rwa (rel 1e-9) and the output
     files written by write_output_files (thm_experiments.out included).
  3. thm_distortion: |M|^2 and M_PW at E = 0.6 MeV against the mpmath values
     of tests/reference/thm_distortion_reference.py (case he3_qf); R_model (the
     interpolated weight) against R on a grid (1e-6); R(E_ref) = 1.
  4. No folding: calculate_rwa(on)/calculate_rwa(off) equals R_model at the data
     energies (1e-12); a table (distortion=table:<file>) multiplies by its w(E);
     an experiment without distortion, an unknown one and an energy that
     leaves the spectator no energy raise.

Needs the compiled engine and an AZURE2 binary for 2-4; skips them cleanly.

Run from anywhere:  python3 tests/pyazr/thm_distortion_test.py
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
KIN = "beam=18O target=3He spectator=d Ebeam=115"
# tests/reference/thm_distortion_reference.py, case he3_qf (E, |M|^2, M_PW)
REF_HE3_QF = (0.6, 3.425574843087e+02, 6.583432872958e+01)

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


with tempfile.TemporaryDirectory() as tmp:
    # -- 1. AzrModel ----------------------------------------------------------
    print("1. AzrModel: the distortion keys")
    proj = os.path.join(tmp, "model")
    fresh_copy(proj)
    src = os.path.join(proj, AZR)
    m = AzrModel.from_file(src)
    ws = [100, 4.2, 0.7, 20, 4.2, 0.7, 0, 0, 0, 4.2]
    m.set_thm_experiment("A", [1, 2], beam="18O", target="3He", spectator="d", Ebeam=115,
                         distortion="optical", opticalAA=ws, opticalSF="plane",
                         spectatorAngle=8, distortionRef=0.7, distortionRatio="dw",
                         boundState="yukawa:2")
    rec = {"segments": [1, 2], "background": "none", "beam": "18O", "target": "3He",
           "spectator": "d", "Ebeam": 115.0, "distortion": "optical",
           "opticalAA": "100,4.2,0.7,20,4.2,0.7,0,0,0,4.2", "opticalSF": "plane",
           "spectatorAngle": "8", "distortionRef": "0.7", "distortionRatio": "dw",
           "boundState": "yukawa:2"}
    check("record carries the distortion keys", m.thm_experiments() == {"A": rec},
          m.thm_experiments())
    one = os.path.join(proj, "one.azr")
    m.write(one)
    line = (f"experiment[A] segments=1,2 {KIN} distortion=optical "
            "opticalAA=100,4.2,0.7,20,4.2,0.7,0,0,0,4.2 opticalSF=plane spectatorAngle=8 "
            "distortionRef=0.7 distortionRatio=dw boundState=yukawa:2\n")
    check("canonical line", line in open(one).read(), open(one).read()[-400:])
    check("reloaded: same record", AzrModel.from_file(one).thm_experiments() == {"A": rec})
    for kwargs, frag in [
            (dict(distortion="coulomb"), "needs the kinematics"),
            (dict(distortion="coulomb", beam="18O", target="3He", spectator="d", Ebeam=115,
                  opticalAA="plane"), "need distortion=optical"),
            (dict(distortion="optical", beam="18O", target="3He", spectator="d", Ebeam=115,
                  opticalSF=[50, 1.2, 0.6]), "ten numbers"),
            (dict(beam="18O", target="3He", spectator="d", Ebeam=115, spectatorAngle=10),
             "needs distortion=coulomb or distortion=optical")]:
        try:
            AzrModel.from_file(src).set_thm_experiment("B", [1], **kwargs)
            check(f"refused: {kwargs}", False, "no ValueError")
        except ValueError as err:
            check(f"refused like AZURE2: {frag}", frag in str(err), str(err))
    for text, frag in [
            (f"experiment[A] segments=1 {KIN} distortion=dwba", "expected none, coulomb, optical or table"),
            (f"experiment[A] segments=1 {KIN} distortion=optical opticalAA=50,0,0.6,0,0,0,0,0,0,0",
             "ten numbers"),
            (f"experiment[A] segments=1 {KIN} distortion=coulomb spectatorAngle=200", "expected qf"),
            (f"experiment[A] segments=1 {KIN} distortion=coulomb distortionRatio=pw", "expected dwpw or dw"),
            (f"experiment[A] segments=1 {KIN} distortion=coulomb boundState=hulthen",
             "expected whittaker or yukawa"),
            (f"experiment[A] segments=1 {KIN} theta=5", "expected all or thmin-thmax"),
            ("experiment[A] segments=1 distortionRef=0.6", "needs distortion=coulomb")]:
        path = os.path.join(proj, "bad.azr")
        with open(path, "w") as f:
            f.write(open(src).read().rstrip("\n") + "\n<thm>\n" + text + "\n</thm>\n")
        try:
            AzrModel.from_file(path).thm_experiments()
            check(f"file refused: {text!r}", False, "no ValueError")
        except ValueError as err:
            check(f"file refused like AZURE2: {frag}", frag in str(err), str(err))

    # Coulomb consistency (CheckThmCoulombConsistency, CheckThmExperiments):
    # C_l with R(E) on a distorted a + A wave counts the x-A Coulomb twice; a
    # ps window with distortionRatio=dw counts |phi|^2 twice.
    print("1b. AzrModel: the Coulomb consistency rules")
    kin = dict(beam="18O", target="3He", spectator="d", Ebeam=115)
    mc = AzrModel.from_file(src)
    mc.set_thm_option("coulombIntegral", True)
    for kwargs in (dict(distortion="coulomb"), dict(distortion="optical", opticalSF="plane"),
                   dict(distortion="optical", opticalAA=ws)):
        try:
            mc.set_thm_experiment("C", [1, 2], **kin, **kwargs)
            check(f"coulombIntegral=1 refused with {kwargs}", False, "no ValueError")
        except ValueError as err:
            check(f"coulombIntegral=1 refused with {kwargs}", "counted twice" in str(err), str(err))
    mc.set_thm_experiment("C", [1, 2], **kin, distortion="optical", opticalAA="plane")
    check("coulombIntegral=1 with opticalAA=plane accepted",
          mc.thm_experiments()["C"]["opticalAA"] == "plane")
    mr = AzrModel.from_file(src)
    mr.set_thm_experiment("C", [1, 2], **kin, distortion="coulomb")
    before = mr.thm_options()
    try:
        mr.set_thm_option("coulombIntegral", True)
        check("set_thm_option(coulombIntegral) refused with R(E)", False, "no ValueError")
    except ValueError as err:
        check("set_thm_option(coulombIntegral) refused with R(E), model unchanged",
              "counted twice" in str(err) and mr.thm_options() == before, str(err))
    try:
        AzrModel.from_file(src).set_thm_experiment("C", [1, 2], **kin, ps="hulthen:0-30",
                                                   distortion="coulomb", distortionRatio="dw")
        check("ps window with distortionRatio=dw refused", False, "no ValueError")
    except ValueError as err:
        check("ps window with distortionRatio=dw refused", "distortionRatio=dwpw" in str(err), str(err))
    AzrModel.from_file(src).set_thm_experiment("C", [1, 2], **kin, ps="hulthen:0-30",
                                               distortion="coulomb", distortionRatio="dwpw")
    check("ps window with distortionRatio=dwpw accepted", True)

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

    def project(name, block, fold=True):
        d = os.path.join(tmp, name)
        fresh_copy(d)
        text = open(os.path.join(d, AZR)).read()
        if not fold:
            text = re.sub(r"<targetInt>.*?</targetInt>", "<targetInt>\n</targetInt>", text, flags=re.S)
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

    block = f"experiment[A] segments=1,2 {KIN} distortion=coulomb spectatorAngle=30"
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
    from_py = output_files(py_dir)
    for name in sorted(n for n in from_cli if n.startswith("AZUREOut_") or n in
                       ("chiSquared.out", "normalizations.out", "thm_experiments.out")):
        check(f"write_output_files: {name} identical to the CLI's", from_py.get(name) == from_cli[name])
    check("thm_experiments.out has the distortion rows",
          from_cli["thm_experiments.out"].count("distortion_point") == 3)

    print("3. thm_distortion")
    qf_dir = project("qf", f"experiment[A] segments=1,2 {KIN} distortion=coulomb")
    with azure2(os.path.join(qf_dir, "run.azr"), cwd=qf_dir) as s:
        r = s.thm_distortion("A", [REF_HE3_QF[0]])
        check(f"kind coulomb, theta_cm 180 (Trojan horse = target, qf): {r['kind']}, {r['theta_cm']}",
              r["kind"] == "coulomb" and abs(r["theta_cm"][0] - 180) < 1e-12 and r["x"][0] == 1.0)
        check(f"|M|^2 at 0.6 MeV vs mpmath ({r['M2'][0]:.10e})", rel(r["M2"][0], REF_HE3_QF[1]) < 2e-5,
              f"{r['M2'][0]} vs {REF_HE3_QF[1]}")
        check(f"M_PW^2 at 0.6 MeV vs mpmath ({r['M2_PW'][0]:.10e})",
              rel(r["M2_PW"][0], REF_HE3_QF[2]**2) < 4e-6, f"{r['M2_PW'][0]} vs {REF_HE3_QF[2]**2}")
        grid = np.linspace(0.40, 1.00, 61) + 0.0037
        r = s.thm_distortion("A", grid)
        worst = float(np.max(np.abs(r["R_model"] / r["R"] - 1)))
        check(f"R_model (interpolated) == R on a grid ({worst:.1e})", worst < 1e-6)
        check("R(E_ref) = 1", abs(s.thm_distortion("A", [r["E_ref"]])["R"][0] - 1) < 1e-14)
        check(f"R = rho/rho(E_ref), rho = |M|^2/|M_PW|^2",
              np.max(np.abs(r["R"] * (r["M2_PW"] / r["M2"]) / np.mean(r["R"] * r["M2_PW"] / r["M2"]) - 1))
              < 1e-12)
        for bad, frag in [(("nope", [0.6]), "no THM experiment"),
                          (("A", [12.0]), "no energy left")]:
            try:
                s.thm_distortion(*bad)
                check(f"raises: {bad}", False)
            except Exception as err:
                check(f"raises ({err})", frag in str(err), str(err))

    print("4. no folding: the model ratio is R_model; a table is its w(E)")
    table = os.path.join(tmp, "w.dat")
    with open(table, "w") as f:
        f.write("# E w\n0.4 2.0\n1.0 0.5\n")
    runs = {"off": f"experiment[A] segments=1,2 {KIN}",
            "on": f"experiment[A] segments=1,2 {KIN} distortion=coulomb distortionRatio=dw",
            "tab": f"experiment[A] segments=1,2 distortion=table:{table}"}
    model = {}
    for name, blk in runs.items():
        d = project("nf_" + name, blk, fold=False)
        with azure2(os.path.join(d, "run.azr"), cwd=d) as s:
            x = np.asarray(s.params_rwa, float)
            model[name] = np.concatenate([np.asarray(v) for v in s.calculate_rwa(x)])
            ecm = np.concatenate([np.asarray(v) for v in s.energies])
            if name == "off":
                try:
                    s.thm_distortion("A", [0.6])
                    check("an experiment without distortion raises", False)
                except Exception as err:
                    check(f"an experiment without distortion raises ({err})", "has no distortion" in str(err))
            else:
                model[name + "_R"] = s.thm_distortion("A", ecm)["R_model"]
    worst = float(np.max(np.abs(model["on"] / model["off"] / model["on_R"] - 1)))
    check(f"coulomb: calculate_rwa(on)/calculate_rwa(off) == R_model at {ecm.size} points "
          f"({model['on_R'].min():.4f}-{model['on_R'].max():.4f}, worst {worst:.1e})", worst < 1e-12)
    w = 2.0 * np.exp((ecm - 0.4) / 0.6 * np.log(0.25))
    worst = float(np.max(np.abs(model["tab"] / model["off"] / w - 1)))
    check(f"table: ratio == w(E) log-linear ({worst:.1e})", worst < 1e-12)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM distortion checks passed")
