#!/usr/bin/env python3
"""Self-test of the azure2-fit scripts on tests/upgrade/T3 (8Be, synthetic, 8 free parameters, 498 data rows) in a
temporary directory: lsq_fit.py, engine_fit.py, convert_sav_names.py.

T3_start.azr has several local minima (objective 755.22, 591.83, 580.33, 509.35); the checks:
  1. lsq_fit.py from the start ends where the engine's Levenberg-Marquardt ends (755.2236, 1e-6 relative);
  2. lsq_fit.py with 8 starts (jitter 3 parabolic errors, seed 1) reaches the engine's MIGRAD with the gradient
     (509.3479) or goes below it by less than 1e-6 relative;
  3. a staged fit (--free '^width_') leaves every other parameter at its start, bit for bit, and lowers the
     objective;
  4. the fitted project written by --save-fit reopens with the same objective (1e-8 relative);
  5. engine_fit.py: lm, migrad-gradient and migrad converge, and an unknown method exits 2;
  6. convert_sav_names.py renames an old-format file, checks it (exit 0), refuses a mismatched model (exit 2) and a
     current-format file (exit 1).
Prints PASS / FAIL lines; exit 0 all passed, 1 a failure, 77 the engine cannot be imported.  About 15 s at one
thread.  Run: python3 selftest.py   (PYTHONPATH: the checkout; PYAZR_ENGINE_DIR where the build put the engine)
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
failed = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail and not ok else ""), flush=True)
    if not ok:
        failed.append(name)


def run(script, *args, cwd):
    p = subprocess.run([sys.executable, os.path.join(HERE, script), *args], cwd=cwd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def main():
    try:
        import numpy as np
        import pyazr
        pyazr.kinds()
    except Exception as e:
        print(f"skip: the engine cannot be imported ({type(e).__name__}: {e})")
        return 77
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix="azure2-fit-selftest-") as tmp:
        shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T3"), os.path.join(tmp, "T3"),
                        ignore=shutil.ignore_patterns("output", "checks"))
        os.makedirs(os.path.join(tmp, "T3", "output"))
        start = os.path.join("T3", "T3_start.azr")
        with pyazr.Session(os.path.join(tmp, start)) as s:
            lm = s.fit_with_engine("lm").objective
            mg = s.fit_with_engine("migrad-gradient").objective

        code, out = run("lsq_fit.py", start, "--out", "one.json", "--save-fit", os.path.join("T3", "fit.azr"), cwd=tmp)
        one = json.load(open(os.path.join(tmp, "one.json"))) if code == 0 else {"objective": float("nan")}
        check(f"lsq_fit from the start = the engine's LM minimum: {one['objective']:.8f} vs {lm:.8f}",
              code == 0 and abs(one["objective"] - lm) <= 1e-6 * lm, out[-300:] if code else "")

        code, out = run("lsq_fit.py", start, "--starts", "8", "--jitter", "3", "--seed", "1", "--out", "ms.json",
                        cwd=tmp)
        ms = json.load(open(os.path.join(tmp, "ms.json"))) if code == 0 else {"objective": float("nan"), "starts": []}
        ends = sorted({round(r["objective"], 2) for r in ms["starts"] if r["objective"]})
        check(f"lsq_fit with 8 starts reaches MIGRAD-gradient's minimum: {ms['objective']:.8f} vs {mg:.8f} "
              f"(the starts ended at {ends[:4]})", code == 0 and ms["objective"] <= mg * (1 + 1e-6))

        code, out = run("lsq_fit.py", start, "--free", "^width_", "--out", "stage.json", cwd=tmp)
        st = json.load(open(os.path.join(tmp, "stage.json")))
        held = [i for i, n in enumerate(st["names"]) if not n.startswith("width_")]
        same = all(st["x"][i] == st["x0"][i] for i in held)
        check(f"staged fit of the widths: {len(held)} other parameters unchanged, objective "
              f"{st['objective_start']:.1f} -> {st['objective']:.1f}", same and st["objective"] < st["objective_start"])

        with pyazr.Session(os.path.join(tmp, "T3", "fit.azr")) as s:
            f = s.objective()
            back = f(f.x0)
        check(f"--save-fit reopens at the fit's objective ({back:.10g})", abs(back - one["objective"]) <=
              1e-8 * one["objective"])

        code, out = run("engine_fit.py", start, "--method", "lm", "--method", "migrad-gradient", "--method",
                        "migrad", "--out", "eng.json", cwd=tmp)
        eng = json.load(open(os.path.join(tmp, "eng.json"))) if code in (0, 1) else {"fits": []}
        check("engine_fit: lm, migrad-gradient, migrad converge: "
              + ", ".join(f"{r['method']} {r['objective']:.4f}" for r in eng["fits"]),
              code == 0 and len(eng["fits"]) == 3 and abs(eng["fits"][1]["objective"] - mg) <= 1e-9 * mg)
        code, out = run("engine_fit.py", start, "--method", "simplex", cwd=tmp)
        check("engine_fit: an unknown method exits 2", code == 2, out[-200:])

        old = os.path.join(tmp, "old.sav")
        with open(old, "w") as fh:
            fh.write("  j=1_la=1_energy 3.03\n  j=1_la=1_ch=1_rwa 0.5\n  j=2_la=1_energy 17.64\n"
                     "  j=2_la=1_ch=2_rwa -0.1\n  segment_1_norm 1.02\n")
        with open(os.path.join(tmp, "model.par"), "w") as fh:
            fh.write("energy_1 3.0\nwidth_1_1 0.4\nenergy_2 17.6\nwidth_2_2 -0.2\n")
        with open(os.path.join(tmp, "other.par"), "w") as fh:
            fh.write("energy_1 3.0\nwidth_1_1 0.4\nenergy_2 17.6\nwidth_2_1 -0.2\n")
        code, out = run("convert_sav_names.py", "old.sav", "new.sav", "--check", "model.par", cwd=tmp)
        names = [l.split()[0] for l in open(os.path.join(tmp, "new.sav")) if l.strip()] if code == 0 else []
        check(f"convert_sav_names renames and checks: {names}", code == 0 and names == [
            "energy_1", "width_1_1", "energy_2", "width_2_2", "segment_1_norm"], out[-200:])
        code, _ = run("convert_sav_names.py", "old.sav", "new2.sav", "--check", "other.par", cwd=tmp)
        check("convert_sav_names refuses a model with other channels (exit 2)", code == 2)
        code, _ = run("convert_sav_names.py", "new.sav", "new3.sav", cwd=tmp)
        check("convert_sav_names: a current-format file exits 1", code == 1)
    print(f"{'FAIL' if failed else 'PASS'}: azure2-fit self-test, {time.time() - t0:.1f} s")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
