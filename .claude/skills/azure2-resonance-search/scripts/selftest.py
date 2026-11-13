#!/usr/bin/env python3
"""Self-test of the azure2-resonance-search scripts on tests/upgrade/T3 (8Be, synthetic) in a temporary directory:
dump_residuals.py, rank_map.py, peakfind.py, eval_peaks.py, blur_map.py, make_broad_only.py, blur_compare.py,
blur_figure.py, level_test.py.

The truth of T3 has a 2- level at E_x = 19.0006 MeV (p + 7Li and n + 7Be, l = 0).  The best fit of T3_start.azr
(azure2-fit's lsq_fit.py, objective 509.348) is the base; the same model with that level removed and refitted
(objective about 89939) is the model that misses it.
  1. dump_residuals.py: sum(res^2) equals the data chi-squared (1e-9).
  2. The locators find the missing level: rank_map.py's first region and blur_map.py's first region on the
     reduced model peak within 200 keV of 19.0006 MeV; peakfind.py's first candidate on the data alone lies within
     20 keV of it.
  3. level_test.py scan adds a 2- level at 18.5, 19.0 and 19.5 MeV to the reduced model and refits each: the best
     is the one started at 19.0, it ends within 5 keV of 19.0006 MeV, the F test's p < 1e-10 and BIC falls.
  4. level_test.py removing the level from the base with n7be_total held out: the held-out chi-squared grows
     (more than tenfold) and the removal's Delta is negative.
  5. make_broad_only.py zeroes the two levels narrower than 200 keV; blur_compare.py and blur_figure.py run on the
     full and broad tables; eval_peaks.py runs.
Prints PASS / FAIL lines; exit 0 all passed, 1 a failure, 77 the engine or the program is missing.  About 25 s at
one thread.  Run: python3 selftest.py   (PYTHONPATH: the checkout; PYAZR_ENGINE_DIR; AZURE2_BIN: the program)
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
FIT = os.path.join(ROOT, ".claude", "skills", "azure2-fit", "scripts", "lsq_fit.py")
TRUE_E = 19.000623
failed = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail and not ok else ""), flush=True)
    if not ok:
        failed.append(name)


def run(script, *args, cwd):
    path = script if os.path.isabs(script) else os.path.join(HERE, script)
    p = subprocess.run([sys.executable, path, *args], cwd=cwd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def first_peak(path, column):
    """The value of `column` (0-based, whitespace split) of the first ranked row of a locator's text output."""
    for line in open(path):
        t = line.split()
        if t and t[0] == "1":
            return float(t[column].rstrip("k"))
    return float("nan")


def main():
    try:
        import numpy as np
        import pyazr
        pyazr.kinds()
    except Exception as e:
        print(f"skip: the engine cannot be imported ({type(e).__name__}: {e})")
        return 77
    if not os.environ.get("AZURE2_BIN"):
        print("skip: AZURE2_BIN (the program, for parameters.out) is not set")
        return 77
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix="azure2-resonance-search-selftest-") as tmp:
        t3 = os.path.join(tmp, "T3")
        shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T3"), t3,
                        ignore=shutil.ignore_patterns("output", "checks"))
        for d in ("output", "checks"):
            os.makedirs(os.path.join(t3, d))
        code, out = run(FIT, "T3_start.azr", "--starts", "8", "--jitter", "3", "--seed", "1", "--save-fit", "best.azr",
                        cwd=t3)
        if code != 0:
            check("the best fit of T3 (azure2-fit's lsq_fit.py)", False, out[-300:])
            return 1
        p = pyazr.Project.open(os.path.join(t3, "best.azr"))
        label = next(lv.label for lv in p.levels() if lv.jpi == "2-")
        p.remove_level(label)
        p.save(os.path.join(t3, "reduced.azr"))
        run(FIT, "reduced.azr", "--save-fit", "reduced_fit.azr", cwd=t3)
        subprocess.run([os.environ["AZURE2_BIN"], "--no-readline", "best.azr"], input="1\n\n\n7\n", text=True,
                       cwd=t3, capture_output=True)
        par = os.path.join("output", "parameters.out")

        code, out = run("dump_residuals.py", "best.azr", "full.npz", cwd=t3)
        with pyazr.Session(os.path.join(t3, "best.azr")) as s:
            chi2 = s.evaluate(s.objective().x0).chi2
        res = np.load(os.path.join(t3, "full.npz"))["res"] if code == 0 else np.array([np.nan])
        check(f"dump_residuals: sum(res^2) {np.sum(res ** 2):.8f} = data chi2 {chi2:.8f}",
              abs(np.sum(res ** 2) - chi2) <= 1e-9 * chi2, out[-300:])

        run("dump_residuals.py", "reduced_fit.azr", "reduced.npz", cwd=t3)
        code, out = run("rank_map.py", "reduced.npz", "regions.txt", cwd=t3)
        e = first_peak(os.path.join(t3, "regions.txt"), 2) if code == 0 else float("nan")
        check(f"rank_map: first region of the reduced model peaks at {e:.3f} MeV", abs(e - TRUE_E) < 0.2, out[-300:])
        code, out = run("blur_map.py", "reduced.npz", "blur.txt", cwd=t3)
        e = first_peak(os.path.join(t3, "blur.txt"), 2) if code == 0 else float("nan")
        check(f"blur_map: first region of the reduced model peaks at {e:.3f} MeV", abs(e - TRUE_E) < 0.2, out[-300:])
        code, out = run("peakfind.py", "full.npz", "peaks.txt", "--levels", par, cwd=t3)
        e = first_peak(os.path.join(t3, "peaks.txt"), 1) if code == 0 else float("nan")
        check(f"peakfind: first candidate from the data alone at {e:.4f} MeV", abs(e - TRUE_E) < 0.02, out[-300:])
        code, out = run("eval_peaks.py", "peaks.txt", par, cwd=t3)
        check("eval_peaks runs on the candidates", code == 0 and "sets >= 1" in out, out[-300:])

        code, out = run("level_test.py", "reduced_fit.azr", "--variant",
                        "scan:2-:18.5,19.0,19.5:2/0/2/100000,3/0/2/50000", "--workdir", "scan", "--out", "scan.json",
                        cwd=t3)
        rows = json.load(open(os.path.join(t3, "scan.json"))) if code == 0 else []
        best = min((r for r in rows[1:] if "objective" in r), key=lambda r: r["objective"], default=None)
        ok = best is not None and best["name"] == "add 2- 19"
        if ok:
            found = [v for n, v in zip(best["names"], best["x"]) if n.startswith("energy_") and abs(v - 19.0) < 0.1]
            ok = len(found) == 1 and abs(found[0] - TRUE_E) < 0.005 and best["p_value"] < 1e-10 and \
                best["bic"] < rows[0]["objective"] + rows[0]["k"] * np.log(rows[0]["N"])
            detail = f"{best['name']}: E = {found[0] if found else None}, Delta {best['delta_objective']:.0f}, " \
                     f"p {best['p_value']:.1e}"
        else:
            detail = out[-300:]
        check(f"level_test scan finds the removed level ({detail})", ok)

        code, out = run("level_test.py", "best.azr", "--variant", f"remove:{label}", "--hold-out", "n7be_total",
                        "--workdir", "remove", "--out", "remove.json", cwd=t3)
        rows = json.load(open(os.path.join(t3, "remove.json"))) if code == 0 else [{}, {}]
        ok = code == 0 and rows[1].get("held_out_chi2", 0) > 10 * rows[0]["held_out_chi2"] and \
            rows[1]["delta_objective"] < 0
        check("level_test remove with n7be_total held out: held-out chi2 "
              f"{rows[0].get('held_out_chi2', 0):.1f} -> {rows[1].get('held_out_chi2', 0):.1f}", ok, out[-300:])

        code, out = run("make_broad_only.py", "best.azr", par, "broad.azr", "--gmax", "200", cwd=t3)
        check("make_broad_only zeroes the two levels narrower than 200 keV", code == 0 and out.startswith("2 levels"),
              out[-300:])
        run("dump_residuals.py", "broad.azr", "broad.npz", cwd=t3)
        code, out = run("blur_compare.py", "full.npz", "broad.npz", "compare.txt", cwd=t3)
        check("blur_compare writes its table", code == 0 and os.path.getsize(os.path.join(t3, "compare.txt")) > 0,
              out[-300:])
        code, out = run("blur_figure.py", "full.npz", "broad.npz", par, "blur.png", "0.1", "--panel",
                        "p7li_pa:-:18.0:20.0:7Li(p,a)", cwd=t3)
        check("blur_figure draws its figure", code == 0 and os.path.exists(os.path.join(t3, "blur.png")), out[-300:])
    print(f"{'FAIL' if failed else 'PASS'}: azure2-resonance-search self-test, {time.time() - t0:.1f} s")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
