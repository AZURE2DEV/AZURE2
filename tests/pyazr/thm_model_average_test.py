#!/usr/bin/env python3
"""End-to-end smoke test of scripts/thm_model_average.py.

On a copy of tests/18O_p_a_thm with its two THM segments made one experiment
and the four widths freed:

  1. --dry-run (no numpy, no engine): the variant grid of 2 radii x pw/dw is
     listed with its priors, and the dw variants are reported as refused
     (no distortion=), as AzrModel refuses them.  The line-shape and R(E)
     axes: labels, R(E) collapsed for dw, and the refusals of a project
     without the kinematics keys.
  2. The run (needs numpy, scipy and the engine; skips cleanly otherwise),
     a few least-squares steps per variant, each in its own subprocess: two
     variants fitted, two skipped with the reason, the chi2 of each fitted
     variant is that of its saved snapshot, the weights are Akaike weights
     of the table, the averaged values are the weighted means of the
     variants, and the averaged .azr loads and carries them.
  3. lineshape=on + distortion=coulomb with a --penalty-hook row of 3 and
     --x-scale 1, no fit steps: the chi2 is that of the same project edited
     by hand plus 9, in a subprocess and with --in-process alike.

Run from anywhere:  python3 tests/pyazr/thm_model_average_test.py
"""
import csv
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "1")   # before numpy (OpenBLAS buffers)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SOURCE = os.path.join(ROOT, "tests", "18O_p_a_thm")
AZR = "18O_p_a_thm.azr"
SCRIPT = os.path.join(ROOT, "scripts", "thm_model_average.py")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def load(name):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(ROOT, "pyazr", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


AzrModel = load("azrfile").AzrModel
ma = load("modelavg")

tmp = tempfile.mkdtemp(prefix="thm_modelavg_")
try:
    proj = os.path.join(tmp, "o18")
    shutil.copytree(SOURCE, proj)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(proj, junk), ignore_errors=True)
        os.makedirs(os.path.join(proj, junk))
    model = AzrModel.from_file(os.path.join(proj, AZR))
    model.set_thm_experiment("A", [1, 2])
    for lv in model.levels:
        for c in lv.channels:
            c.channel_fixed = False
    azr = model.write(os.path.join(proj, AZR))
    common = [sys.executable, SCRIPT, azr, "--radius-pairs", "1",
              "--radii", "4.6", "5.1", "--vertex-model", "pw", "dw",
              "--prior", "radius=4.6:0.5"]
    env = dict(os.environ)

    print("1. --dry-run")
    r = subprocess.run(common + ["--dry-run"], capture_output=True, text=True,
                       env=env, timeout=120)
    lines = [l for l in r.stdout.splitlines() if l.startswith("  ")]
    check("exit 0", r.returncode == 0, r.stderr[-500:])
    check("4 variants listed", r.stdout.startswith("4 variant(s)") and len(lines) == 4,
          r.stdout)
    check("labels", [l.split()[0] for l in lines]
          == ["r4.6_pw", "r4.6_dw", "r5.1_pw", "r5.1_dw"], lines)
    check("priors", ["prior 0.5" in lines[0], "prior 1" in lines[2]] == [True, True])
    check("dw refused without distortion=",
          all("refused" in l and "distortion" in l for l in lines if "_dw" in l)
          and not any("refused" in l for l in lines if "_pw" in l))
    check("dry run writes nothing", not os.path.exists(os.path.join(tmp, "out")))

    kin = AzrModel.from_file(azr)
    kin.set_thm_experiment("A", [1, 2], beam="18O", target="d", spectator="n",
                           Ebeam=54)
    kazr = kin.write(os.path.join(proj, "kin.azr"))
    axes = ["--vertex-model", "pw", "dw", "--optical", "coulomb",
            "--lineshape", "on", "off", "--distortion", "none", "coulomb"]
    labels = ["pw_ls-on_R-none", "pw_ls-on_R-coulomb", "pw_ls-off_R-none",
              "pw_ls-off_R-coulomb", "dw_coulomb_ls-on", "dw_coulomb_ls-off"]
    r = subprocess.run([sys.executable, SCRIPT, kazr] + axes + ["--dry-run"],
                       capture_output=True, text=True, env=env, timeout=120)
    lines = [l for l in r.stdout.splitlines() if l.startswith("  ")]
    check("N_C x R(E) x vertex: 6 variants, R(E) collapsed for dw",
          r.returncode == 0 and [l.split()[0] for l in lines] == labels, r.stdout + r.stderr)
    check("with the kinematics keys none is refused",
          not any("refused" in l for l in lines), lines)
    r = subprocess.run([sys.executable, SCRIPT, azr] + axes + ["--dry-run"],
                       capture_output=True, text=True, env=env, timeout=120)
    ref = {l.split()[0]: "refused" in l for l in r.stdout.splitlines() if l.startswith("  ")}
    check("without them N_C and R(E) are refused, the plain pw variant is not",
          ref.get("pw_ls-off_R-none") is False and ref.get("pw_ls-on_R-none")
          and ref.get("pw_ls-off_R-coulomb"), r.stdout)

    # The script's helpers, without a subprocess (pure Python).
    sspec = importlib.util.spec_from_file_location("thm_model_average", SCRIPT)
    tma = importlib.util.module_from_spec(sspec)
    sspec.loader.exec_module(tma)
    # --optical on a project whose experiment is already dw: one variant per
    # potential pair, not one collapsed (pw) variant
    g = tma.variant_grid({"optical": ["ancai06/kd03", "daehnick80/kd03"]}, project_dw=True)
    check("project already dw: the optical axis is kept",
          [l for l, _ in g] == ["ancai06+kd03", "daehnick80+kd03"], g)
    g = tma.variant_grid({"optical": ["ancai06/kd03", "daehnick80/kd03"]})
    check("project pw: the optical axis collapses", len(g) == 1, g)
    # a ps=table: path does not make the label (a file name) a path
    g = tma.variant_grid({"ps": ["table:tabs/ps.dat", "hulthen:0-40"]})
    check("labels are file names", all(os.sep not in l and "/" not in l for l, _ in g)
          and len({l for l, _ in g}) == 2, g)
    # prior rules as typed: CLI axis spelling, radius as a number; one that
    # matches nothing is refused instead of leaving every prior at 1
    g = tma.variant_grid({"radius_pairs": [1], "radii": [4.6, 6.1], "vertex_model": ["pw", "dw"]})
    pr = tma.priors_for(g, {"radius=6.10": 0.5, "vertex-model=dw": 0.1})
    check("prior rules: radius=6.10 and vertex-model=dw apply",
          pr == {"r4.6_pw": 1.0, "r4.6_dw": 0.1, "r6.1_pw": 0.5, "r6.1_dw": 0.05}, pr)
    try:
        tma.priors_for(g, {"radius=7": 0.5})
        check("a prior rule that matches nothing is refused", False)
    except SystemExit:
        check("a prior rule that matches nothing is refused", True)
    # --out inside the project directory: no copy of the copy
    inner = os.path.join(proj, "modelavg", "work")
    tma.copy_project(proj, inner, AZR)
    check("copy_project with the destination inside the source",
          os.path.isfile(os.path.join(inner, AZR))
          and not os.path.exists(os.path.join(inner, "modelavg")), os.listdir(inner))
    shutil.rmtree(os.path.join(proj, "modelavg"))

    print("\n2. the run")
    try:
        import numpy                                             # noqa: F401
        import scipy.optimize                                    # noqa: F401
        sys.path.insert(0, ROOT)
        from pyazr import azure2
    except Exception as err:
        print(f"skip the run: engine or scipy not available ({type(err).__name__}: {err})")
        sys.exit(1 if failures else 0)
    out = os.path.join(tmp, "out")
    r = subprocess.run(common + ["--out", out, "--max-nfev", "4"],
                       capture_output=True, text=True, env=env, timeout=900)
    check("exit 0", r.returncode == 0, (r.stdout + r.stderr)[-1500:])
    rows = list(csv.DictReader(open(os.path.join(out, "variants.csv"))))
    check("one row per variant", [x["label"] for x in rows]
          == ["r4.6_pw", "r4.6_dw", "r5.1_pw", "r5.1_dw"], rows)
    st = {x["label"]: x for x in rows}
    check("pw fitted, dw skipped with the reason",
          st["r4.6_pw"]["status"] == st["r5.1_pw"]["status"] == "fitted"
          and st["r4.6_dw"]["status"] == st["r5.1_dw"]["status"] == "skipped"
          and "distortion" in st["r4.6_dw"]["reason"], rows)
    variants = json.load(open(os.path.join(out, "variants.json")))
    avg = json.load(open(os.path.join(out, "average.json")))
    check("k = 4 widths + the profiled norm",
          all(v["nfree"] == 5 and v["npoints"] == 118 for v in variants))
    for v in variants:
        fit = os.path.join(out, "work", v["label"] + "_fit.azr")
        import numpy as np
        with azure2(fit, cwd=os.path.join(out, "work")) as s:
            chi = float(np.sum(s.calculate_chi2_rwa(np.asarray(s.params_rwa, float))))
        check(f"{v['label']}: chi2 is that of its snapshot", close(chi, v["chi2"], 1e-5),
              f"{chi} vs {v['chi2']}")
        check(f"{v['label']}: chi2 went down", v["chi2"] <= v["fit"]["start_chi2"])
    aic = [v["chi2"] + 2 * v["nfree"] for v in variants]
    raw = [v["prior"] * math.exp(-0.5 * (a - min(aic))) for v, a in zip(variants, aic)]
    want = [x / sum(raw) for x in raw]
    check("Akaike weights with the priors",
          all(close(w["weight"], x) for w, x in zip(avg["weights"], want)),
          (avg["weights"], want))
    ok = True
    for name, a in avg["parameters"].items():
        m = sum(w * v["values"][name] for w, v in zip(want, variants))
        ok &= close(a["mean"], m, 1e-9)
    check("averaged values = weighted means of the variants", ok)
    back = ma.Variant.from_dict(variants[0])
    check("variants.json carries the covariance", back.covariance is not None
          and len(back.covariance[0]) == 4)
    aazr = os.path.join(out, "averaged", "18O_p_a_thm_avg.azr")
    check("averaged .azr written beside its data", os.path.exists(aazr)
          and os.path.isdir(os.path.join(out, "averaged", "data")))
    slots = {}
    m2 = AzrModel.from_file(aazr)
    for (_, n), lv in m2.engine_level_keys().items():
        for c in lv.channels:
            slots[ma.width_label(lv.J, lv.parity, n, c.pair, c.L, c.S)] = c.gamma
    check("it carries the means", all(
        close(slots[k], a["mean"], 1e-6) for k, a in avg["parameters"].items()),
        {k: (slots.get(k), a["mean"]) for k, a in avg["parameters"].items()})
    with azure2(aazr, cwd=os.path.join(out, "averaged")) as s:
        check("it opens in the engine", len(s.params_rwa) == 4)
    summ = open(os.path.join(out, "summary.txt")).read()
    check("summary lists the skipped variants", "r5.1_dw" in summ and "weight" in summ)
    check("the subprocesses' fit lines are in run.log",
          open(os.path.join(out, "run.log")).read().count("end chi2") == 2)

    print("\n3. line shape, R(E), penalty hook, in a subprocess and in-process")
    hook = os.path.join(tmp, "hook.py")
    with open(hook, "w") as fh:
        fh.write("def penalty(session, x):\n    return [3.0]\n")
    hand = AzrModel.from_file(kazr)
    hand.set_thm_experiment("A", [1, 2], beam="18O", target="d", spectator="n",
                            Ebeam=54, lineshape=True, distortion="coulomb")
    hazr = hand.write(os.path.join(proj, "hand.azr"))
    with azure2(hazr, cwd=proj) as s:
        want = float(np.sum(s.calculate_chi2_rwa(np.asarray(s.params_rwa, float)))) + 9.0
    with azure2(kazr, cwd=proj) as s:
        plain = float(np.sum(s.calculate_chi2_rwa(np.asarray(s.params_rwa, float))))
    got = []
    for extra in ([], ["--in-process"]):
        o3 = os.path.join(tmp, "out3" + "".join(extra))
        r = subprocess.run([sys.executable, SCRIPT, kazr, "--out", o3, "--lineshape", "on",
                            "--distortion", "coulomb", "--max-nfev", "0", "--x-scale", "1",
                            "--penalty-hook", hook] + extra,
                           capture_output=True, text=True, env=env, timeout=600)
        v = json.load(open(os.path.join(o3, "variants.json")))
        got.append(v[0]["chi2"] if r.returncode == 0 and len(v) == 1 else float("nan"))
        check(f"{' '.join(extra) or 'subprocess'}: chi2 = hand-edited project + 9",
              close(got[-1], want, 1e-7) and v[0]["fit"]["npenalty"] == 1,
              f"{got[-1]} vs {want}; {(r.stdout + r.stderr)[-800:]}")
    check("the axes changed the model", not close(want - 9.0, plain, 1e-4), (want, plain))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("all THM model-average driver checks passed")
