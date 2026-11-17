#!/usr/bin/env python3
"""init --sav and the file's FIXED parameters (2026-11-10), with the real engine on 13N.
AZURE2's parameter-file reader sets every parameter it names, fixed level energies and widths included; rmfit's seed took
only the free values and the sharded evaluator took the fixed ones from the .azr, so an .azr + .sav pair whose fixed values
disagree was evaluated as a model nobody fitted (a benchmark initialised at more than twice the archived objective; two dozen fixed
background widths differed by up to 10x) -- while a script-built base means the .azr and uses the .sav as a seed.
Campaign.apply_sav_fixed: default `ignore` keeps the .azr's fixed values and logs the disagreement with the objective both ways;
`--sav-fixed apply` writes <stem>_sav.azr (+ .sav) with the file applied and init uses it as the base.
Model: 13N with the E1 photon width of the 3/2- 3.5032 MeV level FIXED; a .sav written from a fit vector, then its fixed
photon width edited to 3x the value (the disagreeing pair).
  - the pair with equal fixed values is left alone (base unchanged, no _sav.azr);
  - the disagreeing pair, default: baseline = the .azr's fixed values, WARNING with both objectives in campaign.log and readme;
  - the disagreeing pair, apply: init's baseline objective equals pyazr's load_sav objective of (.azr + edited .sav), not the .azr's;
  - the campaign's base_azr is the _sav.azr and a readme note records it.
Run: python3 rmfit/tests/test_init_sav_fixed_13n.py   (about a minute)"""
import os, sys, shutil, tempfile, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3 (its pyazr 3 version is test_p3_init_sav_fixed_13n.py)")
    sys.exit(77)
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.campaign import Campaign
from pyazr import AzrModel, azure2
failures = []
FLAGS = dict(use_brune=True)


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_initsav_")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    lv = next(l for l in mdl.active_levels if l.jpi == "3/2-" and abs(float(l.energy) - 3.5032) < 1e-3)
    ph = next(c for c in lv.channels if c.is_photon)
    ph.channel_fixed = True
    mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))
    # a .sav of the model as it is (fixed values equal) and one with the fixed photon width changed
    with azure2(base, data_mode=True, **FLAGS) as m:
        x = np.asarray(m.params_rwa, float).copy()
        azr_same, sav_same = m.save_fit(os.path.join(tmp, "same.azr"), x, param_sav=True, verify=True)
        obj_azr = m.objective(x)
    lines = open(sav_same).read().splitlines()
    with azure2(base, data_mode=True, **FLAGS) as m:
        name_fixed = next(p.name for p in m.parameters if p.fixed and p.name.startswith("width_"))
    out = []
    for ln in lines:
        t = ln.split()
        if t and t[0] == name_fixed:
            v = float(t[1]) * 3.0 if float(t[1]) != 0 else 1.0e-2
            ln = f"{t[0]:>28s} {v: .7e} {0.0: .7e}"
        out.append(ln)
    sav_diff = os.path.join(tmp, "diff.sav")
    open(sav_diff, "w").write("\n".join(out) + "\n")
    with azure2(base, data_mode=True, **FLAGS) as m:
        xd = m.load_sav(sav_diff)
        obj_sav = m.objective(xd)
    print(f"   fixed parameter edited: {name_fixed}; objective .azr {obj_azr:.4f}, .azr+.sav {obj_sav:.4f}")
    check("the edited pair evaluates to another objective", abs(obj_sav - obj_azr) > 1e-5 * max(1.0, abs(obj_azr)), f"{obj_sav} vs {obj_azr}")
    # 1. equal pair: untouched
    d1 = os.path.join(tmp, "c_same"); os.makedirs(d1)
    c1 = Campaign.init(d1, base, sav=sav_same, flags=FLAGS, nshards=2, policy=dict(bg_energy=15.0, e_window=0.5), note="same")
    check("equal fixed values: base unchanged", os.path.normpath(c1.base_azr) == os.path.normpath(base), c1.base_azr)
    check("equal fixed values: no _sav.azr written", not os.path.exists(os.path.splitext(base)[0] + "_sav.azr"))
    _, _, o1, _ = c1.incumbent(); c1.close()
    check("equal fixed values: baseline = the .azr objective", abs(o1 - obj_azr) < 1e-6 * max(1.0, abs(obj_azr)), f"{o1} vs {obj_azr}")
    # 2a. disagreeing pair, default mode: the .azr's fixed values are used, the disagreement is logged with both objectives
    d0 = os.path.join(tmp, "c_ign"); os.makedirs(d0)
    c0 = Campaign.init(d0, base, sav=sav_diff, flags=FLAGS, nshards=2, policy=dict(bg_energy=15.0, e_window=0.5), note="ignore")
    _, _, o0, _ = c0.incumbent(); c0.close()
    check("default (ignore): base unchanged, no _sav.azr", os.path.normpath(c0.base_azr) == os.path.normpath(base)
          and not os.path.exists(os.path.splitext(base)[0] + "_sav.azr"), c0.base_azr)
    check("default (ignore): baseline = the .azr's fixed values with the .sav's free ones", abs(o0 - obj_azr) < 1e-6 * max(1.0, abs(obj_azr)), f"{o0} vs {obj_azr}")
    lg0 = open(os.path.join(d0, "campaign.log")).read()
    check("default (ignore): WARNING names the count and the objective both ways",
          "WARNING 1 FIXED parameter(s)" in lg0 and "with the .azr's fixed values" in lg0 and f"{obj_sav:.3f} with the .sav's" in lg0
          and "--sav-fixed apply" in lg0, lg0[-600:])
    check("default (ignore): readme carries the warning", "WARNING at init" in open(os.path.join(d0, "readme")).read() and "FIXED parameter" in open(os.path.join(d0, "readme")).read())
    # 2b. disagreeing pair, --sav-fixed apply: base rewritten with the .sav applied
    d2 = os.path.join(tmp, "c_diff"); os.makedirs(d2)
    c2 = Campaign.init(d2, base, sav=sav_diff, flags=FLAGS, nshards=2, policy=dict(bg_energy=15.0, e_window=0.5), note="diff", sav_fixed="apply")
    _, _, o2, _ = c2.incumbent(); c2.close()
    exp = os.path.splitext(base)[0] + "_sav.azr"
    check("disagreeing fixed values: campaign base is <stem>_sav.azr", os.path.normpath(c2.base_azr) == os.path.normpath(exp), c2.base_azr)
    check("disagreeing fixed values: _sav.azr + .sav written", os.path.exists(exp) and os.path.exists(os.path.splitext(exp)[0] + ".sav"))
    check("disagreeing fixed values: baseline = pyazr load_sav objective (.azr + .sav as the CLI reads it)",
          abs(o2 - obj_sav) < 1e-6 * max(1.0, abs(obj_sav)), f"{o2} vs {obj_sav} (the .azr alone gives {obj_azr})")
    cfg = json.load(open(os.path.join(d2, "campaign.json")))
    check("campaign.json base_azr points at the _sav.azr", cfg["base_azr"].endswith("13N_sav.azr"), cfg["base_azr"])
    rd = open(os.path.join(d2, "readme")).read()
    check("readme records the applied .sav", "the .azr with the .sav applied" in rd)
    print(f"   baselines: same {o1:.4f}  diff {o2:.4f}")
    if failures:
        print(f"FAILED {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    shutil.rmtree(tmp, ignore_errors=True)
    print("all init --sav fixed checks passed")


if __name__ == "__main__":
    main()
