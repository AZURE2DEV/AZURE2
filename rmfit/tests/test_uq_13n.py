#!/usr/bin/env python3
"""rmfit uq end to end on the 13N fixture (queue job tests/_queue_uq/run_crc_uq_test, ~10 min on 4 slots, 2 chains): export of a tiny campaign, saturation check,
Laplace, NUTS (2 x (150 + 150)), gates with test thresholds, bands on the fixture's own <segmentsTest> grids, UQ.md.
Run:  cd rmfit/tests/_queue_uq && qsub run_crc_uq_test   (or RMFIT_UQ_LOCAL=1 python3 rmfit/tests/test_uq_13n.py on a login node, < 20 min)"""
import json, os, shutil, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.campaign import Campaign
from rmfit import uq as U

FAILS = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""))
    sys.stdout.flush()
    if not ok:
        FAILS.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_uq13n_")
    src = os.path.join(AZURE2_ROOT, "tests", "13N")
    shutil.copytree(os.path.join(src, "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    from rmfit.azr3 import Model as AzrModel
    mdl = AzrModel.from_file(os.path.join(src, "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    for f in ("artemov", "seagrave", "skowronski_fels", "kettner_0."):     # the target-integrated sets: 2.7 s per gradient
        mdl.set_segment_active(f, False)                                    # with them, 0.02 s without (NUTS needs thousands)
    mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))
    t0 = time.time()
    c = Campaign.init(tmp, base, nshards=3, note="uq fixture", policy=dict(theta2_cap=10.0, bg_theta2_cap=10.0),   # interior mode (the
                      flags=dict(use_brune=True, ignore_externals=True, use_long_wavelength=True))   # fixture's 1/2+ sits at theta2 2.3); externals off
    c.run_polish(max_nfev=80, tighten=True)
    c.close()
    c = Campaign(tmp)
    cfg = dict(enabled=True, chains=2, warmup=150, draws=150, target_accept=0.9, max_tree_depth=6, ess_min=20, rhat_max=1.2,
               divergent_max_frac=0.2, bfmi_min=0.1, max_attempts=1, band_draws=60, bands_max=2,
               bands=[dict(reaction="12C(p,g)13N[total]", quantity="sigma", unit="b",
                           energies="list " + " ".join(f"{0.1 + 0.02 * k:.2f}" for k in range(96)) + " MeV lab",
                           quantity_band="sfactor", label="12C(p,g)13N")])
    # the fixture has no <segmentsTest>, so the band comes from a spec grid (as a campaign's spec does)
    st = U.run(campaign=c, spec_uq=cfg, log=print, nprocs=2, threads=1)
    c.close()
    S = st["stages"]
    print("   stages:", list(S))
    check("EXPORT wrote uq/<name>.azr + .sav", os.path.exists(S["EXPORT"]["azr"]) and os.path.exists(S["EXPORT"]["sav"]))
    check("the uq copy writes under uq/output", "uq/output" in open(S["EXPORT"]["azr"]).read().replace("\\", "/"))
    su = S["SETUP"]
    check("SETUP: Park session reproduces the objective", abs(su["park_objective"] - su["objective"]) < 1e-6 * su["objective"], str((su["park_objective"], su["objective"])))
    check("SETUP: Laplace positive definite", su["eig_min"] > 0, str(su["eig_min"]))
    check("SETUP: sampled set = free minus held", su["n_sampled"] == su["nfree"] - len(su["held"]), str((su["n_sampled"], su["nfree"], su["held"])))
    sa = S["SAMPLE"]
    check("SAMPLE: 2 chains x 150 draws", sa["chains"] == 2 and sa["draws"] == 150, str(sa))
    z = np.load(os.path.join(tmp, "uq", f"chains_{sa['tag']}.npz"))
    check("chains file holds samples (2, 150, n)", z["samples"].shape[:2] == (2, 150) and z["samples"].shape[2] == su["n_sampled"], str(z["samples"].shape))
    ga = S["GATES"]
    check("GATES computed R-hat / ESS / divergences", all(k in ga for k in ("rhat_max", "ess_min", "divergent_frac", "bfmi_min")), str(ga))
    print("   gates:", {k: ga[k] for k in ("rhat_max", "ess_min", "ess_median", "divergent_frac", "bfmi_min", "sd_ratio_median", "passed")})
    check("GATES: sd(sampled)/sd(Laplace) within 0.3-3 at the median", 0.3 < ga["sd_ratio_median"] < 3, str(ga["sd_ratio_median"]))
    ba = S.get("BANDS", {})
    check("BANDS ran on the spec grid", bool(ba.get("grids")), str(ba.get("note")))
    if ba.get("grids"):
        g = ba["grids"][0]
        print("   band:", {k: g[k] for k in ("label", "quantity", "e0", "best0", "sigma_lin0", "median0", "ratio0", "offset0", "ratio_median")})
        check("band table and figure written", os.path.exists(os.path.join(tmp, "uq", g["table"])) and os.path.exists(os.path.join(tmp, "uq", f"uq_band_{g['key']}.png")))
        check("sampled half-width comparable with the linear sigma (0.3-3)", 0.3 < g["ratio_median"] < 3, str(g["ratio_median"]))
        check("linear sigma positive and finite", np.isfinite(g["sigma_lin0"]) and g["sigma_lin0"] > 0, str(g["sigma_lin0"]))
    check("UQ.md written", os.path.exists(os.path.join(tmp, "UQ.md")))
    # resume: a second run must not redo any stage
    t1 = time.time()
    st2 = U.run(campaign=Campaign(tmp), spec_uq=cfg, log=print, nprocs=2, threads=1)
    check("second run resumes (no stage redone, < 20 s)", time.time() - t1 < 20 and st2["stages"]["SAMPLE"]["tag"] == sa["tag"], f"{time.time()-t1:.0f} s")
    print(f"   total {(time.time()-t0)/60:.1f} min; dir {tmp}")
    print(open(os.path.join(tmp, "UQ.md")).read()[:1500])
    if FAILS:
        print(f"FAILED {len(FAILS)}: {', '.join(FAILS)}")
        sys.exit(1)
    shutil.rmtree(tmp, ignore_errors=True)
    print("all uq 13N checks passed")


if __name__ == "__main__":
    main()
