#!/usr/bin/env python3
"""End-to-end campaign on the 13N model: init, tightening polish, a sign
round, a level-add round, report/export.  Run:  python3 rmfit/tests/test_campaign_13n.py
"""
import json, os, shutil, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.campaign import Campaign
from rmfit.structure import LevelCandidate
from pyazr import AzrModel

failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok: failures.append(name)

def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_camp_")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))
    t0 = time.time()
    c = Campaign.init(tmp, base, flags={}, nshards=3, policy=dict(bg_energy=15.0, e_window=0.5),
                      search=dict(top_k=3, pairs_per_group=2, stage1_nfev=30, stage2_nfev=50, union_nfev=50), note="13N test")
    cid, x, obj, keys = c.incumbent()
    check("init produced a verified baseline", np.isfinite(obj))
    fr = c.run_polish(max_nfev=80, tighten=True)
    cid2, x2, obj2, _ = c.incumbent()
    check("polish improved the incumbent", obj2 < obj, f"{obj2} vs {obj}")
    st = c.run_sign_round()
    check("sign round ran", "proposed" in st and st["proposed"] > 0, str(st))
    cand = LevelCandidate(energy=3.0, channels=[dict(pair=1, L=0, S=0.5, gamma=20000.0)], source="test")
    res = c.run_level_add(cand, jpis=["1/2+", "3/2-"], pinned_nfev=25, global_nfev=40)
    check("level-add round produced an outcome", "outcome" in res and res["outcome"] != "scan_failed", str(res.get("outcome")))
    print("   level add:", {k: res[k] for k in ("jpi", "delta", "outcome") if k in res})
    off = c.run_dataset_off("meyer_84", nfev=30)
    check("dataset-off test reports an extra recovery", "extra_recovery" in off and np.isfinite(off["extra_recovery"]), str(off))
    print("   dataset off:", {k: off[k] for k in ("own_chi2", "rest_before", "rest_after", "extra_recovery") if k in off})
    rep = c.report()
    check("report written", os.path.exists(rep))
    azr, sav, sc = c.export(path=os.path.join(tmp, "export.azr"))
    check("export verified", os.path.exists(azr) and os.path.exists(sav))
    c.close()
    print(f"   total {time.time()-t0:.0f}s; STATUS.md:\n" + open(os.path.join(c.dir, "STATUS.md")).read()[:600])
    shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all campaign checks passed")

if __name__ == "__main__":
    main()
