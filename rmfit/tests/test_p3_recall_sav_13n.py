#!/usr/bin/env python3
"""RECALL children start from the incumbent's norms and shifts, not their nominal values (2026-11-20).

bake() writes the levels into the exported .azr and every free norm / shift only into its companion .sav, so a child
campaign initialised from the baked .azr alone starts with norms at nominal (one data set's chi2
rose 37-fold).  run_recall now passes the incumbent's .sav.  Here: 13N with meyer_84's norm free, the incumbent's norm
moved off nominal by a short polish, a one-level recall round; the child's baseline (candidate 1) must carry the parent's norm.
Run:  python3 rmfit/tests/test_p3_recall_sav_13n.py (pyazr 3; the pyazr 2 version is test_recall_sav_13n.py)"""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if not _PYAZR3:
    print("SKIP: a pyazr 3 test; this process runs pyazr 2 (the pyazr 2 version is test_recall_sav_13n.py)")
    sys.exit(77)
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.campaign import Campaign
from rmfit.azr3 import Model as AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_recall_sav_")
    try:
        shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
        os.makedirs(os.path.join(tmp, "output"))
        mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
        mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
        mdl.set_output_dir(os.path.join(tmp, "output"))
        base = mdl.write(os.path.join(tmp, "13N.azr"))
        c = Campaign.init(tmp, base, flags={}, nshards=2, policy=dict(bg_energy=15.0, e_window=0.5), note="recall .sav test")
        c.run_polish(max_nfev=40)                       # the fitted norm moves off its nominal 1
        cid, x, obj, keys = c.incumbent()
        j = next(i for i, k in enumerate(keys) if k.kind == "norm")
        v = float(x[j])
        check("incumbent's fitted norm is off nominal (test setup)", abs(v - 1.0) > 1e-4, f"{v}")
        res = c.run_recall(n=1, min_width_kev=0.0, global_nfev=3, pinned_nfev=3, local_nfev=0)
        check("recall ran one child", res["n"] == 1, str(res))
        r0 = res["results"][0]
        child = Campaign(os.path.join(tmp, r0["dir"]), enforce_engine=False)
        xc, kc = child.ledger.load(1)
        jc = next((i for i, k in enumerate(kc) if str(k) == str(keys[j])), None)
        check("child has the norm as a free parameter", jc is not None, str([str(k) for k in kc]))
        if jc is not None:
            check("child baseline starts at the incumbent's norm, not nominal", abs(float(xc[jc]) - v) < 1e-6 * max(1.0, abs(v)),
                  f"child {xc[jc]} vs incumbent {v}")
        child.close()
        c.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all recall .sav checks passed")


if __name__ == "__main__":
    main()
