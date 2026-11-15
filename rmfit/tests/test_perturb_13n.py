#!/usr/bin/env python3
"""perturb -> sibling campaign -> a sign round recovers (13N)."""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.campaign import Campaign
from pyazr import AzrModel
def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_pert_")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0); mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))
    c = Campaign.init(tmp, base, flags={}, nshards=3, policy=dict(bg_energy=15.0, e_window=0.5),
                      search=dict(top_k=3, pairs_per_group=2, stage1_nfev=30, stage2_nfev=50, union_nfev=50))
    c.run_polish(max_nfev=80)
    cid, x, obj, keys = c.incumbent()
    st = c.perturb(os.path.join(tmp, "pert"), n_flips=2, seed=1)
    print("perturbed:", st["objective"], "target", st["target"], st["flips"])
    ok = st["objective"] > st["target"]
    c.close()
    c2 = Campaign(os.path.join(tmp, "pert"))
    r = c2.run_sign_round()
    print("round:", r)
    _, _, obj2, _ = c2.incumbent()
    print(f"after one round: {obj2:.3f} (target {st['target']:.3f})")
    c2.close()
    shutil.rmtree(tmp, ignore_errors=True)
    print("PERTURB OK" if ok else "PERTURB NOT WORSE")
if __name__ == "__main__":
    main()
