#!/usr/bin/env python3
"""Campaign.init --spec on the 13N model (queue; ~15 min): spec copied, sha recorded, policy/flags derived, DECISIONS.md
usable, a level add with adopt_level: false ends accept-pending (if it would have been accepted) with a queued decision.
Run:  python3 rmfit/tests/test_campaign_spec_13n.py"""
import json, os, shutil, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
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
    tmp = tempfile.mkdtemp(prefix="rmfit_cspec_")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))
    spec = {"reaction": "12C+p", "campaign": "test", "owners": ["test"], "parametrization": "brune",
            "engine": {"flags": {"use_long_wavelength": True, "ignore_externals": False}},
            "nshards": 3, "bounds": {"bg_energy": 15.0, "e_window": 0.5},
            "search": {"top_k": 3, "pairs_per_group": 2, "stage1_nfev": 30, "stage2_nfev": 50, "union_nfev": 50},
            "data_sets": [{"file": "meyer_84", "role": "free", "norm": {"nominal": 1.0, "fixed": False, "prior_sigma": 0.05}}],
            "pre_authorized": {"adopt_level": False, "flip_adopt": True},
            "standing_decisions": [{"date": "2026-11-07", "who": "test", "text": "no level is adopted without the evaluator"}]}
    sp = os.path.join(tmp, "myspec.json")
    json.dump(spec, open(sp, "w"), indent=1)
    t0 = time.time()
    c = Campaign.init(tmp, base, flags=None, nshards=24, note="13N spec test", spec=sp)
    check("spec.json copied into the campaign", os.path.exists(os.path.join(tmp, "spec.json")))
    check("spec sha recorded", c.cfg.get("spec_sha") == c.spec.sha, str(c.cfg.get("spec_sha"))[:12])
    check("policy derived from the spec", c.cfg["policy"].get("bg_energy") == 15.0 and c.policy.e_window == 0.5, str(c.cfg["policy"]))
    check("flags derived from the spec", c.flags.get("use_brune") is True and c.flags.get("ignore_externals") is False, str(c.flags))
    check("nshards derived from the spec", c.nshards == 3, str(c.nshards))
    check("engine pin recorded", bool(c.cfg.get("engine", {}).get("so_md5")))
    cid, x, obj, keys = c.incumbent()
    check("baseline verified", np.isfinite(obj))
    did = c.decisions.add("test-question", "is the queue wired?", options=("yes", "no"))
    check("decisions queue works in the campaign", c.decisions.get(did)["status"] == "pending")
    c.decisions.answer(did, "yes", by="test")
    c.run_polish(max_nfev=60, tighten=True)
    cand = LevelCandidate(energy=3.0, channels=[dict(pair=1, L=0, S=0.5, gamma=20000.0)], source="test")
    res = c.run_level_add(cand, jpis=["1/2+", "3/2-"], pinned_nfev=25, global_nfev=40)
    out = res.get("outcome")
    print("   level add outcome:", out, {k: res[k] for k in ("jpi", "delta") if k in res})
    check("level add ran", out not in (None, "scan_failed"), str(out))
    if out == "accept-pending":
        pend = c.decisions.pending()
        check("accept-pending queued a decision", any("adopt-level" in d["slug"] for d in pend), str([d.id for d in pend]))
        check("not adopted: base model unchanged", os.path.basename(c.cfg["base_azr"]) == os.path.basename(base), c.cfg["base_azr"])
    else:
        check("no 'accept' adoption slipped through the gate", out != "accept", out)
    c.close()
    print(f"   total {time.time()-t0:.0f}s")
    shutil.rmtree(tmp, ignore_errors=True)
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all campaign-spec checks passed")

if __name__ == "__main__":
    main()
