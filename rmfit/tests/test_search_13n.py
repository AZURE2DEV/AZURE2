#!/usr/bin/env python3
"""One sign round on the small 13N model: exercises moves.exposure,
single_flips/pair_flips, search.screen, polish_flip (sign lock), greedy union,
and the ledger.  Run:  python3 rmfit/tests/test_search_13n.py
"""
import os
import shutil
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.engine import Sharded, verify
from rmfit.keys import take
from rmfit.ledger import Ledger
from rmfit.local import BoundsPolicy, default_bounds, polish
from rmfit.moves import exposure, single_flips, pair_flips
from rmfit.search import SearchConfig, sign_round
from pyazr import AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_search_")
    work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    mdl = AzrModel.from_file(os.path.join(work, "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    base = mdl.write(os.path.join(work, "base.azr"))
    policy = BoundsPolicy(bg_energy=15.0, e_window=0.5)
    cfg = SearchConfig(top_k=4, pairs_per_group=2, stage1_nfev=30, stage2_nfev=60, union_nfev=60)

    sh = Sharded(base, os.path.join(tmp, "shards"), nshards=3, flags={})
    x = sh.x0.copy()
    b = default_bounds(sh, x, policy)
    f0 = polish(sh, x, bounds=b, policy=policy, max_nfev=120)
    x0, o0 = f0.x, f0.objective
    print(f"  polished start: {o0:.3f} ({f0.status})")

    print("1. exposure and proposals")
    ex = exposure(sh, x0)
    check("exposure has an entry per width", len(ex) == sum(k.kind == "width" for k in sh.keys))
    singles = single_flips(x0, sh.keys, exposure_map=ex)
    check("single flips proposed", len(singles) > 0, str(len(singles)))
    jpis = sorted({m.jpi for m in singles})
    pairs = pair_flips(x0, sh.keys, jpis[0], ex, top=3)
    print(f"  {len(singles)} single flips over {jpis}; {len(pairs)} pair flips in {jpis[0]}")

    print("2. a sign round")
    t0 = time.time()
    xb, ob, stat, trials = sign_round(sh, x0, o0, cfg, policy, log=print)
    print(f"  round in {time.time()-t0:.0f}s: {o0:.3f} -> {ob:.3f}; screened {stat.n_screened}, "
          f"polished {stat.n_polished}, accepted {stat.n_accepted}, reverted {stat.n_reverted}")
    check("round returns a vector no worse than the start", ob <= o0 + 1e-9)
    check("every polished trial has an outcome", all(t.outcome in
          ("accept", "worse", "marginal", "reverts", "crashed", "unphysical") for t in trials))

    print("3. ledger round trip")
    led = Ledger(os.path.join(tmp, "ledger.sqlite"))
    sid = led.add_structure("base", base, sh.keys)
    cid = led.add_candidate(sid, xb, sh.keys, move=dict(kind="test"))
    sc, keys1 = verify(base, take(xb, sh.keys, keys1 if False else sh.keys, xb), {}, workdir=os.path.join(tmp, "variants")) \
        if False else (None, None)
    led.record(cid, type("S", (), dict(objective=ob, chi2=ob, penalty=0.0, prior=0.0, npts=sh.npts,
                                        segment_chi2={}, nfev=0, wall=0.0, status=None, message=None))(), "verified")
    x2, k2 = led.load(cid)
    check("candidate round-trips", np.allclose(x2, xb) and k2 == sh.keys)
    check("best is the candidate", led.best(sid) == (cid, ob))
    led.write_status_md()
    check("STATUS.md written", os.path.exists(os.path.join(tmp, "STATUS.md")))
    sh.close()
    shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
        sys.exit(1)
    print("all search checks passed")


if __name__ == "__main__":
    main()
