#!/usr/bin/env python3
"""The multi-level channel-pattern screen on the small 13N model: a synthetic
multi_flip move goes through search.screen with backtracking GN iterations,
and a sign round persists its screen table.  Run:
python3 rmfit/tests/test_multiscreen_13n.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.engine import Sharded
from rmfit.local import BoundsPolicy, default_bounds, polish
from rmfit.moves import Move
from rmfit.search import SearchConfig, screen, sign_round, crossover_round
from rmfit.moves import crossover_moves
from rmfit.keys import flip, sign_distance
from pyazr import AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_multiscreen_")
    work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    mdl = AzrModel.from_file(os.path.join(work, "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    base = mdl.write(os.path.join(work, "base.azr"))
    policy = BoundsPolicy(bg_energy=15.0, e_window=0.5)
    cfg = SearchConfig(top_k=2, pairs_per_group=1, stage1_nfev=20, stage2_nfev=30, union_nfev=30,
                       enumerate_channels=True, frozen_top=2, gn_iters_multi=3, multi_top=1)

    sh = Sharded(base, os.path.join(tmp, "shards"), nshards=3, flags={})
    x = sh.x0.copy()
    b = default_bounds(sh, x, policy)
    f0 = polish(sh, x, bounds=b, policy=policy, max_nfev=60)
    x0, o0 = f0.x, f0.objective
    print(f"  polished start: {o0:.3f} ({f0.status})")

    print("1. synthetic multi-level pattern through the screen")
    widths = [(j, k) for j, k in enumerate(sh.keys) if k.kind == "width" and x0[j] != 0]
    by_jpi = {}
    for j, k in widths:
        by_jpi.setdefault(k.jpi, []).append(k)
    jpi, ks = max(by_jpi.items(), key=lambda t: len(t[1]))
    mv = Move("multi_flip", ks[:2], description=f"flip 2 of {jpi}: " + ", ".join(k.level for k in ks[:2]))
    single = Move("flip", ks[:1], description=f"flip {ks[0]}")
    res = screen(sh, x0, [mv, single], b, cfg, log=print)
    check("both moves screened", len(res) == 2)
    for r in res:
        check(f"finite screen numbers for {r.move.kind}", np.isfinite(r.obj_flip) and np.isfinite(r.obj_gn))
        check(f"GN never worse than the flipped point for {r.move.kind}", r.obj_gn <= r.obj_flip + 1e-9,
              f"{r.obj_gn} vs {r.obj_flip}")
        check(f"x_gn keeps the flipped signs for {r.move.kind}",
              all(np.sign(r.x_gn[sh.keys.index(k)]) == -np.sign(x0[sh.keys.index(k)]) for k in r.move.targets))

    print("2. a sign round persists its screen table")
    xb, ob, stat, trials = sign_round(sh, x0, o0, cfg, policy, log=print)
    check("round no worse than the start", ob <= o0 + 1e-9)
    check("screen table recorded", len(stat.screen) == stat.n_screened and stat.n_screened > 0)
    check("screen rows carry move/kind/frozen/gn",
          all(set(r) == {"move", "kind", "frozen", "gn"} for r in stat.screen))
    check("screen table sorted by gn", all(stat.screen[i]["gn"] <= stat.screen[i + 1]["gn"]
                                           for i in range(len(stat.screen) - 1)))
    print("3. basin crossover with a donor vector")
    # a donor with a genuinely different sign pattern (the raw flip: a polish
    # of this small model walks straight back to x0's pattern)
    donor_p = flip(x0, sh.keys, [ks[0], ks[-1]])
    d0, diff = sign_distance(x0, donor_p, sh.keys)
    check("donor differs from the incumbent", d0 >= 1, str(d0))
    mv = crossover_moves(x0, donor_p, sh.keys)
    check("crossover proposes a move per differing channel (+ groups/all)", (d0 == 0) == (len(mv) == 0),
          f"distance {d0}, {len(mv)} moves")
    check("every crossover move flips only differing amplitudes",
          all(set(m.targets) <= set(diff) for m in mv))
    xc, oc, statc, trialsc = crossover_round(sh, x0, o0, donor_p, cfg, policy, log=print, label="test donor")
    check("crossover round no worse than the start", oc <= o0 + 1e-9)
    check("crossover polished every proposal", statc.n_polished == len(mv), f"{statc.n_polished} vs {len(mv)}")
    sh.close()
    shutil.rmtree(tmp, ignore_errors=True)
    print("FAILURES: " + ", ".join(failures) if failures else "ALL OK")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
