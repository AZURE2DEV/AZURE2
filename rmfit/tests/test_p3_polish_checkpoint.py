#!/usr/bin/env python3
"""polish writes its current iterate every 10th Jacobian when RMFIT_POLISH_CHECKPOINT names a file (2026-11-22).

Why: with target integration folded in, a 400-nfev polish can run ~20 h on one shard, and its vector lived
only in memory until least_squares returned; a qdel, node loss or crash lost all of it.

Fixture tests/fixtures_p3/shift_fn_groups with the 0.8 and 1.0 MeV levels freed (rmfit.azr3 Model), start vector perturbed so the polish has work to do, on a
SingleSession and on a Sharded evaluator (3 shards):
  (a) without the variable: no file is written;
  (b) with it: the polish result is bit-identical to (a), the file exists, holds a multiple of 10 Jacobians, the
      evaluator's own keys, and an objective that a fresh evaluation of the stored x reproduces; it lies between the
      start and the final objective;
  (c) load_checkpoint refuses a key list that differs.
    AZURE2_ROOT=<build> python3 rmfit/tests/test_p3_polish_checkpoint.py      (about a minute; pyazr 3 -- the pyazr 2 version is
    test_polish_checkpoint.py)"""
import os
import shutil
import sys
import tempfile

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.pop("RMFIT_POLISH_CHECKPOINT", None)

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from _fixtures import fixture  # noqa: E402
from rmfit.backend import PYAZR3  # noqa: E402
if not PYAZR3:
    print("SKIP: a pyazr 3 test; this process runs pyazr 2 (the pyazr 2 version is test_polish_checkpoint.py)")
    sys.exit(77)
FIX = fixture("shift_fn_groups")  # fixtures_p3/ on pyazr 3
failures = []


def check(cond, what):
    print(("  ok    " if cond else "  FAIL  ") + what, flush=True)
    if not cond:
        failures.append(what)


def run(ev, label, work):
    from rmfit.local import polish, default_bounds, BoundsPolicy, load_checkpoint
    x1 = np.array(ev.x0, float)
    for j, k in enumerate(ev.keys):
        # 1.4 / 0.01 (were 1.15 / 0.004, 2026-11-26): on engine d2a4cbe8 the Sharded polish of the milder start converged in
        # 9 Jacobians, short of the 10 that write a checkpoint, even with max_nfev 80
        if k.kind == "width":
            x1[j] *= 1.4
        elif k.kind == "shift":
            x1[j] += 0.01
    b = default_bounds(ev, x1, BoundsPolicy())
    # max_nfev 80 (was 40): on engine d2a4cbe8 the Sharded polish of this start stopped at njev 9, short of the 10 that
    # write a checkpoint (2026-11-26)
    kw = dict(bounds=b, max_nfev=80, ftol=1e-15, xtol=1e-15, gtol=1e-15)
    ck = os.path.join(work, f"ckpt_{label}.npz")
    os.environ.pop("RMFIT_POLISH_CHECKPOINT", None)
    a = polish(ev, x1, **kw)
    check(not os.path.exists(ck), f"({label}a) no checkpoint without the variable")
    os.environ["RMFIT_POLISH_CHECKPOINT"] = ck
    try:
        b2 = polish(ev, x1, **kw)
    finally:
        os.environ.pop("RMFIT_POLISH_CHECKPOINT", None)
    check(np.array_equal(a.x, b2.x) and a.objective == b2.objective and a.njev == b2.njev,
          f"({label}b) result bit-identical with the checkpoint on ({a.objective:.6f}, njev {a.njev})")
    if a.njev < 10:
        check(False, f"({label}b) the polish reached 10 Jacobians (njev {a.njev}): the test has no teeth")
        return
    check(os.path.exists(ck) and not os.path.exists(ck + ".tmp"), f"({label}b) checkpoint written, no temp file left")
    x, info = load_checkpoint(ck, ev.keys)
    check(info["njev"] % 10 == 0 and 10 <= info["njev"] <= a.njev, f"({label}b) njev {info['njev']} (final {a.njev})")
    o = ev.objective(x)
    check(abs(o - info["objective"]) <= 1e-9 * abs(o), f"({label}b) stored objective {info['objective']:.9f} = fresh {o:.9f}")
    check(a.objective - 1e-9 <= info["objective"] <= a.start_objective + 1e-9,
          f"({label}b) between final {a.objective:.6f} and start {a.start_objective:.6f}")
    try:
        load_checkpoint(ck, list(ev.keys)[:-1])
        check(False, f"({label}c) mismatched keys refused")
    except ValueError:
        check(True, f"({label}c) mismatched keys refused")


def main():
    import rmfit  # noqa: F401  (puts AZURE2_ROOT on sys.path)
    from rmfit.engine import SingleSession, Sharded
    work = tempfile.mkdtemp(prefix="rmfit_ckpt_")
    try:
        fix = os.path.join(work, "fix")
        shutil.copytree(FIX, fix)
        os.makedirs(os.path.join(fix, "output"), exist_ok=True); os.makedirs(os.path.join(fix, "checks"), exist_ok=True)
        azr = os.path.join(fix, "sfg.azr")
        from rmfit.azr3 import Model
        mdl = Model.from_file(azr)
        for lv in mdl.levels:
            if float(lv.energy) in (0.8, 1.0):       # free these levels' energies and widths (as test_p3_jac_skip_columns)
                lv.set_fixed(False)
                for c in lv.channels:
                    c.channel_fixed = False
        mdl.write(azr)
        S = SingleSession(azr, cwd=fix, name="ck")
        try:
            run(S, "single", work)
        finally:
            S.close()
        H = Sharded(azr, os.path.join(work, "sh"), nshards=3)
        try:
            print(f"  info  sharded: {len(H.shards)} shard(s)")
            run(H, "sharded", work)
        finally:
            H.close()
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("ALL PASSED" if not failures else f"FAILED: {failures}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
