#!/usr/bin/env python3
"""verify() maps a vector by key: a Sharded evaluator orders energy-shift parameters differently from a single session (2026-11-21).

On a model with energy-shift / shift-function parameters, Sharded.keys and SingleSession.keys hold the same parameters in
a different order.  run_level_add passed the sharded polish vector to verify() positionally (``take(x, ev2.keys,
ev2.keys, ...)``, an identity), so every level add on such a model was "verified" well above its polish and rejected
(a benchmark with a shared energy-shift group: +221,000 on a 2.28e6 objective); verify_candidate maps onto the
ledger structure's keys, which come from the sharded evaluator for every structure adopted after init, so polish and sign
candidates were mis-scored the same way.  Fixture fixtures/shift_fn_groups/sfg.azr (p+p; two constant shifts and a
shift function, all free):
  (a) precondition: Sharded and SingleSession hold the same keys in a different order, and a positional verify of a
      sharded-order vector scores a different fit (what the old call sites did);
  (b) verify(..., keys=sharded keys) equals the sharded score; a vector lacking a parameter is refused;
  (c) Campaign.verify_candidate on a structure recorded with the sharded key order equals the sharded score.
Run:  python3 rmfit/tests/test_verify_key_order.py"""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3")
    sys.exit(77)
import numpy as np  # noqa: E402
from _fixtures import fixture  # noqa: E402
from rmfit.campaign import Campaign  # noqa: E402
from rmfit.engine import Sharded, verify  # noqa: E402

FIX = fixture("shift_fn_groups")
FLAGS = dict(use_brune=True)
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def close(a, b):
    return abs(a - b) <= 1e-7 * max(1.0, abs(a), abs(b))


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_verify_keys_")
    try:
        d = os.path.join(tmp, "A")
        os.makedirs(d)
        shutil.copytree(os.path.join(FIX, "data"), os.path.join(d, "data"))
        shutil.copy2(os.path.join(FIX, "sfg.azr"), d)
        c = Campaign.init(d, os.path.join(d, "sfg.azr"), flags=FLAGS, nshards=2, policy=dict(bg_energy=15.0))
        # a second file of the same model: the structure an adopted add / removal would be
        shutil.copy2(c.base_azr, os.path.join(d, "sfg_b.azr"))
        path_b = os.path.join(d, "sfg_b.azr")
        ev = Sharded(path_b, os.path.join(d, "shards_b"), nshards=2, flags=FLAGS, objective=c.objective)
        ks = list(ev.keys)
        # distinct values for every free parameter, so any permutation changes the fit
        truth = {"shift[seg1]": -0.004, "shift[seg2]": -0.007018, "shift_fn[seg1 a0]": 0.006, "shift_fn[seg1 a1]": 0.0015}
        x = np.array(ev.x0, dtype=float)
        for j, k in enumerate(ks):
            if str(k) in truth:
                x[j] = truth[str(k)]
        s_sh = ev.score(x).objective
        ev.close()
        sc_pos, kv = verify(path_b, x, FLAGS, workdir=d, name="v_pos", objective=c.objective)
        check("(a) Sharded and single session: same parameters, different order",
              sorted(map(str, ks)) == sorted(map(str, kv)) and list(map(str, ks)) != list(map(str, kv)),
              f"sharded {list(map(str, ks))} single {list(map(str, kv))}")
        check("(a) a positional verify of the sharded-order vector scores a different fit (the old call sites)",
              not close(sc_pos.objective, s_sh), f"{sc_pos.objective} vs sharded {s_sh}")
        sc_key, _ = verify(path_b, x, FLAGS, workdir=d, name="v_key", objective=c.objective, keys=ks)
        check("(b) verify(keys=sharded keys) = the sharded score", close(sc_key.objective, s_sh), f"{sc_key.objective} vs {s_sh}")
        try:
            verify(path_b, x[:-1], FLAGS, workdir=d, name="v_short", objective=c.objective, keys=ks[:-1])
            check("(b) a vector lacking a parameter is refused", False, "no error")
        except ValueError:
            check("(b) a vector lacking a parameter is refused", True)
        sid = c.ledger.add_structure("sfg_b.azr", path_b, ks, parent=c.structure_id, note="sharded key order")
        c.adopt_structure(sid, note="(test)")
        cid, sc_vc = c.verify_candidate(x, ks, move=dict(kind="test"))
        check("(c) verify_candidate on a structure with the sharded key order = the sharded score",
              close(sc_vc.objective, s_sh), f"{sc_vc.objective} vs {s_sh}")
        c.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all verify key-order checks passed")


if __name__ == "__main__":
    main()
