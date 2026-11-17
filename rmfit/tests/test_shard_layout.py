#!/usr/bin/env python3
"""Shard layout (2026-10-10): composite "sum" segments are split by energy with a per-component cost, and every segment with a
free energy shift gets a shard of its own; RMFIT_SPLIT_SUM=0 / RMFIT_SHIFT_SOLO=0 restore the old layout.  Parsing only (no engine).
Uses an archive merged model (a sum segment over 8 exits + 6 free-shift segments); skipped if it is not on disk.
Run: python3 rmfit/tests/test_shard_layout.py"""
import os, sys, tempfile, shutil
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3 (its pyazr 3 version is test_p3_shard_layout.py)")
    sys.exit(77)
from rmfit.model import build_shards, plan_pieces
from pyazr.datasets import SegmentSet
AZR = os.path.join(os.environ.get("RMFIT_ARCHIVE", os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))), "an archived campaign")  # archive model; RMFIT_ARCHIVE


def layout(env):
    old = {k: os.environ.get(k) for k in env}
    os.environ.update(env)
    tmp = tempfile.mkdtemp(prefix="layout_")
    try:
        return build_shards(AZR, 16, tmp)
    finally:
        for k, v in old.items():
            if v is None: os.environ.pop(k, None)
            else: os.environ[k] = v
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    if not os.path.exists(AZR):
        print("SKIPPED (no merged model on disk)"); return 0
    ok = True
    def check(name, cond, detail=""):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}")); ok &= bool(cond)
    segs = SegmentSet.from_file(AZR)
    comp = [s.key for s in segs if s.active and s.operation == "sum"]
    shift = sorted(s.key for s in segs if s.active and s.vary_shift)
    new = layout({}); old = layout({"RMFIT_SPLIT_SUM": "0", "RMFIT_SHIFT_SOLO": "0"})
    npts = lambda L: sum(p.npts for s in L for p in s.pieces)
    check("same number of points in both layouts", npts(new) == npts(old), f"{npts(new)} vs {npts(old)}")
    cnt = lambda L, k: sum(1 for s in L for p in s.pieces if p.base_key == k)
    split_new = [k for k in comp if cnt(new, k) > 1]
    check("the large sum segment (Cierjacks, key 3) is split (new) and whole (old)", 3 in split_new and cnt(old, 3) == 1,
          f"new {cnt(new, 3)}, old {cnt(old, 3)}")
    check("sum segments with a free shift are never split (BandH 9, Brandenburg 81)", all(cnt(new, k) == 1 for k in comp if k in shift),
          str({k: cnt(new, k) for k in comp if k in shift}))
    check("small sum segments stay whole (below the target cost)", all(cnt(new, k) == 1 for k in (2, 79) if k in comp))
    solo = [s for s in new if len(s.pieces) == 1 and s.pieces[0].base_key in shift]
    check("every free-shift segment alone in a shard (new)", sorted(s.pieces[0].base_key for s in solo) == shift, f"{[s.name for s in solo]} vs {shift}")
    mixed = [s.name for s in old if sum(1 for p in s.pieces if p.base_key in shift) and len(s.pieces) > 1]
    check("old layout mixes shift segments with others (the switch works)", bool(mixed), str(mixed))
    pieces = plan_pieces(AZR)
    check("plan_pieces without a target cost keeps one piece per segment", len(pieces) == sum(1 for s in segs if s.active), f"{len(pieces)}")
    print("ALL PASSED" if ok else "FAILED"); return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
