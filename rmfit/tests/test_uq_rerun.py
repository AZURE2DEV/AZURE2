#!/usr/bin/env python3
"""rmfit.uq re-run guards and the per-chain wall estimate (2026-11-13).  Login node, seconds:
    python3 rmfit/tests/test_uq_rerun.py"""
import contextlib, io, os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.uq import _rerun_plan, Checkpoint, config

FAILS = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""))
    if not ok:
        FAILS.append(name)


def main():
    cfg = config({})
    check("defaults: rerun_rhat_max 1.5, rerun_max_hours 12", cfg["rerun_rhat_max"] == 1.5 and cfg["rerun_max_hours"] == 12.0)
    # R-hat 1.013 after 1.39 h -> re-run (~2.8 h)
    ok, why = _rerun_plan({"rhat_max": 1.013}, cfg, 0, 1.39 * 3600)
    check("a near-miss re-runs", ok and "attempt 2 of 2" in why, why)
    # attempts used up
    ok, why = _rerun_plan({"rhat_max": 1.013}, cfg, 1, 100.0)
    check("no re-run when the attempts are used", not ok and "2 of 2" in why, why)
    # R-hat 2.0 -> no re-run, whatever the time
    ok, why = _rerun_plan({"rhat_max": 2.0}, cfg, 0, 600.0)
    check("R-hat 2.0 > 1.5: no automatic re-run", not ok and "rerun_rhat_max" in why, why)
    # a 7 h attempt would re-run for ~14 h > 12 h
    ok, why = _rerun_plan({"rhat_max": 1.05}, cfg, 0, 7 * 3600)
    check("re-run estimate 14 h > 12 h: no automatic re-run", not ok and "rerun_max_hours" in why, why)
    # both guards off -> the old behaviour
    ok, why = _rerun_plan({"rhat_max": 3.0}, dict(cfg, rerun_rhat_max=0, rerun_max_hours=0), 0, 1e6)
    check("guards off (0): re-run as before", ok, why)
    # Checkpoint prints the draw-phase estimate after the first 100 draws
    tmp = tempfile.mkdtemp(prefix="rmfit_uqr_")
    cb = Checkpoint(tmp, 150, "a0", n_draws=150)
    st = dict(accept=0.9, n_leapfrog=7, depth=3, divergent=False, energy=1.0)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        for it in range(300):
            cb(it, [0.0, 1.0], st)
    out = buf.getvalue()
    check("estimate line printed once, after 100 draws", out.count("estimate after 100 draws") == 1 and "h left in attempt a0" in out, out[-300:])
    cb2 = Checkpoint(tmp, 150, "a0")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        for it in range(300):
            cb2(it, [0.0], st)
    check("no estimate without n_draws (old callers)", "estimate" not in buf.getvalue())
    if FAILS:
        print(f"FAILED {len(FAILS)}: {', '.join(FAILS)}")
        sys.exit(1)
    print("all uq re-run checks passed")


if __name__ == "__main__":
    main()
