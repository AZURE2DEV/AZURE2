#!/usr/bin/env python3
"""local.polish_refreeze / Campaign.run_polish(refreeze_every) -- the relvector covariance re-frozen between polish
passes (2026-10-24, the evaluator: "Let's do 2").

Fixture fixtures/relvector_p12c (p+12C, from an archived campaign): ngrel.azr = normgroups + a covariance
with model-scaled (relvector) 2 % terms; ngcov.azr = the same with absolute terms only.  Needs a pyazr with data
covariances (AZURE2 >= 4b71177).  Checks:
  (a) the problem: one polish from a bad start (widths x1.6), frozen there, ends on an objective a fresh session at its
      result does not reproduce;
  (b) polish_refreeze from the same start (tolerance 1e-6): several passes with shrinking re-freeze jumps, ends
      'converged', last jump within tolerance,
      total evaluations <= max_nfev, and its objective = a FRESH SingleSession's score at its x (1e-9);
  (c) every=0 on ngrel, and any 'every' on ngcov (no relvector): identical to polish (same x, same objective);
  (d) Sharded (3 shards): polish_refreeze's objective = a fresh SingleSession's score at its x;
  (e) [2026-10-25] campaigns now fix the relvector model at a reference (test_cov_reference.py), so
      Campaign.run_polish ignores refreeze_every; the old campaign-level check was retired.
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rmfit.engine import SingleSession, Sharded
from rmfit.keys import take
from rmfit.local import polish, polish_refreeze, default_bounds

FIX = os.path.join(HERE, "fixtures", "relvector_p12c")
FLAGS = dict(use_brune=True)
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def fresh(azr, x, keys, tmp, tag):
    """Objective of x in a NEW single session (frozen at x, as rmfit verify sees it)."""
    with SingleSession(azr, cwd=tmp, flags=FLAGS, name="fresh_" + tag) as f:
        return f.score(take(x, keys, f.keys, f.x0)).objective


def bad_start(ev):
    x = np.array(ev.x0, float)
    for j, k in enumerate(ev.keys):
        if k.kind == "width":
            x[j] *= 1.6
    return x


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_refreeze_")
    shutil.copytree(os.path.join(FIX, "data"), os.path.join(tmp, "data"))
    for f in ("ngrel.azr", "ngcov.azr"):
        shutil.copy2(os.path.join(FIX, f), tmp)
    os.makedirs(os.path.join(tmp, "output"), exist_ok=True)
    os.makedirs(os.path.join(tmp, "checks"), exist_ok=True)
    rel, absn = os.path.join(tmp, "ngrel.azr"), os.path.join(tmp, "ngcov.azr")
    try:
        with SingleSession(rel, cwd=tmp, flags=FLAGS, name="rel") as ev:
            check("ngrel has relvector terms", ev.cov is not None and ev.cov.relative)
            xb = bad_start(ev)
            b = default_bounds(ev, xb)
            ev.refreeze_covariance()
            p1 = polish(ev, xb, bounds=b, max_nfev=60)
            f1 = fresh(rel, p1.x, ev.keys, tmp, "p1")
            print(f"  (a) one polish: {p1.start_objective:.4f} -> {p1.objective:.4f} ({p1.status}); fresh at its x {f1:.4f}")
            check("(a) one frozen polish disagrees with a fresh session at its result",
                  abs(f1 - p1.objective) > 1e-6 * abs(f1), f"{p1.objective} vs {f1}")

            # tolerance tightened so the loop is driven to the fixed point (this fixture's relvector effect is
            # small: the default tolerance of 1 stops after one pass, jump 0.83)
            r = polish_refreeze(ev, xb, bounds=b, max_nfev=400, every=20, tol_abs=1e-6, tol_rel=0.0)
            fr_ = fresh(rel, r.x, ev.keys, tmp, "r")
            print(f"  (b) polish_refreeze: {r.start_objective:.4f} -> {r.objective:.4f} ({r.status}: {r.message}); "
                  f"nfev {r.nfev}; fresh {fr_:.4f}")
            for p in r.passes or []:
                print(f"      pass {p['pass_']}: {p['start']:.4f} -> {p['end']:.4f} ({p['status']}, {p['nfev']}), "
                      f"re-frozen {p['refrozen']:.4f} (jump {p['jump']:+.4g})")
            check("(b) several passes", r.passes is not None and len(r.passes) >= 2, f"{r.passes}")
            check("(b) converged", r.status == "converged", r.status)
            last = r.passes[-1] if r.passes else {}
            check("(b) last re-freeze jump within tolerance", bool(last) and abs(last["jump"]) <= 1e-6, f"{last}")
            js = [abs(p["jump"]) for p in r.passes or []]
            check("(b) re-freeze jumps shrink pass by pass (fixed-point iteration)",
                  len(js) >= 2 and all(b2 <= b1 * 1.01 + 1e-9 for b1, b2 in zip(js, js[1:])), f"{js}")
            check("(b) total evaluations <= max_nfev", r.nfev <= 400, f"{r.nfev}")
            check("(b) objective = a fresh session's at its x", np.isclose(r.objective, fr_, rtol=1e-9, atol=1e-9),
                  f"{r.objective} vs {fr_}")
            check("(b) start objective = frozen at the start point", np.isclose(r.start_objective, p1.start_objective,
                                                                               rtol=1e-12), f"{r.start_objective}")

            ev.refreeze_covariance()
            q0 = polish(ev, xb, bounds=b, max_nfev=30)
            ev.refreeze_covariance()
            q1 = polish_refreeze(ev, xb, bounds=b, max_nfev=30, every=0)
            check("(c) every=0 = polish", np.array_equal(q0.x, q1.x) and q0.objective == q1.objective and q1.passes is None,
                  f"{q0.objective} vs {q1.objective}")

        with SingleSession(absn, cwd=tmp, flags=FLAGS, name="abs") as ea:
            check("ngcov has no relvector terms", ea.cov is not None and not ea.cov.relative)
            xa = bad_start(ea)
            ba = default_bounds(ea, xa)
            a0 = polish(ea, xa, bounds=ba, max_nfev=30)
            a1 = polish_refreeze(ea, xa, bounds=ba, max_nfev=30, every=5)
            check("(c) no relvector: polish_refreeze = polish", np.array_equal(a0.x, a1.x) and a0.objective == a1.objective
                  and a1.passes is None, f"{a0.objective} vs {a1.objective}")

        sh = Sharded(rel, os.path.join(tmp, "shards"), nshards=3, flags=FLAGS)
        try:
            with SingleSession(rel, cwd=tmp, flags=FLAGS, name="rel2") as ev:
                xs = take(bad_start(ev), ev.keys, sh.keys, sh.x0)
            rs = polish_refreeze(sh, xs, max_nfev=300, every=20)
            fs = fresh(rel, rs.x, sh.keys, tmp, "s")
            print(f"  (d) Sharded: {rs.start_objective:.4f} -> {rs.objective:.4f} ({rs.status}), {len(rs.passes or [])} "
                  f"passes; fresh {fs:.4f}")
            check("(d) Sharded: objective = a fresh session's at its x", np.isclose(rs.objective, fs, rtol=1e-9, atol=1e-9),
                  f"{rs.objective} vs {fs}")
        finally:
            sh.close()

    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("FAILED: " + ", ".join(failures) if failures else "all passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
