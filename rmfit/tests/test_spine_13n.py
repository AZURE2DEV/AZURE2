#!/usr/bin/env python3
"""Phase-1 gates of the rmfit spine on the small 13N model (seconds per eval):

  * a Sharded evaluator reproduces the SingleSession objective, residuals,
    per-segment chi-squared and Jacobian (columns of split segments summed);
  * the exact sign symmetries hold (level flip; particle-pair flip) and a
    single-channel flip inside one J-group is *not* a symmetry;
  * canonical_signs is invariant under the symmetries;
  * polish() from a jittered start converges back below the start.

Run:  python3 rmfit/tests/test_spine_13n.py [nshards]
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
from rmfit.engine import SingleSession, Sharded
from rmfit.keys import canonical_signs, sign_pattern, flip, take, index_map
from rmfit.local import polish, default_bounds, BoundsPolicy, sanity, free_mask
from pyazr import AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def main(nshards=4):
    tmp = tempfile.mkdtemp(prefix="rmfit_13n_")
    work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(work, junk), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    # free two norms so the penalty rows are exercised
    mdl = AzrModel.from_file(os.path.join(work, "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    mdl.set_segment_norm("kettner_55", vary=True, sys_error=10.0)
    base = mdl.write(os.path.join(work, "base.azr"))
    flags = {}

    print("1. single vs sharded")
    t0 = time.time()
    single = SingleSession(base, cwd=work, flags=flags)
    print(f"   single session up in {time.time()-t0:.1f}s: {len(single.keys)} free, {single.npts} points")
    t0 = time.time()
    sh = Sharded(base, os.path.join(tmp, "shards"), nshards=nshards, flags=flags)
    print(f"   {len(sh.shards)} shards up in {time.time()-t0:.1f}s: "
          + ", ".join(f"{s.name}:{s.npts}pts/{len(s.pieces)}pc" for s in sh.shards))
    check("same key set", set(single.keys) == set(sh.keys),
          f"{set(single.keys) ^ set(sh.keys)}")
    check("same point count", single.npts == sh.npts, f"{single.npts} vs {sh.npts}")
    x = take(single.x0, single.keys, sh.keys, sh.x0)
    check("x0 agrees", np.allclose(x, sh.x0))
    rng = np.random.default_rng(1)
    for trial, xs in enumerate([single.x0, single.x0 * (1 + 0.05 * rng.standard_normal(single.x0.size))]):
        xsh = take(xs, single.keys, sh.keys, sh.x0)
        o1, o2 = single.objective(xs), sh.objective(xsh)
        check(f"objective agrees (trial {trial})", np.isclose(o1, o2, rtol=1e-10), f"{o1} vs {o2}")
        s1, s2 = single.score(xs), sh.score(xsh)
        check(f"per-segment chi2 agrees (trial {trial})",
              all(np.isclose(s1.segment_chi2[k], s2.segment_chi2[k], rtol=1e-10) for k in s1.segment_chi2)
              and set(s1.segment_chi2) == set(s2.segment_chi2))
        r1, J1 = single.residuals(xs, jac=True)
        r2, J2 = sh.residuals(xsh, jac=True)
        check(f"sum of squares agrees (trial {trial})", np.isclose(np.sum(r1**2), np.sum(r2**2), rtol=1e-10))
        # compare Jacobians column-by-column through the key map, using the
        # Gram matrix (row order differs between the two)
        m_ = index_map(single.keys, sh.keys)
        G1 = J1.T @ J1
        G2 = J2[:, m_].T @ J2[:, m_]
        check(f"Jacobian Gram matrices agree (trial {trial})",
              np.allclose(G1, G2, rtol=1e-8, atol=1e-8 * np.max(np.abs(G1))),
              f"max abs diff {np.max(np.abs(G1-G2)):.3e} of {np.max(np.abs(G1)):.3e}")
        g1, g2 = J1.T @ r1, J2[:, m_].T @ r2
        check(f"gradients agree (trial {trial})", np.allclose(g1, g2, rtol=1e-8, atol=1e-8 * np.max(np.abs(g1))))
    single.close()

    print("\n2. sign symmetries (sharded evaluator)")
    x = sh.x0.copy()
    o0 = sh.objective(x)
    widths = [(j, k) for j, k in enumerate(sh.keys) if k.kind == "width" and x[j] != 0]
    by_level = {}
    for j, k in widths:
        by_level.setdefault(k.level, []).append(k)
    lvl, ks = max(by_level.items(), key=lambda kv: len(kv[1]))
    o_lvl = sh.objective(flip(x, sh.keys, ks))
    check(f"level flip is a symmetry ({lvl}, {len(ks)} amplitudes)", np.isclose(o_lvl, o0, rtol=1e-9), f"{o_lvl} vs {o0}")
    pairs = sorted({k.pair for _, k in widths})
    for pair in pairs:
        ks = [k for _, k in widths if k.pair == pair]
        o_p = sh.objective(flip(x, sh.keys, ks))
        print(f"      pair {pair} flip ({len(ks)} amplitudes): {o_p:.6f} vs {o0:.6f}"
              + ("  (symmetry)" if np.isclose(o_p, o0, rtol=1e-9) else "  (NOT a symmetry -- photon pair with external capture?)"))
    if len(ks) > 1:
        o_1 = sh.objective(flip(x, sh.keys, ks[:1]))
        check("a single amplitude flip changes the objective", not np.isclose(o_1, o0, rtol=1e-6), f"{o_1} vs {o0}")
    c0 = canonical_signs(x, sh.keys)
    c1 = canonical_signs(flip(x, sh.keys, by_level[lvl]), sh.keys)
    check("canonical pattern invariant under a level flip", c0 == c1, f"{c0} vs {c1}")
    p1 = [k for _, k in widths if k.pair == pairs[0]]
    c2 = canonical_signs(flip(x, sh.keys, p1), sh.keys)
    check("canonical pattern invariant under a pair flip", c0 == c2, f"{c0} vs {c2}")
    c3 = canonical_signs(flip(x, sh.keys, ks[:1]), sh.keys)
    check("canonical pattern changes under a single flip", c0 != c3)
    print(f"      raw {sign_pattern(x, sh.keys)}  canonical {c0}")

    print("\n3. polish from a jittered start")
    policy = BoundsPolicy(bg_energy=15.0, e_window=0.5)
    b = default_bounds(sh, x, policy)
    rep = sanity(sh, x, policy, b)
    print(f"      sanity at x0: {rep.summary()}; max theta2 {rep.max_theta2:.3g}")
    xj = x * (1 + 0.03 * rng.standard_normal(x.size))
    xj = b.clip(xj)
    oj = sh.objective(xj)
    t0 = time.time()
    fr = polish(sh, xj, bounds=b, policy=policy, max_nfev=200, verbose=1)
    print(f"      jittered {oj:.3f} -> {fr.objective:.3f} ({fr.status}; {fr.message}; "
          f"nfev={fr.nfev} njev={fr.njev}; {time.time()-t0:.0f}s; dead={len(fr.dead)}; at bound={len(fr.bounds_hit)})")
    check("polish improves the jittered start", fr.objective < oj)
    check("polish ends near or below the original", fr.objective <= o0 * 1.05, f"{fr.objective} vs {o0}")

    print("\n4. bake + fresh-session verification")
    from rmfit.engine import bake
    # the rwa -> physical -> rwa round trip through <levels> is not exact
    # (AZURE2's TransformIn is numerical); on this jittered 13N point it is
    # ~4e-5 in chi2, on a converged large fit ~1e-8
    azr, sav, sc = bake(base, fr.x, sh.keys, os.path.join(tmp, "baked.azr"), flags,
                        workdir=os.path.join(tmp, "variants"), verify_rtol=1e-3)
    print(f"      round-trip noise: {sc.roundtrip}")
    check("bake round-trips the polished objective", np.isclose(sc.objective, fr.objective, rtol=1e-3),
          f"{sc.objective} vs {fr.objective}")
    check("companion .sav written", os.path.exists(sav))
    fg = polish(sh, xj, free=free_mask(sh, jpi=["1/2+"]), bounds=b, policy=policy, max_nfev=30)
    check("group polish runs and does not worsen", fg.objective <= oj + 1e-9, fg.message)
    sh.close()
    shutil.rmtree(tmp, ignore_errors=True)

    print()
    if failures:
        print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
        sys.exit(1)
    print("all spine checks passed")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
