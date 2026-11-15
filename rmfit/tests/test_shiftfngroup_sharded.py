#!/usr/bin/env python3
"""rmfit Sharded with shift-function groups (`shiftfngroup <label>`, AZURE2 dev 62392fc, 2026-10-24).

Every member of a shift-function group uses the shift function of the group's first active member, which declares
it; so build_shards keeps a group's pieces in ONE shard and never splits a member (a split owner would declare the
function twice in one shard, an engine error; a shard without the owner would give its first member the function's
role).  Fixture tests/fixtures/shift_fn_groups/sfg.azr (AZURE2 tests/pyazr/data/shift_fn_groups extended):
  F  segments 1, 2 (24 points each): function a0, a1 free, constant shifts free (shift-solo unit)
  G  segments 3, 4 (48 points each): function a0 = 0.002 and the constant shifts FIXED (a packed unit, would be split)
     segments 5, 6 (48 points each): no group (split into pieces)
Per parametrization (standard, Brune), at a perturbed vector:
  (a) SingleSession's chi2 and rows equal the engine's (calculate_chi2_rwa, residual_jacobian);
  (b) Sharded (4 shards): every group in one shard, G's members unsplit, each shard file declares each group's function
      once, on the group's first line; objective, gradient and per-segment chi2 equal SingleSession's;
  (c) Sharded's gradient agrees with central differences on the shift columns.
Needs a pyazr that knows shiftfngroup (AZURE2 dev >= 62392fc):
    AZURE2_ROOT=<build> python3 rmfit/tests/test_shiftfngroup_sharded.py      (about a minute)"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
FIX = os.path.join(HERE, "fixtures", "shift_fn_groups")
failures = []


def check(cond, what):
    print(("  ok    " if cond else "  FAIL  ") + what, flush=True)
    if not cond:
        failures.append(what)


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


def perturbed(keys, x0):
    x = np.array(x0, float)
    for j, k in enumerate(keys):
        if k.kind == "shift":
            x[j] += 0.004 if k.seg == 1 else -0.003
        elif k.kind == "shift_fn":
            x[j] += 0.002 * (k.term + 1)
    return x


def main():
    import rmfit  # noqa: F401  (puts AZURE2_ROOT on sys.path)
    import pyazr.datasets as ds
    if "shiftfngroup" not in getattr(ds, "_KEYWORDS", {}):
        print(f"SKIP: the pyazr at {os.path.dirname(ds.__file__)} does not know shiftfngroup (AZURE2 dev >= 62392fc)")
        return 0
    from pyazr import azure2
    from rmfit.engine import SingleSession, Sharded
    from rmfit.keys import take
    from rmfit.covariance import shift_fn_groups, line_keyword, segment_lines
    work = tempfile.mkdtemp(prefix="rmfit_sfg_")
    try:
        shutil.copytree(FIX, os.path.join(work, "fix"))
        fix = os.path.join(work, "fix")
        os.makedirs(os.path.join(fix, "output"), exist_ok=True); os.makedirs(os.path.join(fix, "checks"), exist_ok=True)
        azr = os.path.join(fix, "sfg.azr")
        groups = shift_fn_groups(azr)
        check(groups.units() == {1: 1, 2: 1, 3: 3, 4: 3}, f"fixture groups F = 1,2 and G = 3,4: {groups.units()}")
        for brune in (False, True):
            flags = dict(use_brune=brune)
            print(f"== use_brune={brune}", flush=True)
            S = SingleSession(azr, cwd=fix, flags=flags, name=f"s_{int(brune)}")
            E = azure2(azr, cwd=fix, **flags)
            try:
                kinds = sorted({k.kind for k in S.keys})
                check(sum(k.kind == "shift_fn" for k in S.keys) == 2 and sum(k.kind == "shift" for k in S.keys) == 2,
                      f"free parameters: F's a0, a1 and two constant shifts ({len(S.keys)} keys, kinds {kinds})")
                x1 = perturbed(S.keys, S.x0)
                # (a)
                rS, JS = S.residuals(x1, jac=True)
                rE, JE = E.residual_jacobian(x1)
                n = S.npts
                check(rel(rS[:n], rE) < 1e-9, f"(a) rows = engine's (rel {rel(rS[:n], rE):.1e})")
                sc = S.score(x1)
                etot = float(np.sum(E.calculate_chi2_rwa(x1)))
                check(abs(sc.chi2 - etot) < 1e-10 * etot, f"(a) chi2 {sc.chi2:.8f} = engine {etot:.8f}")
                sseg = np.array([sc.segment_chi2[k] for k, _ in S.row_ranges])
                # (b)
                H = Sharded(azr, os.path.join(work, f"sh_{int(brune)}"), nshards=4, flags=flags)
                try:
                    lay = [sorted(p.base_key for p in s.pieces) for s in H.shards]
                    print("  shards (base keys of their pieces):", lay, flush=True)
                    for lab, mem in (("F", (1, 2)), ("G", (3, 4))):
                        holders = [i for i, keys in enumerate(lay) if any(k in mem for k in keys)]
                        check(len(holders) == 1, f"(b) group {lab} sits in one shard ({holders})")
                    check(all(sum(k == m for keys in lay for k in keys) == 1 for m in (3, 4)),
                          "(b) G's members are not split")
                    check(any(sum(k == m for keys in lay for k in keys) > 1 for m in (5, 6)),
                          "(b) an ungrouped 48-point segment IS split (the test can split)")
                    for s in H.shards:
                        lines = segment_lines(s.azr)
                        for lab in ("F", "G"):
                            mine = [t for t in lines if line_keyword(t, "shiftfngroup") == lab]
                            if mine:
                                decl = [i for i, t in enumerate(mine) if "shiftfn" in [x.lower() for x in t]]
                                check(decl == [0], f"(b) {os.path.basename(s.azr)}: group {lab} declares its function "
                                                   f"once, on its first line ({decl})")
                    h1 = take(x1, S.keys, H.keys, H.x0)
                    check(set(H.keys) == set(S.keys), "(b) Sharded and Single have the same free parameters")
                    oS, oH = S.objective(x1), H.objective(h1)
                    check(abs(oS - oH) < 1e-9 * oS, f"(b) objective Single {oS:.8f} = Sharded {oH:.8f}")
                    rH, JH = H.residuals(h1, jac=True)
                    gS = 2 * JS.T @ rS
                    gH = take(2 * JH.T @ rH, H.keys, S.keys, np.zeros(len(S.keys)))
                    check(rel(gH, gS) < 1e-7, f"(b) gradient Single = Sharded (rel {rel(gH, gS):.1e})")
                    hs = H.score(h1)
                    hseg = np.array([hs.segment_chi2[k] for k, _ in S.row_ranges])
                    check(rel(hseg, sseg) < 1e-9, f"(b) per-segment chi2 Single = Sharded (rel {rel(hseg, sseg):.1e})")
                    # (c)
                    gh = 2 * JH.T @ rH
                    worst = 0.0
                    cols = [j for j, k in enumerate(H.keys) if k.kind in ("shift", "shift_fn")]
                    for j in cols:
                        h = 1e-6 * max(1.0, abs(h1[j]))
                        xp, xm = h1.copy(), h1.copy()
                        xp[j] += h; xm[j] -= h
                        fd = (H.objective(xp) - H.objective(xm)) / (2 * h)
                        worst = max(worst, abs(fd - gh[j]) / max(abs(gh).max(), 1e-12))
                    check(worst < 1e-5, f"(c) Sharded gradient = central differences on {len(cols)} shift columns "
                                        f"(worst {worst:.1e} of the largest)")
                finally:
                    H.close()
            finally:
                S.close()
                E.close()
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("ALL PASSED" if not failures else f"FAILED: {failures}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
