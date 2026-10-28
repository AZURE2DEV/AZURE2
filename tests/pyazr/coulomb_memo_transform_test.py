#!/usr/bin/env python3
"""The transform to physical parameters leaves nothing in the Coulomb memo (2026-10-28).

CNuc::TransformOut works at the LEVEL energies -- penetrabilities and shift functions there, and the external part
of every radiative width, an external-capture radial integral at the level energy.  A fit or a sampler moves those
energies at every call, so the (l, radius, energy) entries it stored never recurred: each transform_rwa added ~5,000
entries per worker for good in the 12C(a,g) MCMC chains of 2026-10-26 (--fix-gamma-total calls it once or more per
log-probability), ~4 GB per worker after 17,000 evaluations.  The transform now runs inside a
CoulFuncCache::TransientScope: it shares Coulomb functions within the call and keeps none afterwards.

Checks on tests/13N (12C + p with external capture), Brune and Park:
  (1) 20 transform_rwa calls at level energies moved by ~1 keV leave the memo's entry count unchanged;
  (2) the transform is a pure function of the parameters: the first point evaluated again gives the same vector
      bit for bit, and so does a fresh session (nothing depends on what the memo held);
  (3) the objective at those points is unchanged by the transforms in between (bit for bit).
dev faed9a5 fails (1).  Needs the compiled engine; skips cleanly without it.
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "tests", "13N")
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from pyazr import azure2
except ImportError as e:
    print(f"SKIP: {e}")
    sys.exit(0)

work = tempfile.mkdtemp(prefix="coulomb_memo_transform.")
try:
    for mode in ("brune", "park"):
        d = os.path.join(work, mode)
        os.makedirs(os.path.join(d, "output")); os.makedirs(os.path.join(d, "checks"))
        shutil.copytree(os.path.join(SRC, "data"), os.path.join(d, "data"))
        shutil.copy(os.path.join(SRC, "13N.azr"), os.path.join(d, "13N.azr"))
        kw = dict(use_park=True, use_brune=False) if mode == "park" else {}
        M = azure2(os.path.join(d, "13N.azr"), cwd=d, **kw)
        x0 = np.asarray(M.params_rwa, float).copy()
        kinds = {p.free_index: p.kind for p in M.parameters if p.free_index is not None and not p.fixed}
        energy = np.array([kinds.get(i) == "energy" for i in range(x0.size)])
        rng = np.random.default_rng(7)
        X = []
        for _ in range(20):
            x = x0.copy()
            x[energy] += 1e-3 * rng.standard_normal(energy.sum())
            X.append(x)
        o_before = M.objective(X[0])
        first = np.asarray(M.transform_rwa(X[0]), float)
        e0 = M.cache_stats()["entries"]
        for x in X[1:]:
            M.transform_rwa(x)
        grown = M.cache_stats()["entries"] - e0
        check(f"{mode}: 19 transforms at moved level energies add no memo entries ({energy.sum()} free energies)",
              grown == 0, f"+{grown} entries")
        again = np.asarray(M.transform_rwa(X[0]), float)
        check(f"{mode}: the transform repeats bit for bit", np.array_equal(first, again),
              f"max |d| {np.max(np.abs(first - again)):.3e}")
        check(f"{mode}: the objective is unchanged by the transforms in between",
              M.objective(X[0]) == o_before, f"{M.objective(X[0])!r} vs {o_before!r}")
        del M
        F = azure2(os.path.join(d, "13N.azr"), cwd=d, **kw)
        check(f"{mode}: a fresh session gives the same transform", np.array_equal(np.asarray(F.transform_rwa(X[0]), float), first))
        del F
finally:
    shutil.rmtree(work, ignore_errors=True)

if failures:
    print(f"\n{len(failures)} check(s) FAILED")
    sys.exit(1)
print("\nall checks passed")
