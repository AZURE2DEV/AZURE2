#!/usr/bin/env python3
"""Jacobian-cost shard balancing (2026-10-26, the evaluator: "go ahead with the Jacobian-cost shard balancing").

Sharded measures each shard's residual and resid+Jacobian time once, fits R (alpha + beta k) (k = finite-differenced
shift columns of the shard), repacks with model.pack_jac and caches the model in <workdir>/jac_balance.json.
Checks:
  (a) pack_jac on synthetic pieces: shifted units spread one per shard, the heaviest shift unit gets no extra data,
      predicted loads = R (alpha + beta k), and a greedy plan beats objective-cost LPT on the modelled time;
  (b) shift_fn_groups fixture (free constant shifts + shift functions, 4 shards), repack forced
      (RMFIT_BALANCE_MIN_GAIN=-10): balanced Sharded = unbalanced Sharded = SingleSession in objective, gradient
      J^T r and per-segment chi2; groups still in one shard; jac_balance.json written with repacked=true;
  (c) a second Sharded on the same workdir uses the cache (no measurement) and gives the same objective;
  (d) relvector_p12c fixture (covariance + normgroups, 3 shards), forced repack: objective and gradient = unbalanced;
  (e) RMFIT_JAC_BALANCE=0 / balance=False: no measurement, no cache file.
Needs a pyazr that knows shiftfngroup (AZURE2 dev >= 62392fc)."""
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if not _PYAZR3:
    print("SKIP: a pyazr 3 test; this process runs pyazr 2 (the pyazr 2 version is test_jac_balance.py)")
    sys.exit(77)
from _fixtures import fixture  # noqa: E402  (fixtures/ or fixtures_p3/ by backend)
import numpy as np

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


def synthetic():
    from rmfit.model import Piece, pack_jac, pack
    alpha, beta = 1.35, 3.7
    mk = lambda key, c: Piece(key, [], 10, c, "x")
    shift_units = [[mk(1, 3.0)], [mk(2, 0.5)], [mk(3, 0.4)]]
    rest = [[mk(10 + i, c)] for i, c in enumerate([2.0, 1.5, 1.2, 1.0, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3])]
    shifted = {1, 2, 3}
    kc = lambda u: sum(1 for p in u if p.base_key in shifted)
    bins, loads = pack_jac(shift_units + rest, 4, alpha, beta, kc)
    ks = [kc(b) for b in bins]
    check("(a) shifted units spread one per shard", sorted(ks) == [0, 1, 1, 1], f"{ks}")
    heavy = next(b for b in bins if any(p.base_key == 1 for p in b))
    check("(a) the heaviest shifted unit gets no extra data", len(heavy) == 1, f"{[p.base_key for p in heavy]}")
    pred = [sum(p.cost for p in b) * (alpha + beta * kc(b)) for b in bins]
    check("(a) predicted loads = R (alpha + beta k)", np.allclose(pred, loads), f"{pred} vs {loads}")
    lpt = pack([p for u in shift_units + rest for p in u], 4)
    lpt_max = max(sum(p.cost for p in b) * (alpha + beta * kc(b)) for b in lpt)
    check("(a) balanced plan beats objective-cost LPT on the modelled Jacobian time", max(loads) < lpt_max,
          f"{max(loads):.2f} vs {lpt_max:.2f}")
    print(f"  (a) modelled slowest shard: balanced {max(loads):.2f}, objective-cost LPT {lpt_max:.2f}")


def compare(name, azr, cwd, tmp, nshards, flags, perturb):
    from rmfit.engine import SingleSession, Sharded
    from rmfit.keys import take
    with SingleSession(azr, cwd=cwd, flags=flags, name="single") as ss:
        x = perturb(ss.keys, ss.x0)
        o_s = ss.objective(x)
        r_s, J_s = ss.residuals(x, jac=True)
        g_s = J_s.T @ r_s
        seg_s = ss.score(x).segment_chi2
        keys_s = ss.keys
    os.environ["RMFIT_BALANCE_MIN_GAIN"] = "-10"
    try:
        wb = os.path.join(tmp, f"{name}_bal")
        sb = Sharded(azr, wb, nshards=nshards, flags=flags, balance=True)
    finally:
        os.environ.pop("RMFIT_BALANCE_MIN_GAIN", None)
    su = Sharded(azr, os.path.join(tmp, f"{name}_plain"), nshards=nshards, flags=flags, balance=False)
    try:
        info = sb.balance_info or {}
        check(f"({name}) balanced evaluator measured and repacked", info.get("repacked") is True, f"{info}")
        cache = os.path.join(wb, "jac_balance.json")
        check(f"({name}) jac_balance.json written", os.path.exists(cache) and json.load(open(cache)).get("repacked"))
        check(f"({name}) balance=False: no cache file", not os.path.exists(os.path.join(tmp, f"{name}_plain",
                                                                                        "jac_balance.json")))
        res = {}
        for tag, ev in (("bal", sb), ("plain", su)):
            xe = take(x, keys_s, ev.keys, ev.x0)
            o = ev.objective(xe)
            r, J = ev.residuals(xe, jac=True)
            g = take(J.T @ r, ev.keys, keys_s, np.zeros(len(keys_s)))
            seg = ev.score(xe).segment_chi2
            res[tag] = (o, g, seg)
        for tag in ("bal", "plain"):
            o, g, seg = res[tag]
            check(f"({name}) {tag}: objective = single", np.isclose(o, o_s, rtol=1e-9), f"{o} vs {o_s}")
            check(f"({name}) {tag}: gradient = single (rel 1e-7)", rel(g, g_s) < 1e-7, f"{rel(g, g_s):.2e}")
            check(f"({name}) {tag}: per-segment chi2 = single",
                  all(np.isclose(seg[k], seg_s[k], rtol=1e-8, atol=1e-10) for k in seg_s), "")
        print(f"  ({name}) {info.get('measured_max', float('nan')):.3f} s measured slowest, "
              f"{info.get('predicted_max', float('nan')):.3f} s predicted after repack; alpha "
              f"{info.get('model', {}).get('alpha', float('nan')):.2f} beta {info.get('model', {}).get('beta', float('nan')):.2f}")
    finally:
        sb.close()
        su.close()
    # (c) cache reuse
    sc = Sharded(azr, wb, nshards=nshards, flags=flags)
    try:
        check(f"({name}) second evaluator on the same workdir uses the cache",
              (sc.balance_info or {}).get("source") == "cache", f"{sc.balance_info}")
        o = sc.objective(take(x, keys_s, sc.keys, sc.x0))
        check(f"({name}) cached plan: objective = single", np.isclose(o, o_s, rtol=1e-9), f"{o} vs {o_s}")
        return sc.shards
    finally:
        sc.close()


def main():
    synthetic()
    import rmfit  # noqa: F401
    tmp = tempfile.mkdtemp(prefix="rmfit_jacbal_")
    try:
        if True:
            fix = fixture("shift_fn_groups")
            d = os.path.join(tmp, "sfg")
            shutil.copytree(fix, d)

            def perturb(keys, x0):
                x = np.array(x0, float)
                for j, k in enumerate(keys):
                    if k.kind == "shift":
                        x[j] += 0.004 if k.seg == 1 else -0.003
                    elif k.kind == "shift_fn":
                        x[j] += 0.002 * (k.term + 1)
                    elif k.kind == "width":
                        x[j] *= 1.05
                return x
            shards = compare("b", os.path.join(d, "sfg.azr"), d, tmp, 4, dict(use_brune=True), perturb)
            from rmfit.covariance import shift_fn_groups, shift_groups
            units = dict(shift_fn_groups(os.path.join(d, "sfg.azr")).units())
            units.update(shift_groups(os.path.join(d, "sfg.azr")).units())
            where = {}
            for s in shards:
                for p in s.pieces:
                    if p.base_key in units:
                        where.setdefault(units[p.base_key], set()).add(s.name)
            check("(b) every shift / shift-function group in one shard", all(len(v) == 1 for v in where.values()),
                  f"{where}")
        else:
            print("SKIP (b)/(c): this pyazr does not know shiftfngroup")
        fix = fixture("relvector_p12c")
        d = os.path.join(tmp, "p12c")
        shutil.copytree(fix, d)

        def perturb2(keys, x0):
            x = np.array(x0, float)
            for j, k in enumerate(keys):
                if k.kind == "norm":
                    x[j] *= 1.1
                elif k.kind == "width":
                    x[j] *= 1.03
            return x
        compare("d", os.path.join(d, "ngrel.azr"), d, tmp, 3, dict(use_brune=True), perturb2)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("FAILED: " + ", ".join(failures) if failures else "all passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
