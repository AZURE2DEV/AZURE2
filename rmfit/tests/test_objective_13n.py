#!/usr/bin/env python3
"""The configurable objective (ObjectiveSpec): reduced-chi2 data-set weights
and Sivia's loss, on the small 13N model.  Run:
python3 rmfit/tests/test_objective_13n.py
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
from rmfit.engine import SingleSession, Sharded, ObjectiveSpec, rho_sivia, verify
from rmfit.model import single_variant
from rmfit.local import BoundsPolicy, default_bounds, polish
from pyazr import AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    print("1. Sivia loss derivatives")
    z = np.array([0.0, 1e-5, 5e-4, 1e-3, 2e-3, 0.1, 1.0, 4.0, 30.0, 500.0, 2000.0])
    rho, d1, d2 = rho_sivia(z)
    check("rho(0) = 0 and rho'(0) = 1/2", abs(rho[0]) < 1e-12 and abs(d1[0] - 0.5) < 1e-12)
    h = 1e-6
    zz = np.array([2e-3, 0.1, 1.0, 4.0, 30.0, 500.0])
    fd1 = (rho_sivia(zz + h)[0] - rho_sivia(zz - h)[0]) / (2 * h)
    fd2 = (rho_sivia(zz + h)[1] - rho_sivia(zz - h)[1]) / (2 * h)
    check("rho' matches finite differences", np.allclose(fd1, rho_sivia(zz)[1], rtol=1e-5, atol=1e-8),
          str(np.max(np.abs(fd1 - rho_sivia(zz)[1]))))
    check("rho'' matches finite differences", np.allclose(fd2, rho_sivia(zz)[2], rtol=1e-4, atol=1e-8),
          str(np.max(np.abs(fd2 - rho_sivia(zz)[2]))))
    a, b = rho_sivia(np.array([1e-3 - 1e-9]))[0][0], rho_sivia(np.array([1e-3 + 1e-9]))[0][0]
    check("series and closed form agree at the switch", abs(a - b) < 1e-9, f"{a} vs {b}")
    check("rho' positive everywhere", np.all(d1 > 0))
    check("large-z tail ~ 2 log z", abs(rho[-1] - (2 * np.log(2000.0) - 2 * np.log(2))) < 1e-6)

    from rmfit.engine import rho_sivia_scaled
    c = np.array([0.5, 4.0, 25.0, 200.0])
    zz = np.array([0.3, 8.0, 25.0, 900.0])
    h = 1e-6
    fd1 = (rho_sivia_scaled(zz + h, c)[0] - rho_sivia_scaled(zz - h, c)[0]) / (2 * h)
    fd2 = (rho_sivia_scaled(zz + h, c)[1] - rho_sivia_scaled(zz - h, c)[1]) / (2 * h)
    check("scaled rho' matches finite differences", np.allclose(fd1, rho_sivia_scaled(zz, c)[1], rtol=1e-5, atol=1e-8))
    check("scaled rho'' matches finite differences", np.allclose(fd2, rho_sivia_scaled(zz, c)[2], rtol=1e-4, atol=1e-8))
    check("scaled loss ~ quadratic when z << c", abs(rho_sivia_scaled(np.array([1e-3]), np.array([1e3]))[0][0] / 5e-4 - 1) < 1e-3)
    check("scaled loss ~ c * 2 log(z/c) when z >> c", abs(rho_sivia_scaled(np.array([4e5]), np.array([4.0]))[0][0]
                                                        - 4 * (2 * np.log(1e5) - 2 * np.log(2))) < 1e-3)

    tmp = tempfile.mkdtemp(prefix="rmfit_objective_")
    work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True)
    os.makedirs(os.path.join(work, "output"))
    mdl = AzrModel.from_file(os.path.join(work, "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    base = mdl.write(os.path.join(work, "base.azr"))
    # reference scales: the plain baseline's chi2/N per data file (what Campaign.init does)
    ref = SingleSession(single_variant(base, "ref", tmp), cwd=tmp, name="ref")
    sc_ref = ref.score(ref.x0)
    ds_of = ref.dataset_of_key()
    n_of, chi_of = {}, {}
    for key, (a, b) in ref.row_ranges:
        ds = str(ds_of[int(key)]); n_of[ds] = n_of.get(ds, 0) + (b - a)
        chi_of[ds] = chi_of.get(ds, 0.0) + sc_ref.segment_chi2[key]
    scales = {ds: max(chi_of[ds] / n_of[ds], 1.0) for ds in n_of}
    ref.close()
    spec = ObjectiveSpec(dataset_weight="reduced", dataset_key="file", loss="sivia", loss_scale="dataset", scales=scales)
    policy = BoundsPolicy(bg_energy=15.0, e_window=0.5)

    print("2. single session: plain vs weighted objective")
    plain = SingleSession(single_variant(base, "plain", tmp), cwd=tmp, name="plain")
    x0 = plain.x0.copy()
    rng = np.random.default_rng(3)
    x1 = x0 * (1 + 0.05 * rng.standard_normal(x0.size))
    sc_plain = [plain.score(x) for x in (x0, x1)]
    check("plain spec unchanged: objective = chi2 + penalty", all(abs(s.objective - (s.chi2 + s.penalty + s.prior)) < 1e-9 for s in sc_plain))
    check("plain objective() == score().objective",
          all(abs(plain.objective(x) - s.objective) < 1e-9 * abs(s.objective) for x, s in zip((x0, x1), sc_plain)),
          "; ".join(f"{plain.objective(x)} vs {s.objective}" for x, s in zip((x0, x1), sc_plain)))
    plain.close()
    single = SingleSession(single_variant(base, "weighted", tmp), cwd=tmp, name="weighted", objective=spec)
    check("data-row weights sum to N", abs(np.sum(single.row_sw ** 2) - single.npts) < 1e-9,
          f"{np.sum(single.row_sw ** 2)} vs {single.npts}")
    per_ds = {}
    for k2, (a, b) in single.row_ranges:
        ds = os.path.basename(single.segments.by_key(k2).data_file)
        per_ds[ds] = per_ds.get(ds, 0.0) + float(np.sum(single.row_sw[a:b] ** 2))
    check("every data set gets equal total weight",
          np.allclose(list(per_ds.values()), single.npts / len(per_ds)), str(per_ds))
    half = SingleSession(single_variant(base, "half", tmp), cwd=tmp, name="half",
                         objective=ObjectiveSpec(dataset_weight="reduced", dataset_key="file", weight_power=0.5))
    check("tempered weights (p=0.5) sum to N", abs(np.sum(half.row_sw ** 2) - half.npts) < 1e-9)
    ratio_full = np.max(single.row_sw ** 2) / np.min(single.row_sw ** 2)
    ratio_half = np.max(half.row_sw ** 2) / np.min(half.row_sw ** 2)
    check("tempered weights have a smaller dynamic range", ratio_half < ratio_full and abs(ratio_half - np.sqrt(ratio_full)) < 1e-6 * ratio_full,
          f"{ratio_half} vs sqrt({ratio_full})")
    half.close()
    sc_w = [single.score(x) for x in (x0, x1)]
    check("weighted objective() == score().objective", all(abs(single.objective(x) - s.objective) < 1e-9 * abs(s.objective) for x, s in zip((x0, x1), sc_w)))
    check("plain chi2 reported unchanged", all(abs(a.chi2 - b.chi2) < 1e-7 * a.chi2 for a, b in zip(sc_plain, sc_w)),
          "; ".join(f"{a.chi2} vs {b.chi2}" for a, b in zip(sc_plain, sc_w)))
    r, J = single.residuals(x1, jac=True)
    check("rows_objective(residuals) == objective", abs(single.rows_objective(r) - single.objective(x1)) < 1e-9 * single.objective(x1))
    check("objective differs from chi2 (weights and loss act)", abs(sc_w[0].objective - sc_plain[0].objective) > 1e-3)
    check("score carries the objective kind", "reduced" in sc_w[0].objective_kind and "sivia" in sc_w[0].objective_kind)
    # at the reference vector the standardized residuals are O(1): the scaled loss
    # keeps most of the chi-squared (unscaled Sivia would compress it ~10x)
    # c rho(z/c) increases with c (rho concave, rho(0) = 0) and tends to z: with
    # scales >= 1 the objective is ordered  unscaled Sivia <= scaled Sivia <= weighted chi2
    w_chi2 = float(np.sum((single.residuals(x0)[0][:single.npts]) ** 2)) + sc_w[0].penalty
    spec_none = ObjectiveSpec(dataset_weight="reduced", dataset_key="file", loss="sivia", loss_scale="none")
    single_none = SingleSession(single_variant(base, "unscaled", tmp), cwd=tmp, name="unscaled", objective=spec_none)
    o_none = single_none.objective(x0)
    check("unscaled Sivia <= scaled Sivia <= weighted chi2", o_none <= sc_w[0].objective <= w_chi2,
          f"{o_none} <= {sc_w[0].objective} <= {w_chi2}")
    single_none.close()
    spec_glob = ObjectiveSpec(dataset_weight="none", loss="sivia", loss_scale="global",
                              scales=dict(scales, **{"*": sc_ref.chi2 / ref.npts}))
    single_glob = SingleSession(single_variant(base, "global", tmp), cwd=tmp, name="global", objective=spec_glob)
    check("global scale ignores the per-file entries", np.allclose(single_glob.row_c, sc_ref.chi2 / ref.npts))
    single_glob.close()

    print("3. sharded == single with the weighted objective")
    sh = Sharded(base, os.path.join(tmp, "shards"), nshards=3, flags={}, objective=spec)
    for name, x in (("x0", x0), ("x1", x1)):
        o1, o2 = single.objective(x), sh.objective(x)
        check(f"objective agrees at {name}", abs(o1 - o2) <= 1e-7 * max(1.0, abs(o1)), f"{o1} vs {o2}")
    r1, J1 = single.residuals(x1, jac=True)
    r2, J2 = sh.residuals(x1, jac=True)
    # the sharded rows are in shard-packing order: compare per base segment, a split segment's pieces joined in
    # energy order (2026-10-26: Jacobian-cost balancing splits segments more finely)
    seg1 = {k: (a, b) for k, (a, b) in single.row_ranges}
    seg2 = {}
    for entries in sh.layout:
        for p, _, (ga, gb), _, _ in entries:
            seg2.setdefault(int(p.base_key), []).append((p.e_lo if p.e_lo is not None else -np.inf, ga, gb))
    ok_r = ok_J = True
    for k, (a, b) in seg1.items():
        idx = np.concatenate([np.arange(ga, gb) for _, ga, gb in sorted(seg2[int(k)])])
        ok_r &= idx.size == b - a and np.allclose(r1[a:b], r2[idx], rtol=1e-7, atol=1e-9)
        ok_J &= idx.size == b - a and np.allclose(J1[a:b], J2[idx], rtol=1e-6, atol=1e-8)
    check("weighted residual rows agree per segment", ok_r)
    check("weighted Jacobians agree per segment", ok_J)
    check("penalty rows agree", np.allclose(r1[single.npts:], r2[sh.npts:], rtol=1e-9, atol=1e-12))

    print("4. polish with the Sivia loss")
    b = default_bounds(sh, x1, policy)
    f = polish(sh, x1, bounds=b, policy=policy, max_nfev=60)
    check("polish does not worsen the objective", f.objective <= sh.objective(x1) + 1e-9, f"{f.objective} vs {sh.objective(x1)}")
    check("polish result objective == evaluator objective", abs(f.objective - sh.objective(f.x)) < 1e-7 * max(1, f.objective),
          f"{f.objective} vs {sh.objective(f.x)}")
    check("polish improves from a 5% jitter", f.objective < sh.objective(x1) - 1e-3)
    sc_v, _ = verify(base, f.x, {}, workdir=tmp, name="v", objective=spec)
    check("verify with the spec reproduces the sharded objective", abs(sc_v.objective - f.objective) < 1e-6 * f.objective,
          f"{sc_v.objective} vs {f.objective}")
    sh.close(); single.close()
    shutil.rmtree(tmp, ignore_errors=True)
    print("FAILURES: " + ", ".join(failures) if failures else "ALL OK")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
