#!/usr/bin/env python3
"""The joint (mixture) fit on two unrelated test models sharing one data set.

Members: A = ~/AZURE2/tests/13N (12C+p), B = ~/AZURE2/tests/identical_pp_res
(p+p).  Both get the same c.m.-differential 90-degree elastic segment on
1-4 MeV (lab), whose data are made from the models themselves:
y = 0.6 f_A + 0.4 f_B with 2 % noise.  Checks:

  1. the spec validation refuses a shared segment that is not the same data;
  2. with weights (1, 0) the shared chi-squared equals member A's own;
  3. the joint Jacobian agrees with central finite differences;
  4. chi2 + penalty is the objective;
  5. from perturbed parameters and weights (0.5, 0.5) the polish recovers
     chi2/N ~ 1 on the shared set and the weights 0.6 / 0.4;
  6. freeze_and_exchange writes standalone files that reproduce the joint
     fit's shared chi-squared;
  7. member-tagged keys print and compare as expected.

Run from anywhere (takes a few minutes; the 13N model has 16 data segments):
    python3 rmfit/tests/test_mixture_13n_pp.py
"""
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
import numpy as np                                     # noqa: E402
from rmfit import AZURE2_ROOT                          # noqa: E402
from rmfit.keys import ParamKey                        # noqa: E402

fails = []


def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + ("" if cond else f"  -- {detail}"))
    if not cond:
        fails.append(name)


def prepare(work):
    """Member variants with the shared segment, and the shared data file."""
    from pyazr.azrfile import AzrModel
    from pyazr import azure2
    from rmfit.mixture import _run_pyazr
    from rmfit.model import block_lines
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(os.path.join(work, "shared"))
    e_lab = np.round(np.arange(1.0, 4.0 + 1e-9, 0.1), 6)
    shared = os.path.join(work, "shared", "shared.dat")
    with open(shared, "w") as f:
        for e in e_lab:
            f.write(f"{e:.4f} 90.0 1.0 1.0\n")
    srcs = {"A": os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"),
            "B": os.path.join(AZURE2_ROOT, "tests", "identical_pp_res", "identical_pp_res.azr")}
    azrs, keys = {}, {}
    for name, src in srcs.items():
        d = os.path.join(work, name)
        os.makedirs(os.path.join(d, "output"))
        os.symlink(os.path.join(os.path.dirname(src), "data"), os.path.join(d, "data"))
        m = AzrModel.from_file(src)
        m.set_output_dir(os.path.join(d, "output"))
        m.add_data_segment("../shared/shared.dat", entrance=1, exit=1, observable="differential-cm",
                           energy_min=1.0, energy_max=4.0, angle_min=90.0, angle_max=90.0,
                           norm=1.0, vary_norm=True, norm_error=5.0)
        azrs[name] = os.path.join(d, f"{name}.azr")
        m.write(azrs[name])
        keys[name] = sum(1 for ln in block_lines(open(azrs[name]).read(), "segmentsData") if ln.strip())
    # the members' models on the grid, each in its own process
    f = {}
    for name in srcs:
        res = _run_pyazr(f"""
with azure2({azrs[name]!r}, cwd={os.path.dirname(azrs[name])!r}, use_brune=True) as s:
    vals = s.calculate_rwa(s.params_rwa)
    i = [int(d.key) for d in s.active_datasets].index({keys[name]})
    result['f'] = np.asarray(vals[i], float).tolist()
""", cwd=os.path.dirname(azrs[name]))
        f[name] = np.asarray(res["f"], float)
    assert f["A"].size == e_lab.size and f["B"].size == e_lab.size, (f["A"].size, f["B"].size, e_lab.size)
    rng = np.random.default_rng(13)
    y = 0.6 * f["A"] + 0.4 * f["B"]
    err = 0.02 * y
    y = y * (1.0 + 0.02 * rng.standard_normal(y.size))
    with open(shared, "w") as fh:
        for e, v, s_ in zip(e_lab, y, err):
            fh.write(f"{e:.4f} 90.0 {v:.8e} {s_:.8e}\n")
    return azrs, keys, f, y, err


def write_spec(work, azrs, keys, wA=0.6, wB=0.4, free=False, sigma=0.0, sum_to_one=True):
    cfg = dict(members=[dict(name="A", azr=os.path.relpath(azrs["A"], work), flags=dict(use_brune=True),
                             weight=wA, weight_free=free, weight_sigma=sigma),
                        dict(name="B", azr=os.path.relpath(azrs["B"], work), flags=dict(use_brune=True),
                             weight=wB, weight_free=free, weight_sigma=sigma)],
               shared=[dict(name="shared", segments={"A": keys["A"], "B": keys["B"]})],
               weights=dict(sum_to_one=sum_to_one, sum_sigma=1e-3))
    with open(os.path.join(work, "mixture.json"), "w") as fh:
        json.dump(cfg, fh, indent=1)


def main():
    from rmfit.mixture import MixtureSpec, Mixture, freeze_and_exchange, shared_gid
    from rmfit.engine import SingleSession
    from rmfit.local import polish, default_bounds, BoundsPolicy
    work = os.path.join(HERE, "_mixture_work")
    azrs, keys, f, y, err = prepare(work)
    print(f"shared data: {y.size} points, f_A {f['A'].mean():.4g} f_B {f['B'].mean():.4g} (segment keys {keys})")

    # ---- 1. validation ---------------------------------------------------------
    print("\n1. validation")
    write_spec(work, azrs, keys)
    spec = MixtureSpec.load(work)
    check("a consistent spec loads", len(spec.members) == 2 and len(spec.shared) == 1)
    bad = json.load(open(os.path.join(work, "mixture.json")))
    bad["shared"][0]["segments"]["B"] = 1          # B's first segment: a different data file
    with open(os.path.join(work, "mixture_bad.json"), "w") as fh:
        json.dump(bad, fh)
    try:
        MixtureSpec.load(os.path.join(work, "mixture_bad.json"))
        check("a shared segment on different data is refused", False, "no error")
    except ValueError as e:
        check("a shared segment on different data is refused", "disagree" in str(e), str(e))

    # ---- 2. weights (1, 0): the shared chi2 is member A's own -----------------
    print("\n2. weights (1,0)")
    write_spec(work, azrs, keys, wA=1.0, wB=0.0)
    with Mixture(MixtureSpec.load(work), log=print) as ev:
        print(ev.describe())
        sc = ev.score(ev.x0)
        gid = shared_gid(0)
        joint_shared = sc.segment_chi2[gid]
        nA = ev.nR[0]
        check("global R block = A then B", ev.n_rmatrix == ev.nR[0] + ev.nR[1] and ev.keys[0].member == "A"
              and ev.keys[nA].member == "B", f"{ev.keys[0]} {ev.keys[nA]}")
        with SingleSession(azrs["A"], flags=dict(use_brune=True)) as single:
            scA = single.score(single.x0)
            own = scA.segment_chi2[keys["A"]]
        check("shared chi2 == A's own at the seed", abs(joint_shared - own) <= 1e-8 * max(1.0, own),
              f"{joint_shared:.10g} vs {own:.10g}")
        # the private rows of A are A's own rows too
        privA = sum(v for k, v in sc.segment_chi2.items() if 10000 <= k < 20000)
        privA_single = sum(v for k, v in scA.segment_chi2.items() if k != keys["A"])
        # 1e-5: the two sessions differ in their external-capture cache (one reuses
        # output/intEC.dat written by the data-making run, the variant computes
        # its own), which is a few 1e-6 of the 13N chi-squared
        check("A's private chi2 unchanged (1e-5)", abs(privA - privA_single) <= 1e-5 * max(1.0, privA_single),
              f"{privA:.10g} vs {privA_single:.10g}")

    # ---- 3./4. Jacobian and objective with free weights -------------------------
    print("\n3. Jacobian")
    write_spec(work, azrs, keys, wA=0.6, wB=0.4, free=True, sigma=0.2)
    with Mixture(MixtureSpec.load(work), log=print) as ev:
        x = ev.x0.copy()
        r, J = ev.residuals(x, jac=True)
        check("Jacobian shape", J.shape == (r.size, x.size), str(J.shape))
        wcols = [c for c in ev.weight_col if c is not None]
        ncol = ev.shared_norm_col[0]
        check("weights and the shared norm are free columns", len(wcols) == 2 and ncol is not None)
        # columns: two R-matrix of A, two of B, both weights, the shared norm
        test_cols = [0, 1, ev.nR[0], ev.nR[0] + 1] + wcols + [ncol]
        worst = 0.0
        for c in test_cols:
            h = 1e-5 * max(1.0, abs(x[c]))
            xp, xm = x.copy(), x.copy()
            xp[c] += h; xm[c] -= h
            rp, _ = ev.residuals(xp, jac=False)
            rm, _ = ev.residuals(xm, jac=False)
            fd = (rp - rm) / (2 * h)
            scale = max(np.max(np.abs(J[:, c])), 1e-12)
            worst = max(worst, np.max(np.abs(fd - J[:, c])) / scale)
        check("Jacobian columns agree with central differences (1e-4)", worst < 1e-4, f"worst rel {worst:.2e}")
        sc = ev.score(x)
        check("objective = chi2 + penalty", abs(sc.objective - (sc.chi2 + sc.penalty)) < 1e-9,
              f"{sc.objective} vs {sc.chi2 + sc.penalty}")
        check("sum of squares of the rows = objective", abs(float(np.sum(r ** 2)) - sc.objective) < 1e-9)
        # the sum-to-one row: at x0 the weights sum to 1 exactly
        check("sum-to-one row is zero at the seed", abs(ev.penalties.residuals(x)[-1]) < 1e-12)

        # ---- 5. recovery ---------------------------------------------------------
        print("\n5. recovery from perturbed parameters and weights (0.5, 0.5)")
        x1 = x.copy()
        rng = np.random.default_rng(7)
        for k in range(2):
            lo = sum(ev.nR[:k])
            widths = [lo + j for j in range(ev.nR[k]) if ev.keys[lo + j].kind == "width"][:3]
            for j in widths:
                x1[j] *= 1.0 + 0.3 * rng.choice([-1, 1])
            energies = [lo + j for j in range(ev.nR[k]) if ev.keys[lo + j].kind == "energy"]
            if energies:
                x1[energies[0]] += 0.020
        for c in wcols:
            x1[c] = 0.5
        sc1 = ev.score(x1)
        n_sh = y.size
        print(f"  perturbed: objective {sc1.objective:.2f}, shared chi2/N {sc1.segment_chi2[gid] / n_sh:.2f}")
        policy = BoundsPolicy()
        bounds = default_bounds(ev, x1, policy)
        check("weights are bounded in [0, 1]", all(bounds.lb[c] == 0.0 and bounds.ub[c] >= 1.0 for c in wcols))
        res = polish(ev, x1, bounds=bounds, policy=policy, max_nfev=60, verbose=1)
        sc2 = ev.score(res.x)
        w = ev.weights(res.x)
        print(f"  polished: objective {sc2.objective:.2f} ({res.status}), shared chi2/N {sc2.segment_chi2[gid] / n_sh:.2f}, "
              f"weights {w}, norm {ev.shared_norm(res.x, 0):.4f}")
        check("shared chi2/N recovered to < 1.6", sc2.segment_chi2[gid] / n_sh < 1.6, f"{sc2.segment_chi2[gid] / n_sh:.2f}")
        check("weights recovered within 0.06", abs(w["A"] - 0.6) < 0.06 and abs(w["B"] - 0.4) < 0.06, str(w))
        check("weights still sum to one (1e-2)", abs(w["A"] + w["B"] - 1.0) < 1e-2, str(w))

        # ---- 6. freeze and exchange ----------------------------------------------
        print("\n6. freeze and exchange")
        frozen = os.path.join(work, "frozen")
        rep = freeze_and_exchange(ev, res.x, frozen, log=print)
        check("two standalone files written", len(rep["members"]) == 2 and all(os.path.exists(v["azr"]) for v in rep["members"].values()))
        check("two tables written", len(rep["tables"]) == 2)
        check("every standalone reproduces the joint shared chi2", bool(rep["verify"]) and all(v["ok"] for v in rep["verify"].values()),
              json.dumps(rep["verify"]))
        txt = open(rep["members"]["A"]["azr"]).read()
        check("A's shared line carries the mixture block with B's table",
              "mixture" in txt and "shared_B.dat" in txt and f"{w['A']:.6g}"[:5] in txt)

    # ---- 7. keys -----------------------------------------------------------------
    print("\n7. keys")
    k = ParamKey("width", level="22Ne:1/2+#3", pair=1, L=0, S=0.5, member="22Ne")
    check("member-tagged key prints its member", str(k) == "g[22Ne:1/2+#3 p1 L0 S0.5]", str(k))
    check("jpi strips the member prefix", k.jpi == "1/2+" and k.ordinal == 3, f"{k.jpi} {k.ordinal}")
    check("untagged keys unchanged", str(ParamKey("norm", seg=4)) == "norm[seg4]"
          and ParamKey("norm", seg=4) == ParamKey("norm", seg=4, member=""))
    check("a weight key prints", str(ParamKey("weight", level="22Ne")) == "w[22Ne]")

    print()
    if fails:
        print(f"FAILED: {len(fails)} check(s): {', '.join(fails)}")
        return 1
    print("all mixture checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
