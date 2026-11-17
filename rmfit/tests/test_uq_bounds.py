#!/usr/bin/env python3
"""rmfit.uq BoundedTarget (2026-11-12): the tanh map, its Jacobian and gradient, NUTS on a truncated normal, checkpoint
back-mapping, and _box from a BoundsPolicy.  Login node, seconds:  python3 rmfit/tests/test_uq_bounds.py"""
import os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
import numpy as np
from rmfit.uq import BoundedTarget, _box, _chains_from_checkpoints
from rmfit.local import BoundsPolicy

FAILS = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""))
    if not ok:
        FAILS.append(name)


class Gauss:
    """log p(x) = -0.5 (x - mu)^T A (x - mu); factory protocol like AzurePosterior."""
    def __init__(self, mu, A):
        self.mu, self.A = np.asarray(mu, float), np.asarray(A, float)
    def __call__(self, x=None):
        if x is None:
            return self
        d = np.asarray(x, float) - self.mu
        return -0.5 * float(d @ self.A @ d), -self.A @ d


def main():
    rng = np.random.default_rng(3)
    A = np.array([[2.0, 0.3, 0.0], [0.3, 1.0, -0.2], [0.0, -0.2, 4.0]])
    T = Gauss([0.2, -0.1, 0.5], A)
    lo, hi = np.array([-1.0, -np.inf, -0.8]), np.array([1.0, np.inf, 0.8])
    B = BoundedTarget(T, lo, hi)
    check("bounded mask", B.b.tolist() == [True, False, True])
    x = np.array([0.3, 1.7, -0.5])
    u = B.to_u(x)
    check("to_x(to_u(x)) = x", np.allclose(B.to_x(u), x, atol=1e-12), str(B.to_x(u)))
    check("unbounded entry untouched", u[1] == x[1])
    # gradient in u by central differences
    lp, g = B(u)
    num = np.zeros(3)
    for i in range(3):
        e = np.zeros(3); e[i] = 1e-6
        num[i] = (B(u + e)[0] - B(u - e)[0]) / 2e-6
    check("gradient in u matches finite differences", np.allclose(g, num, rtol=1e-6, atol=1e-8), f"{g} vs {num}")
    # log p(u) = log p(x) + log |dx/du|
    d = B.dxdu(u)
    check("Jacobian term", np.isclose(lp, T(x)[0] + np.log(d[0]) + np.log(d[2])))
    # the box is respected for any u, and the edges map to +-large u
    check("to_x inside the box (the edge only where tanh saturates)", np.all(np.abs(B.to_x(rng.normal(size=(50, 3)) * 10)[:, [0, 2]]) <= np.array([1.0, 0.8])))
    lp_far, g_far = B(np.array([40.0, 0.0, -40.0]))
    check("saturated u: finite log p and gradient", np.isfinite(lp_far) and np.all(np.isfinite(g_far)), f"{lp_far} {g_far}")
    # 2-D arrays (chains x draws)
    U = rng.normal(size=(2, 5, 3))
    check("to_x on (chains, draws, dim)", np.allclose(B.to_x(U)[:, :, 1], U[:, :, 1]) and np.all(np.abs(B.to_x(U)[:, :, 0]) < 1))
    # chain starts: clipped inside the box, unbounded entries untouched and finite (a queue job: nan from inf arithmetic)
    st = B.inside(np.array([[1.3, 5.0, -0.79], [0.0, -2.0, 0.0]]))
    check("inside(): bounded entries clipped, unbounded untouched", np.all(np.isfinite(st)) and st[0, 0] < 1.0 and st[0, 1] == 5.0 and st[1, 1] == -2.0
          and st[0, 2] > -0.8 and np.all(np.isfinite(B(B.to_u(st[0]))[0:1])), str(st))
    # metric in u: inverse of H_u = D A D + diag(2 sech^2 u) on the bounded entries (no gradient term without g0)
    C = np.linalg.inv(A)
    Cu = B.metric_u(C, u, sd_max=1e9)
    Hu = A * np.outer(d, d)
    t = np.tanh(u[B.b]); Hu[np.flatnonzero(B.b), np.flatnonzero(B.b)] += 2 * (1 - t * t)
    check("metric_u = inv(D A D + Jacobian curvature)", np.allclose(Cu, np.linalg.inv(Hu)), str(Cu))
    # the background-pole pathology: a mode held at the box edge along a direction the unbounded Laplace leaves wide.
    # log p = -0.5 ((x1 - 3 x2)/0.2)^2 - 0.5 (x2/100)^2 + 5 x2 with x2 in [-0.8, 0.8]: mode x2 = 0.8 (pushed by the linear
    # term), Laplace sd(x2) = 100, sd(x1) = 300 for a box of half-width 0.8 (the ridge is curved in u, so the chains are
    # longer than the truncated-normal check's; the old metric's chains are short: they do not move).  Exact: p(x2) on the box ~ exp(5 x2).
    class Ridge:
        def __call__(self, x=None):
            if x is None:
                return self
            x1, x2 = float(x[0]), float(x[1])
            r = (x1 - 3 * x2) / 0.2
            lp = -0.5 * r * r - 0.5 * (x2 / 100) ** 2 + 5 * x2
            return lp, np.array([-r / 0.2, 3 * r / 0.2 - x2 / 100 ** 2 + 5])
    R = Ridge()
    BR = BoundedTarget(R, [-np.inf, -0.8], [np.inf, 0.8])
    A_r = np.array([[1 / 0.2 ** 2, -3 / 0.2 ** 2], [-3 / 0.2 ** 2, 9 / 0.2 ** 2 + 1 / 100 ** 2]])
    C_r = np.linalg.inv(A_r)
    x_mode = np.array([2.4, 0.8])
    u_mode = BR.to_u(BR.inside(x_mode))
    g_mode = R(BR.to_x(u_mode))[1]
    Cu_old = C_r / np.outer(BR.dxdu(u_mode), BR.dxdu(u_mode))
    Cu_new = BR.metric_u(C_r, u_mode, g0=g_mode)
    check("edge mode: the old C / dd^T metric is huge in u", np.sqrt(Cu_old[1, 1]) > 1e3, f"{np.sqrt(Cu_old[1, 1]):.3g}")
    check("edge mode: the Hessian metric's bounded sd <= 1 in u", np.sqrt(Cu_new[1, 1]) <= 1.0 + 1e-9, f"{np.sqrt(Cu_new[1, 1]):.3g}")
    check("edge mode: the unbounded coordinate's sd in u follows the ridge, not the Laplace", np.sqrt(Cu_new[0, 0]) < 10, f"{np.sqrt(Cu_new[0, 0]):.3g} (Laplace {np.sqrt(C_r[0, 0]):.3g})")
    from rmfit.sampling import nuts, split_rhat
    xg = np.linspace(-0.8, 0.8, 20001); wg = np.exp(5 * xg)
    m_exact = float((xg * wg).sum() / wg.sum())
    res = []
    for lab, M in (("old", Cu_old), ("new", Cu_new)):
        S = []
        for c in range(2):
            nw, nd = (100, 200) if lab == "old" else (500, 2000)
            rr = nuts(BR, u_mode + 0.1 * np.sqrt(np.diag(M)) * np.random.default_rng(c).normal(size=2), n_draws=nd, n_warmup=nw, inv_metric=M, target_accept=0.9, max_tree_depth=8, seed=10 + c)
            S.append(BR.to_x(rr.samples[0]))
        S = np.stack(S)
        rh = float(np.max(split_rhat(S)))
        res.append((lab, rh, float(S[:, :, 1].mean()), float(rr.step_size[0]) if np.ndim(rr.step_size) else float(rr.step_size)))
        print(f"   {lab} metric: R-hat {rh:.3f}, mean x2 {S[:, :, 1].mean():.4f} (exact {m_exact:.4f}), step size {res[-1][3]:.3g}")
    lab, rh, m2, eps = res[1]
    check("edge mode: NUTS with the Hessian metric mixes (R-hat < 1.05)", rh < 1.05, f"{rh:.3f}")
    check("edge mode: NUTS with the Hessian metric gives the exact mean of x2 within 0.03", abs(m2 - m_exact) < 0.03, f"{m2:.4f} vs {m_exact:.4f}")
    # NUTS on a truncated normal N(0, 1) on [-0.5, 2]: mean and sd against the closed form
    from rmfit.sampling import nuts
    from scipy.stats import truncnorm
    T1 = Gauss([0.0], [[1.0]])
    B1 = BoundedTarget(T1, [-0.5], [2.0])
    r = nuts(B1, B1.to_u(np.array([0.3])), n_draws=3000, n_warmup=500, inv_metric=np.array([[1.0]]), target_accept=0.9, seed=5)
    xs = B1.to_x(r.samples[0])[:, 0]
    tn = truncnorm(-0.5, 2.0)
    check("NUTS truncated normal: mean", abs(xs.mean() - tn.mean()) < 0.06, f"{xs.mean():.3f} vs {tn.mean():.3f}")
    check("NUTS truncated normal: sd", abs(xs.std() - tn.std()) < 0.06, f"{xs.std():.3f} vs {tn.std():.3f}")
    check("NUTS truncated normal: every draw inside", xs.min() > -0.5 and xs.max() < 2.0)
    # checkpoints in u are mapped back to x
    tmp = tempfile.mkdtemp(prefix="rmfit_uqb_")
    for pid in (1, 2):
        q = rng.normal(size=(30, 3)) * 5
        np.savez(os.path.join(tmp, f"ckpt_t_{pid}.npz"), q=q, st=np.zeros((30, 5)), nw=10)
    out = _chains_from_checkpoints(tmp, "t", 10, 20, lo, hi)
    check("checkpoint assembly maps u -> x", out is not None and out["samples"].shape == (2, 20, 3) and np.all(np.abs(out["samples"][:, :, 0]) < 1.0))
    # _box from a policy: physical cap 1 below bg_energy, background cap 3 above; no cap for a non-particle name
    pol = BoundsPolicy(theta2_cap=1.0, bg_theta2_cap=3.0, bg_energy=20.0)
    names = ["energy_1", "width_1_2", "width_2_2", "width_1_5", "segment_1_norm"]
    wig = {"width_1_2": 0.64, "width_2_2": 0.64}
    le = {"energy_1": 16.9, "width_1_2": 16.9, "width_2_2": 26.9, "width_1_5": 16.9}
    logs = []
    lo2, hi2, rows = _box(dict(bounds="policy", theta2_cap=None), pol, names, [0, 1, 2, 3, 4], np.array([16.9, 0.5, -0.9, 0.1, 1.0]), wig, le, logs.append)
    check("_box: physical amplitude at gamma_W", np.isclose(hi2[1], 0.8) and np.isclose(lo2[1], -0.8))
    check("_box: background amplitude at sqrt(3) gamma_W", np.isclose(hi2[2], np.sqrt(3 * 0.64)))
    check("_box: energies, photon widths, norms unbounded", np.isinf(lo2[0]) and np.isinf(hi2[3]) and np.isinf(hi2[4]))
    lo3, hi3, rows3 = _box(dict(bounds="policy", theta2_cap=None), pol, names, [1], np.array([16.9, 0.9, 0, 0, 1]), wig, le, logs.append)
    check("_box: a mode outside its bound widens it with a warning", np.isclose(hi3[0], 0.945) and any("WARNING" in l for l in logs))
    lo4, hi4, _ = _box(dict(bounds="none"), pol, names, [1, 2], np.zeros(5), wig, le, logs.append)
    check("_box: bounds none", np.all(np.isinf(lo4)))
    lo5, hi5, _ = _box(dict(bounds="policy", theta2_cap=4.0), None, names, [1, 2], np.zeros(5), wig, le, logs.append)
    check("_box: theta2_cap without a policy", np.isclose(hi5[0], 1.6) and np.isclose(hi5[1], 1.6))
    if FAILS:
        print(f"FAILED {len(FAILS)}: {', '.join(FAILS)}"); sys.exit(1)
    print("all uq bounds checks passed")


if __name__ == "__main__":
    main()
