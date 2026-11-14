#!/usr/bin/env python3
"""No-U-turn sampler (multinomial, Betancourt 2017, arXiv:1701.02434; Stan's extra U-turn checks across subtrees)
with a dense or diagonal metric, dual-averaging step size (Hoffman & Gelman 2014) and Stan's windowed warm-up, and
one chain driver (`run_chain`).

Generic: it needs only ``logp_grad(z) -> (logp, grad)``.  A point where the model fails must return
``(-inf, anything)``: the trajectory then counts as divergent and the sampler stays where it was.  Ported from the
user-level skill rmatrix-posterior (nuts.py and the warm-up of chain.py; the checkpointing is left out).
sample_nuts.py runs it on an AZURE2 objective in coordinates whitened at the MAP; selftest.py checks it on an
ill-conditioned Gaussian and a hard wall.

    python3 nuts.py --help     (a module: run sample_nuts.py)
"""
import argparse
import sys

import numpy as np


# ------------------------------------------------------------------ the tree
class _Leaf:
    __slots__ = ("z", "p", "g", "lp", "ps")

    def __init__(self, z, p, g, lp, ps):
        self.z, self.p, self.g, self.lp, self.ps = z, p, g, lp, ps


class _Tree:
    __slots__ = ("left", "right", "prop", "logw", "rho", "n", "acc", "turn", "div")


def _no_turn(ps_left, ps_right, rho):
    return float(ps_left @ rho) > 0.0 and float(ps_right @ rho) > 0.0


class NUTS:
    def __init__(self, logp_grad, dim, max_depth=10, target_accept=0.85,
                 max_energy_error=1000.0):
        self.f = logp_grad
        self.dim = dim
        self.max_depth = max_depth
        self.delta = target_accept
        self.dmax = max_energy_error

    # -- metric: M^-1 = Minv; momentum p ~ N(0, M), M = Minv^-1 ---------------
    def _set_metric(self, Minv):
        self.Minv = np.asarray(Minv, float)
        self.Mchol = np.linalg.cholesky(np.linalg.inv(self.Minv))

    def _kin(self, p):
        return 0.5 * float(p @ (self.Minv @ p))

    def _leapfrog(self, leaf, eps):
        p = leaf.p + 0.5 * eps * leaf.g
        z = leaf.z + eps * (self.Minv @ p)
        lp, g = self.f(z)
        self.ngrad += 1
        if not np.isfinite(lp):
            return _Leaf(z, p, g, -np.inf, p)
        p = p + 0.5 * eps * g
        return _Leaf(z, p, g, lp, self.Minv @ p)

    def _build(self, leaf, v, depth, eps, H0, rng):
        if depth == 0:
            new = self._leapfrog(leaf, v * eps)
            t = _Tree()
            H = -new.lp + self._kin(new.p) if np.isfinite(new.lp) else np.inf
            dH = H - H0
            t.left = t.right = t.prop = new
            t.logw = -dH if np.isfinite(dH) else -np.inf
            t.rho = new.p.copy()
            t.n = 1
            t.acc = float(np.exp(min(-dH, 0.0))) if np.isfinite(dH) else 0.0
            t.div = (not np.isfinite(dH)) or dH > self.dmax
            t.turn = False
            return t
        a = self._build(leaf, v, depth - 1, eps, H0, rng)
        if a.div or a.turn:
            return a
        b = self._build(a.right if v > 0 else a.left, v, depth - 1, eps, H0, rng)
        t = _Tree()
        t.n, t.acc, t.div = a.n + b.n, a.acc + b.acc, b.div
        if b.div or b.turn:
            t.turn, t.left, t.right, t.prop = b.turn, a.left, a.right, a.prop
            t.logw, t.rho = a.logw, a.rho
            return t
        t.logw = np.logaddexp(a.logw, b.logw)
        t.prop = b.prop if np.log(rng.random()) < b.logw - t.logw else a.prop
        first, second = (a, b) if v > 0 else (b, a)          # time order
        t.left, t.right = first.left, second.right
        t.rho = first.rho + second.rho
        t.turn = not (_no_turn(t.left.ps, t.right.ps, t.rho)
                      and _no_turn(t.left.ps, second.left.ps, first.rho + second.left.p)
                      and _no_turn(first.right.ps, t.right.ps, second.rho + first.right.p))
        return t

    def transition(self, z, lp, g, eps, rng):
        """One NUTS transition from (z, lp, g).  Returns the new point and stats."""
        self.ngrad = 0
        p = self.Mchol @ rng.standard_normal(self.dim)
        start = _Leaf(z, p, g, lp, self.Minv @ p)
        H0 = -lp + self._kin(p)
        tree = _Tree()
        tree.left = tree.right = tree.prop = start
        tree.logw, tree.rho, tree.n, tree.acc, tree.turn, tree.div = 0.0, p.copy(), 0, 0.0, False, False
        depth, divergent = 0, False
        while depth < self.max_depth:
            v = 1 if rng.random() < 0.5 else -1
            sub = self._build(tree.right if v > 0 else tree.left, v, depth, eps, H0, rng)
            tree.n += sub.n
            tree.acc += sub.acc
            depth += 1
            if sub.div:
                divergent = True
                break
            if sub.turn:
                break
            if np.log(rng.random()) < sub.logw - tree.logw:        # biased progressive
                tree.prop = sub.prop
            first, second = (tree, sub) if v > 0 else (sub, tree)
            left, right = first.left, second.right
            rho = first.rho + second.rho
            turn = not (_no_turn(left.ps, right.ps, rho)
                        and _no_turn(left.ps, second.left.ps, first.rho + second.left.p)
                        and _no_turn(first.right.ps, right.ps, second.rho + first.right.p))
            tree.left, tree.right, tree.rho = left, right, rho
            tree.logw = np.logaddexp(tree.logw, sub.logw)
            if turn:
                break
        q = tree.prop
        stats = dict(accept=tree.acc / max(tree.n, 1), depth=depth, n_leapfrog=tree.n,
                     divergent=divergent, energy=-q.lp + self._kin(q.p), ngrad=self.ngrad)
        return q.z, q.lp, q.g, stats


# ------------------------------------------------------------------ warm-up
def warmup_schedule(n_warmup, init=75, term=50, base=25):
    """Stan's windows: (start, end) of every slow window that ends with a metric update."""
    wins, start, size = [], init, base
    end_slow = n_warmup - term
    while start < end_slow:
        end = start + size
        if end + 2 * size > end_slow:          # last window absorbs the remainder
            end = end_slow
        wins.append((start, end))
        start, size = end, 2 * size
    return wins


class DualAveraging:
    def __init__(self, eps0, delta, gamma=0.05, t0=10.0, kappa=0.75):
        self.mu, self.delta, self.gamma, self.t0, self.kappa = np.log(10 * eps0), delta, gamma, t0, kappa
        self.t, self.Hbar, self.logeps_bar = 0, 0.0, 0.0

    def update(self, accept):
        self.t += 1
        w = 1.0 / (self.t + self.t0)
        self.Hbar = (1 - w) * self.Hbar + w * (self.delta - accept)
        logeps = self.mu - np.sqrt(self.t) / self.gamma * self.Hbar
        eta = self.t ** -self.kappa
        self.logeps_bar = eta * logeps + (1 - eta) * self.logeps_bar
        return float(np.exp(logeps))

    def restart(self, eps):
        self.mu, self.t, self.Hbar, self.logeps_bar = np.log(10 * eps), 0, 0.0, 0.0

    def final(self):
        return float(np.exp(self.logeps_bar))

    def state(self):
        return dict(mu=self.mu, t=self.t, Hbar=self.Hbar, logeps_bar=self.logeps_bar)

    def load(self, s):
        self.mu, self.t, self.Hbar, self.logeps_bar = s["mu"], s["t"], s["Hbar"], s["logeps_bar"]


def regularized_covariance(samples, metric="auto", shrink_to=1e-3):
    """Metric estimate (M^-1) from one warm-up window.

    dense: Stan's estimate, shrunk towards shrink_to * I (Stan: 1e-3 for raw coordinates;
    1 in coordinates already whitened by the MAP Hessian, which chain.py passes).  It needs many more draws than
    dimensions: with 200 draws in 202 dimensions it is rank-deficient and its smallest
    eigenvalues collapse to the 1e-3 floor (cond ~3e5 against the true covariance), and
    NUTS then needs ~60x more gradients per draw.  diag: variances only, same shrinkage.
    identity: keep the metric the whitening gave.  auto: dense when the window holds at
    least 10 draws per dimension, diag otherwise.
    """
    n, d = samples.shape
    if metric == "auto":
        metric = "dense" if n >= 10 * d else "diag"
    if metric == "identity":
        return np.eye(d)
    if metric == "diag":
        return np.diag((n / (n + 5.0)) * np.var(samples, axis=0, ddof=1) + shrink_to * (5.0 / (n + 5.0)))
    C = np.cov(samples, rowvar=False)
    return (n / (n + 5.0)) * C + shrink_to * (5.0 / (n + 5.0)) * np.eye(d)


def run_chain(logp_grad, dim, warmup, draws, seed, init=None, init_scale=1.0, eps0=0.5, max_depth=10,
              target_accept=0.85, metric="auto", shrink_to=1.0, log=None):
    """One chain: warm-up (dual averaging; Stan's windows update the metric) then `draws` draws.  Returns a dict of
    arrays over the kept draws: z, logp, accept, depth, n_leapfrog, divergent, ngrad, and the final step size
    "eps" and metric "Minv".  shrink_to = 1 suits coordinates already whitened by the MAP Hessian."""
    rng = np.random.default_rng(seed)
    sampler = NUTS(logp_grad, dim, max_depth=max_depth, target_accept=target_accept)
    z = np.asarray(init, float).copy() if init is not None else init_scale * rng.standard_normal(dim)
    Minv, eps, buf = np.eye(dim), eps0, []
    sampler._set_metric(Minv)
    da = DualAveraging(eps, target_accept)
    wins = warmup_schedule(warmup)
    win_ends = {e for _, e in wins}
    lp, g = logp_grad(z)
    if not np.isfinite(lp):
        raise ValueError("the target fails at the chain's starting point; start nearer the MAP")
    keep = {k: [] for k in ("z", "logp", "accept", "depth", "n_leapfrog", "divergent", "ngrad")}
    for it in range(warmup + draws):
        z, lp, g, st = sampler.transition(z, lp, g, eps, rng)
        if it < warmup:
            eps = da.update(st["accept"])
            if wins and wins[0][0] <= it < wins[-1][1]:
                buf.append(z.copy())
            if it + 1 in win_ends:
                Minv = regularized_covariance(np.array(buf), metric, shrink_to)
                sampler._set_metric(Minv)
                buf = []
                da.restart(eps)
            if it + 1 == warmup:
                eps = da.final()
            continue
        keep["z"].append(z.copy())
        keep["logp"].append(lp)
        for k in ("accept", "depth", "n_leapfrog", "divergent", "ngrad"):
            keep[k].append(st[k])
        if log and (it + 1 - warmup) % max(1, draws // 5) == 0:
            log(f"  draw {it + 1 - warmup}/{draws}  log p {lp:.3f}  eps {eps:.3f}  depth {st['depth']}")
    out = {k: np.array(v) for k, v in keep.items()}
    out["eps"], out["Minv"] = eps, Minv
    return out


if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__.split("\n\n")[0]).parse_args()
    print("nuts.py is a module; sample_nuts.py runs it on a project.")
    sys.exit(0)
