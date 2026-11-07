"""Gradient-based posterior sampling: the No-U-Turn sampler (NUTS) for AZURE2 posteriors.

Ensemble samplers (emcee) need no gradient, but their autocorrelation time grows with the number of parameters, and an
R-matrix posterior with tens to hundreds of strongly correlated parameters mixes slowly.  Hamiltonian Monte Carlo moves
along the gradient of the log posterior, so one trajectory crosses the whole posterior; with a mass matrix taken from the
posterior covariance (the Laplace approximation, or a covariance adapted during warm-up) the strong correlations of
level energies and widths stop mattering.  AZURE2 has an analytic gradient of chi-squared
(:meth:`pyazr.azure2.chi2_and_grad`, reverse-mode adjoint, about two forward evaluations for all level parameters), so
a gradient costs a few evaluations, not ``2 n``.

The sampler is the multinomial NUTS of Betancourt (2017) with the additional U-turn checks across sub-trees used by Stan,
dual-averaging step-size adaptation (Hoffman & Gelman 2014) and windowed covariance adaptation of a dense or diagonal
inverse metric.  It works on any function ``theta -> (log p, d log p / d theta)``; :class:`AzurePosterior` builds that
function for an ordinary AZURE2 fit (chi-squared + AZURE2's own penalties + optional Gaussian priors).

References: M. D. Hoffman and A. Gelman, JMLR 15 (2014) 1593; M. Betancourt, arXiv:1701.02434 (2017);
A. Vehtari et al., Bayesian Analysis 16 (2021) 667 (R-hat and effective sample size).
"""
import math
import os
import time
import warnings
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

__all__ = ["nuts", "sample", "NUTSResult", "ess", "split_rhat", "AzurePosterior"]


# ---------------------------------------------------------------------------------------------------------------------
# metric (mass matrix)

class _Metric:
    """Inverse metric M^-1 (the posterior covariance scale): dense (n x n) or diagonal (n,)."""

    def __init__(self, inv_metric, n):
        if inv_metric is None:
            inv_metric = np.ones(n)
        a = np.asarray(inv_metric, float)
        self.dense = a.ndim == 2
        if self.dense:
            if a.shape != (n, n):
                raise ValueError(f"inv_metric has shape {a.shape}, expected ({n}, {n}) or ({n},)")
            a = 0.5 * (a + a.T)
            try:
                self.chol = np.linalg.cholesky(a)
            except np.linalg.LinAlgError as err:
                raise ValueError("inv_metric is not positive definite") from err
        else:
            if a.shape != (n,):
                raise ValueError(f"inv_metric has shape {a.shape}, expected ({n},) or ({n}, {n})")
            if np.any(a <= 0) or not np.all(np.isfinite(a)):
                raise ValueError("diagonal inv_metric must be positive and finite")
            self.sd = np.sqrt(a)
        self.inv = a

    def velocity(self, p):
        return self.inv @ p if self.dense else self.inv * p

    def kinetic(self, p):
        return 0.5 * float(p @ self.velocity(p))

    def draw_momentum(self, rng):
        z = rng.standard_normal(self.inv.shape[0])
        if self.dense:
            # p ~ N(0, M), M = (L L^T)^-1  ->  p = L^-T z
            from scipy.linalg import solve_triangular
            return solve_triangular(self.chol, z, trans="T", lower=True)
        return z / self.sd


# ---------------------------------------------------------------------------------------------------------------------
# NUTS core

class _Point:
    __slots__ = ("q", "p", "g", "lp")

    def __init__(self, q, p, g, lp):
        self.q, self.p, self.g, self.lp = q, p, g, lp


def _evaluate(f, q):
    try:
        lp, g = f(q)
        lp = float(lp)
        g = np.asarray(g, float)
        if not np.isfinite(lp) or not np.all(np.isfinite(g)):
            return -np.inf, np.zeros_like(q)
        return lp, g
    except (ValueError, RuntimeError, FloatingPointError, ArithmeticError):
        return -np.inf, np.zeros_like(q)


def _leapfrog(f, pt, eps, metric, counter):
    p = pt.p + 0.5 * eps * pt.g
    q = pt.q + eps * metric.velocity(p)
    lp, g = _evaluate(f, q)
    counter[0] += 1
    p = p + 0.5 * eps * g
    return _Point(q, p, g, lp)


def _no_uturn(ps_minus, ps_plus, rho):
    return float(ps_plus @ rho) > 0.0 and float(ps_minus @ rho) > 0.0


class _Tree:
    __slots__ = ("beg", "end", "ps_beg", "ps_end", "rho", "prop", "lsw", "n", "acc", "div", "turn")


def _build(f, pt, direction, depth, eps, metric, H0, rng, counter, max_dH):
    if depth == 0:
        new = _leapfrog(f, pt, direction * eps, metric, counter)
        t = _Tree()
        H = -new.lp + metric.kinetic(new.p) if np.isfinite(new.lp) else np.inf
        dH = H - H0
        t.beg = t.end = new
        t.ps_beg = t.ps_end = metric.velocity(new.p)
        t.rho = new.p.copy()
        t.prop = new
        t.lsw = -dH if np.isfinite(dH) else -np.inf
        t.n = 1
        t.acc = min(1.0, math.exp(-dH)) if np.isfinite(dH) and dH > -700 else (1.0 if np.isfinite(dH) else 0.0)
        t.div = not np.isfinite(dH) or dH > max_dH
        t.turn = False
        return t
    left = _build(f, pt, direction, depth - 1, eps, metric, H0, rng, counter, max_dH)
    if left.div or left.turn:
        return left
    right = _build(f, left.end, direction, depth - 1, eps, metric, H0, rng, counter, max_dH)
    t = _Tree()
    t.n = left.n + right.n
    t.acc = left.acc + right.acc
    t.beg, t.ps_beg = left.beg, left.ps_beg
    t.end, t.ps_end = right.end, right.ps_end
    t.div, t.turn = right.div, right.turn
    if right.div or right.turn:
        t.rho, t.prop, t.lsw = left.rho, left.prop, left.lsw
        return t
    t.lsw = np.logaddexp(left.lsw, right.lsw)
    t.prop = right.prop if math.log(rng.random()) < right.lsw - t.lsw else left.prop
    t.rho = left.rho + right.rho
    ok = _no_uturn(left.ps_beg, right.ps_end, t.rho)
    # the extra checks across the two halves (Stan >= 2.26): they catch U-turns the outer ends miss
    ok = ok and _no_uturn(left.ps_beg, right.ps_beg, left.rho + right.beg.p)
    ok = ok and _no_uturn(left.ps_end, right.ps_end, right.rho + left.end.p)
    t.turn = not ok
    return t


def _transition(f, cur, eps, metric, rng, max_depth, max_dH):
    """One NUTS transition from point ``cur`` (q, g, lp; p is drawn here)."""
    counter = [0]
    p0 = metric.draw_momentum(rng)
    start = _Point(cur.q, p0, cur.g, cur.lp)
    H0 = -cur.lp + metric.kinetic(p0)
    minus = plus = start
    ps_minus = ps_plus = metric.velocity(p0)
    rho = p0.copy()
    prop, lsw = start, 0.0
    depth, acc, div = 0, 0.0, False
    for depth in range(max_depth):
        direction = 1 if rng.random() < 0.5 else -1
        edge = plus if direction == 1 else minus
        sub = _build(f, edge, direction, depth, eps, metric, H0, rng, counter, max_dH)
        acc += sub.acc
        if sub.div:
            div = True
            break
        if sub.turn:
            break
        if math.log(rng.random()) < sub.lsw - lsw:          # biased progressive sampling at the top level
            prop = sub.prop
        lsw = np.logaddexp(lsw, sub.lsw)
        if direction == 1:
            ok = (_no_uturn(ps_minus, sub.ps_beg, rho + sub.beg.p)
                  and _no_uturn(ps_plus, sub.ps_end, sub.rho + plus.p))
            plus, ps_plus = sub.end, sub.ps_end
        else:
            ok = (_no_uturn(sub.ps_beg, ps_plus, rho + sub.beg.p)
                  and _no_uturn(sub.ps_end, ps_minus, sub.rho + minus.p))
            minus, ps_minus = sub.end, sub.ps_end
        rho = rho + sub.rho
        ok = ok and _no_uturn(ps_minus, ps_plus, rho)
        if not ok:
            break
    n = counter[0]
    energy = -prop.lp + metric.kinetic(prop.p)
    return _Point(prop.q, None, prop.g, prop.lp), dict(
        accept=acc / max(n, 1), n_leapfrog=n, depth=depth + 1, divergent=div, energy=energy)


def _find_step(f, cur, metric, rng, eps=1.0):
    """Hoffman & Gelman Alg. 4: double/halve eps until a single leapfrog step's acceptance crosses 1/2."""
    counter = [0]
    p0 = metric.draw_momentum(rng)
    H0 = -cur.lp + metric.kinetic(p0)
    start = _Point(cur.q, p0, cur.g, cur.lp)

    def logacc(e):
        new = _leapfrog(f, start, e, metric, counter)
        if not np.isfinite(new.lp):
            return -np.inf
        return H0 - (-new.lp + metric.kinetic(new.p))
    la = logacc(eps)
    up = 1.0 if la > math.log(0.5) else -1.0
    for _ in range(60):
        if up > 0 and not la > math.log(0.5):
            break
        if up < 0 and la > math.log(0.5):
            break
        eps *= 2.0 ** up
        la = logacc(eps)
    return eps


def _windows(n_warmup, init=75, term=50, base=25):
    """End iterations (exclusive, 0-based) of Stan's slow metric-adaptation windows."""
    if n_warmup < init + term + base:
        init, term = int(0.15 * n_warmup), int(0.1 * n_warmup)
        base = n_warmup - init - term
    ends, start, w = [], init, base
    while start + w < n_warmup - term:
        nxt = start + w
        if nxt + 2 * w > n_warmup - term:       # stretch the last window to the terminal buffer
            nxt = n_warmup - term
        ends.append(nxt)
        start, w = nxt, 2 * w
    if not ends:
        ends = [n_warmup - term]
    return init, ends


@dataclass
class NUTSResult:
    """Draws and diagnostics of one or more NUTS chains.

    Attributes
    ----------
    samples : array (chains, draws, dim) -- post-warm-up positions, in the parameterization of the target.
    logp : array (chains, draws) -- log posterior at each draw.
    stats : dict of arrays (chains, draws) -- ``accept``, ``n_leapfrog``, ``depth``, ``divergent``, ``energy``.
    step_size : array (chains,) -- adapted step size.
    inv_metric : list of arrays -- adapted inverse metric per chain (dense or diagonal).
    warmup : dict -- the same per-iteration stats for the warm-up phase.
    seconds : float -- wall time.
    """
    samples: np.ndarray
    logp: np.ndarray
    stats: dict
    step_size: np.ndarray
    inv_metric: list
    warmup: dict = field(default_factory=dict)
    seconds: float = 0.0

    def summary(self, names=None):
        """Text table: mean, sd, split R-hat and bulk ESS per parameter, plus divergences, tree depth and E-BFMI."""
        s = self.samples
        names = names or [f"theta[{i}]" for i in range(s.shape[2])]
        r, e = split_rhat(s), ess(s)
        lines = [f"{'parameter':28s} {'mean':>13s} {'sd':>11s} {'R-hat':>7s} {'ESS':>8s}"]
        for i, nm in enumerate(names):
            lines.append(f"{nm:28s} {s[:, :, i].mean():13.6g} {s[:, :, i].std():11.4g} {r[i]:7.3f} {e[i]:8.0f}")
        st = self.stats
        en = st["energy"]
        bfmi = [float(np.sum(np.diff(c) ** 2) / np.sum((c - c.mean()) ** 2)) for c in en]
        lines.append(f"divergent {int(st['divergent'].sum())} of {st['divergent'].size}; mean leapfrog steps "
                     f"{st['n_leapfrog'].mean():.1f}; max depth {int(st['depth'].max())}; accept {st['accept'].mean():.3f}; "
                     f"step sizes {np.array2string(self.step_size, precision=3)}; E-BFMI {np.array2string(np.array(bfmi), precision=2)}")
        lines.append(f"gradient evaluations {int(st['n_leapfrog'].sum() + self.warmup.get('n_leapfrog', np.zeros(1)).sum())}"
                     f" ({self.seconds:.0f} s); min ESS {e.min():.0f}")
        return "\n".join(lines)


def nuts(logp_and_grad: Callable, x0, n_draws=1000, n_warmup=1000, inv_metric=None, adapt_metric=True, dense=True,
         target_accept=0.8, max_tree_depth=10, step_size=None, seed=None, max_delta_H=1000.0, callback=None):
    """Run one NUTS chain.

    Parameters
    ----------
    logp_and_grad : callable ``theta -> (log p, grad log p)``.  Returning ``-inf`` (or raising ValueError/RuntimeError)
        marks a forbidden point; the trajectory there counts as divergent.
    x0 : start, shape (dim,).
    n_draws, n_warmup : iterations kept / used for adaptation.
    inv_metric : starting inverse metric: (dim, dim) covariance (e.g. a Laplace covariance) or (dim,) variances;
        default identity.  It sets the scale of every direction, so start from a covariance when one is available.
    adapt_metric : re-estimate the inverse metric in Stan's windows during warm-up (dense if ``dense``, else diagonal).
    target_accept : dual-averaging target of the mean acceptance statistic (0.8; raise to 0.9-0.95 if divergent).
    max_tree_depth : at most 2**max_tree_depth leapfrog steps per draw.
    step_size : starting step size (default: Hoffman-Gelman heuristic).
    seed : int or numpy Generator.
    callback : optional ``callback(iteration, point_q, stats)`` after every iteration (progress, checkpoints).

    Returns
    -------
    NUTSResult with one chain.

    Raises
    ------
    ValueError if the start has no finite log posterior or the metric is invalid.
    """
    t_start = time.time()
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    q = np.asarray(x0, float).copy()
    n = q.size
    lp, g = _evaluate(logp_and_grad, q)
    if not np.isfinite(lp):
        raise ValueError("nuts: the start x0 has no finite log posterior")
    cur = _Point(q, None, g, lp)
    metric = _Metric(inv_metric, n)
    if not dense and metric.dense:
        metric = _Metric(np.diag(metric.inv).copy(), n)
    eps = step_size or _find_step(logp_and_grad, cur, metric, rng)
    init_buf, ends = _windows(n_warmup) if n_warmup > 0 else (0, [])
    # dual averaging
    mu, Hbar, log_eps_bar, m_da = math.log(10.0 * eps), 0.0, 0.0, 0
    gamma, t0, kappa = 0.05, 10.0, 0.75
    window = []
    keys = ("accept", "n_leapfrog", "depth", "divergent", "energy")
    wstats = {k: np.zeros(n_warmup) for k in keys}
    dstats = {k: np.zeros(n_draws) for k in keys}
    draws = np.zeros((n_draws, n))
    lps = np.zeros(n_draws)
    win_start = init_buf
    for it in range(n_warmup + n_draws):
        cur, st = _transition(logp_and_grad, cur, eps, metric, rng, max_tree_depth, max_delta_H)
        if it < n_warmup:
            for k in keys:
                wstats[k][it] = st[k]
            m_da += 1
            Hbar = (1.0 - 1.0 / (m_da + t0)) * Hbar + (target_accept - st["accept"]) / (m_da + t0)
            log_eps = mu - math.sqrt(m_da) / gamma * Hbar
            eta = m_da ** (-kappa)
            log_eps_bar = eta * log_eps + (1.0 - eta) * log_eps_bar
            eps = math.exp(log_eps)
            if adapt_metric and it >= init_buf and ends and it < ends[-1]:
                window.append(cur.q.copy())
            if adapt_metric and ends and it + 1 == ends[0]:
                w = np.array(window)
                k = w.shape[0]
                if k > 2:
                    if dense:
                        cov = np.atleast_2d(np.cov(w, rowvar=False))
                        cov = (k / (k + 5.0)) * cov + 1e-3 * (5.0 / (k + 5.0)) * np.diag(np.diag(cov))
                    else:
                        v = w.var(axis=0, ddof=1)
                        cov = (k / (k + 5.0)) * v + 1e-3 * (5.0 / (k + 5.0)) * v
                    try:
                        metric = _Metric(cov, n)
                    except ValueError:
                        warnings.warn("nuts: adapted covariance not positive definite; metric kept")
                window = []
                ends = ends[1:]
                eps = _find_step(logp_and_grad, cur, metric, rng, eps)
                mu, Hbar, log_eps_bar, m_da = math.log(10.0 * eps), 0.0, 0.0, 0
            if it + 1 == n_warmup:
                eps = math.exp(log_eps_bar)
        else:
            j = it - n_warmup
            draws[j], lps[j] = cur.q, cur.lp
            for k in keys:
                dstats[k][j] = st[k]
        if callback is not None:
            callback(it, cur.q, st)
    return NUTSResult(samples=draws[None], logp=lps[None], stats={k: v[None] for k, v in dstats.items()},
                      step_size=np.array([eps]), inv_metric=[metric.inv], warmup={k: v[None] for k, v in wstats.items()},
                      seconds=time.time() - t_start)


# ---------------------------------------------------------------------------------------------------------------------
# several chains in parallel

def _chain_worker(args):
    factory, x0, kw, seed = args
    return nuts(factory(), x0, seed=seed, **kw)


def sample(target_factory, x0, n_chains=4, nprocs=None, seed=None, **kw):
    """Run ``n_chains`` NUTS chains in separate processes and stack them.

    Parameters
    ----------
    target_factory : picklable callable with no arguments that returns ``theta -> (log p, grad)`` inside the worker
        (an :class:`AzurePosterior` is its own factory: it builds its AZURE2 session on first use in each process).
    x0 : (dim,) start for every chain, or (n_chains, dim) starts.
    n_chains, nprocs : chains and worker processes (default ``n_chains``).
    seed : base seed; chain c uses seed + c.
    **kw : passed to :func:`nuts`.

    Returns
    -------
    NUTSResult with ``n_chains`` chains.
    """
    import multiprocessing as mp
    x0 = np.asarray(x0, float)
    starts = np.repeat(x0[None], n_chains, axis=0) if x0.ndim == 1 else x0
    if starts.shape[0] != n_chains:
        raise ValueError(f"sample: {starts.shape[0]} starts for {n_chains} chains")
    base = 0 if seed is None else int(seed)
    jobs = [(target_factory, starts[c], kw, base + c) for c in range(n_chains)]
    t = time.time()
    if (nprocs or n_chains) == 1:
        res = [_chain_worker(j) for j in jobs]
    else:
        with mp.get_context("spawn").Pool(nprocs or n_chains) as pool:
            res = pool.map(_chain_worker, jobs)
    return NUTSResult(samples=np.concatenate([r.samples for r in res]), logp=np.concatenate([r.logp for r in res]),
                      stats={k: np.concatenate([r.stats[k] for r in res]) for k in res[0].stats},
                      step_size=np.concatenate([r.step_size for r in res]), inv_metric=sum((r.inv_metric for r in res), []),
                      warmup={k: np.concatenate([r.warmup[k] for r in res]) for k in res[0].warmup},
                      seconds=time.time() - t)


# ---------------------------------------------------------------------------------------------------------------------
# diagnostics

def _autocov(x):
    n = x.size
    f = np.fft.rfft(x - x.mean(), n=2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n] / n
    return ac


def ess(samples):
    """Effective sample size per parameter (Vehtari et al. 2021, Geyer's initial monotone sequence, all chains).

    ``samples``: (chains, draws, dim) or (chains, draws).  Returns an array (dim,) (or a float)."""
    s = np.asarray(samples, float)
    squeeze = s.ndim == 2
    if squeeze:
        s = s[:, :, None]
    m, n, d = s.shape
    out = np.zeros(d)
    for i in range(d):
        x = s[:, :, i]
        acov = np.array([_autocov(c) for c in x])
        W = acov[:, 0].mean() * n / (n - 1.0)
        B = n * x.mean(axis=1).var(ddof=1) if m > 1 else 0.0
        var_plus = (n - 1.0) / n * W + B / n
        if var_plus <= 0:
            out[i] = float(m * n)
            continue
        rho = 1.0 - (W - acov.mean(axis=0) * n / (n - 1.0)) / var_plus
        rho[0] = 1.0
        tau, prev, t = -1.0, np.inf, 0
        while t + 1 < n:
            P = rho[t] + rho[t + 1]
            if P < 0:
                break
            P = min(P, prev)
            tau += 2.0 * P
            prev = P
            t += 2
        out[i] = m * n / max(tau, 1.0 / math.log10(m * n + 10))
    return float(out[0]) if squeeze else out


def split_rhat(samples):
    """Split-chain potential scale reduction R-hat per parameter (Gelman et al.; ~1.00 when converged, > 1.01 suspect)."""
    s = np.asarray(samples, float)
    squeeze = s.ndim == 2
    if squeeze:
        s = s[:, :, None]
    m, n, d = s.shape
    h = n // 2
    sp = np.concatenate([s[:, :h], s[:, n - h:]], axis=0)
    W = sp.var(axis=1, ddof=1).mean(axis=0)
    B = h * sp.mean(axis=1).var(axis=0, ddof=1)
    r = np.sqrt(((h - 1.0) / h * W + B / h) / np.where(W > 0, W, np.inf))
    return float(r[0]) if squeeze else r


# ---------------------------------------------------------------------------------------------------------------------
# AZURE2 posterior

class AzurePosterior:
    """log posterior of an AZURE2 project for :func:`nuts` / :func:`sample`.

        log p(theta) = -objective(x) / 2  -  sum_k ((x_k - mean_k) / sigma_k)^2 / 2

    with ``objective`` = :meth:`pyazr.azure2.objective` (data chi-squared, divided by ``temper``, plus AZURE2's own
    normalization/shift/... penalties) and ``x`` = the non-fixed RWA parameter vector (``params_rwa`` order) with the
    sampled entries replaced by ``theta``.  The gradient: data chi-squared from the engine's adjoint
    (:meth:`pyazr.azure2.chi2_and_grad`), penalties by central differences (they are cheap closed forms).

    Parameters
    ----------
    azr : path of the .azr project.
    flags : keyword arguments for :class:`pyazr.azure2` (e.g. ``{"use_brune": True, "ignore_externals": True}``);
        take them from the fit's job script, they are CLI-only.
    sample : indices into ``params_rwa`` or parameter names to sample (default: all non-fixed parameters).
    x_template : the full non-fixed vector the unsampled entries are held at (default: the .azr's ``params_rwa``).
    priors : ``{name or index: (mean, sigma)}`` Gaussian priors in RWA space.
    temper : divide the data chi-squared by this (global error inflation, s^2).
    workdir : if given, each process works on a private copy in ``<workdir>/w<pid>`` (output directory and EC cache of
        its own; every other entry of the project's directory is linked), as parallel chains need.

    The engine session is created on the first call inside each process, so the object pickles cheaply.
    """

    def __init__(self, azr, flags=None, sample=None, x_template=None, priors=None, temper=1.0, workdir=None):
        self.azr = os.path.abspath(azr)
        self.flags = dict(flags or {})
        self._sample = sample
        self._x_template = None if x_template is None else np.asarray(x_template, float)
        self._priors = dict(priors or {})
        self.temper = float(temper)
        self.workdir = None if workdir is None else os.path.abspath(workdir)
        self._m = None

    def __getstate__(self):
        d = dict(self.__dict__)
        d["_m"] = None
        return d

    def __call__(self, theta=None):
        if theta is None:                       # factory protocol for sample(): returns itself
            return self
        return self.logp_and_grad(theta)

    # -- setup ---------------------------------------------------------------------------------------------------
    def _session(self):
        if self._m is not None:
            return self._m
        from .azure2 import azure2
        path, cwd = self.azr, os.path.dirname(self.azr)
        if self.workdir:
            from .azrfile import AzrModel
            d = os.path.join(self.workdir, f"w{os.getpid()}")
            os.makedirs(os.path.join(d, "output"), exist_ok=True)
            os.makedirs(os.path.join(d, "checks"), exist_ok=True)
            for e in os.listdir(cwd):
                if e in ("output", "checks") or e.endswith(".azr"):
                    continue
                if not os.path.lexists(os.path.join(d, e)):
                    os.symlink(os.path.join(cwd, e), os.path.join(d, e))
            mdl = AzrModel.from_file(self.azr)
            mdl.set_output_dir("output/")
            path, cwd = mdl.write(os.path.join(d, os.path.basename(self.azr))), d
        self._m = azure2(path, cwd=cwd, **self.flags)
        m = self._m
        names = [p.name for p in sorted((p for p in m.parameters if not p.fixed), key=lambda p: p.free_index)]
        self.names_all = names
        x = np.asarray(m.params_rwa, float) if self._x_template is None else self._x_template.copy()
        if x.size != len(names):
            raise ValueError(f"AzurePosterior: x_template has {x.size} entries, the model {len(names)} free parameters")
        self.x_template = x
        if self._sample is None:
            idx = np.arange(len(names))
        else:
            idx = np.array([names.index(s) if isinstance(s, str) else int(s) for s in self._sample], int)
        self.index = idx
        self.names = [names[i] for i in idx]
        pi, pm, ps = [], [], []
        for k, (mu, sg) in self._priors.items():
            j = names.index(k) if isinstance(k, str) else int(k)
            pi.append(j); pm.append(float(mu)); ps.append(float(sg))
        self.prior_idx, self.prior_mean, self.prior_sig = np.array(pi, int), np.array(pm), np.array(ps)
        # the penalties are closed forms in the normalization/shift/target parameters; level energies and widths enter
        # only through a priorCovariance block, so only those columns need finite differences
        kinds = {p.free_index: p.kind for p in m.parameters if not p.fixed}
        has_cov = float(m.penalties(x)["prior_covariance"]) != 0.0
        self._pen_cols = np.array([a for a, j in enumerate(idx) if has_cov or kinds.get(int(j)) not in ("energy", "width")], int)
        return m

    @property
    def x0(self):
        """Start vector theta (the template's sampled entries)."""
        self._session()
        return self.x_template[self.index].copy()

    def full(self, theta):
        """The full non-fixed RWA vector for ``theta``."""
        self._session()
        x = self.x_template.copy()
        x[self.index] = theta
        return x

    def _penalty(self, x):
        pen = self._m.penalties(x)
        return (float(np.sum(pen["norm"])) + float(np.sum(pen["shift"])) + float(np.sum(pen["shift_sqrt"]))
                + float(np.sum(pen["shift_fn"])) + pen["screening"] + pen["density"] + pen["doppler"]
                + pen["thickness"] + pen["prior_covariance"])

    def _prior(self, x):
        if self.prior_idx.size == 0:
            return 0.0, np.zeros_like(x)
        z = (x[self.prior_idx] - self.prior_mean) / self.prior_sig
        g = np.zeros_like(x)
        np.add.at(g, self.prior_idx, -z / self.prior_sig)
        return -0.5 * float(z @ z), g

    # -- target ----------------------------------------------------------------------------------------------------
    def logp_and_grad(self, theta):
        """(log p, d log p / d theta) at ``theta``; ``-inf`` where the engine fails."""
        m = self._session()
        x = self.full(np.asarray(theta, float))
        chi2, gchi2 = m.chi2_and_grad(x)
        if not np.isfinite(chi2):
            return -np.inf, np.zeros(self.index.size)
        P = self._penalty(x)
        gP = np.zeros(x.size)
        for j in self.index[self._pen_cols]:
            h = 1e-6 * max(1.0, abs(x[j]))
            xp, xm = x.copy(), x.copy()
            xp[j] += h; xm[j] -= h
            gP[j] = (self._penalty(xp) - self._penalty(xm)) / (2.0 * h)
        lpp, gpp = self._prior(x)
        lp = -0.5 * (chi2 / self.temper + P) + lpp
        g = -0.5 * (gchi2 / self.temper + gP) + gpp
        return lp, g[self.index]

    def laplace_covariance(self, theta):
        """Gauss-Newton Laplace covariance at ``theta``: (J^T J / temper + H_pen/2 + prior)^-1 on the sampled entries.

        J = :meth:`pyazr.azure2.residual_jacobian`; the penalty Hessian by central differences of the penalty gradient.
        A good ``inv_metric`` for :func:`nuts` when ``theta`` is near the mode."""
        m = self._session()
        x = self.full(np.asarray(theta, float))
        _, J = m.residual_jacobian(x)
        Js = J[:, self.index]
        H = Js.T @ Js / self.temper
        k = self.index.size
        Hp = np.zeros((k, k))
        cols = self._pen_cols
        for a in cols:
            j = self.index[a]
            h = 1e-4 * max(1.0, abs(x[j]))
            for s in (+1, -1):
                xs = x.copy(); xs[j] += s * h
                gs = np.zeros(k)
                for b in cols:
                    i = self.index[b]
                    hi = 1e-6 * max(1.0, abs(xs[i]))
                    xp, xm = xs.copy(), xs.copy(); xp[i] += hi; xm[i] -= hi
                    gs[b] = (self._penalty(xp) - self._penalty(xm)) / (2 * hi)
                Hp[a] += s * gs / (2 * h)
        H = H + 0.25 * (Hp + Hp.T)
        if self.prior_idx.size:
            pos = {j: a for a, j in enumerate(self.index)}
            for j, sg in zip(self.prior_idx, self.prior_sig):
                if j in pos:
                    H[pos[j], pos[j]] += 1.0 / sg ** 2
        return np.linalg.inv(H)
