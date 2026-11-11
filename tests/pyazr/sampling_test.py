"""pyazr.sampling: the NUTS sampler and the AZURE2 posterior it samples.

1. A 12-dimensional correlated Gaussian (four decades of scale, rotated): means and variances recovered, R-hat ~ 1,
   no divergences -- with the identity metric adapted during warm-up, and with the exact covariance as the metric.
2. ess / split_rhat on independent draws (ESS ~ N, R-hat ~ 1) and on two chains stuck apart (R-hat >> 1).
3. AzurePosterior on tests/shift_fn: log p = -objective/2 at the .azr values, the gradient against central
   differences, the Laplace covariance positive definite, and a short NUTS run with finite draws.
Run from anywhere:  python3 tests/pyazr/sampling_test.py"""
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


import importlib.util                                   # the sampler itself needs no engine
spec = importlib.util.spec_from_file_location("pyazr_sampling", os.path.join(ROOT, "pyazr", "sampling.py"))
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)

rng = np.random.default_rng(1)
d = 12
Q, _ = np.linalg.qr(rng.standard_normal((d, d)))
C = Q @ np.diag(np.logspace(-2, 2, d)) @ Q.T
P = np.linalg.inv(C)
mu = rng.standard_normal(d)
gauss = lambda x: (-0.5 * (x - mu) @ P @ (x - mu), -P @ (x - mu))
for label, kw in (("adapted metric", {}), ("exact covariance metric", dict(inv_metric=C, adapt_metric=False))):
    r = S.nuts(gauss, mu + 0.3, n_draws=1500, n_warmup=700, seed=3, **kw)
    s = r.samples[0]
    z = np.linalg.solve(np.linalg.cholesky(C), s.mean(0) - mu)
    vr = np.diag(np.cov(s, rowvar=False)) / np.diag(C)
    e = S.ess(r.samples)
    check(f"Gaussian, {label}: whitened mean error rms {np.sqrt(np.mean(z ** 2)):.3f} < 0.1",
          np.sqrt(np.mean(z ** 2)) < 0.1)
    check(f"Gaussian, {label}: variance ratios {vr.min():.2f}..{vr.max():.2f} within 0.8..1.25", 0.8 < vr.min() and vr.max() < 1.25)
    check(f"Gaussian, {label}: no divergences ({int(r.stats['divergent'].sum())}), min ESS {e.min():.0f} > 300",
          r.stats["divergent"].sum() == 0 and e.min() > 300)

# a good starting metric must survive warm-up adaptation in dimension > window size (25 draws, n = 40): before the
# shrinkage toward the current metric, the 7Be run (n = 54) went from ~22 to ~800 leapfrog steps per draw
d2 = 40
Q2, _ = np.linalg.qr(rng.standard_normal((d2, d2)))
C2 = Q2 @ np.diag(np.logspace(-2, 2, d2)) @ Q2.T
P2 = np.linalg.inv(C2)
g2 = lambda x: (-0.5 * x @ P2 @ x, -P2 @ x)
fixed = S.nuts(g2, np.zeros(d2), n_draws=300, n_warmup=300, inv_metric=C2, adapt_metric=False, seed=4)
adapt = S.nuts(g2, np.zeros(d2), n_draws=300, n_warmup=300, inv_metric=C2, adapt_metric=True, seed=4)
lf, la = fixed.stats["n_leapfrog"].mean(), adapt.stats["n_leapfrog"].mean()
check(f"Laplace metric kept through adaptation (n = {d2}): {la:.1f} leapfrog steps per draw vs {lf:.1f} unadapted (< 2.5x; 5.5x before the fix)",
      la < 2.5 * lf)

iid = rng.standard_normal((4, 2000, 2))
check(f"ess of iid draws {S.ess(iid).min():.0f} within 6000..10000", 6000 < S.ess(iid).min() < 10000)
check(f"split R-hat of iid draws {S.split_rhat(iid).max():.3f} < 1.01", S.split_rhat(iid).max() < 1.01)
stuck = np.stack([rng.standard_normal(1000), 5 + rng.standard_normal(1000)])
check(f"split R-hat of two separated chains {S.split_rhat(stuck):.2f} > 1.5", S.split_rhat(stuck) > 1.5)

try:
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    from pyazr import azure2  # noqa: F401
    from pyazr.sampling import AzurePosterior, nuts
except ImportError as err:
    print(f"skip: engine not available ({err}); sampler checks above ran")
    sys.exit(1 if failures else 0)

with tempfile.TemporaryDirectory() as tmp:
    w = os.path.join(tmp, "shift_fn")
    shutil.copytree(os.path.join(ROOT, "tests", "shift_fn"), w)
    for junk in ("output", "checks"):
        shutil.rmtree(os.path.join(w, junk), ignore_errors=True)
        os.makedirs(os.path.join(w, junk))
    T = AzurePosterior(os.path.join(w, "shift_fn.azr"), workdir=os.path.join(tmp, "work"))
    x0 = T.x0
    lp, g = T(x0)
    obj = T._m.objective(T.full(x0))
    check(f"log p {lp:.6f} = -objective/2 {-0.5 * obj:.6f}", abs(lp + 0.5 * obj) < 1e-9 * max(1.0, abs(obj)))
    worst = 0.0
    for j in range(x0.size):
        h = 1e-6 * max(1.0, abs(x0[j]))
        xp, xm = x0.copy(), x0.copy(); xp[j] += h; xm[j] -= h
        fd = (T(xp)[0] - T(xm)[0]) / (2 * h)
        worst = max(worst, abs(fd - g[j]) / max(1.0, abs(fd)))
    check(f"gradient vs central differences on {x0.size} parameters: worst relative {worst:.2e} < 1e-4", worst < 1e-4)
    Cl = T.laplace_covariance(x0)
    check("Laplace covariance positive definite", np.all(np.linalg.eigvalsh(0.5 * (Cl + Cl.T)) > 0))
    r = nuts(T, x0, n_draws=20, n_warmup=20, inv_metric=Cl, seed=1)
    check(f"short NUTS run: {r.samples.shape[1]} finite draws", np.all(np.isfinite(r.samples)))

print("sampling: all checks passed" if not failures else f"sampling: {len(failures)} check(s) FAILED")
sys.exit(1 if failures else 0)
