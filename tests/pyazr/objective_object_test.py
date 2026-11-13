#!/usr/bin/env python3
"""pyazr.Objective (Session.fit_objective()) against today's methods, the engine's sampler and itself.

(a) On T1, T2, T3 (tests/upgrade), 13N_fast (tests/mcmc_physical_norm), a project of each nuisance kind with a prior
    (screening: tests/electron_screening/grad.azr; stopping scale: tests/stopping_yield fit.azr; a kernel's width:
    tests/rectangular_resolution fit.azr) and the 15O gamma_decay record with priors (tests/gamma_decay
    record.azr, the variant of tests/reference/sampler_objective_test.cpp), at the start and at a displaced point:
    f(x) and f.value_and_grad(x) equal Session.objective(x) and objective_and_grad(x) bit for bit; the first rows of
    f.residuals(x) are Session.residuals(x) bit for bit, and the first rows and columns of f.residuals_and_jacobian(x)
    those of residual_jacobian(x); the sum of the squares of the rows is f(x) (1e-10), from the forward rows and from
    the Jacobian's rows; f.log_prob(x) = -0.5 f(x); x0, names and bounds are the free rows of the records.
(b) The gradient against a central difference of f, on T1 and 13N_fast (analytic there), and the Jacobian against a
    central difference of the rows on T1.
(c) log_prob is -inf at the three excluded points of the record (a branching below 0 and above 1, Gamma_gamma below
    0), where Park's J < 0 (T2 under Park, widths scaled up) and at a T2 norm below its lower limit (the fit's
    limits of every kind are walls: PosteriorSupport); and it equals the engine's sampler's ln P
    minus its constant (Session.sampler_log_probability) on the record with priors and on T2.
(d) Two threads calling one Objective at once get what each gets alone, bit for bit.

Run: python3 tests/pyazr/objective_object_test.py
"""
import os
import shutil
import sys
import tempfile
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
TESTS = os.path.join(ROOT, "tests")

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    import pyazr
    from pyazr import Session
    pyazr.kinds()
except Exception as err:  # engine not built
    print(f"skip: engine not available ({type(err).__name__}: {err})")
    sys.exit(77)


def copy_dir(src, dst, project=None, edits=()):
    """A copy of src (outputs dropped) with `project` edited by the (old, new) text pairs; the project's path."""
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("output", "checks"))
    for d in ("output", "checks"):
        os.makedirs(os.path.join(dst, d), exist_ok=True)
    if project is None:
        return None
    path = os.path.join(dst, project)
    text = open(path).read()
    for old, new in edits:
        assert old in text, (project, old)
        text = text.replace(old, new)
    open(path, "w").write(text)
    return path


RECORD_PRIORS = [
    ("    total_width  0.0385027 eV fit\n", "    total_width  0.0385027 eV fit prior 0.040 0.004\n"),
    ("    remainder    g+15O[5.183]\n", "    remainder    g+15O[5.183] prior 0.15 0.01\n"),
    ("    branching    g+15O[6.1763] 0.5843746023006179 fit\n",
     "    branching    g+15O[6.1763] 0.5843746023006179 fit prior 0.6 0.01\n"),
]


def same(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return a.shape == b.shape and bool(np.all((a == b) | (np.isnan(a) & np.isnan(b))))


def displaced(x):
    return np.array([v * (1.0 + 0.02 * (1 if k % 2 else -1)) for k, v in enumerate(x)])


def against_today(label, s):
    f = s.fit_objective()
    free = s.parameters.free
    check(f"{label}: x0, names, bounds are the free records'", same(f.x0, s.params_rwa) and f.names == free.names
          and f.bounds == [tuple(b) for b in free.bounds], f"{f.names[:3]} {free.names[:3]}")
    for name, x in (("start", f.x0), ("displaced", displaced(f.x0))):
        v = f(x)
        check(f"{label} {name}: f(x) = objective(x) bit for bit ({v:.10g})", v == s.objective(x), s.objective(x))
        val, g = f.value_and_grad(x)
        val0, g0 = s.objective_and_grad(x)
        check(f"{label} {name}: value_and_grad = objective_and_grad bit for bit", val == val0 and same(g, g0))
        rows = f.residuals(x)
        old = s.residuals(x)
        check(f"{label} {name}: the first {len(old)} rows are residuals(x) bit for bit, {len(rows)} rows",
              same(rows[:len(old)], old) and len(rows) >= len(old))
        rel = abs(float(np.sum(rows ** 2)) - v) / max(1.0, abs(v))
        check(f"{label} {name}: sum of the rows' squares = f(x) ({rel:.1e})", rel <= 1e-10, f"{np.sum(rows ** 2)} {v}")
        check(f"{label} {name}: log_prob = -f/2", f.log_prob(x) == -0.5 * v, f"{f.log_prob(x)} {-0.5 * v}")
        try:
            r_old, J_old = s.residual_jacobian(x)
        except RuntimeError:
            r_old = None
        if r_old is None:
            check(f"{label} {name}: no analytic Jacobian today, none here either", _refused(f, x))
            continue
        r, J = f.residuals_and_jacobian(x)
        n = len(r_old)
        check(f"{label} {name}: the Jacobian's first {n} rows are residual_jacobian(x)'s bit for bit",
              same(r[:n], r_old) and same(J[:n], J_old) and J.shape == (len(rows), len(x)))
        rel = abs(float(np.sum(r ** 2)) - v) / max(1.0, abs(v))
        check(f"{label} {name}: the Jacobian's rows' squares add up to f(x) ({rel:.1e})", rel <= 1e-10)
        check(f"{label} {name}: forward rows = the Jacobian's rows (1e-9)",
              len(r) == len(rows) and np.allclose(r, rows, rtol=1e-9, atol=1e-12))
    return f


def _refused(f, x):
    try:
        f.jacobian(x)
    except pyazr.EngineError:
        return True
    return False


def central_gradient(label, f, x, rel_tol):
    g = f.grad(x)
    fd = np.zeros_like(x)
    for k in range(len(x)):
        h = 1e-6 * max(abs(x[k]), 1e-3)
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] -= h
        fd[k] = (f(xp) - f(xm)) / (2 * h)
    scale = np.max(np.abs(g)) + 1e-30
    worst = float(np.max(np.abs(g - fd)) / scale)
    check(f"{label}: gradient against a central difference ({worst:.1e} of the largest entry)", worst <= rel_tol,
          f"{g} {fd}")


def central_jacobian(label, f, x, rel_tol):
    _, J = f.residuals_and_jacobian(x)
    fd = np.zeros_like(J)
    for k in range(len(x)):
        h = 1e-6 * max(abs(x[k]), 1e-3)
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] -= h
        fd[:, k] = (f.residuals(xp) - f.residuals(xm)) / (2 * h)
    worst = float(np.max(np.abs(J - fd)) / (np.max(np.abs(J)) + 1e-30))
    check(f"{label}: Jacobian against a central difference of the rows ({worst:.1e})", worst <= rel_tol)


with tempfile.TemporaryDirectory() as tmp:
    def session(name, src, project, edits=(), **options):
        path = copy_dir(src, os.path.join(tmp, name), project, edits)
        return Session(path, **options)

    print("(a) the object against today's methods")
    upgrade = os.path.join(TESTS, "upgrade")
    with session("T1", os.path.join(upgrade, "T1"), "T1.azr") as s:
        f = against_today("T1", s)
        central_gradient("T1", f, displaced(f.x0), 1e-5)
        central_jacobian("T1", f, displaced(f.x0), 1e-5)
    with session("T2", os.path.join(upgrade, "T2"), "T2.azr") as s:
        f = against_today("T2", s)
        # (c) the fit's limits of every kind are walls of the posterior: a norm below its lower limit
        x = f.x0
        k = next(i for i, p in enumerate(s.parameters.free) if p.kind == "norm")
        lower = f.bounds[k][0]
        x[k] = lower - 1e-3
        check(f"T2: a norm below its lower limit ({lower:g}) has log_prob -inf, at the limit a finite one",
              np.isfinite(lower) and f.log_prob(x) == -np.inf and np.isfinite(f.log_prob(f.x0)), f.log_prob(x))
        for name, x in (("start", f.x0), ("displaced", displaced(f.x0))):
            lp, constant = s.sess.sampler_log_probability(f.x0, x)
            got = f.log_prob(x)
            check(f"T2 {name}: log_prob = the engine sampler's ln P - its constant ({got:.10g}, C {constant:.6g})",
                  abs(got - (lp - constant)) <= 1e-10 * max(1.0, abs(got)), f"{lp - constant}")
    with session("T3", os.path.join(upgrade, "T3"), "T3.azr") as s:
        against_today("T3", s)
    fast = os.path.join(tmp, "fast_src")
    os.makedirs(fast)
    os.symlink(os.path.join(TESTS, "13N", "data"), os.path.join(fast, "data"))
    shutil.copy(os.path.join(TESTS, "mcmc_physical_norm", "13N_fast.azr"), fast)
    with session("13N_fast", fast, "13N_fast.azr") as s:
        f = against_today("13N_fast", s)
        central_gradient("13N_fast", f, displaced(f.x0), 1e-5)
    scr = os.path.join(tmp, "scr_src")
    os.makedirs(scr)
    os.symlink(os.path.join(TESTS, "13N", "data"), os.path.join(scr, "data"))
    shutil.copy(os.path.join(TESTS, "electron_screening", "grad.azr"), scr)
    with session("screening", scr, "grad.azr") as s:
        against_today("screening (13N data)", s)
    with session("stopping", os.path.join(TESTS, "stopping_yield", "models"), "fit.azr",
                 [("scale 1 fit prior 1 4%", "scale 1 fit prior 1.04 4%")]) as s:
        against_today("stopping scale", s)
    with session("kernel", os.path.join(TESTS, "rectangular_resolution", "models"), "fit.azr",
                 [("rectangular width 14 keV fit\n", "rectangular width 14 keV fit prior 16 2\n")]) as s:
        against_today("kernel width", s)

    print("(c) log_prob at the excluded points")
    record_dir = os.path.join(TESTS, "gamma_decay", "models")
    with session("record", record_dir, "record.azr", RECORD_PRIORS) as s:
        f = against_today("gamma_decay record with priors", s)
        check("record: the five free numbers have finite bounds, the walls of log_prob",
              all(np.isfinite(lo) and np.isfinite(hi) for lo, hi in f.bounds[:5]), f.bounds[:5])
        for k, value in ((1, -0.01), (2, 1.01), (0, -1e-4)):
            x = f.x0
            x[k] = value
            lp = f.log_prob(x)
            check(f"record: entry {k} = {value} (outside its limits): log_prob = -inf", lp == -np.inf, lp)
            check(f"record: entry {k} = {value}: the evaluation says excluded", s.evaluate(x).excluded)
        for name, x in (("start", f.x0), ("displaced", displaced(f.x0))):
            lp, constant = s.sess.sampler_log_probability(f.x0, x)
            got = f.log_prob(x)
            check(f"record {name}: log_prob = the engine sampler's ln P - its constant ({got:.10g})",
                  abs(got - (lp - constant)) <= 1e-10 * max(1.0, abs(got)), f"{lp} {constant}")
    with session("T2park", os.path.join(upgrade, "T2"), "T2.azr", use_park=True) as s:
        f = s.fit_objective()
        widths = [i for i, p in enumerate(s.parameters.free) if p.kind == "width"]
        x = f.x0
        for scale in (3.0, 10.0, 30.0, 100.0):
            x = f.x0
            x[widths] *= scale
            if np.min(s.park_norms(x)) < 0:
                break
        check(f"T2 under Park: widths x{scale:g} give a level with J < 0", np.min(s.park_norms(x)) < 0,
              np.min(s.park_norms(x)))
        check("T2 under Park, J < 0: log_prob = -inf while f(x) is finite (the fit's quadratic wall)",
              f.log_prob(x) == -np.inf and np.isfinite(f(x)))
        check("T2 under Park at the start (J > 0): log_prob = -f/2", f.log_prob(f.x0) == -0.5 * f(f.x0))

    print("(d) two threads on one Objective")
    with session("T1threads", os.path.join(upgrade, "T1"), "T1.azr") as s:
        f = s.fit_objective()
        points = [f.x0, displaced(f.x0)]
        alone = [(f(x), f.grad(x)) for x in points]
        got = {}

        def work(i):
            got[i] = [(f(points[i]), f.grad(points[i])) for _ in range(3)]

        threads = [threading.Thread(target=work, args=(i,)) for i in range(2)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        check("two threads at once give each thread's values alone, bit for bit",
              all(v == alone[i][0] and same(g, alone[i][1]) for i in range(2) for v, g in got[i]))

if failures:
    print(f"FAIL: {len(failures)} check(s): {failures}")
    sys.exit(1)
print("ok: Objective")
