#!/usr/bin/env python3
"""Self-test of the azure2-uncertainty scripts on tests/upgrade/T3 (8Be, synthetic) in a temporary directory:
fisher.py, write_covariance.py, profile_scan.py, nuts.py, sample_nuts.py, sample_emcee.py, mcmc_stats.py,
draw_bands.py.

The best fit of T3_start.azr (objective 509.348, 8 free parameters) is found first with azure2-fit's lsq_fit.py;
"T3_two" is that fit with every level fixed except the 1+ level at 17.64 MeV (energy_4, width_4_1: two free
parameters, constrained by one data set).
  1. fisher.py: F^-1 equals the engine's Levenberg-Marquardt covariance (J^T J + priors)^-1 at the same point
     (largest difference 1e-3 of sqrt(C_ii C_jj)); dividing the errors of the only data set that sees width_4_1 by 2
     halves its error (0.5 to 1e-6).
  2. write_covariance.py: a covariance of the error column alone leaves chi-squared unchanged (1e-9); a covariance
     joining two data sets with a 5 % absolute term per file and one common term gives the chi-squared r^T R^-1 r of
     an independent numpy calculation (1e-8), the program's and the whitened rows' splits differ set by set and
     agree in sum; the same C written as one full matrix gives the same chi-squared; --check refuses a matrix that
     is not positive definite (exit 1).
  3. profile_scan.py on T3_two: the profile interval of energy_4 at Delta = 1 equals MINOS's within 5 %.
  4. nuts.py on a 10-dimensional Gaussian with condition number 1e6, whitened: the draws' covariance is the
     target's (each variance within 25 %), no divergence; a hard wall (log p = -inf beyond it) is never crossed.
  5. sample_nuts.py on T3_two (exact Hessian at the fit): R-hat < 1.05, no divergence, each posterior sd equal to
     the Laplace sd within 3 Monte-Carlo errors + 3 %.
  6. sample_emcee.py on T3_two: posterior sds within 15 % of the Laplace sds.
  7. mcmc_stats.py: the autocorrelation time of an AR(1) chain (rho = 0.8: tau = 9) within 20 %; R-hat 1 for
     chains of one target and > 1.1 for chains with shifted means; the command's exit codes.
  8. draw_bands.py --laplace on T3_two: the band from draws equals the linear band within 10 %.
Prints PASS / FAIL lines; exit 0 all passed, 1 a failure, 77 the engine cannot be imported.  About 90 s at one
thread.  Run: python3 selftest.py   (PYTHONPATH: the checkout; PYAZR_ENGINE_DIR where the build put the engine)
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
FIT = os.path.join(ROOT, ".claude", "skills", "azure2-fit", "scripts", "lsq_fit.py")
failed = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail and not ok else ""), flush=True)
    if not ok:
        failed.append(name)


def run(script, *args, cwd):
    path = script if os.path.isabs(script) else os.path.join(HERE, script)
    p = subprocess.run([sys.executable, path, *args], cwd=cwd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def fisher_checks(np, pyazr, tmp):
    code, out = run("fisher.py", "T3/two.azr", "--improve", "p7li_pp_lab150=2", "--json", "fisher.json", cwd=tmp)
    js = json.load(open(os.path.join(tmp, "fisher.json"))) if code == 0 else {"names": ["width_4_1"], "improve": {
        "p7li_pp_lab150=2": [0]}, "sigma": [1]}
    sys.path.insert(0, HERE)
    import fisher
    with pyazr.Session(os.path.join(tmp, "T3", "best.azr")) as s:
        f = s.objective()
        x = f.x0
        F, _, _, _ = fisher.fisher(s, f, x)
        C_lm = np.asarray(s.fit_with_engine("lm", x0=x).covariance)
    C = np.linalg.inv(F)
    d = np.sqrt(np.outer(np.diag(C), np.diag(C)))
    dev = np.max(np.abs(C - C_lm) / d)
    check(f"fisher: F^-1 = the engine's LM covariance (largest difference {dev:.1e} of sqrt(C_ii C_jj))",
          code == 0 and dev < 1e-3, out[-300:])
    j = js["names"].index("width_4_1")
    ratio = js["improve"]["p7li_pp_lab150=2"][j] / js["sigma"][j]
    check(f"fisher: on two parameters seen by one data set, its errors / 2 give sigma(width_4_1) x {ratio:.6f}",
          abs(ratio - 0.5) < 1e-6)


def covariance_checks(np, pyazr, tmp):
    import write_covariance as wc
    t3 = os.path.join(tmp, "T3")
    files = ["data/p7li_pa.dat", "data/p7li_pn.dat"]
    code, out = run("write_covariance.py", "cov_stat.txt", "--file", files[0], cwd=t3)
    p = pyazr.Project.open(os.path.join(t3, "best.azr"))
    with p.session() as s:
        base = s.evaluate(s.objective().x0).dataset("p7li_pa").chi2
    p.set("p7li_pa", covariance="cov_stat.txt")
    p.save(os.path.join(t3, "cov_stat.azr"))
    with pyazr.Session(os.path.join(t3, "cov_stat.azr")) as s:
        c = s.evaluate(s.objective().x0).dataset("p7li_pa").chi2
    check(f"write_covariance: the error column alone leaves chi2 unchanged ({base:.8f} vs {c:.8f})",
          code == 0 and abs(c - base) <= 1e-9 * base, out[-300:])

    code, out = run("write_covariance.py", "cov_joint.txt", "--file", files[0], "--file", files[1],
                    "--common-per-file", "0.05", "--common", "0.03", "--absolute", cwd=t3)
    p = pyazr.Project.open(os.path.join(t3, "best.azr"))
    for name in ("p7li_pa", "p7li_pn"):
        p.set(name, covariance="cov_joint.txt")
    p.save(os.path.join(t3, "cov_joint.azr"))
    with pyazr.Session(os.path.join(t3, "cov_joint.azr")) as s:
        f = s.objective()
        ev = s.evaluate(f.x0)
        r = f.residuals(f.x0)
        sets = [ev.dataset("p7li_pa"), ev.dataset("p7li_pn")]
        prog = sum(d.chi2 for d in sets)
        white = sum(float(np.sum(r[d.rows] ** 2)) for d in sets)
        split = [(d.chi2, float(np.sum(r[d.rows] ** 2))) for d in sets]
        res = np.concatenate([d.residuals for d in sets])
    _, terms = wc.read(os.path.join(t3, "cov_joint.txt"))
    vals = np.concatenate([wc.data_rows(os.path.join(t3, f))[0] for f in files])
    _, Cm = wc.dense(terms, len(vals), vals)
    D = np.sqrt(np.diag(Cm))
    R = Cm / np.outer(D, D)
    ref = float(res @ np.linalg.solve(R, res))
    check(f"write_covariance: joined sets, the program's chi2 {prog:.8f} = numpy r^T R^-1 r {ref:.8f}",
          code == 0 and abs(prog - ref) <= 1e-8 * ref, out[-300:])
    check("the two splits differ per set and agree in sum: " + ", ".join(f"{a:.5f}/{b:.5f}" for a, b in split),
          abs(prog - white) <= 1e-8 * prog and abs(split[0][0] - split[0][1]) > 1e-4)

    np.savetxt(os.path.join(t3, "C.txt"), Cm)
    code, out = run("write_covariance.py", "cov_matrix.txt", "--file", files[0], "--file", files[1], "--no-stat",
                    "--matrix", "C.txt", cwd=t3)
    p = pyazr.Project.open(os.path.join(t3, "cov_joint.azr"))
    for name in ("p7li_pa", "p7li_pn"):
        p.set(name, covariance="cov_matrix.txt")
    p.save(os.path.join(t3, "cov_matrix.azr"))
    with pyazr.Session(os.path.join(t3, "cov_matrix.azr")) as s:
        ev = s.evaluate(s.objective().x0)
        cm = ev.dataset("p7li_pa").chi2 + ev.dataset("p7li_pn").chi2
    check(f"write_covariance: the same C as one matrix term gives {cm:.8f}",
          code == 0 and abs(cm - prog) <= 1e-8 * prog, out[-300:])

    bad = -np.eye(len(vals)) * 1e-6 + np.diag(Cm.diagonal()) * 0
    np.savetxt(os.path.join(t3, "bad.txt"), bad)
    code, out = run("write_covariance.py", "cov_bad.txt", "--file", files[0], "--file", files[1], "--no-stat",
                    "--matrix", "bad.txt", cwd=t3)
    check("write_covariance --check refuses a matrix that is not positive definite (exit 1)", code == 1, out[-200:])


def two_parameter_project(pyazr, tmp):
    p = pyazr.Project.open(os.path.join(tmp, "T3", "best.azr"))
    for lv in p.levels():
        if lv.jpi in ("2+", "2-") and not lv.fixed:
            p.set_level(lv.label, fixed=True)
            for c in lv.channels:
                p.set_channel(lv.label, c.pair, c.l, c.s, fixed=True)
    return p.save(os.path.join(tmp, "T3", "two.azr"))


def profile_checks(tmp):
    code, out = run("profile_scan.py", "T3/two.azr", "--parameter", "energy_4", "--points", "7", "--span", "2",
                    "--minos", "--out", "profile.json", cwd=tmp)
    js = json.load(open(os.path.join(tmp, "profile.json"))) if code == 0 else {}
    ok = code == 0 and js.get("minos")
    if ok:
        lo, hi = js["interval"][0] - js["best"], js["interval"][1] - js["best"]
        mlo, mhi = js["minos"]
        ok = abs(lo / mlo - 1) < 0.05 and abs(hi / mhi - 1) < 0.05
        out = f"profile {lo:.3g}/{hi:.3g}, MINOS {mlo:.3g}/{mhi:.3g}"
    check(f"profile_scan: the Delta = 1 interval of energy_4 equals MINOS within 5 % ({out[-120:]})", ok)


def nuts_numpy_checks(np):
    import nuts
    d = 10
    rng = np.random.default_rng(1)
    Q, _ = np.linalg.qr(rng.standard_normal((d, d)))
    lam = np.logspace(0, -6, d)
    C = (Q * lam) @ Q.T
    P = np.linalg.inv(C)
    L = np.linalg.cholesky(C)                          # whitening x = L z
    lp = lambda z: (-0.5 * float((L @ z) @ P @ (L @ z)), -(L.T @ P @ (L @ z)))
    out = nuts.run_chain(lp, d, warmup=200, draws=1000, seed=[1, 0])
    X = out["z"] @ L.T
    ratio = np.diag(np.cov(X.T)) / lam.max() / (np.diag(C) / lam.max())
    check(f"nuts: whitened Gaussian, cond 1e6, variances / target in [{ratio.min():.2f}, {ratio.max():.2f}], "
          f"{int(out['divergent'].sum())} divergences",
          np.all(np.abs(ratio - 1) < 0.25) and out["divergent"].sum() == 0)

    def wall(z):
        if z[0] > 1.0:
            return -np.inf, np.zeros(2)
        return -0.5 * float(z @ z), -z
    out = nuts.run_chain(wall, 2, warmup=100, draws=400, seed=[2, 0])
    check(f"nuts: a hard wall at z0 = 1 is never crossed (max {out['z'][:, 0].max():.3f}, "
          f"{int(out['divergent'].sum())} divergent transitions)", out["z"][:, 0].max() <= 1.0)


def sampler_checks(np, tmp):
    code, out = run("sample_nuts.py", "T3/two.azr", "--chains", "2", "--warmup", "200", "--draws", "400", "--out",
                    "nuts.npz", cwd=tmp)
    ok = code in (0, 1) and os.path.exists(os.path.join(tmp, "nuts.npz"))
    if ok:
        import mcmc_stats
        d = np.load(os.path.join(tmp, "nuts.npz"))
        x, lap = d["x"], d["laplace_sigma"]
        rh = [mcmc_stats.split_rhat(x[:, :, j]) for j in range(2)]
        es = [mcmc_stats.ess(x[:, :, j]) for j in range(2)]
        sd = x.reshape(-1, 2).std(axis=0)
        tol = [3 / np.sqrt(2 * e) + 0.03 for e in es]
        ok = max(rh) < 1.05 and d["divergent"].sum() == 0 and all(abs(sd[j] / lap[j] - 1) < tol[j] for j in range(2))
        out = (f"R-hat {max(rh):.3f}, ESS {min(es):.0f}, sd / Laplace sd {sd[0] / lap[0]:.3f}, {sd[1] / lap[1]:.3f}, "
               f"{int(d['divergent'].sum())} divergences, {d['ngrad'].mean():.1f} gradients per draw")
    check(f"sample_nuts on two parameters recovers the Laplace covariance ({out[-160:]})", ok)
    lap = np.load(os.path.join(tmp, "nuts.npz"))["laplace_sigma"] if ok else None
    code, out = run("sample_emcee.py", "T3/two.azr", "--walkers", "16", "--steps", "700", "--burn", "200", "--out",
                    "emcee.npz", cwd=tmp)
    ok = code in (0, 1) and lap is not None
    if ok:
        d = np.load(os.path.join(tmp, "emcee.npz"))
        sd = d["x"].reshape(-1, 2).std(axis=0)
        ok = all(abs(sd / lap - 1) < 0.15)
        out = f"sd / Laplace sd {sd[0] / lap[0]:.3f}, {sd[1] / lap[1]:.3f}, acceptance {float(d['acceptance']):.2f}"
    check(f"sample_emcee on two parameters: {out[-120:]}", ok)


def stats_checks(np, tmp):
    import mcmc_stats
    rng = np.random.default_rng(3)
    rho, n = 0.8, 20000
    x = np.zeros((4, n))
    e = rng.standard_normal((4, n))
    for t in range(1, n):
        x[:, t] = rho * x[:, t - 1] + np.sqrt(1 - rho ** 2) * e[:, t]
    tau = mcmc_stats.tau(x)
    check(f"mcmc_stats: AR(1) rho 0.8, tau {tau:.2f} (exact 9)", abs(tau / 9 - 1) < 0.2)
    shifted = x + np.arange(4)[:, None] * 1.0
    r1, r2 = mcmc_stats.split_rhat(x), mcmc_stats.split_rhat(shifted)
    check(f"mcmc_stats: R-hat {r1:.4f} for one target, {r2:.2f} for shifted chains", r1 < 1.01 and r2 > 1.1)
    np.savez(os.path.join(tmp, "good.npz"), x=x[:, :, None], names=np.array(["a"]))
    np.savez(os.path.join(tmp, "bad.npz"), x=shifted[:, :, None])
    c1, _ = run("mcmc_stats.py", "--npz", "good.npz", cwd=tmp)
    c2, _ = run("mcmc_stats.py", "--npz", "bad.npz", cwd=tmp)
    c3, _ = run("mcmc_stats.py", "--npz", "good.npz", "--key", "none", cwd=tmp)
    check(f"mcmc_stats.py exit codes: converged {c1}, not converged {c2}, no such array {c3}",
          (c1, c2, c3) == (0, 1, 2))


def band_checks(np, tmp):
    code, out = run("draw_bands.py", "T3/two.azr", "--laplace", "300", "--dataset", "p7li_pp_lab150", "--out",
                    "bands.npz", cwd=tmp)
    ok = code == 0
    if ok:
        d = np.load(os.path.join(tmp, "bands.npz"))
        half = 0.5 * (d["p7li_pp_lab150__q84"] - d["p7li_pp_lab150__q16"])
        lin = d["p7li_pp_lab150__linear"]
        ratio = float(np.median(half[lin > 0] / lin[lin > 0]))
        ok = abs(ratio - 1) < 0.1
        out = f"median ratio {ratio:.3f}"
    check(f"draw_bands: band from Laplace draws = linear band ({out[-120:]})", ok)


def main():
    try:
        import numpy as np
        import pyazr
        pyazr.kinds()
    except Exception as e:
        print(f"skip: the engine cannot be imported ({type(e).__name__}: {e})")
        return 77
    sys.path.insert(0, HERE)
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix="azure2-uncertainty-selftest-") as tmp:
        shutil.copytree(os.path.join(ROOT, "tests", "upgrade", "T3"), os.path.join(tmp, "T3"),
                        ignore=shutil.ignore_patterns("output", "checks"))
        os.makedirs(os.path.join(tmp, "T3", "output"))
        code, out = run(FIT, os.path.join("T3", "T3_start.azr"), "--starts", "8", "--jitter", "3", "--seed", "1",
                        "--save-fit", os.path.join("T3", "best.azr"), cwd=tmp)
        if code != 0:
            check("the best fit of T3 (azure2-fit's lsq_fit.py)", False, out[-300:])
            return 1
        two_parameter_project(pyazr, tmp)
        for part in (lambda: fisher_checks(np, pyazr, tmp), lambda: covariance_checks(np, pyazr, tmp),
                     lambda: profile_checks(tmp), lambda: nuts_numpy_checks(np), lambda: sampler_checks(np, tmp),
                     lambda: stats_checks(np, tmp), lambda: band_checks(np, tmp)):
            t = time.time()
            try:
                part()
            except Exception as e:  # a crash of one part is a failure of that part, the others still run
                check(f"a part of the self-test raised {type(e).__name__}: {e}", False)
            print(f"      ({time.time() - t:.1f} s)", flush=True)
    print(f"{'FAIL' if failed else 'PASS'}: azure2-uncertainty self-test, {time.time() - t0:.1f} s")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
