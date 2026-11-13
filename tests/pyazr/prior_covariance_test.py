#!/usr/bin/env python3
"""The prior parameter covariance (``priorCovariance=<file>`` in ``<config>``) through pyazr.

Pure Python (no engine): ``pyazr.priors`` -- ``write_prior`` / ``read_prior`` round trip and refusals,
``read_covariance_dat`` (the ``#columns`` / ``#scaled`` header), ``prior_from_fit`` (means from param.sav, the
covariance sub-block for a ``names`` subset, the basis tag, ``columns=`` for a header-less file, refusals);
the ``<config>`` key through ``Project.config`` / ``set_config(priorCovariance=)`` (get / set / remove).

With the engine, on the p+p model of tests/norm_groups (every level parameter free; each data set its own free
norm, the group D removed):
``terms["prior covariance"]`` at moved parameters = the hand-computed quadratic form; a diagonal prior on the
two norms gives the same ``objective()`` as their normError priors (1e-10); the last rows of ``objective().residuals_and_jacobian``
are L^-1 (p - mu) (their squared sum = the prior term, the whole sum = evaluate().chi2 + the prior) and their
columns match central differences of ``prior_covariance_residuals``; a fixed listed parameter is marginalized
(``prior_covariance_names``); a refused file makes ``Session()`` raise.

Run from anywhere:  python3 tests/pyazr/prior_covariance_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


try:
    import numpy as np                          # noqa: E402
    from pyazr import Project                   # noqa: E402
    from pyazr import priors                    # noqa: E402
except Exception as err:
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)

PROJECT = os.path.join(ROOT, "tests", "norm_groups")
AZR = os.path.join(PROJECT, "norm_groups.azr")
N1, N2 = "xs_40_s08.norm", "xs_90_s09.norm"   # the two data sets' norm parameters
tmp = tempfile.mkdtemp(prefix="pyazr_prior_cov_")


def raises(fn, text):
    try:
        fn()
    except ValueError as e:  # FormatError is a ValueError
        return text in str(e)
    return False


try:
    print("1. pyazr.priors without the engine")
    C = np.array([[4e-4, 1.2e-4, 0.0], [1.2e-4, 9e-4, -1e-5], [0.0, -1e-5, 1e-4]])
    p = os.path.join(tmp, "p.txt")
    priors.write_prior(p, ["segment_1_norm", "width_1_1", "energy_2"], [1.02, 0.6, 0.81], C, basis=1, comment="test")
    r = priors.read_prior(p)
    check("write_prior / read_prior round trip (names, means, covariance, tag)",
          r["names"] == ["segment_1_norm", "width_1_1", "energy_2"] and np.array_equal(r["mean"], [1.02, 0.6, 0.81])
          and np.array_equal(r["cov"], C) and r["basis"] == 1, r)
    check("write_prior refuses a non-symmetric matrix",
          raises(lambda: priors.write_prior(p + "x", ["a", "b"], [0, 0], [[1, 0.1], [0.2, 1]]), "not symmetric"))
    check("write_prior refuses a non-positive-definite matrix",
          raises(lambda: priors.write_prior(p + "x", ["a", "b"], [0, 0], [[1, 2], [2, 1]]), "not positive definite"))
    check("write_prior refuses a size mismatch",
          raises(lambda: priors.write_prior(p + "x", ["a", "b"], [0], np.eye(2)), "size mismatch"))

    # a fit's output: param.sav (tagged) and covariance.dat with the #columns / #scaled header
    out = os.path.join(tmp, "fit")
    os.makedirs(out)
    with open(os.path.join(out, "param.sav"), "w") as fh:
        fh.write("    #parametrization       1.0000000e+00       0.0000000e+00\n")
        for name, v in [("energy_1", 0.0), ("width_1_1", 0.61595294), ("energy_2", 0.8), ("width_2_1", 0.82062684),
                        ("width_3_1", 0.80801316), ("segment_1_norm", 1.25)]:
            fh.write(f"{name:>20s}{v:20.7e}{0.1:20.7e}\n")
    cov = np.array([[6e-8, -1.5e-8, 1.5e-7], [-1.5e-8, 5e-5, -2.1e-5], [1.5e-7, -2.1e-5, 2.1e-5]])
    with open(os.path.join(out, "covariance.dat"), "w") as fh:
        fh.write("#columns width_1_1 width_2_1 width_3_1\n#scaled 1.5\n")
        for row in cov:
            fh.write(" ".join(f"{v:.17e}" for v in row) + "\n")
    Cr, cols, scale = priors.read_covariance_dat(os.path.join(out, "covariance.dat"))
    check("read_covariance_dat: matrix, #columns, #scaled",
          np.array_equal(Cr, cov) and cols == ["width_1_1", "width_2_1", "width_3_1"] and scale == 1.5, (cols, scale))
    pf = priors.prior_from_fit(os.path.join(out, "param.sav"), os.path.join(out, "covariance.dat"),
                               os.path.join(tmp, "pf.txt"))
    r = priors.read_prior(pf)
    check("prior_from_fit: names = the columns, means = param.sav, covariance = covariance.dat, tag 1",
          r["names"] == cols and np.allclose(r["mean"], [0.61595294, 0.82062684, 0.80801316], rtol=0, atol=0)
          and np.array_equal(r["cov"], cov) and r["basis"] == 1, r)
    check("prior_from_fit notes the scale factor", "(scaled by 1.5)" in open(pf).read())
    r = priors.read_prior(priors.prior_from_fit(os.path.join(out, "param.sav"), os.path.join(out, "covariance.dat"),
                                                os.path.join(tmp, "pf2.txt"), names=["width_3_1", "width_1_1"]))
    check("prior_from_fit names=: the marginal sub-block, in the order asked",
          r["names"] == ["width_3_1", "width_1_1"] and np.array_equal(r["cov"], cov[np.ix_([2, 0], [2, 0])]), r)
    check("prior_from_fit refuses a name the covariance does not hold",
          raises(lambda: priors.prior_from_fit(os.path.join(out, "param.sav"), os.path.join(out, "covariance.dat"),
                                               os.path.join(tmp, "x.txt"), names=["segment_1_norm"]), "not in"))
    plain = os.path.join(out, "plain.dat")
    np.savetxt(plain, cov)
    check("a covariance.dat without #columns is refused without columns=",
          raises(lambda: priors.prior_from_fit(os.path.join(out, "param.sav"), plain, os.path.join(tmp, "x.txt")),
                 "no '#columns' header"))
    r = priors.read_prior(priors.prior_from_fit(os.path.join(out, "param.sav"), plain, os.path.join(tmp, "pf3.txt"),
                                                columns=["width_1_1", "width_2_1", "width_3_1"]))
    check("... and accepted with columns=", r["names"] == ["width_1_1", "width_2_1", "width_3_1"], r["names"])

    m = Project.open(AZR)
    check("Project.config: no priorCovariance without the key", "priorCovariance" not in m.config())
    m.set_config(priorCovariance="prior.txt")
    text = m.text
    check("set: priorCovariance=prior.txt before </config>", "priorCovariance=prior.txt\n</config>" in text)
    m.save(os.path.join(tmp, "k.azr"))
    check("reads back from the file", Project.open(os.path.join(tmp, "k.azr")).config()["priorCovariance"] == "prior.txt")
    m.set_config(priorCovariance=None)
    check("None removes the key, the file as before", m.text == open(AZR).read())
    check("a path with '#' is refused", raises(lambda: m.set_config(priorCovariance="a#b"), "without '#'"))

    print("\n2. with the engine")
    try:
        from pyazr import Session
        have_engine = True
    except Exception as e:                       # pragma: no cover
        print(f"  skip  engine not built ({e})")
        have_engine = False
    if have_engine:
        def project(name, norms, prior_text=None, fixed_levels=False):
            """The norm_groups project with the group D removed, each data set its own norm (NORMS: two texts),
            and the prior file when given; FIXED_LEVELS fixes every level energy and width."""
            work = os.path.join(tmp, name)
            shutil.copytree(os.path.join(PROJECT, "data"), os.path.join(work, "data"))
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))
            pr = Project.open(AZR)
            pr.set("xs_40_s08", norm=norms[0])
            pr.set("xs_90_s09", norm=norms[1])
            pr.remove("D")
            if fixed_levels:
                for lv in pr.levels():
                    pr.set_level(lv, fixed=True)
                    for c in lv.channels:
                        pr.set_channel(lv, c.pair, c.l, c.s, fixed=True)
            if prior_text is not None:
                pr.set_config(priorCovariance="prior.txt")
                open(os.path.join(work, "prior.txt"), "w").write(prior_text)
            pr.save(os.path.join(work, "t.azr"))
            return os.path.join(work, "t.azr"), work

        # norm 1, free; with a 10 % / 5 % prior in SE
        S0 = ["1 fit", "1 fit"]
        SE = ["1 fit prior 1 10%", "1 fit prior 1 5%"]

        Cb = np.array([[4e-4, 1.2e-4], [1.2e-4, 9e-4]])
        azr, work = project("corr", S0, f"{N1} width_1_1\n1.02 0.6\n4e-4 1.2e-4\n1.2e-4 9e-4\n")
        with Session(azr, cwd=work) as mdl:
            fi = {q.name: q.free_index for q in mdl.parameters if not q.fixed}
            x = np.asarray(mdl.params_rwa, float).copy()
            x[fi[N1]], x[fi["width_1_1"]] = 1.07, 0.63
            d = np.array([1.07 - 1.02, 0.63 - 0.6])
            want = float(d @ np.linalg.solve(Cb, d))
            pen = mdl.evaluate(x)
            check("terms['prior covariance'] = d^T C^-1 d (1e-12)",
                  abs(pen.terms["prior covariance"] - want) <= 1e-12 * want, (pen.terms["prior covariance"], want))
            check("objective() includes it",
                  abs(mdl.objective()(x) - (mdl.evaluate(x).chi2 + want)) <= 1e-10 * mdl.objective()(x))
            check("prior_covariance_names", mdl.prior_covariance_names == [N1, "width_1_1"],
                  mdl.prior_covariance_names)
            r, J = mdl.objective().residuals_and_jacobian(x)
            rp = mdl.prior_covariance_residuals(x)
            check("objective().residuals_and_jacobian: the last two rows are L^-1 (p - mu)", np.allclose(r[-2:], rp, rtol=0, atol=1e-14),
                  (r[-2:], rp))
            check("sum r^2 = evaluate().chi2 + the prior (1e-10)",
                  abs(float(r @ r) - (mdl.evaluate(x).chi2 + want)) <= 1e-10 * float(r @ r))
            worst = 0.0
            for col in range(x.size):
                h = 1e-6 * max(1.0, abs(x[col]))
                xp, xm = x.copy(), x.copy()
                xp[col] += h
                xm[col] -= h
                fd = (mdl.prior_covariance_residuals(xp) - mdl.prior_covariance_residuals(xm)) / (2 * h)
                worst = max(worst, float(np.max(np.abs(J[-2:, col] - fd))))
            Linv = np.linalg.inv(np.linalg.cholesky(Cb))
            check("their Jacobian rows = central differences (every column) and L^-1 on the prior's columns",
                  worst < 1e-6 * np.abs(Linv).max() and
                  np.allclose(J[-2:, [fi[N1], fi["width_1_1"]]], Linv, rtol=1e-12, atol=0), worst)

        azr, work = project("diag", S0, f"{N1} {N2}\n1 1\n0.01 0\n0 0.0025\n")
        azr2, work2 = project("normerr", SE)
        with Session(azr, cwd=work) as a, Session(azr2, cwd=work2) as b:
            x = np.asarray(a.params_rwa, float).copy()
            fa = {q.name: q.free_index for q in a.parameters if not q.fixed}
            fb = {q.name: q.free_index for q in b.parameters if not q.fixed}
            check("same free-parameter layout", fa == fb)
            x[fa[N1]], x[fa[N2]] = 1.07, 0.96
            oa, ob = a.objective()(x), b.objective()(x)
            check("diagonal prior on the norms: objective() = the normError priors' (1e-10)",
                  abs(oa - ob) <= 1e-10 * ob, (oa, ob))
            check("... prior term 1.13 = the norm penalties",
                  abs(a.evaluate(x).terms["prior covariance"] - 1.13) < 1e-12 and
                  abs(float(np.sum(np.array([d.norm_chi2 for d in b.evaluate(x).data]))) - 1.13) < 1e-12)
            check("no prior term without the key", b.evaluate(x).terms["prior covariance"] == 0.0 and
                  b.prior_covariance_residuals(x).size == 0)

        # every level fixed: energy_2 is a fixed parameter -> marginalized
        azr, work = project("fixed", S0, f"{N1} energy_2\n1.02 0.81\n4e-4 1.5e-4\n1.5e-4 1e-4\n", fixed_levels=True)
        with Session(azr, cwd=work) as mdl:
            fi = {q.name: q.free_index for q in mdl.parameters if not q.fixed}
            x = np.asarray(mdl.params_rwa, float).copy()
            x[fi[N1]] = 1.07
            check("a fixed listed parameter is marginalized", mdl.prior_covariance_names == [N1],
                  mdl.prior_covariance_names)
            v = mdl.evaluate(x).terms["prior covariance"]
            check("... the term is the marginal (n - mu)^2 / C_nn", abs(v - 0.05 ** 2 / 4e-4) <= 1e-12 * v, v)

        azr, work = project("bad", S0, f"{N1} segment_9_norm\n1 1\n0.01 0\n0 0.01\n")
        try:
            with Session(azr, cwd=work):
                pass
            check("an unknown name makes Session() raise", False, "no error")
        except Exception as e:
            check("an unknown name makes Session() raise", True)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): " + ", ".join(failures))
    sys.exit(1)
print("all checks passed")
