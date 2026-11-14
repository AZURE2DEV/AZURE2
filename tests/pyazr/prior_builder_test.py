#!/usr/bin/env python3
"""The prior builder (``pyazr.priors.build``) and the one reading of the prior file's grammar.

1. The two readers: a table of small prior files, each with the engine's verdict (measured: ``Session`` on the p+p
   model of tests/norm_groups with ``priorCovariance`` naming the file); ``pyazr.priors.read_prior`` gives the same
   verdict on every file (the engine's reading is the definition).  With the engine, the verdicts are measured again.
2. With the engine, on the p+p model (each data set its own free norm, the group D removed):
   the term of a built prior through a session = d^T C^-1 d by hand (1e-10); the "fit" part and ``prior_from_fit``
   give the same names, means, covariance and number lines; the "chain" part = the draws' mean and covariance; every
   refusal (what, allowed, fix); the notes (fixed, a prior of its own, a centre outside the limits, an untagged or
   scaled fit); a parameter with a prior of its own: the engine adds the multivariate term only; write -> read ->
   build round trip; ``use`` sets ``priorCovariance``.

Run from anywhere:  python3 tests/pyazr/prior_builder_test.py
"""
import json
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
    from pyazr import Project, priors           # noqa: E402
    from pyazr.errors import Error, FormatError, ParameterError, ParameterKeyError  # noqa: E402
except Exception as err:
    print(f"skip: pyazr not importable ({type(err).__name__}: {err})")
    sys.exit(77)

PROJECT = os.path.join(ROOT, "tests", "norm_groups")
AZR = os.path.join(PROJECT, "norm_groups.azr")
N1, N2 = "xs_40_s08.norm", "xs_90_s09.norm"
tmp = tempfile.mkdtemp(prefix="pyazr_prior_builder_")

OK = f"{N1} {N2}\n1 1\n0.01 0.002\n0.002 0.0025\n"
# (label, text, the engine's verdict in a Brune run of the p+p model: True accepted).  "run" marks a refusal that
# needs the run (a name the model lacks, the tag of another basis): the reader alone accepts those files.
READERS = [
    ("plain", OK, True),
    ("comments and blank lines", f"# c\n\n{N1} {N2}   # names\n\n1 1 # means\n0.01 0.002\n\n0.002 0.0025\n# e\n", True),
    ("CRLF and tabs", OK.replace("\n", "\r\n").replace(" ", "\t"), True),
    ("tag 1", "#parametrization 1\n" + OK, True),
    ("tag '# parametrization 1'", "# parametrization 1\n" + OK, True),
    ("tag as param.sav writes it", "    #parametrization       1.0000000e+00       0.0000000e+00\n" + OK, True),
    ("tag 0 in a Brune run", "#parametrization 0\n" + OK, "run"),
    ("tag 3", "#parametrization 3\n" + OK, False),
    ("tag 1.5", "#parametrization 1.5\n" + OK, False),
    ("tag without a value", "#parametrization\n" + OK, False),
    ("tag 'one'", "#parametrization one\n" + OK, False),
    ("tag after the names (a comment)", f"{N1} {N2}\n#parametrization 0\n1 1\n0.01 0.002\n0.002 0.0025\n", True),
    ("two tags: the last counts", "#parametrization 0\n#parametrization 1\n" + OK, True),
    ("'#parametrization=0' (a comment)", "#parametrization=0\n" + OK, True),
    ("a name twice", f"{N1} {N1}\n1 1\n0.01 0.002\n0.002 0.0025\n", False),
    ("nan", f"{N1} {N2}\nnan 1\n0.01 0.002\n0.002 0.0025\n", False),
    ("inf", f"{N1} {N2}\n1 1\ninf 0.002\n0.002 0.0025\n", False),
    ("1_0 (Python's float reads it)", f"{N1} {N2}\n1 1\n0.01 0.002\n0.002 0.002_5\n", False),
    ("hexadecimal 0x1p0 (strtod reads it)", f"{N1} {N2}\n0x1p0 1\n0.01 0.002\n0.002 0.0025\n", True),
    ("+1. and .1e1", f"{N1} {N2}\n+1. .1e1\n0.01 0.002\n0.002 0.0025\n", True),
    ("Fortran 1d0", f"{N1} {N2}\n1d0 1\n0.01 0.002\n0.002 0.0025\n", False),
    ("not symmetric", f"{N1} {N2}\n1 1\n0.01 0.002\n0.003 0.0025\n", False),
    ("asymmetric below 1e-10", f"{N1} {N2}\n1 1\n0.01 0.002\n0.0020000000000001 0.0025\n", True),
    ("not positive definite", f"{N1} {N2}\n1 1\n0.01 0.006\n0.006 0.0025\n", False),
    ("pivot below 1e-13 of its diagonal", f"{N1} {N2}\n1 1\n0.01 0.005\n0.005 0.0025000000000001\n", False),
    ("a row too many", OK + "0 0\n", False),
    ("a short row", f"{N1} {N2}\n1 1\n0.01\n0.002 0.0025\n", False),
    ("a mean short", f"{N1} {N2}\n1\n0.01 0.002\n0.002 0.0025\n", False),
    ("names only", f"{N1} {N2}\n", False),
    ("empty", "", False),
    ("a name the model lacks", f"{N1} segment_9_norm\n1 1\n0.01 0\n0 0.01\n", "run"),
    ("'1#c': a comment after a number", f"{N1} {N2}\n1 1#c\n0.01 0.002\n0.002 0.0025\n", True),
    ("a no-break space between names", f"{N1} {N2}\n1 1\n0.01 0.002\n0.002 0.0025\n", False),
    ("a form feed between means", f"{N1} {N2}\n1\f1\n0.01 0.002\n0.002 0.0025\n", True),
]


def parts3(e):
    return bool(e.what and e.allowed and e.fix)


def refused(fn, text, cls=ParameterError):
    """True when fn() raises ``cls`` with its three parts and ``text`` in its message."""
    try:
        fn()
    except cls as e:
        return text in str(e) and parts3(e)
    except Exception as e:  # noqa: BLE001
        print(f"        (raised {type(e).__name__}: {e})")
        return False
    return False


try:
    print("1. the two readers (the engine's verdicts as measured)")
    for label, text, verdict in READERS:
        path = os.path.join(tmp, "reader.txt")
        with open(path, "w", newline="") as fh:
            fh.write(text)
        try:
            priors.read_prior(path)
            python = True
        except FormatError as e:
            python = False
            check(f"   the refusal of '{label}' has its three parts", parts3(e), str(e))
        check(f"read_prior on '{label}': {'accepted' if python else 'refused'}", python == (verdict is not False),
              f"the engine: {verdict}")
    r = priors.read_prior(os.path.join(tmp, "reader.txt"))
    check("read_prior: names and means as written", r["names"] == [N1, N2] and r["mean"].tolist() == [1.0, 1.0], r)

    print("\n2. with the engine")
    try:
        from pyazr import Session
    except Exception as e:                       # pragma: no cover
        print(f"  skip  engine not built ({e})")
        Session = None
    if Session is not None:
        def project(name, norms=("1 fit", "1 fit"), prior_text=None):
            """The norm_groups project with the group D removed, each data set its own norm, and the prior file."""
            work = os.path.join(tmp, name)
            shutil.copytree(os.path.join(PROJECT, "data"), os.path.join(work, "data"))
            os.makedirs(os.path.join(work, "output"))
            os.makedirs(os.path.join(work, "checks"))
            pr = Project.open(AZR)
            pr.set("xs_40_s08", norm=norms[0])
            pr.set("xs_90_s09", norm=norms[1])
            pr.remove("D")
            if prior_text is not None:
                pr.set_config(priorCovariance="prior.txt")
                with open(os.path.join(work, "prior.txt"), "w", newline="") as fh:
                    fh.write(prior_text)
            pr.save(os.path.join(work, "t.azr"))
            return os.path.join(work, "t.azr"), work

        def engine_accepts(text):
            azr, work = project(f"r{engine_accepts.k}", prior_text=text)
            engine_accepts.k += 1
            try:
                with Session(azr, cwd=work):
                    return True
            except Exception:  # noqa: BLE001
                return False
        engine_accepts.k = 0
        devnull = os.open(os.devnull, os.O_WRONLY)
        saved = os.dup(1)
        sys.stdout.flush()
        os.dup2(devnull, 1)       # the engine's ERROR lines
        try:
            measured = [(label, engine_accepts(text), verdict) for label, text, verdict in READERS]
        finally:
            sys.stdout.flush()
            os.dup2(saved, 1)
        for label, got, verdict in measured:
            check(f"the engine on '{label}': {'accepted' if got else 'refused'}", got == (verdict is True),
                  f"table: {verdict}")

        # the session of the builder: norms free with a 10 % prior of their own on the first
        azr, work = project("base", ("1 fit prior 1 10%", "1 fit"))
        with Session(azr, cwd=work) as s:
            P = s.parameters
            check("Parameters.basis = Session.basis (1, Brune)", P.basis == s.basis == 1, (P.basis, s.basis))
            parts = [{"kind": "width", "names": [N1, N2], "percent": 5},
                     {"kind": "width", "names": ["width_1_1"], "sigma": 0.05, "centre": 0.6},
                     {"kind": "correlation", "pairs": [[N1, "width_1_1", 0.4], [N2, N1, -0.2]]}]
            b = priors.build(P, parts, directory=work)
            C = np.array([[0.05 ** 2, -0.2 * 0.05 * 0.05, 0.4 * 0.05 * 0.05],
                          [-0.2 * 0.05 * 0.05, 0.05 ** 2, 0.0],
                          [0.4 * 0.05 * 0.05, 0.0, 0.05 ** 2]])
            check("build: names, centres (the start values, a given centre), covariance",
                  b.names == [N1, N2, "width_1_1"] and b.mean.tolist() == [1.0, 1.0, 0.6]
                  and np.allclose(b.cov, C, rtol=1e-15, atol=0), (b.names, b.mean, b.cov))
            check("build: positive definite, eigenvalues and condition number",
                  b.refusal is None and b.checks["positive_definite"] and b.checks["symmetric"]
                  and abs(b.checks["condition_number"] - np.linalg.cond(C)) < 1e-9 * np.linalg.cond(C), b.checks)
            check("note: a prior of its own is replaced", any(N1 in n and "replaces" in n for n in b.notes), b.notes)
            check("the record prints itself (names, sigmas, correlation, checks, notes)",
                  all(w in str(b) for w in ("width_1_1", "0.05", "correlation:", "-0.200", "condition number",
                                            "note: ")), str(b))
            # round trip: write -> read -> build
            path = b.write(os.path.join(work, "built.txt"))
            r = priors.read_prior(path)
            check("write -> read_prior: names, means, covariance, the run's tag",
                  r["names"] == b.names and np.array_equal(r["mean"], b.mean) and np.array_equal(r["cov"], b.cov)
                  and r["basis"] == 1, r)
            again = priors.build(P, priors.read_parts(path)["parts"], directory=work)
            check("read_parts -> build: the same text", again.text() == b.text())
            # the parameter table: fixed shift -> marginalized note; a centre outside the limits
            n2 = priors.build(P, [{"kind": "width", "names": ["xs_40_s08.shift.offset"], "sigma": 0.01},
                                  {"kind": "width", "names": [N2], "sigma": 0.1, "centre": -1}])
            check("note: a fixed parameter is marginalized by the engine",
                  any("xs_40_s08.shift.offset is fixed" in n and "marginalizes" in n for n in n2.notes), n2.notes)
            check("note: a centre outside the parameter's limits",
                  any(N2 in n and "outside the parameter's limits" in n for n in n2.notes), n2.notes)
            # not positive definite: refusal set, write refused
            bad = priors.build(P, [{"kind": "width", "names": [N1, N2, "width_1_1"], "sigma": 0.1},
                                   {"kind": "correlation", "pairs": [[N1, N2, 0.9], [N1, "width_1_1", 0.9],
                                                                     [N2, "width_1_1", -0.9]]}])
            check("not positive definite: refusal with its three parts, the parts concerned, eigenvalue < 0",
                  bad.refusal and "not positive definite" in bad.refusal["what"] and bad.refusal["allowed"]
                  and bad.refusal["fix"] and bad.refusal_parts == [0, 1] and bad.checks["smallest_eigenvalue"] < 0,
                  (bad.refusal, bad.refusal_parts, bad.checks))
            check("... and write refuses it", refused(lambda: bad.write(os.path.join(work, "x.txt")),
                                                      "the prior is not written"))
            check("... the summary has no text", bad.summary()["text"] is None and b.summary()["text"] == b.text())

            # refusals
            def bld(prts):
                return lambda: priors.build(P, prts, directory=work)
            check("refused: a name the table lacks, the close names said",
                  refused(bld([{"kind": "width", "names": ["width_1_2"], "sigma": 1}]),
                          "did you mean 'width_4_2', 'width_1_1'", ParameterKeyError))
            check("refused: no part", refused(bld([]), "no part"))
            check("refused: a kind that does not exist", refused(bld([{"kind": "widht", "names": [N1]}]),
                                                                  "did you mean 'width'"))
            check("refused: a field the kind lacks", refused(bld([{"kind": "width", "names": [N1], "sigma": 1,
                                                                   "rho": 1}]), "has no field rho"))
            check("refused: a required field missing", refused(bld([{"kind": "width", "sigma": 1}]),
                                                               "needs names"))
            check("refused: sigma and percent both", refused(bld([{"kind": "width", "names": [N1], "sigma": 1,
                                                                   "percent": 2}]), "both sigma and percent"))
            check("refused: neither sigma nor percent", refused(bld([{"kind": "width", "names": [N1]}]),
                                                                "neither sigma nor percent"))
            check("refused: percent of a centre 0", refused(bld([{"kind": "width", "names": ["energy_1"],
                                                                  "percent": 5}]), "percent of a centre 0"))
            check("refused: a number that is not one", refused(bld([{"kind": "width", "names": [N1],
                                                                     "sigma": "abc"}]), "sigma 'abc'"))
            check("refused: a parameter in two parts",
                  refused(bld([{"kind": "width", "names": [N1], "sigma": 1},
                               {"kind": "width", "names": [N1], "sigma": 2}]), "is in part 1 (width) and in part 2"))
            check("refused: a correlation of a parameter in no part",
                  refused(bld([{"kind": "width", "names": [N1], "sigma": 1},
                               {"kind": "correlation", "pairs": [[N1, N2, 0.1]]}]), f"{N2} is in no other part"))
            check("refused: rho of magnitude 1", refused(bld([{"kind": "width", "names": [N1, N2], "sigma": 1},
                                                              {"kind": "correlation", "pairs": [[N1, N2, 1]]}]),
                                                         "rho 1 between"))
            check("refused: a pair twice", refused(bld([{"kind": "width", "names": [N1, N2], "sigma": 1},
                                                        {"kind": "correlation", "pairs": [[N1, N2, .1],
                                                                                          [N2, N1, .2]]}]),
                                                   "twice"))
            check("refused: only correlation parts", refused(bld([{"kind": "correlation", "pairs": [[N1, N2, .1]]}]),
                                                             "only correlation parts"))
            check("refused: a correlation matrix that is not one",
                  refused(bld([{"kind": "width", "names": [N1, N2], "sigma": 1},
                               {"kind": "correlation", "names": [N1, N2], "matrix": [[1, .3], [.2, 1]]}]),
                          "not a correlation matrix"))
            m = priors.build(P, [{"kind": "width", "names": [N1, N2], "sigma": 1},
                                 {"kind": "correlation", "names": [N1, N2], "matrix": [[1, .3], [.3, 1]]}])
            check("a correlation matrix = its pairs", m.cov[0, 1] == m.cov[1, 0] == 0.3, m.cov)
            check("use: an absolute path refused", refused(lambda: b.use(Project.open(azr), "/tmp/x.txt"),
                                                           "an absolute path"))

            # an earlier fit: the part and prior_from_fit
            fit = os.path.join(work, "fitA")
            os.makedirs(fit)
            cols = ["width_1_1", "width_2_1", "width_3_1"]
            with open(os.path.join(fit, "param.sav"), "w") as fh:
                fh.write("    #parametrization       1.0000000e+00       0.0000000e+00\n")
                for name in P.names:
                    v = {"width_1_1": 0.61595294, "width_2_1": 0.82062684, "width_3_1": 0.80801316}.get(
                        name, P.by_name(name).start)
                    fh.write(f"{name:>28s}{v:20.7e}{0.1:20.7e}\n")
            covA = np.array([[6e-8, -1.5e-8, 1.5e-7], [-1.5e-8, 5e-5, -2.1e-5], [1.5e-7, -2.1e-5, 2.1e-5]])
            with open(os.path.join(fit, "covariance.dat"), "w") as fh:
                fh.write("#columns " + " ".join(cols) + "\n#scaled 1.5\n")
                for row in covA:
                    fh.write(" ".join(f"{v:.17e}" for v in row) + "\n")
            pf = priors.prior_from_fit(os.path.join(fit, "param.sav"), os.path.join(fit, "covariance.dat"),
                                       os.path.join(work, "pf.txt"), names=["width_3_1", "width_1_1"])
            fb = priors.build(P, [{"kind": "fit", "sav": "fitA/param.sav", "covariance": "fitA/covariance.dat",
                                   "names": ["width_3_1", "width_1_1"]}], directory=work)
            fb.write(os.path.join(work, "fb.txt"))
            a, c = priors.read_prior(pf), priors.read_prior(os.path.join(work, "fb.txt"))
            check("the part fit = prior_from_fit: names, means, covariance, tag (bit for bit)",
                  a["names"] == c["names"] and np.array_equal(a["mean"], c["mean"]) and
                  np.array_equal(a["cov"], c["cov"]) and a["basis"] == c["basis"] == 1, (a, c))

            def numbers(p):
                return [ln for ln in open(p).read().splitlines() if not ln.startswith("#")]
            check("... and the same number lines", numbers(pf) == numbers(os.path.join(work, "fb.txt")))
            check("refused: a correlation inside a part fit",
                  refused(bld([{"kind": "fit", "sav": "fitA/param.sav", "covariance": "fitA/covariance.dat"},
                               {"kind": "correlation", "pairs": [["width_1_1", "width_2_1", 0.1]]}]),
                          "already correlated by part 1 (fit)"))
            check("note: the fit's covariance was scaled", any("scaled by 1.5" in n for n in fb.notes), fb.notes)
            check("refused: a name the fit's covariance lacks",
                  refused(bld([{"kind": "fit", "sav": "fitA/param.sav", "covariance": "fitA/covariance.dat",
                                "names": [N1]}]), "not in"))
            with open(os.path.join(fit, "param0.sav"), "w") as fh:
                fh.write(open(os.path.join(fit, "param.sav")).read().replace("1.0000000e+00", "0.0000000e+00", 1))
            check("refused: a fit of another basis",
                  refused(bld([{"kind": "fit", "sav": "fitA/param0.sav", "covariance": "fitA/covariance.dat"}]),
                          "holds parametrization 0 (standard) amplitudes and this run uses parametrization 1"))
            with open(os.path.join(fit, "untagged.sav"), "w") as fh:
                fh.write("".join(open(os.path.join(fit, "param.sav")).readlines()[1:]))
            u = priors.build(P, [{"kind": "fit", "sav": "fitA/untagged.sav", "covariance": "fitA/covariance.dat"}],
                             directory=work)
            check("note: a fit without a tag is taken in the run's basis",
                  any("no #parametrization tag" in n for n in u.notes), u.notes)

            # a chain as the sampler job writes it
            rng = np.random.default_rng(3)
            names = P.free.names
            chain = 1.0 + 0.01 * rng.standard_normal((40, 6, len(names)))
            np.save(os.path.join(work, "chain_7.npy"), chain)
            with open(os.path.join(work, "chain_7.json"), "w") as fh:
                json.dump({"names": names, "discard": 10}, fh)
            ch = priors.build(P, [{"kind": "chain", "file": "chain_7.npy", "names": [N2, N1]}], directory=work)
            flat = chain[10:].reshape(-1, len(names))[:, [names.index(N2), names.index(N1)]]
            check("the part chain: the draws' mean and covariance after the JSON's burn-in",
                  np.array_equal(ch.mean, flat.mean(axis=0)) and np.array_equal(ch.cov, np.cov(flat, rowvar=False)),
                  (ch.mean, ch.cov))
            ch5 = priors.build(P, [{"kind": "chain", "file": "chain_7.npy", "names": [N1], "discard": 5}],
                               directory=work)
            check("... with a burn-in given", np.array_equal(ch5.mean, chain[5:, :, names.index(N1)].reshape(-1)
                                                             .mean(keepdims=True)), ch5.mean)
            check("refused: a chain without its JSON",
                  refused(lambda: (os.rename(os.path.join(work, "chain_7.json"), os.path.join(work, "c.json")),
                                   bld([{"kind": "chain", "file": "chain_7.npy"}])())[1], "missing or unreadable"))
            os.rename(os.path.join(work, "c.json"), os.path.join(work, "chain_7.json"))
            check("refused: a name the chain lacks",
                  refused(bld([{"kind": "chain", "file": "chain_7.npy", "names": ["xs_40_s08.shift.offset"]}]),
                          "not in the chain"))
            check("refused: fewer draws than parameters",
                  refused(bld([{"kind": "chain", "file": "chain_7.npy", "discard": 39}]), "more draws than"))

        # the built prior through a session: its term = d^T C^-1 d, the prior of its own replaced
        proj = Project.open(azr)
        b.use(proj, "built_prior.txt")
        check("use: priorCovariance set", proj.config().get("priorCovariance") == "built_prior.txt", proj.config())
        proj.save()
        with Session(azr, cwd=work) as s:
            fi = {q.name: q.free_index for q in s.parameters if not q.fixed}
            x = s.parameters.x0.copy()
            x[fi[N1]], x[fi[N2]], x[fi["width_1_1"]] = 1.07, 0.96, 0.63
            d = np.array([1.07 - 1.0, 0.96 - 1.0, 0.63 - 0.6])
            want = float(d @ np.linalg.solve(C, d))
            e = s.evaluate(x)
            got = e.terms["prior covariance"]
            check("a built prior through a session: the term = d^T C^-1 d (1e-10)", abs(got - want) <= 1e-10 * want,
                  (got, want))
            check("the prior of its own (10 % on the first norm) is not added: objective = data + the term",
                  abs(s.objective()(x) - (e.chi2 + want)) <= 1e-10 * s.objective()(x)
                  and all(dd.norm_chi2 == 0.0 for dd in e.data), (s.objective()(x), e.chi2, want))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
if failures:
    print(f"FAILED: {len(failures)} check(s): " + ", ".join(failures))
    sys.exit(1)
print("all checks passed")
