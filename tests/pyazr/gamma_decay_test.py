#!/usr/bin/env python3
"""A level's radiative decay as Gamma_gamma, branchings and mixing ratios in pyazr (the gamma_decay record,
docs/upgrade/gamma_parametrisation.md), on the projects of tests/gamma_decay/models:

  1. the record built in Python (Project.add on the project with level 5's radiative widths at 0) is the record of
     models/record.azr: the same parameters and the same chi^2;
  2. Gauss-Newton fits of the five data sets in both parametrisations (scipy's least_squares on residual_jacobian):
     the same chi^2 minimum, and the errors of the branchings from the record's covariance (J^T J)^-1 equal the errors
     propagated from the partial widths' covariance to b_f = Gamma_f / Gamma_gamma, to the linearisation;
  3. an agent's use: read the level's record, change one branching, see the chi^2 rise, refit, and the branching
     comes back;
  4. external capture into the ground state and 6.1763 MeV (E1) and 6.7931 MeV (M1) on the record's channels
     (the projects of tests/gamma_decay/check.sh (h)): the record's parameters read through Session.parameters, the
     same chi^2 as the twin with partial widths, and a branching set in the free vector gives the chi^2 of the
     project written with that branching;
  5. the record's numbers set in the physical vector (calculate_chi2, update_segments) reach the engine: the chi^2 of
     the project written with them, and of the twin with its radiative widths scaled alike.

Needs the compiled engine and scipy; skips cleanly without them.  Run from anywhere:
python3 tests/pyazr/gamma_decay_test.py
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
MODELS = os.path.join(ROOT, "tests", "gamma_decay", "models")
failures = []


def check(name, ok, detail=""):
    print(("  ok    " if ok else "  FAIL  ") + name + ("" if ok else f"  ({detail})"))
    if not ok:
        failures.append(name)


try:
    sys.path.insert(0, ROOT)
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    import numpy as np
    from scipy.optimize import least_squares
    from pyazr import azure2, Project, Nuisance
except Exception as err:  # pragma: no cover
    print(f"skip: engine or scipy not available ({type(err).__name__}: {err})")
    sys.exit(77)

# the final states of level 5 in its channel order, and the remainder (models/make_record.py)
FINALS = ["g+15O", "g+15O[6.7931]", "g+15O[5.183]", "g+15O[6.1763]", "g+15O[6.8594]"]
REMAINDER = "g+15O[5.183]"


def fit(s):
    """A Gauss-Newton fit of the free parameters (bounds from the kinds): the minimum, chi^2 and (J^T J)^-1."""
    params = [p for p in s.parameters if not p.fixed]
    lo = [0.0 if p.kind == "radiative" and p.decay_part in ("branching", "total_width") else -np.inf for p in params]
    hi = [1.0 if p.kind == "radiative" and p.decay_part == "branching" else np.inf for p in params]
    x0 = np.clip(np.asarray(s.params_rwa, float), np.asarray(lo) + 1e-12, np.asarray(hi) - 1e-12)
    r = least_squares(lambda x: s.residuals(x), x0, jac=lambda x: s.residual_jacobian(x)[1], bounds=(lo, hi),
                      x_scale="jac", xtol=1e-14, ftol=1e-14, gtol=1e-14)
    _, J = s.residual_jacobian(r.x)
    return r.x, float(np.sum(s.residuals(r.x) ** 2)), np.linalg.inv(J.T @ J), [p.name for p in params]


with tempfile.TemporaryDirectory() as tmp:
    d = os.path.join(tmp, "m")
    shutil.copytree(MODELS, d)
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)
    record = os.path.join(d, "record.azr")

    print("1. the record built in Python")
    built = os.path.join(d, "built.azr")
    proj = Project.open(record)
    proj.remove("level5")
    proj.save(built)
    proj = Project.open(built)
    rec = Project.open(record)["level5"]
    proj.add("gamma_decay", "level5", level=5, J="1/2+",
             total_width=Nuisance(0.0385027, "fit", unit="eV"),
             branching=[f"{t}" for t in rec["branching"]], remainder=REMAINDER)
    proj.save(built)
    with azure2(record, cwd=d) as a, azure2(built, cwd=d) as b:
        na = [p.name for p in a.parameters]
        nb = [p.name for p in b.parameters]
        check("the same parameters, by name and kind", na == nb and
              [p.kind for p in a.parameters] == [p.kind for p in b.parameters])
        ca, cb = a.objective(), b.objective()
        check("the same chi^2", abs(ca - cb) <= 1e-12 * abs(ca), f"{ca} {cb}")
        radiative = [p for p in a.parameters if p.kind == "radiative"]
        check("8 radiative numbers: Gamma_gamma, five branchings (one fixed at 0), two mixing ratios",
              [p.decay_part for p in radiative].count("branching") == 5 and len(radiative) == 8,
              [p.name for p in radiative])
        tw = [p for p in radiative if p.decay_part == "total_width"][0]
        check("Gamma_gamma's parameter: eV, the record's value", abs(tw.value - 0.0385027) < 1e-15, tw.value)

    print("2. fits in both parametrisations, and the errors of the branchings")
    with azure2(os.path.join(d, "partial.azr"), cwd=d) as s:
        xp, chi_p, cov_p, names_p = fit(s)
        # b_f of the partial widths (eV) of the five free photon channels, in channel order = FINALS
        def branchings(x):
            w = np.abs(np.asarray(s.transform_rwa(x), float))
            photon = [k for k, n in enumerate(names_p) if n.startswith("width_")]
            g = w[photon]
            return g / g.sum()
        b_p = branchings(xp)
        G = np.zeros((len(b_p), len(xp)))
        for k in range(len(xp)):
            h = 1e-6 * max(abs(xp[k]), 1e-3)
            up, dn = xp.copy(), xp.copy()
            up[k] += h
            dn[k] -= h
            G[:, k] = (branchings(up) - branchings(dn)) / (2 * h)
        sig_p = np.sqrt(np.diag(G @ cov_p @ G.T))
    with azure2(record, cwd=d) as s:
        xr, chi_r, cov_r, names_r = fit(s)
        idx = {n: k for k, n in enumerate(names_r)}
        explicit = [f for f in FINALS if f != REMAINDER]
        cols = [idx[f"level5.branching.{f}"] for f in explicit]
        b_r = {f: xr[c] for f, c in zip(explicit, cols)}
        b_r[REMAINDER] = 1.0 - sum(b_r.values())
        sig_r = {f: np.sqrt(cov_r[c, c]) for f, c in zip(explicit, cols)}
        sig_r[REMAINDER] = np.sqrt(np.sum(cov_r[np.ix_(cols, cols)]))  # 1 - sum b: the sum of the block
    check("the same chi^2 minimum", abs(chi_p - chi_r) <= 1e-8 * chi_p, f"{chi_p} {chi_r}")
    for k, f in enumerate(FINALS):
        check(f"b({f}): the same value", abs(b_p[k] - b_r[f]) <= 1e-6 * b_p[k] + 1e-12, f"{b_p[k]} {b_r[f]}")
        check(f"b({f}): its error {sig_r[f]:.4e} = the propagated {sig_p[k]:.4e}",
              abs(sig_p[k] - sig_r[f]) <= 1e-4 * sig_p[k], f"{sig_p[k]} {sig_r[f]}")

    print("3. an agent's use: read the record, change one branching, refit")
    proj = Project.open(record)
    r = proj["level5"]
    lines = list(r["branching"])
    check("the record reads back as written", any(l.startswith("g+15O[6.7931] 0.2303") for l in lines), lines)
    changed = [("g+15O[6.7931] 0.35 fit" if l.startswith("g+15O[6.7931]") else l) for l in lines]
    proj.set("level5", branching=changed)
    proj.save(os.path.join(d, "agent.azr"))
    with azure2(os.path.join(d, "agent.azr"), cwd=d) as s:
        p = [q for q in s.parameters if q.name == "level5.branching.g+15O[6.7931]"][0]
        check("the session starts from the new branching", abs(p.value - 0.35) < 1e-15, p.value)
        chi0 = s.objective()
        check("the chi^2 rises", chi0 > chi_r + 100, f"{chi0} vs {chi_r}")
        xa, chi_a, _, names_a = fit(s)
        ba = xa[names_a.index("level5.branching.g+15O[6.7931]")]
        check("the refit brings it back", abs(ba - b_r["g+15O[6.7931]"]) <= 1e-6, f"{ba} {b_r['g+15O[6.7931]']}")
        check("the refit's chi^2 is the minimum", abs(chi_a - chi_r) <= 1e-8 * chi_r, f"{chi_a} {chi_r}")

    print("4. external capture on the record's channels")

    def with_capture(src, dst):
        """Column 31 of the photon lines: capture into pairs 2 and 6 by E1, into pair 3 by M1 (check.sh (h))."""
        out = []
        for line in open(src).read().split("\n"):
            t = line.split()
            if len(t) > 30 and t[26] == "10" and t[5] in ("2", "3", "6"):
                t[30] = "2" if t[5] == "3" else "1"
                line = " " + " ".join(t)
            out.append(line)
        open(dst, "w").write("\n".join(out))

    import subprocess
    twin, rec = os.path.join(d, "partial_ec.azr"), os.path.join(d, "record_ec.azr")
    with_capture(os.path.join(d, "partial.azr"), twin)
    subprocess.run([sys.executable, os.path.join(MODELS, "make_record.py"), twin, rec, "5", "5.183", "fit"], check=True)
    with azure2(twin, cwd=d) as t, azure2(rec, cwd=d) as s:
        radiative = {p.name: p for p in s.parameters if p.kind == "radiative"}
        check("the record's parameters by name, Gamma_gamma in eV",
              "level5.branching.g+15O[6.1763]" in radiative
              and abs(radiative["level5.total_width"].value - 0.0385027) < 1e-15, sorted(radiative))
        ct, cr = t.objective(), s.objective()
        check("the same chi^2 as the twin with partial widths", abs(ct - cr) <= 1e-12 * abs(ct), f"{ct} {cr}")
        p = radiative["level5.branching.g+15O[6.1763]"]
        x = np.asarray(s.params_rwa, float).copy()
        x[p.free_index] = 0.55
        c_set = s.objective(x)
        check("a branching set in the free vector reads back", abs(s.transform_rwa(x)[p.free_index] - 0.55) < 1e-15)
    proj = Project.open(rec)
    lines = [("g+15O[6.1763] 0.55 fit" if l.startswith("g+15O[6.1763]") else l) for l in proj["level5"]["branching"]]
    proj.set("level5", branching=lines)
    proj.save(os.path.join(d, "record_ec_055.azr"))
    with azure2(os.path.join(d, "record_ec_055.azr"), cwd=d) as s:
        c_file = s.objective()
    check("its chi^2 is that of the project written with it (gamma^ext at the set point)",
          abs(c_set - c_file) <= 1e-10 * abs(c_file) and abs(c_set - cr) > 1.0, f"{c_set} {c_file} {cr}")

    print("5. the record's numbers set in the physical vector (calculate_chi2, update_segments)")
    # docs/upgrade/gamma_parametrisation.md section 2.8: the physical vector carries the record's numbers; set there,
    # they must reach the engine (before, FillMnParams took the file's numbers and the vector's were ignored)
    for twin_f, rec_f in ((os.path.join(d, "partial.azr"), record), (twin, rec)):
        label = "with external capture" if rec_f == rec else "without external capture"
        scaled = os.path.join(d, "x12_" + os.path.basename(rec_f))
        proj = Project.open(rec_f)
        proj.set("level5", total_width=Nuisance(0.0385027 * 1.2, "fit", unit="eV"))
        proj.save(scaled)
        with azure2(twin_f, cwd=d) as t, azure2(rec_f, cwd=d) as s, azure2(scaled, cwd=d) as f:
            names = [p.name for p in s.parameters if not p.fixed]
            i = names.index("level5.total_width")
            phys = np.asarray(s.sess.params_values(), float)
            check(f"{label}: the physical vector holds Gamma_gamma (eV)", phys[i] == 0.0385027, phys[i])
            x = phys.copy()
            x[i] *= 1.2
            c_vec, c_file = s.calculate_chi2(x)[0], f.calculate_chi2_rwa(f.params_rwa)[0]
            # the twin: level 5's radiative widths x 1.2 in its physical vector
            tn = [p.name for p in t.parameters if not p.fixed]
            xt = np.asarray(t.sess.params_values(), float).copy()
            k5 = [k for k, n in enumerate(tn) if n.startswith("width_5_") and abs(xt[k]) < 1.0]
            xt[k5] *= 1.2
            c_twin = t.calculate_chi2(xt)[0]
            check(f"{label}: calculate_chi2 with Gamma_gamma x 1.2 is the chi^2 of the project written with it",
                  abs(c_vec - c_file) <= 1e-12 * c_file and abs(c_vec - s.calculate_chi2(phys)[0]) > 1.0,
                  f"{c_vec} {c_file}")
            check(f"{label}: ... and the twin's with its {len(k5)} radiative widths x 1.2",
                  len(k5) == 5 and abs(c_vec - c_twin) <= 1e-12 * c_twin, f"{c_vec} {c_twin}")
            n = int(s.sess.update_segments(list(x)))
            y = np.concatenate([np.asarray(s.sess.calculated_segments(k), float).ravel() for k in range(n)])
            nf = int(f.sess.update_segments(list(np.asarray(f.sess.params_values(), float))))
            yf = np.concatenate([np.asarray(f.sess.calculated_segments(k), float).ravel() for k in range(nf)])
            check(f"{label}: update_segments with Gamma_gamma x 1.2 gives the written project's cross sections",
                  n == nf and np.array_equal(y, yf), np.max(np.abs(y - yf) / np.abs(yf)) if n == nf else (n, nf))

if failures:
    print(f"FAILED: {len(failures)} check(s)")
    sys.exit(1)
print("all checks passed")
