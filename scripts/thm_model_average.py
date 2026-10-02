#!/usr/bin/env python3
"""Fit a grid of THM model variants one after the other and average them.

A THM result depends on choices the data cannot make -- the channel radius,
the entrance vertex (plane or distorted waves), the vertex representation,
the optical potentials, the spectator-momentum window.  This driver takes a
project and a small variant grid, writes one project per variant (through
AzrModel, so every refusal AzrModel or the engine knows applies), refits
each with scipy's least_squares on pyazr's residuals and Jacobian, and hands
the results to pyazr's model averaging (``pyazr/modelavg.py``), which keeps
the statistical and the model spread apart.  See
docs/source/theory/thm_implementation.rst, "Model averaging".

One engine session at a time, by design: a THM session of a realistic model
takes 0.4-2 GB, and the variants are fitted sequentially in this process.

Variant axes (a product over every axis given; one value = not varied):

  --radius-pairs K [K ...] --radii R [R ...]
        set the channel radius R (fm) of the file pair keys K (all of them to
        the same value)
  --vertex-model pw dw           vertexModel= of the THM experiment(s)
  --vertex constant perlevel onshell     the <thm> vertex= option
  --optical SPEC [SPEC ...]      the distortion of a dw vertex: "coulomb", or
        "AA/SF" with AA, SF global potential names (ancai06, daehnick80, kd03,
        bg71, liang09, mcfadden66, avrigeanu94, ":extrapolate" allowed), or
        "plane"/"coulomb"; ten numbers V,R,a,W,RW,aW,WD,RD,aD,RC also work.
        Applies to vertexModel=dw variants; pw variants keep the project's
        distortion unless --pw-distortion (then R(E) is applied to them).
  --ps W [W ...]                 the spectator-momentum window: "delta" (none)
        or e.g. "hulthen:0-50"
  --experiment NAME              the THM experiment(s) edited (default all)

or a JSON spec (``--spec file.json``) with the same keys: {"radius_pairs":
[1], "radii": [...], "vertex_model": [...], "vertex": [...], "optical":
[...], "ps": [...], "experiment": [...], "priors": {"radius=6.1": 0.5}}.
CLI flags override the spec.

Derived quantities (optional): ``--strength NAME=JPI@E`` (e.g.
``213=2-@13.057``: the level of that J^pi nearest E, excitation energy in
MeV) gives omega*gamma = (2J+1)/((2j1+1)(2j2+1)) G_in G_out / G in eV, with
``--strength-in`` / ``--strength-out`` the file pair keys of the entrance and
exit channels (G sums the open channels of the level; physical widths, so
channels entered as reduced-width amplitudes are refused).  Their
statistical errors are propagated from the fit covariance.

Outputs in --out: variants.csv (one row per variant, fitted or skipped, with
the reason), variants.json, average.json/.csv, summary.txt, and
averaged/<project>_avg.azr -- the averaged parameter means written into a copy
of the project (its own radii and THM settings), with the project's data
beside it so it opens in the GUI.  work/ keeps each variant's project and its
fitted snapshot (<label>_fit.azr / .sav).

Usage:
    python3 scripts/thm_model_average.py examples/f19_pag_thm/f19_pag_thm.azr \\
        --out modelavg --radius-pairs 1 4 5 --radii 4.1 5.1 6.1 \\
        --vertex-model pw dw --optical ancai06/kd03:extrapolate \\
        --strength 213=2-@13.057 --strength-in 1 --strength-out 6
    python3 scripts/thm_model_average.py ... --dry-run   # list the variants

--dry-run needs neither numpy nor the engine.
"""

import argparse
import csv
import importlib.util
import itertools
import json
import math
import os
import shutil
import sys
import time
import traceback
import warnings

os.environ.setdefault("OMP_NUM_THREADS", "1")    # before any numpy import

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, ROOT)                          # this checkout's pyazr


def _load(name):
    """pyazr/<name>.py without the package (no numpy needed)."""
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(ROOT, "pyazr", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


AzrModel = _load("azrfile").AzrModel
ma = _load("modelavg")


# ---------------------------------------------------------------------------
#  The variant grid
# ---------------------------------------------------------------------------

AXES = ("radius", "vertex_model", "vertex", "optical", "ps")


def _fmt_radius(r):
    return f"{float(r):g}"


def build_spec(args):
    spec = {}
    if args.spec:
        with open(args.spec) as fh:
            spec = json.load(fh)
    for key in ("radius_pairs", "radii", "vertex_model", "vertex", "optical",
                "ps", "experiment"):
        v = getattr(args, key)
        if v:
            spec[key] = v
    if spec.get("radii") and not spec.get("radius_pairs"):
        raise SystemExit("--radii needs --radius-pairs (the file pair keys).")
    return spec


def variant_grid(spec, pw_distortion=False):
    """[(label, settings)] over the product of the axes given."""
    axes = []
    if spec.get("radii"):
        axes.append(("radius", [_fmt_radius(r) for r in spec["radii"]]))
    for key in ("vertex_model", "vertex", "optical", "ps"):
        if spec.get(key):
            axes.append((key, [str(v) for v in spec[key]]))
    if not axes:
        return [("project", {})]
    out, seen = [], set()
    for combo in itertools.product(*[vals for _, vals in axes]):
        s = dict(zip([k for k, _ in axes], combo))
        # the optical axis only means something for a dw vertex (or with
        # --pw-distortion); collapse it otherwise
        if ("optical" in s and s.get("vertex_model", "pw") != "dw"
                and not pw_distortion):
            s["optical"] = "-"
        key = tuple(sorted(s.items()))
        if key in seen:
            continue
        seen.add(key)
        out.append((label_of(s), s))
    return out


def label_of(s):
    bits = []
    if "radius" in s:
        bits.append(f"r{s['radius']}")
    if "vertex_model" in s:
        bits.append(s["vertex_model"])
    if "vertex" in s:
        bits.append(s["vertex"])
    if s.get("optical", "-") != "-":
        bits.append(s["optical"].replace("/", "+").replace(":", "~")
                    .replace(",", "_"))
    if "ps" in s:
        bits.append("ps-" + s["ps"].replace(":", "~").replace(",", "_"))
    return "_".join(bits) or "project"


def apply_variant(model, s, spec):
    """Edit an AzrModel into the variant.  Raises ValueError/KeyError for a
    variant AzrModel refuses (the model is then not usable)."""
    if "radius" in s:
        for k in spec["radius_pairs"]:
            model.set_channel_radius(int(k), float(s["radius"]))
    if "vertex" in s:
        model.set_thm_option("vertex", s["vertex"])
    exps = model.thm_experiments()
    want = spec.get("experiment") or list(exps)
    if any(k in s for k in ("vertex_model", "optical", "ps")):
        if not exps:
            raise ValueError("the project has no THM experiment[...] line to "
                             "carry vertexModel/optical/ps.")
        for name in want:
            if name not in exps:
                raise KeyError(f"no THM experiment {name!r} (have {list(exps)})")
            rec = dict(exps[name])
            if "vertex_model" in s:
                rec["vertexModel"] = None if s["vertex_model"] == "pw" else s["vertex_model"]
            if s.get("optical", "-") != "-":
                o = s["optical"]
                for k in ("distortionRef", "distortionRatio"):
                    rec.pop(k, None)
                if o == "coulomb":
                    rec["distortion"] = "coulomb"
                    rec.pop("opticalAA", None)
                    rec.pop("opticalSF", None)
                else:
                    if "/" not in o:
                        raise ValueError(f"--optical {o!r}: give 'coulomb' or "
                                         "'AA/SF' (a+A and s+F potentials).")
                    aa, sf = o.split("/", 1)
                    rec.update(distortion="optical", opticalAA=aa, opticalSF=sf)
            if "ps" in s:
                if s["ps"] in ("delta", "none"):
                    rec.pop("ps", None)
                    rec.pop("psNodes", None)
                else:
                    rec["ps"] = s["ps"]
            seg = rec.pop("segments")
            model.set_thm_experiment(name, seg, **rec)
    return model


# ---------------------------------------------------------------------------
#  Fitting one variant (numpy, scipy and the engine from here on)
# ---------------------------------------------------------------------------

def penalty_rows(np, m, x):
    """AZURE2's norm and shift priors as signed residual rows and their
    Jacobian (pyazr.azure2.penalties gives them squared)."""
    rows, jac = [], []
    n = len(x)
    norms = {p.segment_key: p for p in m.parameters.norms}
    shifts = {p.segment_key: p for p in m.parameters.shifts}
    for d in m.active_datasets:
        p = norms.get(d.key)
        if p is not None and not p.fixed and p.free_index is not None:
            sig = d.norm / 100.0 * d.norm_error
            if sig:
                rows.append((x[p.free_index] - d.norm) / sig)
                j = np.zeros(n)
                j[p.free_index] = 1.0 / sig
                jac.append(j)
        q = shifts.get(d.key)
        if (d.vary_shift and d.energy_shift_error and q is not None
                and not q.fixed and q.free_index is not None):
            rows.append((x[q.free_index] - d.energy_shift) / d.energy_shift_error)
            j = np.zeros(n)
            j[q.free_index] = 1.0 / d.energy_shift_error
            jac.append(j)
    return np.asarray(rows, float), (np.asarray(jac) if jac else np.zeros((0, n)))


def profiled_count(m, x):
    """Parameters eliminated in closed form (THM norms and backgrounds): they
    are free parameters for an information criterion, though not in x."""
    k = 0
    in_exp = set()
    for r in m.thm_experiments(x).values():
        k += 1 + len(r["b"])
        in_exp.update(r["segments"])
    for d in m.active_datasets:
        if d.thm and d.vary_norm and d.key not in in_exp:
            k += 1
    return k


def parse_strengths(items):
    out = []
    for it in items or []:
        try:
            name, rest = it.split("=", 1)
            jpi, e = rest.split("@", 1)
            J = jpi[:-1]
            J = float(J.split("/")[0]) / float(J.split("/")[1]) if "/" in J else float(J)
            parity = +1 if jpi.endswith("+") else -1 if jpi.endswith("-") else None
            if parity is None:
                raise ValueError
            out.append((name, J, parity, float(e)))
        except ValueError:
            raise SystemExit(f"--strength {it!r}: expected NAME=JPI@E, e.g. 213=2-@13.057")
    return out


class StrengthCalc:
    """omega*gamma of chosen levels from the full physical parameter vector."""

    def __init__(self, m, strengths, pin, pout, tol=0.01):
        self.items = []
        pairs = {p.number: p for p in m.pairs}
        key_of = {p.number: p.key for p in m.pairs}
        params = list(m.parameters)
        levels = {}
        for p in params:
            if p.kind in ("energy", "width"):
                levels.setdefault((p.jgroup, p.level), []).append(p)
        for name, J, parity, E in strengths:
            best = None
            for (jg, lv), ps in levels.items():
                p0 = ps[0]
                if abs(p0.J - J) > 1e-6 or p0.parity != parity:
                    continue
                e = next((q.value for q in ps if q.kind == "energy"), None)
                if e is None:
                    continue
                if abs(e - E) <= tol and (best is None or abs(e - E) < best[0]):
                    best = (abs(e - E), (jg, lv))
            if best is None:
                raise ValueError(f"strength {name}: no {J:g}{'+' if parity > 0 else '-'} "
                                 f"level within {tol * 1e3:g} keV of {E} MeV")
            ps = levels[best[1]]
            eidx = next(q.index for q in ps if q.kind == "energy")
            widths = []
            for q in ps:
                if q.kind != "width":
                    continue
                pr = pairs[q.pair]
                if q.input_is_rwa:
                    raise ValueError(f"strength {name}: channel {q.name} is entered "
                                     "as a reduced-width amplitude; omega*gamma "
                                     "needs partial widths.")
                widths.append((q.index, key_of[q.pair], pr))
            jin = next((pr for _, k, pr in widths if k in pin), None)
            if jin is None:
                raise ValueError(f"strength {name}: the level has no channel of "
                                 f"pair(s) {pin}")
            omega = (2 * J + 1) * jin.i1i2factor
            self.items.append((name, eidx, widths, omega, set(pin), set(pout)))

    def __call__(self, allphys):
        out = {}
        for name, eidx, widths, omega, pin, pout in self.items:
            ex = allphys[eidx]
            gin = gout = tot = 0.0
            for idx, key, pr in widths:
                is_open = pr.is_photon or ex > pr.sep_energy + pr.excitation
                if not is_open:
                    continue
                g = abs(float(allphys[idx]))
                tot += g
                if key in pin:
                    gin += g
                if key in pout:
                    gout += g
            out[f"wg({name})"] = omega * gin * gout / tot if tot > 0 else 0.0
        return out


def fit_variant(path, cwd, args, strengths, log):
    import numpy as np
    from scipy.optimize import least_squares
    from pyazr import azure2

    t0 = time.time()
    with azure2(path, cwd=cwd) as m:
        x0 = np.asarray(m.params_rwa, float)
        nres = len(m.residuals(x0))
        calc = StrengthCalc(m, strengths, set(args.strength_in or []),
                            set(args.strength_out or [])) if strengths else None
        nev = [0]

        def resid(x):
            try:
                r = np.asarray(m.residuals(x), float)
            except Exception as err:           # keep the optimiser alive
                log(f"    residuals failed: {err}")
                r = np.full(nres, 1e3)
            pr, _ = penalty_rows(np, m, x)
            out = np.concatenate([r, pr])
            nev[0] += 1
            return np.where(np.isfinite(out), out, 1e3)

        def jac(x):
            r, J = m.residual_jacobian(x)
            _, pj = penalty_rows(np, m, x)
            return np.vstack([np.asarray(J, float), pj])

        c0 = float(np.sum(resid(x0) ** 2))
        log(f"    start chi2 {c0:.3f}, {x0.size} free, {nres} points")
        if args.max_nfev > 0 and x0.size:
            sol = least_squares(resid, x0, jac=jac, method="trf", x_scale="jac",
                                max_nfev=args.max_nfev, xtol=1e-10, ftol=args.ftol,
                                gtol=1e-10)
            x, fun, J, status = sol.x, sol.fun, sol.jac, sol.status
        else:
            x = x0
            fun = resid(x)
            J = jac(x)
            status = 0
        chi2 = float(np.sum(fun ** 2))
        kprof = profiled_count(m, x)
        log(f"    end chi2 {chi2:.3f} (status {status}, {nev[0]} evaluations, "
            f"{time.time() - t0:.0f} s)")

        # statistical covariance of x (unscaled): (J^T J)^-1, through the SVD
        # of the column-scaled J -- the widths of one model differ by orders
        # of magnitude in their columns, and (J^T J)^-1 formed directly loses
        # the poorly constrained directions to round-off
        if x.size:
            cn = np.sqrt(np.sum(J ** 2, axis=0))
            cn[cn == 0] = 1.0
            P = np.linalg.pinv(J / cn, rcond=1e-12)
            cov_x = (P @ P.T) / np.outer(cn, cn)
        else:
            cov_x = np.zeros((0, 0))
        allx = lambda xx: np.asarray(m.transform_all_rwa(m._all_rwa(xx),
                                                         include_fixed=True), float)
        params = list(m.parameters)
        free = [p for p in params if not p.fixed and p.free_index is not None]
        names = [ma.parameter_label(p, m.pairs) for p in free]

        def quantities(xx):
            a = allx(xx)
            # norms/shifts/cbkg: the free vector itself; R-matrix: physical
            q = [float(a[p.index]) if p.kind in ("energy", "width")
                 else float(xx[p.free_index]) for p in free]
            d = calc(a) if calc else {}
            return np.asarray(q + list(d.values()), float), list(d)

        q, dnames = quantities(x)
        # dq/dx by central differences of the transform (no data evaluation)
        D = np.zeros((q.size, x.size))
        for i in range(x.size):
            h = 1e-6 * max(abs(x[i]), 1e-12)
            xp = x.copy(); xp[i] += h
            xm = x.copy(); xm[i] -= h
            D[:, i] = (quantities(xp)[0] - quantities(xm)[0]) / (2 * h)
        cov_q = D @ cov_x @ D.T
        sig = np.sqrt(np.maximum(np.diag(cov_q), 0.0))
        npar = len(names)
        values = dict(zip(names, map(float, q[:npar])))
        cov = (names, cov_q[:npar, :npar].tolist())
        derived = {n: (float(q[npar + i]), float(sig[npar + i]))
                   for i, n in enumerate(dnames)}
        if args.derived_hook:
            hook = load_hook(args.derived_hook)
            for k, v in hook(m, x).items():
                derived[str(k)] = v
        fitpath = os.path.splitext(path)[0] + "_fit.azr"
        try:
            m.save_fit(fitpath, x)
        except Exception as err:
            log(f"    save_fit failed: {err}")
        thm = {n: dict(norm=r["norm"], b=[float(v) for v in r["b"]], chi2=r["chi2"])
               for n, r in m.thm_experiments(x).items()}
        seg = [float(c) for c in m.segment_chi2(x)]
    return dict(chi2=chi2, npoints=nres, nfree=int(x.size) + kprof,
                values=values, covariance=cov, derived=derived,
                status=int(status), nev=nev[0], seconds=time.time() - t0,
                start_chi2=c0, segment_chi2=seg, thm=thm)


def load_hook(spec):
    path, _, func = spec.partition(":")
    sp = importlib.util.spec_from_file_location("derived_hook", path)
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return getattr(mod, func or "derived")


# ---------------------------------------------------------------------------

def copy_project(src_dir, dst_dir, azr_name):
    """The project directory without its outputs and other .azr files."""
    def ignore(d, names):
        skip = set()
        if os.path.abspath(d) == os.path.abspath(src_dir):
            for n in names:
                if n in ("output", "checks") or (n.endswith(".azr") and n != azr_name):
                    skip.add(n)
        return skip
    shutil.copytree(src_dir, dst_dir, ignore=ignore)
    for sub in ("output", "checks"):
        os.makedirs(os.path.join(dst_dir, sub), exist_ok=True)


def priors_for(grid, rules):
    """{label: prior} from rules {"axis=value": prior} (multiplied)."""
    out = {}
    for label, s in grid:
        p = 1.0
        for rule, val in rules.items():
            k, _, v = rule.partition("=")
            if k in s and str(s[k]) == v:
                p *= float(val)
        out[label] = p
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="See the module docstring for details.")
    ap.add_argument("project", help="the .azr of the adopted model")
    ap.add_argument("--out", default="modelavg", help="output directory")
    ap.add_argument("--spec", help="JSON variant spec")
    ap.add_argument("--radius-pairs", dest="radius_pairs", type=int, nargs="+")
    ap.add_argument("--radii", type=float, nargs="+")
    ap.add_argument("--vertex-model", dest="vertex_model", nargs="+",
                    choices=("pw", "dw"))
    ap.add_argument("--vertex", nargs="+", choices=("constant", "perlevel", "onshell"))
    ap.add_argument("--optical", nargs="+")
    ap.add_argument("--ps", nargs="+")
    ap.add_argument("--experiment", nargs="+")
    ap.add_argument("--pw-distortion", action="store_true",
                    help="apply the --optical axis to pw variants too (R(E))")
    ap.add_argument("--max-nfev", dest="max_nfev", type=int, default=30,
                    help="least_squares evaluation limit per variant (0: no fit)")
    ap.add_argument("--ftol", type=float, default=1e-8)
    ap.add_argument("--weights", default="aic", choices=ma.WEIGHT_METHODS)
    ap.add_argument("--rescale", default=None,
                    help="'best' (chi2/nu of the best variant) or a number")
    ap.add_argument("--prior", action="append", default=[],
                    help="AXIS=VALUE:WEIGHT, e.g. radius=6.1:0.5 (repeatable)")
    ap.add_argument("--strength", action="append", default=[],
                    help="NAME=JPI@E (excitation MeV), e.g. 213=2-@13.057")
    ap.add_argument("--strength-in", dest="strength_in", type=int, nargs="+")
    ap.add_argument("--strength-out", dest="strength_out", type=int, nargs="+")
    ap.add_argument("--derived-hook", dest="derived_hook",
                    help="file.py[:func]; func(session, x) -> {name: value or "
                         "(value, sigma)}")
    ap.add_argument("--dry-run", action="store_true", help="list the variants and stop")
    args = ap.parse_args(argv)

    project = os.path.abspath(args.project)
    src_dir, azr_name = os.path.split(project)
    spec = build_spec(args)
    grid = variant_grid(spec, args.pw_distortion)
    rules = dict(spec.get("priors", {}))
    for p in args.prior:
        rule, _, w = p.rpartition(":")
        if not rule or "=" not in rule:
            raise SystemExit(f"--prior {p!r}: expected AXIS=VALUE:WEIGHT")
        rules[rule] = float(w)
    priors = priors_for(grid, rules)
    strengths = parse_strengths(args.strength)
    if strengths and not (args.strength_in and args.strength_out):
        raise SystemExit("--strength needs --strength-in and --strength-out (file pair keys).")
    if args.rescale not in (None, "best"):
        args.rescale = float(args.rescale)

    print(f"{len(grid)} variant(s) of {azr_name}:")
    checked = []
    for label, s in grid:
        try:
            apply_variant(AzrModel.from_file(project), s, spec)
            note = ""
        except (ValueError, KeyError) as err:
            note = f"  -- refused: {err}"
        checked.append((label, s, note))
        print(f"  {label:<36} {json.dumps(s)}  prior {priors[label]:g}{note}")
    if args.dry_run:
        return 0

    out = os.path.abspath(args.out)
    work = os.path.join(out, "work")
    if os.path.exists(work):
        shutil.rmtree(work)
    os.makedirs(out, exist_ok=True)
    copy_project(src_dir, work, azr_name)
    logf = open(os.path.join(out, "run.log"), "w")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    rows, variants = [], []
    for i, (label, s, _) in enumerate(checked):
        log(f"[{i + 1}/{len(checked)}] {label}")
        row = dict(label=label, **{k: s.get(k, "") for k in AXES}, status="",
                   reason="", chi2="", npoints="", nfree="", chi2_nu="",
                   nev="", seconds="")
        try:
            mdl = apply_variant(AzrModel.from_file(os.path.join(work, azr_name)), s, spec)
            vdir = os.path.join(work, f"output_{label}")
            mdl.set_output_dir(vdir)          # own caches: the radius changes them
            path = mdl.write(os.path.join(work, f"{label}.azr"))
            r = fit_variant(path, work, args, strengths, log)
        except Exception as err:               # refused or failed: record, go on
            reason = f"{type(err).__name__}: {err}".replace("\n", " ")
            log(f"    skipped -- {reason}")
            if not isinstance(err, (ValueError, KeyError, RuntimeError)):
                log("    " + traceback.format_exc().replace("\n", "\n    "))
            row.update(status="skipped", reason=reason[:500])
            rows.append(row)
            continue
        v = ma.Variant(label, s, r["chi2"], r["npoints"], r["nfree"], r["values"],
                       covariance=r["covariance"], derived=r["derived"],
                       prior=priors[label])
        variants.append((v, r))
        row.update(status="fitted", chi2=f"{r['chi2']:.4f}", npoints=r["npoints"],
                   nfree=r["nfree"],
                   chi2_nu=(f"{r['chi2'] / v.dof:.4f}" if v.dof > 0 else ""),
                   nev=r["nev"], seconds=f"{r['seconds']:.0f}")
        for k in v.derived:
            row[k] = f"{v.derived[k]:.6g}"
            if k in v.derived_errors:
                row[k + "_stat"] = f"{v.derived_errors[k]:.3g}"
        rows.append(row)

    cols = list(dict.fromkeys(k for row in rows for k in row))
    with open(os.path.join(out, "variants.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols)
        wr.writeheader()
        for row in rows:
            wr.writerow(row)
    with open(os.path.join(out, "variants.json"), "w") as fh:
        json.dump([dict(v.to_dict(), fit={k: r[k] for k in
                   ("status", "nev", "seconds", "start_chi2", "segment_chi2", "thm")})
                   for v, r in variants], fh, indent=1)
    if not variants:
        log("no variant could be fitted; nothing to average.")
        return 1

    avg = ma.model_average([v for v, _ in variants], weights=args.weights,
                           rescale=args.rescale)
    avg.write_json(os.path.join(out, "average.json"))
    avg.write_csv(os.path.join(out, "average.csv"))
    adir = os.path.join(out, "averaged")
    if os.path.exists(adir):
        shutil.rmtree(adir)
    copy_project(src_dir, adir, azr_name)
    rm = {k: val for k, val in avg.values.items() if k.startswith(("E[", "G["))}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        apath, written, skipped = ma.write_averaged_azr(
            os.path.join(adir, azr_name),
            os.path.join(adir, os.path.splitext(azr_name)[0] + "_avg.azr"), rm)
    os.remove(os.path.join(adir, azr_name))
    skipped_rows = [r for r in rows if r["status"] == "skipped"]
    text = [avg.summary(), "",
            f"averaged .azr: {apath} ({len(written)} values written"
            + (f", {len(skipped)} not in the template" if skipped else "") + ")"]
    def sign_mixed(a):
        vals = [v[2] for v in a.values if v[1] > 0.01] or [a.mean]
        return min(vals) < 0 < max(vals)
    mixed = [a.name for a in list(avg.parameters.values()) + list(avg.derived.values())
             if sign_mixed(a)]
    if mixed:
        text.append("sign differs between variants of weight > 1 % (the mean of "
                    "such a quantity is not meaningful): " + ", ".join(mixed))
    if skipped_rows:
        text.append("skipped variants:")
        text += [f"  {r['label']}: {r['reason']}" for r in skipped_rows]
    text.append("stat = weighted mean of the variances; model = weighted variance "
                "of the means; total = their quadrature sum. The averaged .azr is "
                "a representative model (the project's radii and THM settings), "
                "not a fit.")
    with open(os.path.join(out, "summary.txt"), "w") as fh:
        fh.write("\n".join(text) + "\n")
    log("\n".join(text))
    logf.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
