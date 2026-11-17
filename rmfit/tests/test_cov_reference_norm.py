#!/usr/bin/env python3
"""A fixed covariance reference is a CROSS SECTION and its relvector terms follow the data normalization
(2026-10-26, the evaluator: "build the fix and test it"), after one large model ran a norm from 0.63 to 3.08
against a data-unit reference (its objective a quarter below the self-consistent chi2).

Fixture fixtures/relvector_p12c/ngrel.azr (normgroups A, B with free leader norms; covariance with 2 % relvector
terms; one segment with a fixed norm).  Checks:
  (a) at a point whose norms differ from the reference's, the objective's data chi2 = a DENSE numpy chi2 with
      C = B + vectors + sum_k (f_k * ref_xs / n)(...)^T, n = each row's current norm (ref_xs = the reference
      session's model x its norms) -- and differs from the data-unit (old) reference there (teeth);
  (b) at the reference norms the objective equals the old behaviour (a session frozen at the seed);
  (c) Jacobian: the norm columns (and two width columns) of the whitened rows = central differences of the
      residual vector, rel 1e-5 -- the d(W)/dn term included (and dropping it fails: teeth);
  (d) Sharded = SingleSession at the moved point, objective and the norm columns of the Jacobian;
  (e) a reference without normalizations (format before 2026-10-26) is refused with a clear message;
  (f) Campaign.init writes norm_i arrays into cov_reference.npz and the baseline = a reference-free fresh session;
  (g) no runaway: raising group A's norm from its reference value by x1.5, x3 (fit fixed) costs at least as much
      under the new reference as the self-consistent chi2 would let it, i.e. obj(new) >= obj(old data-unit)
      at every step, and old/new differ (the old reference is the cheap direction).
"""
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rmfit.engine import SingleSession, Sharded, ObjectiveSpec
from rmfit.keys import take
from rmfit.covariance import read_covariance, load_reference, save_reference

FIX = os.path.join(HERE, "fixtures", "relvector_p12c")
FLAGS = dict(use_brune=True)
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def norm_cols(keys):
    return {k.seg: j for j, k in enumerate(keys) if k.kind == "norm"}


def dense_chi2(ev, r, ref_xs_by_group):
    """Plain data chi2 of residual rows ``r`` (from ev.residuals) with a dense covariance per group, the relvector
    terms scaled by ref_xs / current norm."""
    plain = ev._cov_memo_get(r)[2]
    norms, _ = ev._cov_norm_rows(ev._last_x)
    tot = 0.0
    mask = np.ones(ev.npts, bool)
    for g in ev.cov.groups:
        mask[g.pos] = False
        cf = read_covariance(os.path.relpath(g.path, ev.cwd) if os.path.isabs(g.path) else g.path, ev.cwd)
        rows = g.rows
        C = np.zeros((cf.N, cf.N))
        C[np.diag_indices(cf.N)] += cf.diag
        for brows, vals in cf.blocks:
            C[np.ix_(brows, brows)] += vals
        for v in cf.vectors:
            C += np.outer(v, v)
        C = C[np.ix_(rows, rows)]
        meff = ref_xs_by_group[os.path.basename(g.path)] / norms[g.pos]
        for f in cf.rel:
            u = f[rows] * meff
            C += np.outer(u, u)
        y = plain[g.pos] * g.scale * g.sigma          # m - d in data units
        tot += float(y @ np.linalg.solve(C, y))
    return tot + float(np.sum(plain[mask] ** 2))


class Rec:
    """Wrap an evaluator so the dense check knows the x of the last residuals call."""


def run_resid(ev, x, jac=False):
    ev._last_x = np.asarray(x, float)
    return ev.residuals(x, jac=jac)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_covref_norm_")
    try:
        from rmfit.campaign import Campaign
        d = os.path.join(tmp, "A")
        os.makedirs(d)
        shutil.copytree(os.path.join(FIX, "data"), os.path.join(d, "data"))
        shutil.copy2(os.path.join(FIX, "ngrel.azr"), d)
        c = Campaign.init(d, os.path.join(d, "ngrel.azr"), flags=FLAGS, nshards=2, policy=dict(bg_energy=15.0),
                          objective=dict(cov_free_norm_ok=["cov_rel.txt"]))    # free norms on relvector rows on purpose
        ref = os.path.join(c.dir, "cov_reference.npz")
        z = np.load(ref)
        check("(f) cov_reference.npz carries norm_i arrays", any(n.startswith("norm_") for n in z.files), f"{z.files}")
        cid0, x0, obj0, keys0 = c.incumbent()
        with SingleSession(c.base_azr, cwd=c.vdir, flags=FLAGS, name="fresh0") as f:
            fresh0 = f.score(take(x0, keys0, f.keys, f.x0)).objective
        check("(f) baseline = reference-free fresh session at the seed", np.isclose(obj0, fresh0, rtol=1e-10),
              f"{obj0} vs {fresh0}")
        loaded = load_reference(ref)
        ref_xs = {}
        for name, (rows, model, norm) in loaded.items():
            ref_xs[name] = (rows, model * norm)

        spec = ObjectiveSpec(cov_reference=ref)
        with SingleSession(c.base_azr, cwd=c.vdir, flags=FLAGS, name="ref", objective=spec) as ev:
            xs = take(x0, keys0, ev.keys, ev.x0)
            nc = norm_cols(ev.keys)
            print(f"  free norm columns: {nc}; values at the seed {[round(float(xs[j]), 4) for j in nc.values()]}")
            check("(a) fixture has >= 2 free norms on covariance segments", len(nc) >= 2, f"{nc}")
            xm = xs.copy()
            jA, jB = list(nc.values())[:2]
            xm[jA] *= 2.5
            xm[jB] *= 0.6

            def xs_rows(g):
                rows, xsv = ref_xs[os.path.basename(g.path)]
                at = {int(rr): i for i, rr in enumerate(rows)}
                return xsv[[at[int(rr)] for rr in g.rows]]

            # (b) at the reference norms = old behaviour (frozen at the seed)
            sb = ev.score(xs)
            with SingleSession(c.base_azr, cwd=c.vdir, flags=FLAGS, name="old") as old:
                old.score(xs)
                o_seed = old.score(xs).objective
                o_old_m = old.score(xm).objective         # data-unit model frozen at the seed, moved norms
            check("(b) at the reference norms = session frozen at the seed", np.isclose(sb.objective, o_seed, rtol=1e-10),
                  f"{sb.objective} vs {o_seed}")

            # (a) dense check at moved norms
            r, _ = run_resid(ev, xm)
            sc = ev.score(xm)
            dense = dense_chi2(ev, r, {os.path.basename(g.path): xs_rows(g) for g in ev.cov.groups})
            check("(a) chi2 at moved norms = dense chi2 with ref_xs / n", np.isclose(sc.chi2, dense, rtol=1e-9),
                  f"{sc.chi2} vs {dense}")
            check("(a) differs from the data-unit reference there (teeth)",
                  abs(sc.objective - o_old_m) > 1e-4 * abs(o_old_m), f"{sc.objective} vs {o_old_m}")
            print(f"  (a) moved norms: new {sc.objective:.6f}, old data-unit reference {o_old_m:.6f}, "
                  f"self-consistent chi2 {sc.chi2_self:.6f}")

            # (c) Jacobian vs central differences
            r0, J = ev.residuals(xm, jac=True)
            wcols = [j for j, k in enumerate(ev.keys) if k.kind == "width" and abs(xm[j]) > 1e-3][:2]
            worst = 0.0
            worst_drop = 0.0
            for j in list(nc.values()) + wcols:
                h = 1e-6 * max(abs(xm[j]), 1e-3)
                xp, xq = xm.copy(), xm.copy()
                xp[j] += h
                xq[j] -= h
                fd = (ev.residuals(xp)[0] - ev.residuals(xq)[0]) / (2 * h)
                err = np.max(np.abs(J[:, j] - fd)) / max(np.max(np.abs(fd)), 1e-300)
                worst = max(worst, err)
                print(f"  (c) column {j} {ev.keys[j]}: max |J - FD| / max |FD| = {err:.2e}")
            # teeth: the same J without the d(W)/dn term
            norms, cols = ev._cov_norm_rows(xm)
            r_plain, Jp = ev._call("resid", (xm, True))
            _, Jno, _, _ = ev.cov.apply(r_plain, Jp, ev.row_sw, norms=norms, norm_cols=None)
            for j in nc.values():
                h = 1e-6 * max(abs(xm[j]), 1e-3)
                xp, xq = xm.copy(), xm.copy()
                xp[j] += h
                xq[j] -= h
                fd = (ev.residuals(xp)[0] - ev.residuals(xq)[0]) / (2 * h)
                worst_drop = max(worst_drop, np.max(np.abs(Jno[:, j] - fd[:ev.npts])) / max(np.max(np.abs(fd)), 1e-300))
            check("(c) Jacobian norm and width columns = central differences (rel 1e-5)", worst < 1e-5, f"{worst:.2e}")
            check("(c) without the d(W)/dn term the norm columns are wrong (teeth)", worst_drop > 1e-3,
                  f"{worst_drop:.2e}")

            # (g) raising group A's norm: old reference cheaper than new along the ray
            prev_ok = True
            for lam in (1.5, 3.0):
                xg = xs.copy()
                xg[jA] *= lam
                on = ev.score(xg).objective
                with SingleSession(c.base_azr, cwd=c.vdir, flags=FLAGS, name="oldg") as old:
                    old.score(xs)
                    oo = old.score(xg).objective
                print(f"  (g) norm A x{lam}: new {on:.4f}, old data-unit reference {oo:.4f}")
                prev_ok &= on >= oo - 1e-9 * abs(oo)
            check("(g) raising a norm is never cheaper under the new reference than under the old", prev_ok)

            # (d) Sharded = Single at the moved point
            sh = Sharded(c.base_azr, os.path.join(tmp, "shards"), nshards=3, flags=FLAGS, objective=spec)
            try:
                xsh = take(xm, ev.keys, sh.keys, sh.x0)
                osh = sh.objective(xsh)
                check("(d) Sharded objective = single at moved norms", np.isclose(osh, sc.objective, rtol=1e-9),
                      f"{osh} vs {sc.objective}")
                rsh, Jsh = sh.residuals(xsh, jac=True)
                ncs = norm_cols(sh.keys)
                ok = True
                for seg, j in nc.items():
                    a, b = np.sort(J[:ev.npts, j]), np.sort(Jsh[:sh.npts, ncs[seg]])
                    ok &= a.size == b.size and np.allclose(a, b, rtol=1e-7, atol=1e-10 * np.max(np.abs(a)))
                check("(d) Sharded norm Jacobian columns = single (as sorted values; row order differs)", ok)
            finally:
                sh.close()

        # (e) legacy reference refused
        legacy = os.path.join(tmp, "legacy.npz")
        save_reference(legacy, {n: (rows, model) for n, (rows, model, norm) in loaded.items()})
        try:
            with SingleSession(c.base_azr, cwd=c.vdir, flags=FLAGS, name="legacy",
                               objective=ObjectiveSpec(cov_reference=legacy)) as ev:
                ev.score(xs)
            check("(e) legacy reference refused", False, "no error")
        except ValueError as err:
            check("(e) legacy reference refused", "no data normalizations" in str(err), str(err))
        c.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("FAILED: " + ", ".join(failures) if failures else "all passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
