#!/usr/bin/env python3
"""A data covariance's model-scaled (relvector) terms fixed at ONE reference model (2026-10-25, the evaluator: "revert the
default and use a fixed reference model"), after the refreeze iteration oscillated on a large model.

Fixture fixtures/relvector_p12c (p+12C): ngrel.azr = normgroups + a covariance with 2 % relvector terms; ngcov.azr =
the same with absolute terms only.  Checks:
  (a) Campaign.init on ngrel writes cov_reference.npz (the seed's frozen model) into campaign.json's objective, and the
      baseline equals a reference-free fresh session at the seed;
  (b) a SingleSession with that reference, FIRST evaluated far from the seed (the engine covariance check runs there
      and passes), equals an old-style session frozen at the seed and then moved there; its chi2_self equals a
      reference-free fresh session's chi2 at that point; refreeze_covariance() does not move it;
  (c) Sharded with the reference = SingleSession with the reference;
  (d) Campaign polish from a bad start: the log says the model is fixed, --refreeze-every is ignored, the polish's
      objective = the verified candidate's = the export's (one objective throughout);
  (e) init with cov_reference = another campaign's file gives that campaign's objective at the same point;
  (f) no relvector terms (ngcov): init writes no reference and nothing changes;
  (g) local.polish_refreeze returns the BEST self-consistent pass, not the last (scripted evaluator);
  (h) Campaign.run_polish and the CLI default to refreeze_every = 0.
"""
import inspect
import json
import os
import shutil
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rmfit.engine import SingleSession, Sharded, ObjectiveSpec
from rmfit.keys import take
from rmfit import local as L

FIX = os.path.join(HERE, "fixtures", "relvector_p12c")
FLAGS = dict(use_brune=True)
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def bad_start(keys, x0):
    x = np.array(x0, float)
    for j, k in enumerate(keys):
        if k.kind == "width":
            x[j] *= 1.6
    return x


def scripted_polish_refreeze():
    """(g): passes whose self-consistent objectives are 100, 90, 120, 95 -> pass 2's point is returned."""
    sc_vals = {1.0: 100.0, 2.0: 90.0, 3.0: 120.0, 4.0: 95.0}

    class Cov:
        relative, has_reference = True, False

    class Ev:
        cov = Cov()

        def refreeze_covariance(self):
            return 1

        def score(self, x):
            v = sc_vals[float(x[0])]
            return types.SimpleNamespace(objective=v, chi2=v, penalty=0.0, prior=0.0)

    calls = {"n": 0}

    def fake_polish(ev, x, bounds=None, policy=None, max_nfev=0, log=None, **kw):
        calls["n"] += 1
        k = float(calls["n"])
        return types.SimpleNamespace(x=np.array([k]), objective=sc_vals[k] - 5.0, start_objective=200.0,
                                     status="max_nfev", nfev=10, njev=10, message="max_nfev", free=None, dead=None,
                                     bounds_hit=[])
    real = L.polish
    L.polish = fake_polish
    try:
        r = L.polish_refreeze(Ev(), np.array([0.0]), bounds=object(), max_nfev=40, every=10)
    finally:
        L.polish = real
    check("(g) polish_refreeze returns the best self-consistent pass (2 of 4), not the last",
          float(r.x[0]) == 2.0 and r.objective == 90.0 and "pass 2 of 4" in r.message, f"{r.x} {r.objective} {r.message}")


def main():
    scripted_polish_refreeze()
    from rmfit.campaign import Campaign
    from rmfit import cli
    check("(h) Campaign.run_polish defaults to refreeze_every=0",
          inspect.signature(Campaign.run_polish).parameters["refreeze_every"].default == 0)
    check("(h) CLI --refreeze-every defaults to 0", 'add_argument("--refreeze-every", type=int, default=0' in
          open(cli.__file__).read())

    tmp = tempfile.mkdtemp(prefix="rmfit_covref_")
    try:
        def camp(name, azr, **kw):
            d = os.path.join(tmp, name)
            os.makedirs(d)
            shutil.copytree(os.path.join(FIX, "data"), os.path.join(d, "data"))
            shutil.copy2(os.path.join(FIX, azr), d)
            # ngrel's normgroup leaders are free norms on relvector rows: allowed explicitly since the 2026-10-26 guard
            kw["objective"] = dict(kw.get("objective") or {}, cov_free_norm_ok=["cov_rel.txt"])
            return Campaign.init(d, os.path.join(d, azr), flags=FLAGS, nshards=2, policy=dict(bg_energy=15.0), **kw)

        c = camp("A", "ngrel.azr")
        cfg = json.load(open(os.path.join(c.dir, "campaign.json")))
        ref = os.path.join(c.dir, "cov_reference.npz")
        check("(a) init wrote cov_reference.npz and campaign.json names it",
              os.path.exists(ref) and cfg["objective"].get("cov_reference") == "cov_reference.npz", f"{cfg['objective']}")
        cid0, x0, obj0, keys0 = c.incumbent()
        with SingleSession(c.base_azr, cwd=c.vdir, flags=FLAGS, name="fresh0") as f:
            fresh0 = f.score(take(x0, keys0, f.keys, f.x0)).objective
        check("(a) baseline = reference-free fresh session at the seed", np.isclose(obj0, fresh0, rtol=1e-10),
              f"{obj0} vs {fresh0}")

        spec = ObjectiveSpec(cov_reference=ref)
        with SingleSession(c.base_azr, cwd=c.vdir, flags=FLAGS, name="old") as old:
            xs, xb = take(x0, keys0, old.keys, old.x0), bad_start(old.keys, take(x0, keys0, old.keys, old.x0))
            old.score(xs)                                # frozen at the seed (old behaviour)
            o_old = old.score(xb)
        with SingleSession(c.base_azr, cwd=c.vdir, flags=FLAGS, name="free") as fr_:
            s_free = fr_.score(xb)                       # frozen at xb itself
        with SingleSession(c.base_azr, cwd=c.vdir, flags=FLAGS, name="ref", objective=spec) as ev:
            s_ref = ev.score(xb)                         # first evaluation at xb: engine check runs here
            check("(b) reference session at a far point = session frozen at the seed then moved",
                  np.isclose(s_ref.objective, o_old.objective, rtol=1e-10), f"{s_ref.objective} vs {o_old.objective}")
            check("(b) chi2_self = reference-free fresh session's chi2", s_ref.chi2_self is not None and
                  np.isclose(s_ref.chi2_self, s_free.chi2, rtol=1e-10), f"{s_ref.chi2_self} vs {s_free.chi2}")
            check("(b) the far point really differs from the reference (test has teeth)",
                  abs(s_free.chi2 - s_ref.chi2) > 1e-6 * s_free.chi2, f"{s_free.chi2} vs {s_ref.chi2}")
            ev.refreeze_covariance()
            s2 = ev.score(xb)
            check("(b) refreeze_covariance keeps the reference", s2.objective == s_ref.objective,
                  f"{s2.objective} vs {s_ref.objective}")
            print(f"  (b) at the bad start: fixed reference {s_ref.objective:.4f}, self-consistent chi2 "
                  f"{s_ref.chi2_self:.4f}, seed-frozen {o_old.objective:.4f}")
        sh = Sharded(c.base_azr, os.path.join(tmp, "shards"), nshards=3, flags=FLAGS, objective=spec)
        try:
            os_ = sh.objective(take(xb, ev.keys, sh.keys, sh.x0))
            check("(c) Sharded with the reference = single", np.isclose(os_, s_ref.objective, rtol=1e-9),
                  f"{os_} vs {s_ref.objective}")
        finally:
            sh.close()

        # (d) polish from a bad start, one objective throughout
        bcid, _ = c.verify_candidate(xb, ev.keys, move=dict(kind="bad_start"))
        frc = c.run_polish(max_nfev=60, cid=bcid, refreeze_every=15)
        cid1, x1, obj1, keys1 = c.incumbent()
        log = open(os.path.join(c.dir, "campaign.log")).read()
        check("(d) log: model fixed, --refreeze-every ignored",
              "model fixed at the campaign reference" in log and "--refreeze-every 15 ignored" in log)
        check("(d) polish objective = verified candidate", np.isclose(frc.objective, obj1, rtol=1e-9),
              f"{frc.objective} vs {obj1}")
        print(f"  (d) polish {frc.start_objective:.4f} -> {frc.objective:.4f} ({frc.status}); verified {obj1:.4f}")
        _, _, sce = c.export()
        check("(d) export objective = verified", np.isclose(sce.objective, obj1, rtol=1e-9), f"{sce.objective}")
        check("(d) verified and export logs carry the self-consistent chi2",
              "self-consistent chi2" in open(os.path.join(c.dir, "campaign.log")).read())

        # (e) another campaign given A's reference: same objective at the same point
        cB = camp("B", "ngrel.azr", cov_reference=ref)
        cfgB = json.load(open(os.path.join(cB.dir, "campaign.json")))
        with SingleSession(cB.base_azr, cwd=cB.vdir, flags=FLAGS, name="refB", objective=cB.objective) as eB:
            oB = eB.score(take(x1, keys1, eB.keys, eB.x0)).objective
        check("(e) init --cov-reference: copied, named, same objective as campaign A",
              cfgB["objective"].get("cov_reference") == "cov_reference.npz" and np.isclose(oB, obj1, rtol=1e-10),
              f"{cfgB['objective']} {oB} vs {obj1}")

        # (f) no relvector terms
        cC = camp("C", "ngcov.azr")
        cfgC = json.load(open(os.path.join(cC.dir, "campaign.json")))
        check("(f) no relvector terms: no reference written",
              not os.path.exists(os.path.join(cC.dir, "cov_reference.npz")) and
              not cfgC["objective"].get("cov_reference"), f"{cfgC['objective']}")
        for cc in (c, cB, cC):
            cc.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("FAILED: " + ", ".join(failures) if failures else "all passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
