#!/usr/bin/env python3
"""perturb carries the fixed relvector covariance reference into the sibling campaign (2026-11-21).

A campaign whose data covariance has model-scaled (relvector) terms scores against ONE frozen reference model,
cov_reference.npz, named relative in campaign.json's objective.  perturb copied campaign.json to the sibling but not the
file, so the sibling's first session stopped with FileNotFoundError.  Fixture fixtures/relvector_p12c/ngrel.azr (p+12C, normgroups + 2 % relvector terms):
  (a) perturb returns and the sibling holds cov_reference.npz, byte-identical to the parent's;
  (b) the sibling's campaign.json objective names it, and Campaign(sibling).objective resolves to that file;
  (c) the sibling's verified baseline equals a SingleSession at the same vector with the PARENT's objective (one objective).
Run:  python3 rmfit/tests/test_perturb_cov_reference.py"""
import filecmp, json, os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402  (shared by both backends since 2026-11-26)
from _fixtures import fixture  # noqa: E402
from rmfit.campaign import Campaign
from rmfit.engine import SingleSession
from rmfit.keys import take
from rmfit.model import single_variant

FIX = fixture("relvector_p12c")
FLAGS = dict(use_brune=True)
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_perturb_cov_")
    try:
        d = os.path.join(tmp, "A")
        os.makedirs(d)
        shutil.copytree(os.path.join(FIX, "data"), os.path.join(d, "data"))
        shutil.copy2(os.path.join(FIX, "ngrel.azr"), d)
        c = Campaign.init(d, os.path.join(d, "ngrel.azr"), flags=FLAGS, nshards=2, policy=dict(bg_energy=15.0),
                          objective=dict(cov_free_norm_ok=["cov_rel.txt"]))
        ref = os.path.join(d, "cov_reference.npz")
        check("parent wrote cov_reference.npz (relvector fixture)", os.path.exists(ref))
        out = os.path.join(tmp, "P")
        try:
            st = c.perturb(out, n_flips=2, seed=1, jitter=0.1, e_kick=0.0, polish_nfev=8)
            ok = True
        except Exception as err:
            st, ok = {}, False
            check("(a) perturb returns", False, f"{type(err).__name__}: {err}")
        if ok:
            check("(a) perturb returns", True)
            sib = os.path.join(out, "cov_reference.npz")
            check("(a) sibling holds cov_reference.npz, identical to the parent's", os.path.exists(sib) and filecmp.cmp(ref, sib, shallow=False))
            cfg = json.load(open(os.path.join(out, "campaign.json")))
            check("(b) sibling campaign.json names it", (cfg.get("objective") or {}).get("cov_reference") == "cov_reference.npz", str(cfg.get("objective")))
            c2 = Campaign(out)
            check("(b) sibling objective resolves to its own copy", os.path.normpath(c2.objective.cov_reference) == os.path.normpath(sib),
                  c2.objective.cov_reference)
            cid, x, obj, keys = c2.incumbent()
            src = single_variant(c2.base_azr, "chk", c2.vdir)
            with SingleSession(src, cwd=c2.vdir, flags=c.flags, name="chk", objective=c.objective) as ev:
                o = ev.objective(take(x, keys, ev.keys, ev.x0))
            check("(c) sibling baseline = the parent's objective at the same vector", abs(o - obj) <= 1e-6 * max(1.0, abs(obj)), f"{obj} vs {o}")
            c2.close()
        c.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all perturb cov-reference checks passed")


if __name__ == "__main__":
    main()
