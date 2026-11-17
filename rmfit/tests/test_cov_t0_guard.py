#!/usr/bin/env python3
"""t0 guard (2026-10-26, the evaluator: "build the t0 version"): Campaign.init refuses a FREE data normalization on rows of a
relvector covariance unless the covariance file is named in ObjectiveSpec.cov_free_norm_ok.  Reason: a norm once ran
up fivefold with a data-unit reference and down threefold with a norm-following one; the t0 method (NNPDF, JHEP 05 (2010)
075) fixes the norm, puts its prior in C as a relvector term and holds C at a fixed model t0.

Fixture fixtures/relvector_p12c (ngrel.azr: normgroups A, B with free leader norms on cov_rel.txt rows; ngcov.azr: the
same segments on an absolute covariance without relvector terms).  Checks:
  (a) Evaluator.free_norms_on_relvector() on ngrel lists cov_rel.txt with the free-norm segments (leaders AND followers);
  (b) Campaign.init on ngrel refuses with a ValueError naming the file and those segments;
  (c) with objective cov_free_norm_ok = ["cov_rel.txt"] it initialises and logs a WARNING;
  (d) t0 form (ngrel with varyNorm 0 on every covariance segment): nothing listed, init passes without a warning, and the
      relvector factors are built once and stay fixed while the fit parameters move (constant C: chi2 quadratic);
  (e) ngcov (no relvector terms): nothing listed although its norms are free.
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3 (its pyazr 3 version is test_p3_cov_t0_guard.py)")
    sys.exit(77)
from _fixtures import fixture  # noqa: E402  (fixtures/ or fixtures_p3/ by backend)
import numpy as np
from rmfit.covariance import line_keyword, norm_field
from rmfit.engine import SingleSession, ObjectiveSpec

FIX = fixture("relvector_p12c")
FLAGS = dict(use_brune=True)
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


def fix_norms(src, dst):
    """``src`` with varyNorm 0 on every segment line that names a covariance file."""
    text = open(src).read()
    a = text.index("\n", text.index("<segmentsData>")) + 1
    b = text.index("</segmentsData>")
    out = []
    for l in text[a:b].splitlines():
        t = l.split()
        if t and line_keyword(t, "covariance") is not None:
            t[norm_field(t) + 1] = "0"
            l = " ".join(t)
        out.append(l)
    open(dst, "w").write(text[:a] + "\n".join(out) + "\n" + text[b:])


def main():
    from rmfit.campaign import Campaign
    tmp = tempfile.mkdtemp(prefix="rmfit_t0guard_")
    try:
        def fresh(name, azr, src=None):
            d = os.path.join(tmp, name)
            os.makedirs(d)
            shutil.copytree(os.path.join(FIX, "data"), os.path.join(d, "data"))
            if src is None:
                shutil.copy2(os.path.join(FIX, azr), d)
            else:
                fix_norms(os.path.join(FIX, src), os.path.join(d, azr))
            return d

        d = fresh("probe", "ngrel.azr")
        with SingleSession(os.path.join(d, "ngrel.azr"), cwd=d, flags=FLAGS, name="probe") as ev:
            got = ev.free_norms_on_relvector()
            from rmfit.covariance import NormGroups
            ng = NormGroups.from_azr(os.path.join(d, "ngrel.azr"))
            followers = sorted(ng.leader)
        segs = got.get("cov_rel.txt", [])
        check("(a) ngrel: cov_rel.txt listed with free-norm segments", bool(segs), f"{got}")
        check("(a) normgroup followers of a free leader are listed too", all(f in segs for f in followers),
              f"followers {followers} vs {segs}")

        d = fresh("refuse", "ngrel.azr")
        try:
            Campaign.init(d, os.path.join(d, "ngrel.azr"), flags=FLAGS, nshards=2, policy=dict(bg_energy=15.0))
            check("(b) init refuses free norms on relvector rows", False, "no exception")
        except ValueError as e:
            msg = str(e)
            check("(b) init refuses free norms on relvector rows",
                  "cov_rel.txt" in msg and all(str(s) in msg for s in segs) and "t0" in msg, msg)

        d = fresh("allow", "ngrel.azr")
        c = Campaign.init(d, os.path.join(d, "ngrel.azr"), flags=FLAGS, nshards=2, policy=dict(bg_energy=15.0),
                          objective=dict(cov_free_norm_ok=["cov_rel.txt"]))
        log = open(os.path.join(c.dir, "campaign.log")).read()
        check("(c) cov_free_norm_ok: init passes and logs a WARNING", "WARNING free norm on relvector rows" in log, log[-500:])

        d = fresh("t0", "ngrel_t0.azr", src="ngrel.azr")
        with SingleSession(os.path.join(d, "ngrel_t0.azr"), cwd=d, flags=FLAGS, name="t0probe") as ev:
            check("(d) t0 form: nothing listed", ev.free_norms_on_relvector() == {}, f"{ev.free_norms_on_relvector()}")
        c = Campaign.init(d, os.path.join(d, "ngrel_t0.azr"), flags=FLAGS, nshards=2, policy=dict(bg_energy=15.0))
        log = open(os.path.join(c.dir, "campaign.log")).read()
        check("(d) t0 form: init passes without a warning", "WARNING free norm" not in log, log[-500:])
        cid0, x0, obj0, keys0 = c.incumbent()
        spec = ObjectiveSpec(cov_reference=os.path.join(c.dir, "cov_reference.npz"))
        with SingleSession(c.base_azr, cwd=c.vdir, flags=FLAGS, name="t0", objective=spec) as ev:
            from rmfit.keys import take
            x = np.array(take(x0, keys0, ev.keys, ev.x0), float)
            ev.score(x)
            g = next(g for g in ev.cov.groups if g.relative)
            V0 = g.V
            rng = np.random.default_rng(3)
            xb = x * (1.0 + 0.05 * rng.standard_normal(x.size))
            s_b = ev.score(xb)
            check("(d) t0 form: relvector factors fixed while the parameters move", g.V is V0, "refactored")
            check("(d) t0 form: objective at a moved point is finite", np.isfinite(s_b.objective), f"{s_b.objective}")

        d = fresh("abs", "ngcov.azr")
        with SingleSession(os.path.join(d, "ngcov.azr"), cwd=d, flags=FLAGS, name="abs") as ev:
            check("(e) absolute covariance only: nothing listed", ev.free_norms_on_relvector() == {},
                  f"{ev.free_norms_on_relvector()}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("all passed" if not failures else f"FAILED: {failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
