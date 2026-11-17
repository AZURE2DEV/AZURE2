#!/usr/bin/env python3
"""Close same-J^pi level pairs in bake and sanity, with the real engine on 13N (2026-10-10; rule changed 2026-11-08, the evaluator).
An engine >= dev e239a12 keeps levels with different IDs apart at any gap (enginepin.keeps_close_levels), so:
  exact    a free 1/2+ put 0.3 keV from a broad and from a narrow 1/2+: bake writes the vector as it is, nothing moved; the file
           reads back with every key and the same objective; sanity passes with a close_pairs warning
  sharded  the baked file evaluated by Sharded (3 shards) gives SingleSession's objective
  scan     the objective stays finite as the gap shrinks to 1e-6 keV and follows gap^-4 below 1e-3 keV: no round-off in the Brune
           divided differences.  The growth is the formalism, not a bug: the added level has one channel, so the pair's amplitudes are
           parallel and one combination decouples into a state of width ~gap^2 that capture feeds, with no photon width in the level
           matrix to limit it.  A pair with non-parallel amplitudes stays finite (two-level calculation, 2026-11-08).
  same E   an identical energy is still one level: sanity fails it
  old      a build before e239a12 (keeps_close_levels patched False): the broad pair is separated to 2.5 keV, the fixed level does
           not move, the cost cap raises, the narrow pair is refused
Run: python3 rmfit/tests/test_bake_separate_13n.py   (a few minutes; needs an engine >= e239a12)"""
import os, sys, shutil, tempfile, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, os.environ.get("AZURE2_ROOT", os.path.expanduser("~/AZURE2")))
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3 (its pyazr 3 version is test_p3_bake_separate_13n.py)")
    sys.exit(77)
import numpy as np
from pyazr import AzrModel
import rmfit.enginepin as enginepin
from rmfit.engine import SingleSession, Sharded, bake
from rmfit.local import sanity, SAME_ENERGY_KEV
from rmfit.keys import take
from rmfit.model import single_variant
AZ = os.environ.get("AZURE2_ROOT", os.path.expanduser("~/AZURE2"))
failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"), flush=True)
    if not ok:
        failures.append(name)


E_REF = dict(broad=20.0, narrow=2.3682)
# the narrow level is put 0.3 keV away IN THE FILE: 13N's capture segments keep the sub-point grids of the session's start (no
# re-placement with external capture), and moving a 3 keV wide level by 4.7 keV in the vector changed the objective by 2 % against a
# session that starts from the written energy.
DE_SEED = dict(broad=0.005, narrow=0.0003)


def build(tmp, tag, gamma_ev):
    mdl = AzrModel.from_file(os.path.join(AZ, "tests", "13N", "13N.azr"))
    mdl.add_level(0.5, 1, E_REF[tag] + DE_SEED[tag], [dict(pair=1, L=0, S=0.5, gamma=gamma_ev, fixed=False)], level_fixed=False)
    mdl.set_output_dir(os.path.join(tmp, "output"))
    return mdl.write(os.path.join(tmp, f"13N_{tag}.azr"))


def free_index(ev, tag):
    x = np.array(ev.x0, float)
    j = [i for i, k in enumerate(ev.keys) if k.kind == "energy" and k.level.startswith("1/2+") and abs(x[i] - E_REF[tag] - DE_SEED[tag]) < 1e-7]
    assert len(j) == 1, ([str(k) for k in ev.keys], x)
    return x, list(ev.keys), j[0]


def vector(base, tmp, tag, gap_kev=0.3):
    with SingleSession(single_variant(base, f"v_{tag}", tmp), cwd=tmp, flags={}, name=f"v_{tag}") as ev:
        x, keys, j = free_index(ev, tag)
        x[j] = E_REF[tag] + gap_kev * 1e-3
        return x, keys, j


def run_bake(*a, **kw):
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            return bake(*a, **kw), buf.getvalue()
    finally:
        if buf.getvalue().strip():
            print("   " + buf.getvalue().strip().replace("\n", "\n   "))


def written(azr, tag):
    return sorted(float(l.energy) for l in AzrModel.from_file(azr).levels if l.active and l.jpi == "1/2+" and abs(l.energy - E_REF[tag]) < 0.1)


def main():
    if not enginepin.keeps_close_levels():
        print("SKIP: this engine is older than dev e239a12"); return 0
    tmp = tempfile.mkdtemp(prefix="rmfit_sep_")
    shutil.copytree(os.path.join(AZ, "tests", "13N", "data"), os.path.join(tmp, "data")); os.makedirs(os.path.join(tmp, "output"))
    bases = {}
    for tag, gamma in (("broad", 1.0e6), ("narrow", 3.0e3)):
        print(f"== exact: {tag}", flush=True)
        base = bases[tag] = build(tmp, tag, gamma); x, keys, j = vector(base, tmp, tag)
        (azr, sav, sc), out = run_bake(base, x, keys, os.path.join(tmp, f"baked_{tag}.azr"), {}, workdir=tmp)
        check(f"{tag}: nothing separated", "separated" not in out, out)
        es = written(azr, tag)
        check(f"{tag}: written gap is the fitted 0.3 keV", len(es) == 2 and abs((es[1] - es[0]) * 1e3 - 0.3) < 1e-7, str(es))
        check(f"{tag}: bake verified the read-back", np.isfinite(sc.objective), str(sc))
        with SingleSession(single_variant(base, f"s_{tag}", tmp), cwd=tmp, flags={}, name=f"s_{tag}") as ev:
            rep = sanity(ev, x)
            check(f"{tag}: sanity has no merge_pairs and warns", not rep.merge_pairs and len(rep.close_pairs) == 1, rep.summary())
            print("   " + rep.summary()[:200], flush=True)
        # the Sharded fixture: the file that bake wrote, as a campaign restart would read it
        chk = single_variant(azr, f"chk_{tag}", tmp)
        with SingleSession(chk, cwd=tmp, flags={}, name=f"chk_{tag}") as ev:
            sk = list(ev.keys); xs = take(x, keys, sk, ev.x0); oS = ev.score(xs).objective
        H = Sharded(azr, os.path.join(tmp, f"sh_{tag}"), nshards=3, flags={})
        try:
            check(f"{tag}: Sharded has the same free parameters", set(H.keys) == set(sk))
            oH = H.score(take(x, keys, H.keys, H.x0)).objective
            check(f"{tag}: Sharded objective = Single ({oH:.8f} vs {oS:.8f})", abs(oH - oS) < 1e-7 * abs(oS))
            check(f"{tag}: same objective as the baked score", abs(oS - sc.objective) < 1e-6 * abs(sc.objective), f"{oS} {sc.objective}")
        finally:
            H.close()
    print("== scan and identical energy", flush=True)
    tag = "narrow"; base = bases[tag]
    with SingleSession(single_variant(base, "scan", tmp), cwd=tmp, flags={}, name="scan") as ev:
        x, keys, j = free_index(ev, tag); obj = {}
        for gap in (1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
            x[j] = E_REF[tag] + gap * 1e-3
            obj[gap] = ev.score(x).objective; print(f"   gap {gap:g} keV: objective {obj[gap]:.10f}", flush=True)
        check("scan: finite down to 1e-6 keV", all(np.isfinite(v) for v in obj.values()), str(obj))
        check("scan: gap^-4 below 1e-3 keV", all(abs(obj[g * 0.1] / obj[g] / 1e4 - 1) < 0.05 for g in (1e-3, 1e-4)), str(obj))
        x[j] = E_REF[tag]
        rep = sanity(ev, x, transform=False)
        check(f"same E: a gap below {SAME_ENERGY_KEV:g} keV fails sanity", len(rep.merge_pairs) == 1 and not rep.ok, rep.summary())
    print("== build before e239a12 (patched)", flush=True)
    real = enginepin.keeps_close_levels
    enginepin.keeps_close_levels = lambda *a, **k: False
    try:
        tag = "broad"; base = bases[tag]; x, keys, j = vector(base, tmp, tag)
        (azr, sav, sc), out = run_bake(base, x, keys, os.path.join(tmp, "baked_old.azr"), {}, workdir=tmp)
        check("old: bake reports the separation", "separated near-degenerate" in out, out)
        es = written(azr, tag)
        check("old: written levels 2.5 keV apart", len(es) == 2 and abs((es[1] - es[0]) * 1e3 - 2.5) < 1e-3, str(es))
        check("old: the fixed-energy level did not move", abs(es[0] - 20.0) < 1e-6, str(es))
        try:
            run_bake(base, x, keys, os.path.join(tmp, "baked_cap.azr"), {}, workdir=tmp, separate_max_cost=1e-9)
            check("old: cost cap raises", False, "no exception")
        except ValueError as e:
            check("old: cost cap raises", "would cost" in str(e), str(e))
        tag = "narrow"; base = bases[tag]; x, keys, j = vector(base, tmp, tag)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            bake(base, x, keys, os.path.join(tmp, "baked_old_narrow.azr"), {}, workdir=tmp)
        check("old: the narrow pair is refused", "NOT separated" in buf.getvalue(), buf.getvalue())
    finally:
        enginepin.keeps_close_levels = real
    shutil.rmtree(tmp, ignore_errors=True)
    print("ALL PASSED" if not failures else f"FAILED: {failures}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
