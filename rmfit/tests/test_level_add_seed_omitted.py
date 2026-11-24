#!/usr/bin/env python3
"""Level adds seed the J^pi group's omitted open channels: drive default, run_level_add forwarding (2026-11-22).

An ENSDF candidate names no channels, so add_level_variant gave the new level only the lowest wave of its group's first
two open pairs and wrote every other channel gamma 0 FIXED -- for good, in every later structure.  In a
benchmark the recovered 2+ at 13.0 MeV lacked p2 L3 and its three alpha1 channels (the reference's 2+ has all six);
freeing them in the final model was worth 1,305 of chi2, and the reference's own vector on that freed structure is the
reference (22,634.219).  The seeding existed (add_level_variant(seed_omitted=...), structure.effective_k) but only
behind RMFIT_SEED_OMITTED, and run_level_add could not pass it.
  (a) Driver STRUCTURE_ADD passes seed_omitted=True by default and False when spec structure.seed_omitted is false
      (engine-free stub).
  (b) Campaign.run_level_add forwards seed_omitted to the full-model add_level_variant (None when not given): 13N
      campaign, window scan and variant builder stubbed.
  (c) add_level_variant(seed_omitted=True) frees the group's open, omitted channels: shift_fn_groups (one pair; the
      default seeds 2- L1 and leaves L3 out) always, and the many-pair jpi_integer fixture when it is present
      (local only).
Run:  python3 rmfit/tests/test_level_add_seed_omitted.py"""
import json, os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.pop("JOB_ID", None)
os.environ.pop("RMFIT_SEED_OMITTED", None)
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3")
    sys.exit(77)
import numpy as np  # noqa: E402

FAILS = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""), flush=True)
    if not ok:
        FAILS.append(name)


class _Stop(Exception):
    pass


def part_a(tmp):
    import test_drive_dry as tdd
    from rmfit import drive as dv
    base = tdd.make_dir(tmp)
    d = dv.Driver(tmp, base=base, local=True)
    calls = []

    class Rec(tdd.StubCampaign):
        def run_level_add(self, cand, jpis=None, global_nfev=150, pinned_nfev=60, note="", **kw):
            calls.append(kw)
            return dict(outcome="worse", jpi=cand.jpi or "1/2+", delta=-1.0)
    d._campaign = Rec(tmp, d.decisions)
    cpath = os.path.join(tmp, "candidates.json")

    def one(seed_setting):
        json.dump(dict(candidates=[dict(energy=3.3, jpi="1/2+", width_kev=200.0, source="ensdf", channels=[],
                                        hypothesis="ENSDF level without channels")]), open(cpath, "w"))
        d.state["tested_candidates"] = []
        d.counters["level_tests"] = 0
        d.counters["marginal_in_a_row"] = 0
        d.budget["level_tests_max"] = 5
        if seed_setting is None:
            d.structure.pop("seed_omitted", None)
            d.structure.update({k: v for k, v in dv.DEFAULT_STRUCTURE.items() if k == "seed_omitted"})
        else:
            d.structure["seed_omitted"] = seed_setting
        n = len(calls)
        d._s_structure_add()
        return calls[n] if len(calls) > n else None
    kw = one(None)
    check("(a) STRUCTURE_ADD passes seed_omitted=True by default", kw is not None and kw.get("seed_omitted") is True, str(kw))
    kw = one(False)
    check("(a) spec structure.seed_omitted false -> seed_omitted=False", kw is not None and kw.get("seed_omitted") is False, str(kw))


def part_b(tmp):
    from rmfit import AZURE2_ROOT
    from rmfit import campaign as cm
    from rmfit.structure import LevelCandidate, ScanResult
    from pyazr import AzrModel
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
    mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))
    c = cm.Campaign.init(tmp, base, flags={}, nshards=2, policy=dict(bg_energy=15.0, e_window=0.5), note="seed_omitted forwarding test")
    seen = []
    orig_scan, orig_add = cm.scan_jpi_on_window, cm.add_level_variant
    cm.scan_jpi_on_window = lambda *a, **k: [ScanResult("1/2+", "+", 0.0, -1.0, np.zeros(1), "ok", "", None)]

    def rec(*a, **k):
        seen.append(k.get("seed_omitted", "absent"))
        raise _Stop()
    cm.add_level_variant = rec
    try:
        for arg in ({}, {"seed_omitted": True}, {"seed_omitted": False}):
            try:
                c.run_level_add(LevelCandidate(energy=3.3, jpi="1/2+", channels=[], source="test", note="forwarding"), **arg)
            except _Stop:
                pass
    finally:
        cm.scan_jpi_on_window, cm.add_level_variant = orig_scan, orig_add
        c.close()
    check("(b) run_level_add forwards seed_omitted to add_level_variant (None / True / False)", seen == [None, True, False], str(seen))


def part_c(tmp):
    from _fixtures import fixture
    from rmfit.structure import add_level_variant, LevelCandidate, channel_threshold
    from pyazr import AzrModel
    ran = 0
    for name, need_pairs in (("shift_fn_groups", 1), ("jpi_integer", 3)):
        fx = fixture(name)
        if not os.path.isdir(fx):
            print(f"   ({name} fixture not present: skipped)")
            continue
        src = os.path.join(fx, [f for f in os.listdir(fx) if f.endswith(".azr")][0])
        mdl = AzrModel.from_file(src)
        # a J^pi group whose existing level has >= need_pairs particle pairs and >= 2 particle channels
        best = None
        for lv in mdl.levels:
            ch = [c for c in lv.channels if not getattr(c, "is_photon", False)]
            if len(ch) >= 2 and len({c.pair for c in ch}) >= need_pairs:
                best = (lv, max(channel_threshold(c) for c in ch))
                break
        check(f"(c) {name}: a group with >= {need_pairs} pair(s) and >= 2 channels", best is not None)
        if best is None:
            continue
        lv, thr = best
        jpi = f"{lv.J:g}{'+' if lv.parity > 0 else '-'}" if not isinstance(lv.J, str) else lv.J
        e = thr + 1.0                                   # above every threshold of the group
        cand = LevelCandidate(energy=e, jpi=jpi, channels=[], source="test")
        a = add_level_variant(src, os.path.join(tmp, f"{name}_a.azr"), cand, jpi)
        b = add_level_variant(src, os.path.join(tmp, f"{name}_b.azr"), cand, jpi, seed_omitted=True)
        na, nb = len(a["channels"]), len(b["channels"])
        check(f"(c) {name} {jpi} @ {e:.2f} MeV: default seeds the lowest waves of two pairs only",
              na <= 2 and a["seeded_omitted"] == [], str(a["channels"]))
        check(f"(c) {name}: seed_omitted=True adds the group's omitted open channels, free, with a tiny seed",
              nb > na and b["seeded_omitted"] and all(not q["fixed"] and abs(q["gamma"]) <= 1.0 for q in b["seeded_omitted"]),
              str(b["seeded_omitted"]))
        newb = [l for l in AzrModel.from_file(b["path"]).find(jpi=jpi) if abs(l.energy - e) < 1e-6][-1]
        free_b = sum(1 for c in newb.channels if not c.channel_fixed)
        check(f"(c) {name}: in the written file those channels are free", free_b == nb, f"{free_b} free vs {nb} written")
        ran += 1
    check("(c) at least one fixture ran", ran > 0)


def main():
    for part in (part_a, part_b, part_c):
        tmp = tempfile.mkdtemp(prefix=f"rmfit_seedom_{part.__name__}_")
        try:
            part(tmp)
        except Exception as err:
            import traceback
            traceback.print_exc()
            check(f"{part.__name__} ran", False, repr(err))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    if FAILS:
        print(f"{len(FAILS)} FAILED"); sys.exit(1)
    print("all level-add seed_omitted checks passed")


if __name__ == "__main__":
    main()
