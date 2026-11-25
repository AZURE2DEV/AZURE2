#!/usr/bin/env python3
"""run_remove(only=..., adopt=True): backward elimination restricted to named levels, the cheapest REMOVE adopted (2026-11-21).

The drive's PRUNE step passes the levels the same run adopted: in a benchmark a narrow level adopted before the broad
level it stood in for stayed in the final model, while elimination over the whole model would also have removed a known
level the data do not constrain.  Here: 13N with a spurious 3/2- at Ex 8.0 MeV, far
above the data (Ex <= 4.3 MeV), so switching it off costs nothing.
  (a) only=[3/2- @ 8.0, 3/2- @ 6.0 (no such level)]: exactly one level refitted (the 8.0 one); the 6.0 entry is skipped;
      no other level (e.g. the 1/2+ 2.3682 the data constrain) is refitted.
  (b) its verdict is REMOVE and, with adopt=True, the removal becomes the base model: a new ledger structure, the 3/2-
      at 8.0 switched off (widths 0, fixed: deactivate_level_variant keeps the numbering), the incumbent objective within
      the threshold of the control.
  (c) a second run_remove does not count the switched-off 3/2- as a member of its group (last-of-J^pi guard).
Run:  python3 rmfit/tests/test_remove_only_adopt_13n.py"""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3 (its pyazr 3 version is test_p3_remove_only_adopt_13n.py)")
    sys.exit(77)
from rmfit import AZURE2_ROOT
from rmfit.campaign import Campaign
from rmfit.structure import add_level_variant, LevelCandidate
from pyazr import AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_remove_only_")
    try:
        shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
        os.makedirs(os.path.join(tmp, "output"))
        mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
        mdl.set_output_dir(os.path.join(tmp, "output"))
        b0 = mdl.write(os.path.join(tmp, "13N_orig.azr"))
        info = add_level_variant(b0, os.path.join(tmp, "13N.azr"), LevelCandidate(energy=8.0, jpi="3/2-", channels=[], source="test"), "3/2-")
        base = os.path.join(tmp, "13N.azr")
        n32 = len([lv for lv in AzrModel.from_file(base).active_levels if lv.jpi == "3/2-"])
        print(f"  spurious level {info.get('level')} added; 3/2- active levels: {n32}")
        c = Campaign.init(tmp, base, flags={}, nshards=2, policy=dict(bg_energy=15.0, e_window=0.5), note="remove only/adopt test")
        lines = []
        c.log = lambda s, *a, **k: lines.append(str(s))
        sid0 = c.structure_id
        st = c.run_remove(top=5, nfev=10, only=[dict(jpi="3/2-", energy=8.0), dict(jpi="3/2-", energy=6.0)], adopt=True)
        print("\n".join("    | " + l for l in lines if "removal" in l or "->" in l or "skipped" in l))
        res = st.get("results", [])
        check("exactly one level refitted (the 3/2- at 8.0 MeV)", len(res) == 1 and res[0]["jpi"] == "3/2-" and abs(res[0]["energy"] - 8.0) < 1e-6,
              str([(r["level"], r["energy"]) for r in res]))
        check("the 3/2- @ 6.0 entry (no such level) is skipped", any("6.0" in l and "skipped" in l for l in lines), "; ".join(lines[-8:]))
        check("no data-constrained level refitted (1/2+ 2.3682)", not any(r["jpi"] == "1/2+" for r in res))
        check("verdict REMOVE", res and res[0]["verdict"] == "REMOVE", str(res[:1]))
        a = st.get("adopted")
        check("adopt=True: the removal is adopted", bool(a) and c.structure_id != sid0 and c.structure_id == a["structure"],
              f"adopted {a}, structure {sid0} -> {c.structure_id}")
        # a removal keeps the level in the file (numbering unchanged) with every width 0 and fixed
        off = [lv for lv in AzrModel.from_file(c.base_azr).active_levels if lv.jpi == "3/2-" and abs(float(lv.energy) - 8.0) < 1e-6]
        check("new base model: the 3/2- at 8.0 MeV switched off (widths 0, fixed)", len(off) == 1 and all(
            bool(ch.channel_fixed) and float(ch.gamma) == 0.0 for ch in off[0].channels), str([(ch.gamma, ch.channel_fixed) for lv in off for ch in lv.channels]))
        cid, x, obj, keys = c.incumbent()
        check("incumbent = the adopted removal, within the threshold of the control",
              a and abs(obj - a["objective"]) < 1e-6 and obj - st["control"] <= res[0]["threshold"] + 1e-6,
              f"incumbent {obj}, control {st['control']}, threshold {res[0]['threshold'] if res else None}")
        lines.clear()
        st3 = c.run_remove(top=12, nfev=5, adopt=False)
        n_rows = [l for l in lines if "levels ranked" in l]
        check("the switched-off level has left the removal screen (no free widths)", not any(abs(r["energy"] - 8.0) < 1e-6 for r in st3.get("results", [])),
              str([(r["level"], r["energy"]) for r in st3.get("results", [])]))
        # adopt=False leaves the structure alone
        sid1 = c.structure_id
        st2 = c.run_remove(top=5, nfev=10, only=[dict(jpi="1/2+", energy=2.3682)], adopt=False)
        check("adopt=False never changes the structure", c.structure_id == sid1 and st2.get("adopted") is None, str(st2.get("adopted")))
        c.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all remove only/adopt checks passed")


if __name__ == "__main__":
    main()
