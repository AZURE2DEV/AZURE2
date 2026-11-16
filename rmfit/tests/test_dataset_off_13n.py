#!/usr/bin/env python3
"""Dataset-off (leave-one-dataset-out) on the 13N model: Campaign.run_dataset_off.

Checks the totals' bookkeeping (own chi2 + rest = incumbent; extra recovery >= 0 up to
noise, since a polish only descends) and -- added 2026-09-18 -- the per-dataset breakdown
of the REST (frozen vs refitted), which is what turns a leave-one-out table into a
pairwise tension map: the removed dataset must not appear in it, and its deltas must sum
to minus the extra recovery measured on the same (sharded) evaluator.
"""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit import AZURE2_ROOT
from rmfit.campaign import Campaign
from rmfit.engine import base_segments
from pyazr import AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_off_")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
    mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))
    c = Campaign.init(tmp, base, flags={}, nshards=3,
                      policy=dict(bg_energy=15.0, e_window=0.5), note="dataset-off test")
    cid, x, obj, keys = c.incumbent()
    files = sorted({os.path.basename(s.data_file) for s in base_segments(base) if s.active})
    print(f"  baseline objective {obj:.1f}; active data files: {files}")
    check("the fixture has at least two datasets", len(files) >= 2, f"{files}")
    target = os.path.splitext(files[0])[0]

    st = c.run_dataset_off(target, nfev=15, note="test")
    print(f"  off {target}: own {st['own_chi2']:.1f}, rest {st['rest_before']:.1f} -> {st['rest_after']:.1f}, "
          f"extra {st['extra_recovery']:+.2f}")
    check("own chi2 + rest = incumbent objective", abs(st["own_chi2"] + st["rest_before"] - obj) < 1e-6 * max(obj, 1))
    check("a polish only descends: extra recovery >= 0 (to noise)", st["extra_recovery"] > -1e-3 * max(obj, 1),
          f"{st['extra_recovery']}")
    rows = st.get("by_dataset")
    check("the per-dataset breakdown is recorded", bool(rows), f"{rows}")
    if rows:
        names = [r["name"] for r in rows]
        check("the removed dataset is not in the breakdown", target not in names, f"{names}")
        check("every other dataset is", set(names) == {os.path.splitext(f)[0] for f in files} - {target}, f"{names}")
        tot = sum(r["delta"] for r in rows)
        # chi2-only rows; the objective also carries norm penalties, so allow their change
        check("the deltas sum to minus the extra recovery (up to penalty terms)",
              abs(tot + st["extra_recovery"]) <= max(5.0, 0.05 * abs(st["extra_recovery"])),
              f"sum of deltas {tot:+.3f} vs extra {st['extra_recovery']:+.3f}")
        for r in rows:
            print(f"    {r['name']:30s} {r['before']:10.2f} -> {r['after']:10.2f}  {r['delta']:+9.2f}  (N {r['n']})")
    stats = [r for r in c.ledger.rounds() if r["kind"] == "dataset_off"]
    check("the round's ledger entry carries the breakdown", bool(stats) and "by_dataset" in (stats[-1].get("stats") or {}),
          f"{stats[-1].get('stats') if stats else None}")
    shutil.rmtree(tmp, ignore_errors=True)
    print("\n%s" % ("ALL PASSED" if not failures else f"{len(failures)} FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
