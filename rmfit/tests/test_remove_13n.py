#!/usr/bin/env python3
"""Backward elimination on the 13N model: the plan's Phase 4(a), wired as
Campaign.run_remove and `rmfit structure remove`.

Checks that the frozen screen ranks every level, that an honest refit follows
for the cheapest ones, that a removal decision is taken against the plan's
threshold, and -- the property that matters -- that the last level of a J^pi
group is never removed, since it carries that group's hard-sphere term.
"""
import os, shutil, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit import AZURE2_ROOT
from rmfit.campaign import Campaign
from pyazr import AzrModel

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_rm_")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
    mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))
    t0 = time.time()
    c = Campaign.init(tmp, base, flags={}, nshards=3,
                      policy=dict(bg_energy=15.0, e_window=0.5), note="backward-elimination test")
    cid, x, obj, keys = c.incumbent()
    print(f"  baseline objective {obj:.1f}")

    st = c.run_remove(top=2, nfev=20, note="test")
    print(f"  screened {st['screened']} levels, refitted {st['refitted']}, "
          f"removable {st['removable']}, {st['wall_min']:.1f} min")

    check("the frozen screen ranked several levels", st["screened"] >= 4, f"got {st['screened']}")
    check("at least one candidate got an honest refit", st["refitted"] >= 1, f"got {st['refitted']}")
    for r in st["results"]:
        print(f"    {r['level']} at {r['energy']:.4f} MeV: frozen {r['frozen_delta']:+.1f}, "
              f"refit cost {r['cost']:+.1f} vs threshold {r['threshold']:.1f} -> {r['verdict']}")
    if st["results"]:
        r = st["results"][0]
        # A polish only descends, so the removed model ends at or below obj + frozen_delta
        # and the control ends at or below obj; the cost is measured against the control, so
        # the bound is frozen_delta + (obj - control), NOT frozen_delta.  Comparing the cost
        # directly with the frozen estimate mixes two references and is ill-posed -- the
        # first version of this test did exactly that and failed once the control was added.
        bound = r["frozen_delta"] + (obj - st["control"])
        check("the refit cost respects the descent bound against the control",
              r["cost"] <= bound + max(1.0, 1e-6 * abs(bound)),
              f"cost {r['cost']:.1f} exceeds frozen {r['frozen_delta']:.1f} + (obj - control) "
              f"{obj - st['control']:.1f} = {bound:.1f}")
        check("the control polish is at least as good as every removal variant",
              all(x["refit"] >= st["control"] - max(1.0, 1e-9 * abs(st["control"])) for x in st["results"])
              or st["control"] <= min(x["refit"] for x in st["results"]),
              f"control {st['control']:.1f} vs refits "
              + ", ".join(f"{x['refit']:.1f}" for x in st["results"]))
        check("a verdict was taken against a positive threshold",
              r["threshold"] > 0 and r["verdict"] in ("REMOVE", "keep"))

    # the guard that matters: a J^pi group with a single level must never be offered
    m = AzrModel.from_file(base)
    singles = set()
    for lv in m.levels:
        jpi = f"{lv.J:g}{'+' if lv.parity > 0 else '-'}".replace("0.5", "1/2")
        singles.add(jpi)
    offered = {r["jpi"] for r in st["results"]}
    lone = set()
    for jpi in offered:
        if len(list(m.find(jpi=jpi))) <= 1:
            lone.add(jpi)
    check("no lone-level J^pi group was refitted for removal", not lone, f"offered {lone}")
    print(f"  total {time.time() - t0:.0f} s")
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
    print(("\nFAILED: " + ", ".join(failures)) if failures else "\nbackward elimination works on 13N")
    raise SystemExit(1 if failures else 0)
