#!/usr/bin/env python3
"""rmfit drive: PRUNE over the run's own levels, and a clean residual scan before DONE (the evaluator 2026-11-21).

In a benchmark a narrow level adopted before the broad level it stood in for stayed in the final
model.  The rules: "prune-own" (backward elimination over the levels this run adopted, nothing else) and "scan"
(DONE only after a residual scan of the final structure).
  (a) STRUCTURE_ADD records an accepted level in drive_state["drive_added"]; PRUNE passes only those to run_remove
      (adopt=True); an adopted removal marks the level pruned, counts a structure change and re-runs (stay); a kept
      level is not re-tested on the same structure (next).
  (b) LEVELSCAN does not re-append a level this run pruned, nor a J^pi the bank entry was already tested with.
  (c) CYCLE with no clean scan -> LEVELSCAN (forced: runs even with the level-test budget used up); after a clean scan
      on the same structure -> DONE; a structure change after the scan -> another final scan; final_scans_max caps it.
Engine-free, seconds.  Run:  python3 rmfit/tests/test_drive_prune_scan.py"""
import json, os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
os.environ.pop("JOB_ID", None)
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3")
    sys.exit(77)

FAILS = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""))
    sys.stdout.flush()
    if not ok:
        FAILS.append(name)


def main():
    import test_drive_dry as tdd
    from rmfit import drive as dv
    tmp = tempfile.mkdtemp(prefix="rmfit_prune_scan_")
    try:
        base = tdd.make_dir(tmp)
        d = dv.Driver(tmp, base=base, local=True)
        stub = tdd.StubCampaign(tmp, d.decisions)
        stub.structure_id = 5
        removed = []

        def run_remove(top=3, nfev=100, only=None, adopt=False, note="", **kw):
            stub.calls.append(("remove", tuple((o["jpi"], round(o["energy"], 3)) for o in (only or ())), adopt))
            res = []
            for o in only or ():
                rm = o["jpi"] == "6+" and not removed
                res.append(dict(level=f"{o['jpi']}#1", jpi=o["jpi"], energy=o["energy"], cost=5.0 if rm else 500.0, threshold=20.0,
                                verdict="REMOVE" if rm else "keep"))
            a = None
            if adopt and any(r["verdict"] == "REMOVE" for r in res):
                r = next(r for r in res if r["verdict"] == "REMOVE")
                removed.append(r["jpi"])
                stub.structure_id += 1
                a = dict(level=r["level"], jpi=r["jpi"], energy=r["energy"], cost=r["cost"], threshold=r["threshold"],
                         structure=stub.structure_id, candidate=99, objective=stub.obj + 5.0)
            return dict(results=res, adopted=a, control=stub.obj)
        stub.run_remove = run_remove
        d._campaign = stub
        # (a) two accepted adds (a 6+ first, then the 3-), one marginal
        json.dump(dict(candidates=[dict(energy=11.439, jpi="3-", jpis=["3-", "6+"], hypothesis="scan cluster", width_kev=50.0),
                                   dict(energy=11.48, jpi="3-", hypothesis="ENSDF 3-", width_kev=1000.0),
                                   dict(energy=10.9, jpi="2+", hypothesis="ENSDF 2+", width_kev=20.0)]),
                  open(os.path.join(tmp, "candidates.json"), "w"))
        stub.add_outcomes = ["accept", "accept", "marginal"]
        orig_add = stub.run_level_add

        def run_level_add(cand, jpis=None, **kw):
            r = orig_add(cand, jpis=jpis, **kw)
            if r["outcome"] == "accept":
                r["jpi"] = "6+" if abs(cand.energy - 11.439) < 1e-6 else cand.jpi
                r["energy"] = cand.energy
                r["energy_fit"] = cand.energy + 0.008
                stub.structure_id += 1
            return r
        stub.run_level_add = run_level_add
        for _ in range(3):
            d._s_structure_add()
        own = d.state.get("drive_added", [])
        check("STRUCTURE_ADD recorded the two accepted levels as the run's own", [(a["jpi"], a["energy"]) for a in own] == [("3-", 11.48), ("6+", 11.439)],
              str(own))
        r1 = d._s_prune()
        rm = [c for c in stub.calls if c[0] == "remove"]
        check("PRUNE passed only the run's own levels, adopt=True", rm and set(rm[-1][1]) == {("3-", 11.48), ("6+", 11.439)} and rm[-1][2] is True, str(rm))
        check("adopted removal of the 6+ -> stay, marked pruned", r1[0] == "stay" and any(a["jpi"] == "6+" and a["pruned"] for a in d.state["drive_added"]),
              str(r1) + str(d.state["drive_added"]))
        bank = json.load(open(os.path.join(tmp, "candidates.json")))["candidates"]
        check("its bank entry marked pruned", bank[0].get("pruned") is True, str(bank[0]))
        check("SIGN_POST armed after the removal", int(d.counters.get("adopted_since_sign", 0)) == 1)
        r2 = d._s_prune()
        rm = [c for c in stub.calls if c[0] == "remove"]
        check("re-run on the new structure tests only the remaining own level (3-), kept -> next", r2[0] == "next" and rm[-1][1] == (("3-", 11.48),),
              str(r2) + str(rm[-1]))
        n = len(rm)
        r3 = d._s_prune()
        check("not re-tested on the same structure", r3[0] == "next" and len([c for c in stub.calls if c[0] == "remove"]) == n, str(r3))
        # (b) the scan does not re-append the pruned 6+ nor a J^pi the bank entry was tested with
        stub.levelscan_results = [[dict(ex=11.445, score=50.0, run_width=0.01, gamma_kev=0.05, best_global=[
            dict(jpi="6+", desc="a L6", signs="+", gamma_total_ev=50.0, ex=11.447, frozen=-5.0, gn=-10.0, gain=12.0,
                 widths_ev=[50.0], channels=[dict(pair=1, L=6, S=0.0)])])],
            [dict(ex=12.0, score=80.0, run_width=0.05, gamma_kev=50.0, best_global=[
                dict(jpi="4+", desc="a L4", signs="+", gamma_total_ev=5e4, ex=12.0, frozen=-5.0, gn=-10.0, gain=30.0,
                     widths_ev=[5e4], channels=[dict(pair=1, L=4, S=0.0)])])]]
        d.state["force_scan"] = True
        r4 = d._s_levelscan()
        bank = json.load(open(os.path.join(tmp, "candidates.json")))["candidates"]
        check("LEVELSCAN: the pruned 6+ near 11.447 is not appended again", len(bank) == 3 and "0 candidate(s) appended" in r4[1], r4[1])
        check("forced scan clears force_scan and records the structure it ran on", not d.state.get("force_scan")
              and d.state["scan_on"]["structure"] == stub.structure_id, str(d.state.get("scan_on")))
        # (c) CYCLE: a clean scan on the current structure -> DONE
        d.state["cycle_start_objective"] = stub.obj
        d.state["cycle"] = int(d.budget["cycles_max"])
        r5 = d._s_cycle()
        check("CYCLE after a clean scan on the final structure -> DONE", r5[0] == "DONE" and "clean" in r5[1], str(r5))
        # a structure change after the scan -> final scan (forced even with the level-test budget used up)
        d.counters["structure_changes"] = int(d.counters.get("structure_changes", 0)) + 1
        stub.structure_id += 1
        d.counters["level_tests"] = int(d.budget["level_tests_max"])
        r6 = d._s_cycle()
        check("structure changed after the scan -> LEVELSCAN (final scan 1)", r6[0] == "LEVELSCAN" and d.state.get("force_scan") and d.state["final_scans"] == 1, str(r6))
        r7 = d._s_levelscan()
        check("final scan runs although the level-test budget was used up", "appended" in r7[1] and any(c[0] == "levelscan" for c in stub.calls[-2:]), r7[1])
        bank = json.load(open(os.path.join(tmp, "candidates.json")))["candidates"]
        check("...and appended the new 4+ candidate", bank[-1]["jpi"] == "4+", str(bank[-1]))
        r8 = d._s_cycle()
        check("its candidate untested -> another final scan (2)", r8[0] == "LEVELSCAN" and d.state["final_scans"] == 2 and "untested" in r8[1], str(r8))
        d.state["force_scan"] = False
        d.state["tested_candidates"].append(len(bank) - 1)
        r9 = d._s_cycle()
        check("all of the last scan's candidates tested, same structure -> DONE", r9[0] == "DONE", str(r9))
        d.counters["structure_changes"] += 1
        r10 = d._s_cycle()
        check("final_scans_max 2 reached -> DONE, with the reason in the message", r10[0] == "DONE" and "NOT clean" in r10[1], str(r10))
        check("PRUNE and the stop rule are in the state list / defaults", "PRUNE" in dv.LOOP_STATES and dv.LOOP_STATES.index("PRUNE") == dv.LOOP_STATES.index("STRUCTURE_ADD") + 1
              and dv.DEFAULT_STOP["scan_before_done"] == 1 and dv.DEFAULT_STRUCTURE["prune_own"] is True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if FAILS:
        print(f"{len(FAILS)} FAILED"); sys.exit(1)
    print("all prune / final-scan checks passed")


if __name__ == "__main__":
    main()
