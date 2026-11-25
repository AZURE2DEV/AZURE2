#!/usr/bin/env python3
"""Drive defaults of 2026-11-24 (a method proposal from benchmark diagnostics, approved by the evaluator): engine-free, on the stub campaign
of test_drive_dry.

  (a) DEFAULTS_REV: a new drive gets the new defaults; a drive whose drive_state.json predates them keeps the old
      behaviour (add_first off, no norm box, no re-tests, marginal rule from cycle 1) unless its spec sets a key.
  (b) Cycle 1 with structure.add_first: TIGHTEN (norms boxed to structure.norm_range_first) -> STRUCTURE_ADD over the
      broad queue only (marginal results do not end it in cycle 1) -> TIGHTEN with the campaign's own norm box -> SIGN.
  (c) structure.retest_at_cap: an add unphysical only for an amplitude at the theta^2 cap is re-tested once with that
      channel closed and the level's other channels named (seed_omitted False); a total-width floor is not re-tested.
  (d) structure.split_retest_tau: a split with gain >= 10 tau is re-tested with the local release widened to hold the
      partner; still a split -> a non-blocking decision (no WAITING); a small split is left alone.
  (e) a spec with no level adds (budget.level_tests_max 0: the tier-1/2 bench specs) keeps the old cycle-1
      order and the policy's norm box.
Run:  python3 rmfit/tests/test_drive_add_first.py"""
import json, os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.pop("JOB_ID", None)
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if _PYAZR3:
    print("SKIP: a pyazr 2 test; this process runs pyazr 3")
    sys.exit(77)

FAILS = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""), flush=True)
    if not ok:
        FAILS.append(name)


def stub_class():
    import test_drive_dry as tdd
    from rmfit.local import BoundsPolicy

    class Stub(tdd.StubCampaign):
        def __init__(self, d, decisions, outcomes=None):
            super().__init__(d, decisions)
            self.policy = BoundsPolicy(theta2_cap=11.616)
            self.cfg["policy"] = {"norm_range": [0.2, 5.0], "theta2_cap": 11.616}
            self.polish_boxes, self.adds = [], []
            self.scripted = list(outcomes or [])

        def run_polish(self, *a, **kw):
            self.polish_boxes.append(tuple(self.policy.norm_range))
            return super().run_polish(*a, **kw)

        def run_level_add(self, cand, jpis=None, global_nfev=150, pinned_nfev=60, note="", **kw):
            self.adds.append(dict(energy=cand.energy, jpi=cand.jpi, channels=list(cand.channels), box=tuple(self.policy.norm_range), **kw))
            st = dict(self.scripted.pop(0)) if self.scripted else dict(outcome="worse", delta=-1.0)
            st.setdefault("jpi", cand.jpi or "1/2+")
            st.setdefault("energy", cand.energy)
            if st.get("outcome") == "accept":
                self.obj -= 10.0
            return st
    return Stub


def new_driver(tmp, structure=None, stop=None, old_state=False):
    import test_drive_dry as tdd
    from rmfit import drive as dv
    base = tdd.make_dir(tmp)
    sp = json.load(open(os.path.join(tmp, "spec.json")))
    for k in ("add_first", "retest_at_cap", "split_retest_tau"):
        sp["structure"].pop(k, None)
    sp["stop_rules"].pop("marginal_rule_from_cycle", None)
    sp["structure"].update(structure or {})
    sp["stop_rules"].update(stop or {})
    json.dump(sp, open(os.path.join(tmp, "spec.json"), "w"), indent=1)
    if old_state:
        json.dump(dict(state="TIGHTEN", cycle=1), open(os.path.join(tmp, dv.STATE_FILE), "w"))
    return dv, dv.Driver(tmp, base=base, local=True)


def part_a(tmp):
    dv, d = new_driver(os.path.join(tmp, "new"))
    check("(a) a new drive records defaults_rev 2", d.state.get("defaults_rev") == 2, str(d.state.get("defaults_rev")))
    check("(a) new drive: add_first, norm box [0.5, 2], re-tests on, marginal rule from cycle 2",
          d.structure["add_first"] is True and list(d.structure["norm_range_first"]) == [0.5, 2.0]
          and d.structure["retest_at_cap"] is True and d.structure["split_retest_tau"] == 10.0
          and d.stop["marginal_rule_from_cycle"] == 2, f"{d.structure} {d.stop}")
    dv, d = new_driver(os.path.join(tmp, "old"), old_state=True)
    check("(a) a drive started before 2026-11-24 keeps the old behaviour",
          d.structure["add_first"] is False and d.structure["norm_range_first"] is None and d.structure["retest_at_cap"] is False
          and d.structure["split_retest_tau"] == 0.0 and d.stop["marginal_rule_from_cycle"] == 1, f"{d.structure} {d.stop}")
    dv, d = new_driver(os.path.join(tmp, "old_spec"), structure=dict(add_first=True, norm_range_first=[0.6, 1.8]), old_state=True)
    check("(a) ... unless its spec sets the keys (a list value is read)",
          d.structure["add_first"] is True and list(d.structure["norm_range_first"]) == [0.6, 1.8], str(d.structure))


def part_b(tmp):
    from rmfit.local import BoundsPolicy  # noqa: F401
    dv, d = new_driver(tmp, structure=dict(local_release_mev=0.5))
    Stub = stub_class()
    c = Stub(tmp, d.decisions, outcomes=[dict(outcome="marginal", delta=1.0), dict(outcome="marginal", delta=1.0),
                                         dict(outcome="accept", delta=100.0)])
    d._campaign = c
    json.dump(dict(candidates=[
        dict(energy=3.0, jpi="1/2+", width_kev=500.0, source="ensdf", hypothesis="broad A"),
        dict(energy=3.2, jpi="3/2-", width_kev=300.0, source="ensdf", hypothesis="broad B"),
        dict(energy=3.4, jpi="5/2+", width_kev=100.0, source="ensdf", hypothesis="broad C"),
        dict(energy=3.6, jpi="1/2-", width_kev=1.0, source="ensdf", hypothesis="narrow D")]), open(os.path.join(tmp, "candidates.json"), "w"))
    d.state["state"] = "TIGHTEN"
    seen = []
    for _ in range(20):
        st = d.state["state"]
        seen.append(st)
        if st == "SIGN":
            break
        d.step(st)
    check("(b) cycle 1 order TIGHTEN.. -> STRUCTURE_ADD x4 -> TIGHTEN -> SIGN",
          [s for s in seen if s != "TIGHTEN" or True] and seen[-1] == "SIGN" and "POPULATION" not in seen
          and seen.index("STRUCTURE_ADD") < len(seen) - 1 and seen[seen.index("STRUCTURE_ADD"):].count("STRUCTURE_ADD") == 4
          and seen[-2] == "TIGHTEN", str(seen))
    first_add = seen.index("STRUCTURE_ADD")
    check("(b) the add-first polish ran with the norms boxed to 0.5-2.0",
          c.polish_boxes[:first_add] and all(b == (0.5, 2.0) for b in c.polish_boxes[:first_add]), str(c.polish_boxes))
    check("(b) the three broad adds ran inside the box; the narrow candidate was not tested in the phase",
          [a["energy"] for a in c.adds] == [3.0, 3.2, 3.4] and all(a["box"] == (0.5, 2.0) for a in c.adds), str(c.adds))
    check("(b) two marginal results did not end the cycle-1 broad queue", len(c.adds) == 3, str(len(c.adds)))
    check("(b) after the phase: state add_first done and the policy's own box (0.2-5.0) again",
          d.state.get("add_first") == "done" and tuple(c.policy.norm_range) == (0.2, 5.0) and c.polish_boxes[-1] == (0.2, 5.0),
          f"{d.state.get('add_first')} {c.policy.norm_range} {c.polish_boxes}")
    # cycle 2: the rule applies to the broad queue again
    d.state["cycle"] = 2
    d.counters["marginal_in_a_row"] = 2
    d.counters["level_tests"] = 0
    bank = json.load(open(os.path.join(tmp, "candidates.json")))
    bank["candidates"].append(dict(energy=3.8, jpi="3/2+", width_kev=200.0, source="ensdf", hypothesis="broad E"))
    json.dump(bank, open(os.path.join(tmp, "candidates.json"), "w"))
    r = d._s_structure_add()
    check("(b) cycle 2: two marginal adds in a row end the broad queue", r[0] == "next" and "marginal" in r[1], str(r))


def part_c(tmp):
    from rmfit.structure import LevelCandidate
    dv, d = new_driver(tmp)
    Stub = stub_class()
    th = {"g[2+#9 p1 L2 S0]": 0.29, "g[2+#9 p2 L3 S1]": 0.41, "g[2+#9 p3 L4 S2]": 11.6153}
    c = Stub(tmp, d.decisions, outcomes=[dict(outcome="unphysical", delta=4.9e5, theta2=th, at_bound="", energy_fit=13.002, jpi="2+"),
                                         dict(outcome="accept", delta=4.8e5, jpi="2+")])
    d._campaign = c
    cand = LevelCandidate(energy=13.02, jpi="2+", channels=[], source="ensdf", note="2+ ENSDF")
    st = d._retest(c, cand, ["2+"], d._level_add(c, cand, ["2+"]))
    check("(c) unphysical at the cap -> one re-test, accepted", st.get("outcome") == "accept" and len(c.adds) == 2, str(st))
    a2 = c.adds[1] if len(c.adds) > 1 else {}
    chans = sorted((q["pair"], q["L"], q["S"]) for q in a2.get("channels", []))
    check("(c) the re-test names the level's other channels and closes the one at the cap (seed_omitted False)",
          chans == [(1, 2, 0.0), (2, 3, 1.0)] and a2.get("seed_omitted") is False, str(a2))
    check("(c) the result says what was re-tested", "p3 L4" in str(st.get("retest_of")), str(st.get("retest_of")))
    c.scripted = [dict(outcome="unphysical", delta=1e3, theta2={"g[2+#9 p1 L2 S0]": 0.3}, at_bound="Gamma_total[2+#9] at the 10 keV floor")]
    n = len(c.adds)
    st = d._retest(c, cand, ["2+"], d._level_add(c, cand, ["2+"]))
    check("(c) a total-width floor is not re-tested", st.get("outcome") == "unphysical" and len(c.adds) == n + 1, str(st))


def part_d(tmp):
    from rmfit.structure import LevelCandidate
    dv, d = new_driver(tmp, structure=dict(local_release_mev=0.5))
    Stub = stub_class()
    split = json.dumps(dict(label="0+#8@14.5000", sep_kev=500.0, lost_ev=1e5, lost_frac=0.5))
    big = dict(outcome="splits", delta=91458.0, tau=3000.0, split=split, energy_fit=13.93, jpi="0+", structure=21, candidate=88)
    c = Stub(tmp, d.decisions, outcomes=[dict(big), dict(big)])
    d._campaign = c
    cand = LevelCandidate(energy=14.032, jpi="0+", channels=[], source="ensdf", note="0+ ENSDF")
    st = d._retest(c, cand, ["0+"], d._level_add(c, cand, ["0+"]))
    lr = c.adds[1].get("local_release") if len(c.adds) > 1 else None
    check("(d) a large split is re-tested with the local release widened to hold the partner (0.568 MeV)",
          lr is not None and abs(lr - 0.568) < 1e-6, str(c.adds))
    pend = [x for x in d.decisions.pending() if x["slug"].startswith("split-level-0+")]
    check("(d) still a split -> one decision, default 'leave', not a WAITING one",
          len(pend) == 1 and not d.waiting_on() and st.get("outcome") == "splits", f"{[x['slug'] for x in pend]} {d.waiting_on()}")
    c.scripted = [dict(outcome="splits", delta=353.0, tau=284.0, split=split, jpi="0+")]
    n = len(c.adds)
    st = d._retest(c, cand, ["0+"], d._level_add(c, cand, ["0+"]))
    check("(d) a small split is left alone (no re-test, no decision)",
          len(c.adds) == n + 1 and len([x for x in d.decisions.pending() if x["slug"].startswith("split-level")]) == 1, str(st))


def part_e(tmp):
    dv, d = new_driver(tmp)
    d.budget["level_tests_max"] = 0
    Stub = stub_class()
    c = Stub(tmp, d.decisions)
    d._campaign = c
    d.state["state"] = "TIGHTEN"
    seen = []
    for _ in range(6):
        st = d.state["state"]
        seen.append(st)
        if st != "TIGHTEN":
            break
        d.step(st)
    check("(e) level_tests_max 0: TIGHTEN -> SIGN, no add-first phase, the policy's norm box",
          seen[-1] == "SIGN" and d.state.get("add_first") is None and all(b == (0.2, 5.0) for b in c.polish_boxes),
          f"{seen} {d.state.get('add_first')} {c.polish_boxes}")


def main():
    for part in (part_a, part_b, part_c, part_d, part_e):
        tmp = tempfile.mkdtemp(prefix=f"rmfit_addfirst_{part.__name__}_")
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
    print("all drive add-first / re-test checks passed")


if __name__ == "__main__":
    main()
