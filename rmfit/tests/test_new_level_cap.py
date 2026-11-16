#!/usr/bin/env python3
"""new_level_theta2_cap (2026-09-27, the evaluator: "increase the cap for the narrow resonances to 3").

The new level gets its own theta^2 cap through a +-0.5 keV theta2_cap_windows entry at its seed energy;
the window is in force during the test and stays in the campaign policy only if the level is accepted.
Pure: BoundsPolicy.cap_for plus Campaign.run_level_add around a stubbed _run_level_add.
Run:  python3 rmfit/tests/test_new_level_cap.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.local import BoundsPolicy  # noqa: E402
from rmfit.campaign import Campaign  # noqa: E402
from rmfit.structure import LevelCandidate  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok:
        failures.append(name)


def fake(policy, outcome=None, raise_=False):
    c = Campaign.__new__(Campaign)
    c.policy = policy
    c.cfg = {"policy": {"theta2_cap": policy.theta2_cap, "theta2_cap_windows": list(policy.theta2_cap_windows)}}
    c.log = lambda *a, **k: None
    seen = {}

    def impl(cand, *a, **k):
        seen["cap_during"] = c.policy.cap_for(cand.energy)
        seen["cap_neighbour"] = c.policy.cap_for(cand.energy + 0.002)
        seen["cfg_windows_during"] = list(c.cfg["policy"]["theta2_cap_windows"])
        if raise_:
            raise RuntimeError("boom")
        return {"outcome": outcome}
    c._run_level_add = impl
    return c, seen


def main():
    p = BoundsPolicy(theta2_cap=0.5, bg_theta2_cap=3.0, bg_energy=13.3)
    w = BoundsPolicy.seed_window(12.986, 3.0)
    check("seed window spans +-0.5 keV", abs(w[0] - 12.9855) < 1e-9 and abs(w[1] - 12.9865) < 1e-9 and w[2] == 3.0, w)
    p2 = BoundsPolicy(theta2_cap=0.5, bg_theta2_cap=3.0, bg_energy=13.3, theta2_cap_windows=[w])
    check("window level gets 3", p2.cap_for(12.986) == 3.0)
    check("level 2 keV away keeps 0.5", p2.cap_for(12.988) == 0.5)
    check("background pole unchanged", p2.cap_for(13.7) == 3.0 and p2.cap_for(11.0) == 0.5)
    check("field defaults to off", BoundsPolicy().new_level_theta2_cap == 0.0)

    cand = LevelCandidate(energy=12.28, jpi="5/2-")
    # accepted: window kept in policy and in cfg
    c, seen = fake(BoundsPolicy(theta2_cap=0.5, bg_energy=13.3, new_level_theta2_cap=3.0), outcome="accept")
    c.run_level_add(cand, jpis=["5/2-"])
    check("cap 3 in force during the test", seen["cap_during"] == 3.0, seen)
    check("neighbour not affected during the test", seen["cap_neighbour"] == 0.5, seen)
    check("cfg carries the window during the test (adopt_structure writes it)", len(seen["cfg_windows_during"]) == 1, seen)
    check("accept: window kept in policy", c.policy.cap_for(12.28) == 3.0)
    check("accept: window kept in cfg", len(c.cfg["policy"]["theta2_cap_windows"]) == 1)
    # rejected: window removed from both
    c, seen = fake(BoundsPolicy(theta2_cap=0.5, bg_energy=13.3, new_level_theta2_cap=3.0), outcome="unphysical")
    c.run_level_add(cand)
    check("reject: window removed from policy", c.policy.cap_for(12.28) == 0.5)
    check("reject: window removed from cfg", c.cfg["policy"]["theta2_cap_windows"] == [])
    # exception: removed too, and the exception propagates
    c, seen = fake(BoundsPolicy(theta2_cap=0.5, bg_energy=13.3, new_level_theta2_cap=3.0), raise_=True)
    try:
        c.run_level_add(cand); raised = False
    except RuntimeError:
        raised = True
    check("exception propagates", raised)
    check("exception: window removed", c.policy.cap_for(12.28) == 0.5 and c.cfg["policy"]["theta2_cap_windows"] == [])
    # off: no window at all
    c, seen = fake(BoundsPolicy(theta2_cap=0.5, bg_energy=13.3), outcome="accept")
    c.run_level_add(cand)
    check("off: ordinary cap during the test", seen["cap_during"] == 0.5)
    check("off: no window added", c.cfg["policy"]["theta2_cap_windows"] == [])
    # existing windows are kept and still work
    old = [11.0, 11.3, 1.0]
    c, seen = fake(BoundsPolicy(theta2_cap=0.5, bg_energy=13.3, new_level_theta2_cap=3.0, theta2_cap_windows=[old]), outcome="accept")
    c.run_level_add(cand)
    check("existing windows kept", old in c.policy.theta2_cap_windows and c.policy.cap_for(11.1) == 1.0)

    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}")
        sys.exit(1)
    print("all new-level-cap checks passed")


if __name__ == "__main__":
    main()
