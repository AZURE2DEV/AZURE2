#!/usr/bin/env python3
"""add_level_variant(seed_omitted=...) -- 2026-09-20.

Default behaviour must be unchanged (only the candidate's named channels, the rest inert).
With seed_omitted=True (or RMFIT_SEED_OMITTED=1) the group's OTHER particle channels that
are open and comfortably alive (Gamma_W >= 10 x the 1 eV dead rule) are added FREE with a
seed of min(1 eV, 1% of Gamma_W), i.e. theta^2 <= 0.01; photon channels and closed channels
never.  No engine needed.  Run:  python3 rmfit/tests/test_seed_omitted.py
"""
import os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from rmfit.backend import PYAZR3 as _PYAZR3  # noqa: E402
if not _PYAZR3:
    print("SKIP: a pyazr 3 test; this process runs pyazr 2 (the pyazr 2 version is test_seed_omitted.py)")
    sys.exit(77)
from rmfit import AZURE2_ROOT                                   # noqa: E402  (puts pyazr on the path)
from rmfit.structure import add_level_variant, LevelCandidate   # noqa: E402
from rmfit.azr3 import Model as AzrModel  # noqa: E402

failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    if not ok: failures.append(name)

def new_level(path, jpi, energy):
    lv = [l for l in AzrModel.from_file(path).find(jpi=jpi) if abs(l.energy - energy) < 1e-6]
    return lv[-1]

def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_seedom_")
    base = os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr")
    os.environ.pop("RMFIT_SEED_OMITTED", None)

    # 13N 1/2+ group: channels (pair 1 p+12C, L0) and (pair 2 photon, L1).  Name the photon only.
    cand = LevelCandidate(energy=3.0, channels=[dict(pair=2, L=1, S=0.5, gamma=0.5)], source="test")
    a = add_level_variant(base, os.path.join(tmp, "a.azr"), cand, "1/2+")
    check("default: only the named channel, nothing seeded", len(a["channels"]) == 1 and a["seeded_omitted"] == [], str(a["channels"]))
    lv = new_level(a["path"], "1/2+", 3.0)
    check("default: the omitted particle channel is inert and fixed",
          all((c.gamma == 0.0 and c.channel_fixed) for c in lv.channels if c.pair == 1))

    b = add_level_variant(base, os.path.join(tmp, "b.azr"), cand, "1/2+", seed_omitted=True)
    ex = b["seeded_omitted"]
    check("seed_omitted: the open particle channel is added", [(c["pair"], c["L"]) for c in ex] == [(1, 0)], str(ex))
    check("seed is min(1 eV, 1% Gamma_W): theta^2 <= 0.01 and Gamma_W >= 10 eV",
          all(c["gamma"] <= 0.01 * c["wigner_ev"] + 1e-12 and c["gamma"] <= 1.0 and c["wigner_ev"] >= 10.0 for c in ex), str(ex))
    lv = new_level(b["path"], "1/2+", 3.0)
    free = [(c.pair, c.gamma) for c in lv.channels if not c.channel_fixed]
    check("written file: seeded channel is FREE and nonzero", any(p == 1 and g > 0 for p, g in free), str(free))
    check("no photon channel is ever seeded", all(c["pair"] != 2 for c in ex))

    # below the proton threshold (S_p = 1.944 MeV) the particle channel is closed
    low = LevelCandidate(energy=1.0, channels=[dict(pair=2, L=1, S=0.5, gamma=0.5)], source="test")
    c_ = add_level_variant(base, os.path.join(tmp, "c.azr"), low, "1/2+", seed_omitted=True)
    check("closed channel is not seeded", c_["seeded_omitted"] == [], str(c_["seeded_omitted"]))

    # a named channel is never duplicated
    both = LevelCandidate(energy=3.0, channels=[dict(pair=1, L=0, S=0.5, gamma=2000.0)], source="test")
    d = add_level_variant(base, os.path.join(tmp, "d.azr"), both, "1/2+", seed_omitted=True)
    check("named channel not duplicated", d["seeded_omitted"] == [] and len(d["channels"]) == 1, str(d["channels"]))

    # environment switch
    os.environ["RMFIT_SEED_OMITTED"] = "1"
    e = add_level_variant(base, os.path.join(tmp, "e.azr"), cand, "1/2+")
    check("RMFIT_SEED_OMITTED=1 opts in", len(e["seeded_omitted"]) == 1)
    os.environ["RMFIT_SEED_OMITTED"] = "0"
    f = add_level_variant(base, os.path.join(tmp, "f.azr"), cand, "1/2+")
    check("RMFIT_SEED_OMITTED=0 stays off", f["seeded_omitted"] == [])
    os.environ.pop("RMFIT_SEED_OMITTED", None)

    # the real many-channel case, when the archive file is present
    big = os.path.join(os.environ.get("RMFIT_ARCHIVE", os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))), "an archived campaign")  # archive model; RMFIT_ARCHIVE
    if os.path.exists(big):
        one = LevelCandidate(energy=12.0, channels=[dict(pair=1, L=2, S=0.5, gamma=1000.0)], source="test")
        g = add_level_variant(big, os.path.join(tmp, "g.azr"), one, "5/2+", seed_omitted=True)
        ex = g["seeded_omitted"]; pairs = sorted({c["pair"] for c in ex})
        print("   large model 5/2+ @ 12.0: %d channels seeded, pairs %s, Gamma_W %.3g .. %.3g eV, theta2 max %.2g"
              % (len(ex), pairs, min(c["wigner_ev"] for c in ex), max(c["wigner_ev"] for c in ex), max(c["gamma"] / c["wigner_ev"] for c in ex)))
        check("large model: inelastic alpha channels (pairs 4, 6, 8) are among the seeded", {4, 6, 8} <= set(pairs), str(pairs))
        check("large model: every seed has theta^2 <= 0.01 and none is near the dead rule",
              all(c["gamma"] <= 0.01 * c["wigner_ev"] + 1e-12 and c["wigner_ev"] >= 10.0 for c in ex))
        h = add_level_variant(big, os.path.join(tmp, "h.azr"), one, "5/2+")
        check("large model: default still leaves one free channel", len(h["channels"]) == 1)
    else:
        print("   (archive model not present: many-channel block skipped)")

    print()
    if failures:
        print(f"FAILED: {len(failures)}: {', '.join(failures)}"); sys.exit(1)
    print("all seed-omitted checks passed")

if __name__ == "__main__":
    main()
