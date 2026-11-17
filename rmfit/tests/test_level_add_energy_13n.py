#!/usr/bin/env python3
"""run_level_add must FIT the new level's energy (2026-09-19).

Until this date the added level was written with levelFix = 1 (LevelCandidate.energy_fixed
defaults to True and levelscan passes True), so the "released" polish never released its
energy and every level the searches added sat at a round number; freeing four of them
by hand gained clearly.  Here: remove the 3/2- 3.5032 MeV level of 13N cleanly, offer it
back 30 keV too high with its true widths, and require that (a) the adopted structure has a
FREE energy key for the new level, (b) the released stage moved it toward 3.5032, (c) the
move is reported, and the verdict is not "relocated" for a move inside the window.
Slow (13N recomputes capture integrals per evaluation): ~20-40 min with the small budgets.
"""
import os, shutil, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit import AZURE2_ROOT
from rmfit.campaign import Campaign
from rmfit.structure import LevelCandidate
from rmfit.azr3 import Model as AzrModel

failures = []
E_TRUE, OFFSET = 3.5032, 0.030


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}"))
    sys.stdout.flush()
    if not ok:
        failures.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_addE_")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    mdl = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
    true = mdl.find(jpi="3/2-", energy=E_TRUE)[0]
    chans = [dict(pair=c.pair, L=c.L, S=c.S, gamma=float(c.gamma)) for c in true.channels if float(c.gamma) != 0.0]
    print(f"  true level: 3/2- {E_TRUE} MeV, channels {[(c['pair'], c['L'], round(c['gamma'], 1)) for c in chans]}")
    mdl.remove_level(jpi="3/2-", energy=E_TRUE)
    mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))
    t0 = time.time()
    c = Campaign.init(tmp, base, flags={}, nshards=3, policy=dict(bg_energy=15.0, e_window=0.5), note="level-add energy test")
    cid, x, obj, keys = c.incumbent()
    print(f"  objective without the level: {obj:.1f}")
    cand = LevelCandidate(energy=E_TRUE + OFFSET, jpi="3/2-", channels=chans, source="test")   # energy_fixed left at its default (True)
    st = c.run_level_add(cand, jpis=["3/2-"], pinned_nfev=8, global_nfev=25, note="test")
    print(f"  outcome {st.get('outcome')}, delta {st.get('delta')}, {(time.time() - t0) / 60:.1f} min")
    check("the add ran to a verdict", st.get("outcome") not in (None, "scan_failed", "crashed"), str(st.get("outcome")))
    log = open(os.path.join(tmp, "campaign.log")).read()
    check("the released stage announced a free energy window", "released stage:" in log and "energy free within" in log)
    # the add variant on disk: the new level must be written with a FREE energy
    import glob
    adds = sorted(glob.glob(os.path.join(tmp, "variants", "add_3_2-_*.azr")))
    check("the add variant exists", bool(adds), str(os.listdir(os.path.join(tmp, "variants"))))
    if adds:
        lv = [l for l in AzrModel.from_file(adds[-1]).find(jpi="3/2-") if abs(l.energy - (E_TRUE + OFFSET)) < 0.2]
        check("the new level is in the variant with a free energy (levelFix 0)", bool(lv) and not lv[0].fixed,
              f"{[(l.energy, l.fixed) for l in lv]}")
    e_fit = st.get("energy_fit", None)
    moved = st.get("e_moved", None)
    print(f"  reported: e_moved {moved}, energy_fit {e_fit}")
    check("the move of the energy is reported", moved is not None, str(sorted(st.keys())))
    if moved is not None:
        check("the energy moved TOWARD the true value", moved < -0.005, f"e_moved {moved:+.4f} MeV (offset was +{OFFSET})")
        check("and landed within 15 keV of it", abs((E_TRUE + OFFSET + moved) - E_TRUE) < 0.015,
              f"fitted {E_TRUE + OFFSET + moved:.4f} vs true {E_TRUE}")
        check("a move inside the window is not 'relocated'", st.get("outcome") != "relocated", str(st.get("outcome")))
    shutil.rmtree(tmp, ignore_errors=True)
    print("\n%s" % ("ALL PASSED" if not failures else f"{len(failures)} FAILED: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
