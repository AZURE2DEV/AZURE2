#!/usr/bin/env python3
"""Local release (2026-10-02).  the evaluator: a sign flip or a new level should be judged after "a limited polish of the nearby
other levels and background levels" -- until this date the sign round's screen and stage-1 polish freed only the flipped
level's own J^pi group, and the level-add test only the new level's widths before a short global polish.
13N has a 3/2- at 3.5032 and a 5/2+ at 3.5453 MeV (42 keV apart, different J^pi) and 20 MeV background levels:
 1. local_mask(centre 3.5032, +-0.1) frees BOTH levels (all their energy/width keys), every background width (>= bg_energy),
    the norms, and nothing of the 1/2+ 2.3682; exclude_levels drops a level;
 2. the sign round's columns for a 3/2- flip: group-only (local_release 0) has no 5/2+ key; local_release 0.1 has the 5/2+
    3.5453 keys and the background widths of other groups;
 3. polish_flip's stage-1 free mask equals those local columns (polish intercepted);
 4. a sign round with local_release runs to the end and logs the local release;
 5. run_level_add with local_release logs the local stage for the candidate AND for the control, and reaches a verdict.
Slow part (4-5): 13N recomputes capture integrals per evaluation; ~20-40 min with the small budgets."""
import os, shutil, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.engine import Sharded
from rmfit.local import BoundsPolicy, default_bounds, polish, local_mask, level_energies_now
from rmfit.moves import single_flips
from rmfit import search as S
from rmfit.search import SearchConfig, sign_round, move_cols, polish_flip
from rmfit.campaign import Campaign
from rmfit.structure import LevelCandidate
from rmfit.azr3 import Model as AzrModel
failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok: failures.append(name)
E32, E52, HW = 3.5032, 3.5453, 0.1

def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_local_"); work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True); os.makedirs(os.path.join(work, "output"))
    mdl = AzrModel.from_file(os.path.join(work, "13N.azr")); mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    base = mdl.write(os.path.join(work, "base.azr"))
    policy = BoundsPolicy(bg_energy=15.0, e_window=0.5)
    sh = Sharded(base, os.path.join(tmp, "shards"), nshards=3, flags={})
    x = np.array(sh.x0, float); E = level_energies_now(sh, x)
    lab32 = [l for l, e in E.items() if l.startswith("3/2-") and abs(e - E32) < 1e-4][0]
    lab52 = [l for l, e in E.items() if l.startswith("5/2+") and abs(e - E52) < 1e-4][0]
    lab12 = [l for l, e in E.items() if l.startswith("1/2+") and abs(e - 2.3682) < 1e-4][0]
    print(f"  levels: {lab32} {E[lab32]:.4f}, {lab52} {E[lab52]:.4f}, {lab12} {E[lab12]:.4f}")
    # 1. local_mask
    m = local_mask(sh, x, [E32], HW, bg_energy=15.0)
    sel = {str(sh.keys[j]) for j in np.where(m)[0]}
    lv_keys = lambda lab: {str(k) for k in sh.keys if k.level == lab}
    bgw = {str(k) for k in sh.keys if k.level is not None and k.kind == "width" and E[k.level] >= 15.0}
    norms = {str(k) for k in sh.keys if k.kind == "norm"}
    check("local mask holds every key of the 3/2- and of the 5/2+ 42 keV away", lv_keys(lab32) | lv_keys(lab52) <= sel,
          str((lv_keys(lab32) | lv_keys(lab52)) - sel))
    check("local mask holds every background width and every norm", bgw <= sel and norms <= sel, str((bgw | norms) - sel))
    check("local mask holds nothing of the 1/2+ 2.368", not (lv_keys(lab12) & sel), str(lv_keys(lab12) & sel))
    check("local mask holds nothing else", sel == lv_keys(lab32) | lv_keys(lab52) | bgw | norms,
          str(sel - (lv_keys(lab32) | lv_keys(lab52) | bgw | norms)))
    mx = local_mask(sh, x, [E32], HW, bg_energy=15.0, exclude_levels=[lab32])
    check("exclude_levels drops the level", not ({str(sh.keys[j]) for j in np.where(mx)[0]} & lv_keys(lab32)))
    # 2. sign-round columns
    flips = [mv for mv in single_flips(x, sh.keys, jpi=["3/2-"]) if any(k.level == lab32 for k in mv.targets)]
    check("a 3/2- single flip exists", bool(flips))
    mv = flips[0]
    c0 = {str(sh.keys[j]) for j in move_cols(sh, x, mv, SearchConfig(), policy)}
    c1 = {str(sh.keys[j]) for j in move_cols(sh, x, mv, SearchConfig(local_release=HW), policy)}
    check("group-only columns have no 5/2+ key", not (c0 & lv_keys(lab52)), str(c0 & lv_keys(lab52)))
    check("local columns have the 5/2+ 3.5453 keys", lv_keys(lab52) <= c1, str(lv_keys(lab52) - c1))
    other_bg = {s for s in bgw if not s.split("[")[1].startswith("3/2-")}
    check("local columns have background widths of other groups", bool(other_bg) and other_bg <= c1, str(other_bg - c1))
    # 3. polish_flip stage 1
    seen = {}
    class Stop(Exception): pass
    real = S.polish
    def spy(ev, x0, free=None, **kw):
        seen["free"] = None if free is None else {str(ev.keys[j]) for j in np.where(free)[0]}; raise Stop()
    S.polish = spy
    try:
        polish_flip(sh, x, mv, policy, SearchConfig(local_release=HW), log=lambda *a: None)
    except Stop:
        pass
    finally:
        S.polish = real
    check("stage-1 polish frees exactly the local columns", seen.get("free") == c1,
          f"{len((seen.get('free') or set()) ^ c1)} keys differ")
    # 4. a sign round with the local release
    t0 = time.time(); logs = []
    b = default_bounds(sh, x, policy); f0 = polish(sh, x, bounds=b, policy=policy, max_nfev=60); x0, o0 = f0.x, f0.objective
    cfg = SearchConfig(top_k=2, pairs_per_group=1, stage1_nfev=10, stage2_nfev=15, union_nfev=15, local_release=HW,
                       jpi=["3/2-", "5/2+"], channel_pairs=False)
    xb, ob, stat, trials = sign_round(sh, x0, o0, cfg, policy, log=lambda s: (logs.append(str(s)), print("    " + str(s))))
    check("sign round with local release finished", np.isfinite(ob) and ob <= o0 * (1 + 1e-9), f"{o0} -> {ob}")
    check("sign round logged the local release", any("local release" in s for s in logs))
    print(f"  sign round: {(time.time() - t0) / 60:.1f} min")
    sh.close() if hasattr(sh, "close") else None
    # 5. level add with the local release: remove the 3/2- 3.5032 cleanly, offer it back
    tmp2 = tempfile.mkdtemp(prefix="rmfit_local_add_")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N", "data"), os.path.join(tmp2, "data")); os.makedirs(os.path.join(tmp2, "output"))
    m2 = AzrModel.from_file(os.path.join(AZURE2_ROOT, "tests", "13N", "13N.azr"))
    true = m2.find(jpi="3/2-", energy=E32)[0]
    chans = [dict(pair=c.pair, L=c.L, S=c.S, gamma=float(c.gamma)) for c in true.channels if float(c.gamma) != 0.0]
    m2.remove_level(jpi="3/2-", energy=E32); m2.set_output_dir(os.path.join(tmp2, "output"))
    base2 = m2.write(os.path.join(tmp2, "13N.azr"))
    c = Campaign.init(tmp2, base2, flags={}, nshards=3, policy=dict(bg_energy=15.0, e_window=0.5), note="local-release add test")
    t0 = time.time()
    st = c.run_level_add(LevelCandidate(energy=E32, jpi="3/2-", channels=chans, source="test"), jpis=["3/2-"],
                         pinned_nfev=6, global_nfev=15, local_release=HW, local_nfev=10, note="test")
    print(f"  level add: outcome {st.get('outcome')}, delta {st.get('delta')}, {(time.time() - t0) / 60:.1f} min")
    log = open(os.path.join(tmp2, "campaign.log")).read()
    check("level add reached a verdict", st.get("outcome") not in (None, "scan_failed", "crashed"), str(st.get("outcome")))
    check("level add logged the local-release polish", "local-release polish (+-100 keV" in log)
    check("level add gave the control the same local stage", "control local-release polish (level off" in log)
    print("RESULT:", "PASS" if not failures else f"FAIL ({len(failures)}): {failures}")
    sys.exit(1 if failures else 0)

if __name__ == "__main__":
    main()
