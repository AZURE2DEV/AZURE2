#!/usr/bin/env python3
"""Residual-driven level search on 13N: remove the 3/2- 3.5032 MeV level, polish
the rest, find the residual peak, and check that the template screen ranks 3/2-
first at that energy.  Run: python3 rmfit/tests/test_levelscan_13n.py
"""
import os, shutil, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from rmfit import AZURE2_ROOT
from rmfit.engine import SingleSession, Sharded
from rmfit.keys import take
from rmfit.local import BoundsPolicy, default_bounds, polish
from rmfit.levelscan import candidate_energies, template_variant, template_screen, best_per_jpi, to_level_candidate, incumbent_on_variant, dummy_labels
from rmfit.model import single_variant
from rmfit.structure import JPI_GROUPS
from pyazr import AzrModel
# Brune (observed) parameters: without them the boundary condition of a J^pi group is set by its
# first level, so inserting a dummy level below the others redefines their formal parameters
FLAGS = {"use_brune": True}

failures = []
def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + ("" if ok else f"  -- {detail}")); sys.stdout.flush()
    if not ok: failures.append(name)

def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_levelscan_")
    work = os.path.join(tmp, "13N")
    shutil.copytree(os.path.join(AZURE2_ROOT, "tests", "13N"), work)
    shutil.rmtree(os.path.join(work, "output"), ignore_errors=True); os.makedirs(os.path.join(work, "output"))
    mdl = AzrModel.from_file(os.path.join(work, "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    full = mdl.write(os.path.join(work, "full.azr"))
    mdl.remove_level(jpi="3/2-", energy=3.5032)          # a real removal (deactivate_level only zeroes the widths)
    base = mdl.write(os.path.join(work, "no32m.azr"))
    check("the 3/2- 3.5032 level is gone from the file", not AzrModel.from_file(base).find(jpi="3/2-", energy=3.5032))
    policy = BoundsPolicy(bg_energy=15.0, e_window=0.5)
    def scan_and_screen(x_inc, obj_inc, keys_inc, title, must_rank_first):
        print(title)
        with SingleSession(single_variant(base, "scan", tmp), cwd=tmp, name="scan", flags=FLAGS) as ev1:
            xs = take(x_inc, keys_inc, ev1.keys, ev1.x0)
            cands = candidate_energies(ev1, xs, z=2.5, min_run=3, top=5)
        for c in cands:
            print(f"   Ex {c.ex:.3f} MeV  chi2 in runs {c.score:.0f}  {c.npts} pts  segments {c.segments}")
        near = [c for c in cands if abs(c.ex - 3.5032) < 0.08]
        check(f"{title}: a candidate within 80 keV of the removed level", bool(near), str([round(c.ex, 3) for c in cands]))
        ex = near[0].ex if near else 3.5032
        var = template_variant(base, os.path.join(work, "templates.azr"), ex, jpis=JPI_GROUPS,
                               output_dir=os.path.join(work, "out_templates"))
        # the screening grid the campaign applies (5/20): at the converged grid the 13N
        # variant costs ~7 s per evaluation, at 5/20 ~0.2 s (measured 2026-09-13/14)
        from rmfit.campaign import Campaign
        print(f"   screening grid 5/20 applied to {Campaign.set_grid_setting(var['path'], 5.0, 20.0)} <targetInt> line(s)")
        ev2 = Sharded(var["path"], os.path.join(tmp, "shards_t"), nshards=3, flags=FLAGS)
        x_map = incumbent_on_variant(ev2, x_inc, keys_inc, var)
        labs = dummy_labels(ev2, var)
        dummy_w = [(str(k), float(x_map[j])) for j, k in enumerate(ev2.keys) if k.kind == "width" and k.level in set(labs.values())]
        unmapped = [k for k in keys_inc if k not in set(ev2.keys)]
        print(f"   variant ex {var['ex']:.6f}; dummy labels {labs}; nonzero dummy widths in x_map: {[(k, v) for k, v in dummy_w if v != 0]}; "
              f"incumbent keys absent from the variant: {[str(k) for k in unmapped][:5]}; objective at x_map {ev2.objective(x_map):.1f} vs incumbent {obj_inc:.1f}")
        dl = set(labs.values())
        diff = [(str(k), float(x_map[j]), float(ev2.x0[j])) for j, k in enumerate(ev2.keys) if k.level not in dl and abs(x_map[j] - ev2.x0[j]) > 1e-9]
        print(f"   variant at its own file values: {ev2.objective(ev2.x0):.1f}; non-dummy keys where x_map != file value: {len(diff)} e.g. {diff[:4]}")
        check(f"{title}: zero-amplitude dummy levels change the objective by < 0.2%",
              abs(ev2.objective(x_map) - obj_inc) < 2e-3 * obj_inc, f"{ev2.objective(x_map)} vs {obj_inc}")
        rows = template_screen(ev2, x_map, var, gamma_total_ev=60e3, scales=(0.5, 1.0, 2.0), max_signs=16, gn_top=12,
                               policy=policy, ex_offsets=(-0.06, -0.03, 0.0, 0.03, 0.06), log=print)
        ev2.close()
        best = best_per_jpi(rows)
        print("   ranking: " + "; ".join(f"{r['jpi']} {r['desc']} @{r['ex']:.3f} gain {r['gain']:+.0f}" for r in best[:5]))
        if must_rank_first:
            check(f"{title}: some pattern improves on the incumbent", best[0]["gain"] > 0, str(best[0]["gain"]))
            check(f"{title}: 3/2- ranks first", best[0]["jpi"] == "3/2-", best[0]["jpi"])
            cand = to_level_candidate(best[0])
            check("candidate carries signed widths", cand.jpi == best[0]["jpi"] and len(cand.channels) >= 1 and cand.energy_fixed)
        return best

    # A. the level removed, every other parameter untouched: the pure signature of a
    #    missing level -- the screen must identify it
    sh = Sharded(base, os.path.join(tmp, "shards"), nshards=3, flags=FLAGS)
    x0 = sh.x0.copy(); obj0 = sh.objective(x0); keys0 = list(sh.keys)
    print(f"  model without the 3/2- level at the original parameters: {obj0:.1f}")
    # sanity of the template amplitudes: the removed level's own parameters through the dummy
    var = template_variant(base, os.path.join(work, "templates_exact.azr"), 3.5032, jpis=["3/2-"],
                           output_dir=os.path.join(work, "out_exact"))
    ev2 = Sharded(var["path"], os.path.join(tmp, "shards_e"), nshards=3, flags=FLAGS)
    lab = dummy_labels(ev2, var)["3/2-"]
    x_map = incumbent_on_variant(ev2, x0, keys0, var); ref0 = ev2.objective(x_map); xe = x_map.copy()
    # the particle amplitude goes through the SAME conversion template_screen uses since
    # 2026-09-14 -- gamma = gamma_W sqrt(Gamma / Gamma_W) from the engine's Wigner quantities --
    # not the file-derived x0, which on a dense model came out ~100x too small (a failed Brune
    # transform of the 1 eV dummy) and made every template inert while this test passed
    gW_amp = ev2.wigner_amplitudes(); gW_wid = ev2.wigner_widths(x_map)
    for j, k in enumerate(ev2.keys):
        if k.level != lab: continue
        if k.kind == "energy": xe[j] = 3.5032
        elif k.kind == "width" and k.pair == 1 and k.L == 1:
            check("the dummy's open particle channel has a Wigner width from the engine", j in gW_wid and gW_wid[j] > 0)
            xe[j] = gW_amp[j] * np.sqrt(55200.0 / gW_wid[j])
            r = abs(ev2.x0[j]) / np.sqrt(var["gamma_ref_ev"]) / (gW_amp[j] / np.sqrt(gW_wid[j]))
            print(f"   x0-based / Wigner-based reference amplitude for {k}: {r:.3g}")
        elif k.kind == "width" and k.pair == 2 and k.L == 1: xe[j] = -abs(ev2.x0[j]) * np.sqrt(0.49 / var["gamma_ref_photon_ev"])
    oe = ev2.objective(xe); ev2.close()
    with SingleSession(single_variant(full, "full", tmp), cwd=tmp, name="full", flags=FLAGS) as evf:
        objf = evf.objective(evf.x0)
    print(f"   exact-parameter template: {ref0:.1f} -> {oe:.1f}; the full model at its own parameters: {objf:.1f}")
    check("the removed level's own parameters through the dummy recover most of the loss",
          oe < 0.3 * ref0 and oe < 1.5 * objf, f"{oe} vs loss {ref0} / full model {objf}")
    best_a = scan_and_screen(x0, obj0, keys0, "A. pure removal", must_rank_first=True)
    # B. the model re-polished without the level (the other levels compensate): informational
    f = polish(sh, x0, bounds=default_bounds(sh, x0, policy), policy=policy, max_nfev=100)
    print(f"  re-polished without the level: {f.objective:.1f} ({f.status})")
    sh.close()
    best_b = scan_and_screen(f.x, f.objective, keys0, "B. re-polished removal (informational)", must_rank_first=False) if os.environ.get('LEVELSCAN_PART_B', '1') == '1' else []
    rank_b = [r["jpi"] for r in best_b]
    print(f"   B: 3/2- ranks {rank_b.index('3/2-') + 1 if '3/2-' in rank_b else 'outside the top'} of {len(rank_b)}")
    shutil.rmtree(tmp, ignore_errors=True)
    print("FAILURES: " + ", ".join(failures) if failures else "ALL OK")
    return 1 if failures else 0

if __name__ == "__main__":
    sys.exit(main())
