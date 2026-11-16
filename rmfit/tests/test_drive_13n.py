#!/usr/bin/env python3
"""rmfit drive end to end on the 13N model (queue, ~30-60 min, 3 shards): build phase (audit, reproduce against a copied
chiSquared.out, spec, anchor, seeds, pairs, candidates, init) then the loop with tiny budgets, legs repeated until DONE /
WAITING / STOP; pending wait decisions are answered by the test ('accept' / 'continue') so the whole machine is exercised.
Run:  python3 rmfit/tests/test_drive_13n.py   (tests/_queue_drive/run_crc_drive_test submits it)"""
import json, os, shutil, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
os.environ.setdefault("OMP_NUM_THREADS", "1")
from rmfit import AZURE2_ROOT
from rmfit import drive as dv
from rmfit.azr3 import Model as AzrModel

FAILS = []


def check(name, ok, info=""):
    print(("  ok   " if ok else "  FAIL ") + name + (f"  -- {info}" if info and not ok else ""))
    sys.stdout.flush()
    if not ok:
        FAILS.append(name)


def main():
    tmp = tempfile.mkdtemp(prefix="rmfit_drive13n_")
    src = os.path.join(AZURE2_ROOT, "tests", "13N")
    shutil.copytree(os.path.join(src, "data"), os.path.join(tmp, "data"))
    os.makedirs(os.path.join(tmp, "output"))
    if os.path.exists(os.path.join(src, "output", "chiSquared.out")):
        shutil.copy2(os.path.join(src, "output", "chiSquared.out"), os.path.join(tmp, "output", "chiSquared.out"))
    mdl = AzrModel.from_file(os.path.join(src, "13N.azr"))
    mdl.set_segment_norm("meyer_84", vary=True, sys_error=5.0)
    # 2026-11-10: one background channel free at exactly zero -- SEED_CHANNELS must give it a placeholder and seed it (a zero
    # width is never an engine parameter; a from-scratch demo ran with no background amplitudes because of that)
    bg = [lv for lv in mdl.active_levels if float(lv.energy) >= 15.0]
    ch0 = bg[0].channels[0]; ch0.gamma = 0.0; ch0.channel_fixed = False
    print("   zeroed free channel:", bg[0].jpi, float(bg[0].energy), "pair", ch0.pair, "L", ch0.L)
    mdl.set_output_dir(os.path.join(tmp, "output"))
    base = mdl.write(os.path.join(tmp, "13N.azr"))
    spec = {"reaction": "12C+p", "compound_nucleus": "13N", "campaign": "drive-13n", "owners": ["test"], "parametrization": "brune",
            "engine": {"flags": {"use_long_wavelength": True, "ignore_externals": False}}, "nshards": 3,
            "model": {"fit_window": {"Ex_MeV": [1.9, 4.5]}},
            "bounds": {"bg_energy": 15.0, "e_window": 0.5},
            "search": {"top_k": 3, "pairs_per_group": 2, "stage1_nfev": 20, "stage2_nfev": 30, "union_nfev": 30},
            "data_sets": [{"file": "meyer_84", "role": "free", "norm": {"nominal": 1.0, "fixed": False, "prior_sigma": 0.05}}],
            "budget": {"polish_nfev": 30, "polish_rounds_max": 2, "sign_rounds_max": 2, "level_tests_max": 2, "level_nfev": 40,
                       "pinned_nfev": 25, "dataset_off_max": 1, "dataset_off_nfev": 30, "cycles_max": 1, "grid_check": False, "figures": 1,
                       "restarts_per_cycle": 2, "restart_nfev": 30, "crossover_passes": 1, "population_rounds_max": 1},
            "stop_rules": {"polish_gain_rel_min": 5e-4, "sign_rounds_without_accept": 2, "level_adds_marginal_in_a_row": 2},
            "structure": {"broad_first_kev": 10.0, "narrow_phase": True, "levelscan_per_cycle": 1, "levelscan_top": 1,
                          "levelscan_gamma_kev": None, "local_release_mev": 0.3, "local_nfev": 20, "recall_levels": 1, "recall_min_pass": 0.5},
            "pre_authorized": {"adopt_level": False, "flip_adopt": True, "shifts_free": False},
            "standing_decisions": [{"date": "2026-11-08", "who": "test", "text": "fixture"}]}
    json.dump(spec, open(os.path.join(tmp, "spec.json"), "w"), indent=1)
    json.dump(dict(candidates=[dict(energy=3.0, jpi=None, energy_fixed=True, channels=[dict(pair=1, L=0, S=0.5, gamma=20000.0)],
                                    source="test", hypothesis="fixture: a level at 3.0 MeV, J^pi scanned")]),
              open(os.path.join(tmp, "candidates.json"), "w"), indent=1)
    t0 = time.time()
    answers = {"reproduce-mismatch": "accept", "audit-high": "accepted", "audit-subagents": "skip", "spec-high": "accepted",
               "review-high": "continue", "adopt-level": "reject",
               # non-blocking (default adopt-manually), but pending() lists it: a queue job stalled on it (2026-11-10)
               "recall-low": "adopt-manually"}
    status, legs = None, 0
    while legs < 12:
        legs += 1
        d = dv.Driver(tmp, base=base, local=True)
        status = d.run()
        print(f"   leg {legs}: {status} at {d.state['state']} ({(time.time()-t0)/60:.1f} min)")
        if status in ("DONE", "STOP", "ERROR"):
            break
        if status == "WAITING":
            pend = d.decisions.pending()
            for p in pend:
                key = next((k for k in answers if p["slug"].startswith(k)), None)
                if key is None:
                    print("   unexpected decision", p.id, p["slug"], p.get("question"))
                    break
                d.decisions.answer(p.id, answers[key], by="test")
                print(f"   answered {p.id} {p['slug']}: {answers[key]}")
            else:
                continue
            break
    st = json.load(open(os.path.join(tmp, "drive_state.json")))
    done = [h["state"] for h in st["history"]]
    print("   states run:", done)
    print("   message:", st.get("message"))
    check("the drive ended DONE", status == "DONE", f"{status} at {st['state']}: {st.get('message')}")
    check("build phase reached INIT (campaign.json)", os.path.exists(os.path.join(tmp, "campaign.json")))
    check("working copy and job script written", os.path.exists(os.path.join(tmp, "13N_drive.azr")) and os.path.exists(os.path.join(tmp, "run_crc_drive")))
    for s in ("RECALL", "TIGHTEN", "SIGN", "POPULATION", "LEVELSCAN", "STRUCTURE_ADD", "DATASET_OFF", "REVIEW", "EXPORT", "CYCLE"):
        check(f"state {s} ran", s in done)
    check("REVIEW.md and export exist", os.path.exists(os.path.join(tmp, "REVIEW.md")) and any(f.endswith(".sav") for f in os.listdir(tmp)))
    check("DRIVE.md written", os.path.exists(os.path.join(tmp, "DRIVE.md")))
    seed = [h for h in st["history"] if h["state"] == "SEED_CHANNELS"]
    check("SEED_CHANNELS seeded the zero background channel (placeholder path)", seed and "seeds" in seed[0].get("note", "") and "nothing to seed" not in seed[0]["note"],
          str([h.get("note") for h in seed]))
    rep = st.get("reproduce")
    check("REPRODUCE compared against chiSquared.out", rep is not None and rep.get("nref", 0) > 0, str(rep))
    pop = [h for h in st["history"] if h["state"] == "POPULATION"]
    check("POPULATION ran the restarts (not skipped)", pop and "restarts" in pop[0].get("note", "") and "skipped" not in pop[0]["note"],
          str([h.get("note") for h in pop]))
    rec = [h for h in st["history"] if h["state"] == "RECALL"]
    check("RECALL removed and searched for a known level (not skipped)", rec and "recall" in rec[0].get("note", "") and "skipped" not in rec[0]["note"],
          str([h.get("note") for h in rec]))
    rtab = " ".join(st.get("findings", {}).get("RECALL", []))
    check("RECALL child campaign ran to a classification (no error)", rtab and "Error" not in rtab and (st.get("recall") or {}).get("n") == 1, rtab[:300])
    check("RECALL located the removed 1/2+ 2.368 MeV level (J^pi, energy, width)", (st.get("recall") or {}).get("n_located") == 1
          and "as 1/2+ @ 2.36" in rtab, rtab[:300])
    print("   recall table:", rtab[:400])
    print("   recall:", st.get("recall"), [h.get("note") for h in rec])
    ls = [h for h in st["history"] if h["state"] == "LEVELSCAN"]
    check("LEVELSCAN screened the residuals once the bank was exhausted", any("screened" in h.get("note", "") for h in ls),
          str([h.get("note") for h in ls]))
    print("   levelscan:", [h.get("note") for h in ls])
    print(f"   total {(time.time()-t0)/60:.1f} min in {legs} legs; dir {tmp}")
    if FAILS:
        print(f"FAILED {len(FAILS)}: {', '.join(FAILS)}")
        sys.exit(1)
    shutil.rmtree(tmp, ignore_errors=True)
    print("all drive 13N checks passed")


if __name__ == "__main__":
    main()
