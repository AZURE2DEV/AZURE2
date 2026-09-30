---
name: rmfit-campaign
description: The rmfit process for global R-matrix fitting with AZURE2 -- start a campaign on any reaction, drive or resume one (sharded pyazr evaluation, bound-tightening polish, interference-sign rounds, level add/remove rounds, dataset-off diagnostics, ledger + STATUS.md), judge results, export for GUI review, and refine the method. Use whenever the task mentions rmfit, a campaign directory, STATUS.md, sign rounds, a tightening polish, a level-add test, fitting a new reaction "the rmfit way", or continuing a fit that a previous session left running on the cluster.
---

# The rmfit process

`rmfit` (`/groups/rdeboer1/user/rdeboer1/R-matrix/rmfit/`, README there) does what an
evaluator does by hand -- polish, flip interference signs, add or remove levels, switch
suspect data sets off -- as resumable steps that a script runs on one cluster node and a
fresh session can pick up from `STATUS.md`. It is not installed as a package:

```bash
export PYTHONPATH=/groups/rdeboer1/user/rdeboer1/R-matrix
python3 -m rmfit.cli <campaign_dir> <command> [options]
```

Read this file top to bottom before the first command on a new reaction; for resuming a
campaign, "Read the state" and "The loop" are enough. Refinements go in the log at the end.

## 1. Read the state, never recompute it

A campaign directory carries everything a fresh session needs:

- `STATUS.md` -- best verified objective per structure, recent rounds and log lines.
  Regenerate with `status`; `status --json` for a compact dump.
- `readme` -- the lab notebook (rmfit appends one paragraph per round; add your own
  with `log "..."`). Archive convention: read it before touching anything.
- `campaign.log` -- every line the driver printed; `rmfit_sge.log` / `rmfit_*.log` --
  cluster job logs.
- `ledger.sqlite` + `candidates/*.npz` -- every candidate vector and score. The
  incumbent is the best candidate with a `verified` score (a fresh single-session
  evaluation of the full model). Sharded numbers rank; verified numbers are reported.
- `campaign.json` -- base model, engine flags, shard count, bounds policy, search config.
- `<reaction>.azr` + `.sav` -- the GUI-reviewable export of the incumbent.
- `variants/` -- shard files, window models, with/without-level files (disposable).

`qstat -u $USER` first: one campaign command at a time per directory (they share the
ledger and `variants/`). The login node is for `status`, `report`, `log`, `export`,
short tests (<= 6 workers, <= 30 min). Everything else runs as a job.

## 2. Starting a campaign on a new reaction

Prerequisites (do these in the reaction's fit directory, not in the campaign):

1. A `.azr` that opens in the GUI, with `data/`, an `output/` directory, and a `readme`.
   Know the CLI-mode flags this reaction uses (`--use-brune`, `--ignore-externals`,
   `--no-long-wavelength`, ...): they are not in the `.azr`, they vary per reaction
   (`~/R-matrix/CLAUDE.md`), and rmfit needs them as `--flags` JSON.
2. Reproduce the last finished fit before building on it: load its `param.sav` by name
   (`pyazr`), confirm the objective matches `chiSquared.out` (chi2 + norm penalty), then
   bake it into a clean file (`rmfit.engine.bake`) and check it with the CLI in mode 1
   with a blank external file. Never trust a bake that was "verified" with the `.sav`
   as the external file (that overrides `<levels>`), and never use parameters from a
   running or crashed fit -- only the last completed one. Purge inactive level lines.
3. Merge any new data sets onto the baked file (`AzrModel.add_data_segment`), capped
   at the energy range the level scheme covers (`set_segment_energy_range`), with
   normalizations free and a systematic error that reflects the paper. Look at every
   new file: a single negative or spiked cross section can dominate the objective.
4. Record all of that in the reaction directory's `readme`.

Then, in a new dated campaign directory (`<reaction>/<M>-<D>-<YY>_rmfit_campaign/`):

```bash
mkdir <dir> && cp <base>.azr <dir>/ && cp -r data <dir>/data && mkdir <dir>/output
python3 -m rmfit.cli <dir> init --base <dir>/<base>.azr [--sav <base>.sav] --nshards 24 \
    --flags '{"use_long_wavelength": false, "ignore_externals": true, "use_brune": true}' \
    --policy '{"theta2_cap": 1.0, "bg_theta2_cap": 3.0, "bg_energy": <Ex above the data>, "e_window": 0.3, "norm_range": [0.2, 5.0]}' \
    --search '{"top_k": 6, "pairs_per_group": 4, "stage1_nfev": 30, "stage2_nfev": 60, "union_nfev": 100}'
```

Optional `--objective JSON` (campaign.json `objective`, deBoer et al. 2017 Sec. on fitting):
`{"dataset_weight": "reduced", "dataset_key": "file", "loss": "sivia"}`.
`dataset_weight: reduced` is the reduced-chi-squared method (every data set's chi-squared
divided by its number of points, rescaled by N/K so the objective keeps a chi-squared
magnitude; penalty rows unweighted) -- it stops a large set with small errors from
dominating; `loss: sivia` replaces chi-squared by Sivia's broader PDF,
`-2 log[(1 - e^{-R^2/2})/R^2] - 2 log 2` (~R^2/2 for small residuals, ~2 log R^2 for
outliers: outliers are down-weighted as if their error bar were inflated). The Sivia loss
is applied to residuals *standardized by each data set's reference chi2/N* (`loss_scale:
dataset`, scales computed by `init` from the plain baseline and stored in campaign.json):
without that, on a model with chi2/N ~ 20-300 every point sits in the logarithmic tail,
the data lose ~10x their weight against the quadratic norm penalties and the fit drifts
(13C+a, 2026-09-08: plain chi2 701k -> 1,190k in one polish, MANA abandoned, norms snapped
to nominal). Both methods are "statistically incorrect but useful when error bars are
underestimated". The objective
is a property of a campaign: candidates verified under different objectives are not
comparable, so a change of objective means a new campaign directory (seed it with the
previous export). `REPORT.md` and `segment_scores` keep the plain per-dataset chi-squared
whatever the objective. Verified in `rmfit/tests/test_objective_13n.py`.

`init` verifies the baseline (one full-model evaluation), registers the structure,
writes `STATUS.md` and the job script `run_crc_rmfit` (24-core `long` node, mails
`$USER@nd.edu`, `-m abe`). Re-running `init` keeps an existing config unless new values
are given. Policy meanings: `theta2_cap` is the Wigner-limit box for physical levels,
`bg_theta2_cap` for background poles (levels with Ex >= `bg_energy`, whose energies stay
fixed) -- 3.0 is a starting guess, not a physics default: check where the baseline's
background amplitudes already sit (`init`'s cap audit) before trusting it; 12C+alpha
elastic needed ~100 (see the 2026-09-30 log entry); `e_window` is the energy freedom of a level per polish (MeV); `norm_range`
multiplies the nominal normalization. `narrow_window` (optional, e.g.
`{"gamma_kev": 5, "factor": 3, "min_kev": 1}`) keeps a level whose free partial widths sum
to less than `gamma_kev` within ±max(`min_kev`, `factor`·Γ) of its current energy instead of
±`e_window`, so a polish cannot carry a narrow level away from its peak (the 13C+a 9/2+
drifted 180 keV in one tightening polish). Shard count = cores; sharding is exact
(chi-squared is additive over segments), so more shards only buy speed.

Cost check on a new model (`campaign.log` after `init`): one sharded objective should
be a few seconds and a Jacobian tens of seconds. If a segment with target-effect
convolution dominates (as Cierjacks did on 13C+a), rmfit folds it itself (`model.ConvSpec`);
composite "sum" segments are allowed as convolution pieces.

## 3. The loop

Every step is a job (`qsub -v RMFIT_CMD="..." run_crc_rmfit`, or `-hold_jid <id>` to queue
behind the previous one). Every step ends by rewriting `STATUS.md` and appending to the
readme; read both before the next step.

| step | command | purpose | when done |
|---|---|---|---|
| 1 | `polish --nfev 150 --tighten` | walk amplitudes beyond the theta^2 caps into the box (x10 stages), then converge | all theta^2 inside the caps; gain < 0.05% per stage |
| 2 | `round --rounds 2` (repeat) | interference-sign rounds: screen every canonical flip (singles, exposure pairs, same-channel pairs, optional channel patterns), polish the top-k sign-locked then released, greedy union, verify | improving-on-screen and accepted counts fall round over round; stop after two rounds with nothing accepted |
| 3 | `polish --nfev 150 --shifts-free` | warm-start continuation with energy shifts free | gain < 0.05% |
| 4 | `structure add --candidate i` (per entry of `candidates.json`) | forward addition: window J^pi scan at pinned energy, pinned then released full-model polish, control polish without the level, classify, adopt if accepted | all candidates classified |
| 5 | `structure dataset-off --file <substr>` for the worst chi2/N sets | deactivate + refit; reports the extra recovery of the *other* sets beyond the set's own chi2 (the Botek test) | a table for the evaluator |
| 6 | `round --rounds 2` after any structure change | signs may change with the structure | as in 2 |
| 6b | `crossover --donor <cid> [--donor-dir <other campaign>]` | basin crossover: apply another good fit's sign pattern per channel / J^pi group / all and polish every proposal; use it whenever the ledger (or a sibling campaign, e.g. a recovery benchmark) holds a second verified basin | every donor pattern polished; verified union |
| 7 | `export`, then GUI review | bake the incumbent to `<reaction>.azr` + `.sav`, verified in a fresh process | CLI mode 1 with the `.sav` as external file reproduces the totals |

`candidates.json` (write it before step 4): `{"candidates": [{"energy": Ex_MeV, "jpi":
"5/2+" or null to scan, "energy_fixed": true, "channels": [{"pair": n, "L": l, "S": s,
"gamma": width_eV_with_sign}, ...], "source": "...", "note": "..."}]}`. Sources: ENSDF/TUNL
(`pyazr.nds.fetch_levels`), the reaction's literature tables, and residual runs
(`diagnose.residual_runs`) in the fit's worst segments. Omitted channels are seeded at a
small fraction of the Wigner limit.

Residual-driven level search (Phase 7, `rmfit levelscan [--ex E ...] [--top K] [--gamma keV]
[--scales 0.5,1,2] [--ex-offsets -0.04,0,0.04] [--jpi ...] [--refit N]`): candidate energies
from the incumbent's residual runs clustered in Ex across segments (plus dispersion-shaped
pairs: two adjacent clusters of opposite sign with the level between them); at each, one
dummy level per J^pi (all of the group's channels, energy free) is appended and every
(J^pi, width partition, sign pattern, energy offset, scale) is evaluated frozen in the
full model, the best `--gn-top` get a GN step on the level's amplitudes + norms, and the
best per J^pi is reported (ledger round kind `levelscan`); `--refit N` sends the top N to
the ordinary level-add refit. Photon channels get a small trial width with both signs.
Validated on 13N: a level removed with the other parameters untouched is found (its J^pi
ranks first, gain ~90% of the loss); after the rest of the model has re-polished without it
the single-level signature is gone (the compensation limit in the plan) and only the refit
of the top candidates can recover it. Template amplitudes are the engine's own conversion
of a tiny reference width scaled by sqrt(Gamma); candidate energies are rounded to 1 eV
before writing (see the AZURE2 skill's add_level gotcha).

Every level test and relocation writes a review figure to `figures/` (data with the
incumbent and the candidate calculation around the energy: excitation functions per
segment, angular distributions for single-energy segments) -- DeBoer's visual check.

Relocation of an existing level (`rmfit relocate --level 9/2+#1 --ex 8.4654 [--window 5]`):
energy kick to `--ex`, then a full re-polish with that level's energy confined to
+-`--window` keV, verified and adopted if it beats the incumbent by tau. Use it for a
narrow level that has drifted from its peak (see `narrow_window` above for keeping it
there afterwards). 13N test: a level shifted by 60 keV comes back and recovers the full
model's objective.

Ordinary additions: `touch <dir>/STOP` ends a `round --rounds N` loop after the current
round; `verify <cid>` re-scores a candidate; `log "..."` appends to the readme;
`report` writes `REPORT.md` (per-dataset chi2/N, theta^2 table, structure decisions).

## 4. Judging results

- Acceptance thresholds scale with s^2 = objective / N because errors are
  underestimated (s^2 ~ 15-20 on 13C+a). A flip needs a verified gain
  >= max(0.1% obj, 5 s^2 ln N) *beyond the control polish of the incumbent* (the
  incumbent keeps improving on its own by hundreds per polish; flips are credited only
  with what they gain beyond that). A level needs >= max(s^2 k ln N, 0.2% obj), theta^2
  inside the caps, an energy that stays in its window, and a redistribution ratio
  rho <= 0.5 (rho = chi2 lost by some data sets / gained by others).
- Trial outcomes: accept / marginal / worse / reverts (the sign walked back in the
  released polish; never retried on that structure) / drift (the sign walked back but
  the polish found a vector better than the control by >= tau: kept, joins the union)
  / unphysical / collapse (widths -> 0) / blowup / relocated (energy left its window)
  / live-harmful / crashed.
- Convergence of the sign search: improving-on-screen counts decay geometrically
  (archive history 16 -> 4 -> 3 -> 1). Two rounds without an acceptance = converged for
  the present structure.
- `sanity` flags theta^2 over the cap, |E| > 60 MeV, norms outside (0, 100), a
  non-finite Brune transform, and parameters at a bound (listed in the readme entry).
  Any |value| > 1e6 or |E| > 1000 MeV means catastrophic cancellation in
  total-cross-section chi2, i.e. a fake improvement.
- Read `REPORT.md`'s per-dataset table after every structure step: a data set whose
  chi2/N is far above the rest, or whose removal recovers more than its own chi2 in
  the others, is an evaluation question (normalization, energy calibration, resolution),
  not a fitting question.

## 5. What the human decides

rmfit hands these to the evaluator rather than deciding them:

- normalization scale and systematic error of a new data set (and any rows dropped);
- an angle- or energy-dependent trend in fitted normalizations (the norms absorb it
  silently);
- whether a data set that distorts the rest stays in the fit;
- the energy range of the model (data above the level scheme's range are capped, not
  fitted);
- accepting a level whose gain is statistically clear but whose J^pi or energy
  contradicts the literature.

The GUI review: `printf '1\n<reaction>.sav\n\n7\n' | AZURE2 <flags> --no-gui --no-readline
<reaction>.azr` fills `output/` with the fit (the `.sav` carries the fitted norms), then
`AZURE2 <reaction>.azr &` on an X display.

## 6. Diagnostics and benchmarks

- `keys.canonical_signs(x, keys, fixed)`: the sign pattern modulo the exact symmetries
  (all amplitudes of one level; all amplitudes of one particle pair in every level;
  fixed nonzero amplitudes tie a level's gauge to a pair's). Only valid for particle
  channels -- see limits.
- `keys.sign_distance(x1, x2, keys, fixed)`: the minimal number of flips between two
  patterns modulo those symmetries.
- `perturb --out <dir> --flips 8 --jitter 0 --ekick 0 --polish-nfev 150`: the recovery
  benchmark. Flips high-exposure signs of the incumbent, polishes to a genuine worse
  minimum in a sibling campaign directory, and records the target; every later
  verification there logs the sign distance and the gap. Success = back to the target
  within 1e-4. Run it once per reaction before trusting the method there.
- The oracle script (`13C+a/9-6-26_rmfit_recovery_signs/oracle_patterns.py`) applies a
  known target's sign pattern per channel / per group / all at once and polishes: it
  separates "finding the pattern" from "reaching the basin from the pattern".
- `diagnose.dataset_table`, `theta2_table`, `residual_runs`, `compare` feed `report`.

## 7. Known limits (September 2026)

- Validated only on 13C+a (particle channels only, `ignore_externals`). With photon
  channels: the pair-flip symmetry and the shard machinery (`intEC.*` per shard output
  directory) are untested -- verify sharded == single-session at three vectors before
  trusting a round (`init` does it at the baseline).
- The recovery benchmark on 13C+a passed once (signs-only perturbation, 6.4% worse
  start: four rounds + polish -> 0.25% *below* the target, in a different basin, ~15 h on
  one node). Recovery from a start that also has jittered magnitudes and kicked energies
  is untested since the escalations; run `perturb` on every new reaction before trusting
  the method there.
- Energy shifts are finite-differenced (kept fixed during exploration, freed in the
  continuation polish). Norms of split segments are shared across shards correctly;
  segments with target effects are never split.
- Cost on 13C+a (35k points, 330 free, 24 shards): objective 0.6-1.5 s, Jacobian ~20 s,
  a 60-iteration polish ~9 min, a sign round 2-4 h, a level test ~1.5 h.

## 8. Gotchas

- Never open two pyazr sessions in one process (parameter tables desync); rmfit runs
  one process per session and `bake()` verifies in a fresh process. Do not use
  `save_fit` for campaign files. `pyazr` bakes skip inactive levels since 2026-09-05;
  older baked files may have mis-placed values in J-groups with inactive lines.
- Delete `output/intEC.*` whenever segments change (stale-cache rule of the AZURE2 skill).
- New levels are appended at the end of the file and levels are removed by
  deactivation, so every existing parameter key stays valid across structures.
- A released polish lets amplitudes cross zero: sign patterns are not preserved by
  polishing, and a "basin" is only defined up to the amplitudes that are near zero.
- Frozen screens (flip, no refit) are cheap but favour flips of small amplitudes;
  the GN screen (one damped Gauss-Newton step on the group's columns) is the ranking
  that matters; coordinated multi-level flips need several backtracking GN iterations
  (`search.gn_iters_multi`) and their own polish quota (`search.multi_top`).
- Job mail: `-m abe` (start, end, abort) so a silent crash is noticed.

## 9. Refinement log (append dated entries; code changes get a test in rmfit/tests)

- 2026-09-05 -- pyazr `engine_level_keys` counted inactive levels: the live 13C+a 5/2+
  block was corrupted by every bake. Fixed with a regression test; `active_datasets`
  fixed `penalties()` / `dataset_chi2()` mis-indexing with inactive segments.
- 2026-09-05 -- sharding: the convolved Cierjacks segment (17 of 23 s per evaluation)
  is folded by rmfit on a fine grid; agreement with the engine 4e-7.
- 2026-09-06 -- flips were credited with the incumbent's own continuation gain: added
  the control polish. Frozen sign screens are not predictive; GN screens are.
- 2026-09-06 -- level tests: seed only the new level's keys; window scans sharded
  (33 min per level); the 5/2+ 10.55 MeV Goldberg level accepted on 13C+a (rho 0.32),
  10.52 / 10.53 5/2+ and 9.862 5/2- rejected; BandH_no_narrow distorts the fit (+7.9%).
- 2026-09-07 -- recovery benchmark: an 8-flip perturbation polished to a minimum sits
  29 signs from its origin (`sign_distance`), clustered by channel; single/pair rounds
  gain a few thousand per round without reducing the distance. Added same-channel pair
  flips, per-channel pattern enumeration (`enumerate_channels`), a per-group frozen
  quota split by objective and exposure, backtracking GN iterations and a polish quota
  for multi-level patterns, persisted screen tables, `STOP` files, the oracle
  diagnostic. Outcome pending (jobs 1420848, 1420849).
- 2026-09-07 -- MANA 13C(a,a): Heil's shapes at ~0.67x scale with 1% errors; three bad
  rows removed; norms free at 20% come out 1.24-1.53 rising with angle -- an evaluation
  question, not a fit question.
- 2026-09-07 -- round 3 of the signs benchmark (multi-level screen): 17 of 404 improve on
  the screen; a single flip (5/2-#7 p1 L3) gains 22,400 that rounds 1-2 had not found
  (its screen value changed with the incumbent), while a reverted 4-level pattern still
  ended 1,579 better than the control -- added the `drift` outcome so such vectors are
  not thrown away.
- 2026-09-07 -- Validation A PASSED on the signs-only benchmark: 777,349 -> 771,193 ->
  769,848 -> 747,364 -> 730,814 -> 728,879 (target 730,697). What made the difference:
  the control-polish continuation (the round-3 basin needed 120 more iterations), the
  multi-level screen with backtracking GN iterations plus a polish quota (the winning
  4-level 5/2+ p5 L1 pattern screened at 748,229 against an incumbent of 747,364, i.e.
  it would never have made a top-k list), and single flips whose screen value only
  became favourable after the incumbent moved. Lesson: judge a sign-round result by
  the redistribution ratio too -- the final 3,148 gain over the target was a MANA-vs-Heil
  trade (rho 0.90), an evaluation decision rather than a better fit.
- 2026-09-07 -- `verify_candidate` from a heredoc script fails (`spawn` needs a real
  `__main__` file): write helper scripts to a file with an `if __name__ == "__main__"` guard.
- 2026-09-07 -- oracle on the 728,879 incumbent: the *old* incumbent's 7/2+ p5 L1 S2.5
  pattern (3 flips) polishes to 716,659 (+11,930 beyond the control) although its GN
  screen was 791,423 -- no screen ranks such patterns; only a polish tells. Added the
  basin-crossover round (`crossover --donor`, `moves.crossover_moves`,
  `search.crossover_round`): every differing channel/group pattern of a donor basin is
  polished. Two good basins in the ledger are worth more than any screen.
- 2026-09-08 -- first crossover round on 13C+a: 33 donor patterns, 8.1 h, 6 accepted/drift,
  union kept 2 -> 707,332, polish -> 706,809 (3.2% below the campaign's sign-converged
  730,697; rho 0.60, broad gains, Heil aa the one loser). Recipe that produced the best fit
  of this campaign: sign-converge -> `perturb` benchmark (finds a second basin) -> transplant
  -> `crossover --donor <old incumbent>` -> polish. Consider running `perturb` + crossover
  routinely, not only as a benchmark: a second basin is the most valuable proposal source.
- 2026-09-08 -- DeBoer's guidance on underestimated error bars (BandH stays in): added the
  configurable objective (`ObjectiveSpec`: reduced-chi2 per data set, Sivia loss) from
  deBoer et al. 2017; the MANA row cuts stay; MANA judged more accurate than Heil in shape.
  New campaign `13C+a/9-8-26_rmfit_robust/` seeded from the 706,809 export runs it, then
  a perturb + crossover cycle (donors: the perturbed sibling and the chi2 campaign's
  candidates 49 / 86 via `--donor-dir`).
- 2026-09-08 -- the unscaled Sivia loss is degenerate when chi2/N >> 1 (see section 2):
  job 1426076 stopped after its first polish (objective 53k -> 40k while plain chi2
  701k -> 1,190k). Fix: `loss_scale: dataset` (c rho(z/c), c = w_i s_i^2 with s_i^2 the
  baseline's chi2/N of the data set); `init` computes the scales. Lesson: any robust loss
  needs the residual scale of the data it is applied to; test on the real model's
  chi2/N before launching a campaign.
- 2026-09-08 -- objective scan on 13C+a (`9-8-26_rmfit_robust/objective_scan.md`): the
  reduced-chi2-per-data-set weighting is unusable when data-set sizes span orders of
  magnitude (it hands the fit to the small sets; tempering with `weight_power` 0.5 does
  not save it); Sivia scaled per data set abandons whole sets; Sivia with ONE global scale
  (the overall chi2/N, `loss_scale: global`) is the coherent "poorly described sets have
  underestimated errors" objective and was chosen. Rule: run `objective_scan.py`-style
  150-iteration polishes and read the per-file table before committing a campaign to a
  non-chi2 objective; and remember that basins found under any objective are reusable as
  crossover donors under plain chi2.
- 2026-09-08 -- freeing normalizations is a structure change: edit a copy of the base
  file (`AzrModel.set_segment_norm(<file substring>, vary=True, sys_error=<percent>)`
  hits every segment of that file), copy the `.sav` alongside (new norms start at
  nominal), re-`init` (fresh ledger). Crossover donors from a campaign without those
  norms still work: `run_crossover` fills norm/shift keys the donor lacks from the
  incumbent and only refuses a donor missing R-matrix parameters. 13C+a: only 1 of the
  28 Heil (a,a) angles had a free normalization; DeBoer saw the scaling problem in the GUI.
- 2026-09-09 -- PLANNED (DeBoer): residual-driven level search, `rmfit levelscan`: residual
  peaks -> candidate energies; per energy a template screen of every J^pi x (l,s) width
  partition x sign pattern with the level pinned at an estimated total width inside the
  current model (frozen evaluation + one GN step on the level's amplitudes and norms;
  Legendre-coefficient shape comparison per angular segment), then refit of the top few
  with the existing level-add machinery, and a data/fit/candidate angular-distribution
  figure at the peak. Note the physics: the response to a new level is quadratic in its
  amplitudes and saturates at resonance, so templates need a finite trial width -- a
  linear matched filter would be identically zero. Design in the plan file
  (`~/.claude/plans/cached-snacking-bubble.md`, Phase 7). ML surrogate deferred.
- 2026-09-09 -- narrow levels drift: the 13C+a 9/2+ (Gamma_n ~0.4 keV, seen only by the
  111-point Cierjacks window with a resolution function) sat 10 keV off its peak in the
  baseline and slid 180 keV during the first tightening polish, when the objective was
  dominated by the not-yet-normalized MANA set; once more than ~1 keV from the peak a
  narrow level has no gradient pull and never comes back (rmfit's trust-region steps have
  no "initial step" problem, but the same flat surface). Remedies to build (Phase 7 with
  the level search): per-level energy windows of a few times max(Gamma, resolution) for
  narrow levels in `default_bounds`, and a relocation move to the nearest residual peak;
  narrow levels can also be mapped by direct (E, Gamma) scans
  (`9-9-26_seg92_9halfplus/grid2d.py`).
- 2026-09-10 -- robust campaign (global-scaled Sivia, free Heil norms) finished at 183,309
  (plain chi2 931,838): sign rounds empty under this objective too, the crossover with the
  chi2 campaign's old basin gave -1.2%, the perturb + crossover cycle nothing. Objective G
  behaves as designed (ND 2021 chi2/N 16.8 -> 11.2, Heil 18.7 -> 12.3; BandH x4.6, MANA
  x1.65). Lesson: with an already sign-converged seed, a perturbation only 1% worse does not
  produce a genuinely different basin; make the perturbation larger (more flips, jitter)
  when the aim is a donor basin rather than a recovery test.
- 2026-09-10 -- with the error bars corrected by hand (BandH x2, MANA 3% systematic in
  quadrature, Heil norms free) the plain chi2 polish gains in every data set (510,640 ->
  412,816, losses < 2k) while the Sivia loss keeps trading the corrected sets away: once
  the evaluator has fixed the uncertainties, use plain chi2; the robust loss is only a
  stand-in for uncertainties nobody has looked at yet. Record every error-bar change in
  the readme with the formula, before/after statistics and file checksums.
- 2026-09-10 -- `rmfit levelscan` built (rmfit/levelscan.py, test_levelscan_13n.py). Lessons:
  the engine numbers levels inside a J^pi group in its own order (resolve dummy labels from
  `ev.level_energy`, never from file order); pyazr's `deactivate_level` keeps a zero-width
  level with a free energy in the file (use `remove_level` for a real removal); new J^pi
  groups add hard-sphere terms; photon channels need their own tiny trial widths; and a
  16-digit energy in an add_level line silently corrupts the model. Not yet run on 13C+a.
- 2026-09-10 -- `levelscan.review_figure` + `Campaign.review_figure`: written by `relocate`
  and `structure add` (figures/). 13N check: analyzing-power segments plotted against angle
  at their energy, excitation functions against Ex; the relocated level visibly matches
  the Baumann A_y distributions the shifted one misses.
- 2026-09-11 -- the exported `.azr` of a campaign carries the fitted `<levels>` but NOMINAL
  normalizations: anything that re-evaluates an export (a gate script, a side fit) must seed
  the norms from the companion `.sav` by parameter name, or the objective is wrong. Queue a
  follow-up analysis with `qsub -hold_jid <campaign job>` so it runs the moment the node frees.
- 2026-09-11/12 -- PERTURBATION SIZE, both ends measured on the same sign-converged 13C+a
  fit, and neither worked. 8 flips landed 1.0% above the incumbent; its rounds and crossover
  found nothing (12 h). 12 flips landed 22% above it -- a genuinely different basin, but a
  dead one: two sign rounds plus two polishes moved it only 496,184 -> 496,139 (0.009% over
  9.5 h, sign distance to the target 32 -> 31), and the crossover of its best candidate back
  onto the incumbent screened 0 of 24 patterns improving and accepted 0 of 24 after the
  polishes (406,878.3 -> 406,873.9, i.e. control drift). So a 12-flip basin is both
  unrecoverable by local search and useless as a donor: the flips destroy the shared
  structure that makes a donor pattern transplantable. Read the geometry, not just the
  objective gap: inside a wrong sign basin the surface is nearly flat (0.009% over two
  rounds), so the gap tells you nothing about whether search can cross back. Practical rule
  until a better one is measured: for a donor basin, perturb 2-4 flips inside ONE J^pi group
  and keep the rest of the pattern intact, rather than scattering flips across groups; use a
  large scattered perturbation only as a from-scratch recovery benchmark, and expect it to
  need a structure move or a fresh seed, not sign rounds. Always verify the `perturb` line's
  value AND, after the first round, whether the sign distance to the target is falling at
  all -- if it moves by 1 flip in a round, stop the arm.
- 2026-09-11 -- RETRACTION (DeBoer challenged the claim, and he was right): an earlier entry
  here said MINUIT (CLI mode 2) "cannot fit such an energy at all (10% initial step)".
  The hard-coded initial steps are real (`src/CNuc.cpp:1562` uses 0.1*E for a level energy,
  `:1579` 0.1*gamma for a width), but they do not prevent the fit. A basin test on the
  isolated 9/2+ level in Cierjacks segment 92 (energy free, one free neutron width,
  everything else fixed; `9-9-26_seg92_9halfplus/minuit_basin.py`) converges to
  E = 8.4653 MeV, Gamma_n = 479 eV from every start in 8.4620-8.4700 MeV, including a
  start with Gamma_n = 18.77 eV. It fails only from 8.4600 MeV (collapses to 53 eV) and
  runs off from 8.4800 MeV. The original failing case reproduces as a success
  (`minuit_repro.py`). What actually moved that level 8.4554 -> 8.6376 MeV was rmfit's own
  bound-tightening polish while the objective was dominated by the un-normalized MANA set.
  Lesson for the method, not for MINUIT: a narrow level needs a per-level energy window,
  and a single anomalous observation about a widely used code is not evidence -- test it
  across starting points before writing it down.
- 2026-09-12 -- a model variant written into a SUBDIRECTORY needs the base model's `data/`
  linked beside it, because `<segmentsData>` paths in an .azr are relative. `structure.py`
  and `windows.py` already did this; `levelscan.template_variant` did not, so the first
  13C+a gate run died on `scan_A/data/Fow_Joh_Fre_EXFOR_ntotal.dat` after 3 minutes of
  session build. Fixed in `rmfit/levelscan.py` (same 8-line idiom). The 13N test never
  caught it because it writes its variant BESIDE the base model, where the link already
  exists -- when a helper takes a path, test it with a path in a fresh subdirectory, not
  just a sibling name. Related: `qsub` a csh job script whose python call passes a value
  starting with a minus (`--offsets -0.05,0,0.05`) and argparse reads it as an option; use
  the `--opt=value` form. Both failures cost a node handover each, and both would have been
  caught by running the job script's exact command line once on the login node first.
- 2026-09-12 -- THE TEMPLATE SCREEN AND THE NARROW-RESONANCE INTEGRATION FIX INTERACT
  BADLY, and this is now the binding cost problem. Since commit a4095a6 the adaptive
  target-effect grid actually resolves a narrow resonance, so an evaluation whose trial
  level is narrow INSIDE a target-integration segment subdivides enormously: on tests/13N
  (two active `<targetInt>` entries, evaluations normally sub-second) one screen worker
  burned 69 minutes of solid CPU on a single frozen pass and had to be killed. On 13C+a it
  is why the gate spent 77 min on 4788 patterns. Widening the width grid downwards makes it
  worse: a prior-partition template can give a channel 0.7% of the total, and at the 0.3x
  scale that is a few hundred eV.
  CORRECTION (measured the same day, `9-9-26_seg92_9halfplus/width_cost_cli.py`): the cost
  does NOT diverge, it SATURATES, and calling it pathological was wrong. On segment 92
  (111 points, kernel sigma 0.675 keV) the wall time per mode-1 evaluation against the
  trial width is 14.9 s at 50 keV, 14.9 s at 10 keV, 34.5 s at 2 keV, then flat at
  36-40 s from 500 eV down to 10 eV -- a factor 2.6, not orders of magnitude. So the 13N
  screen's 69 CPU-minutes was ~2160 evaluations at ~2 s each, arithmetic I had not done,
  not a hang and not a divergence. A width floor is still right on INFORMATION grounds (a
  template narrower than the resolution is not identifiable from the data) but it buys
  ~2.6x on the affected patterns, not a rescue. THE DOMINANT COST IS PATTERN COUNT.
  THE REAL LEVER IS THE GRID SETTING, and it is nearly free. The last two optional numbers
  on a `<targetInt>` line are resonanceWidthMultiplier and pointsPerWidth
  (`TargetEffect.cpp:92-102`; defaults 20 and 50 since a4095a6). Measured at the fitted
  width, segment 92 (`grid_setting_cost.py`):
       mult  ppw   wall     chi2      deviation
         20   50  113.0 s   929.33    reference (the converged default)
         10   50   61.8 s   929.35    +0.002%
          5   50   36.8 s   929.38    +0.005%
          5   20   16.9 s   929.63    +0.032%     <- 6.7x cheaper, use for SCREENING
          2   20   10.1 s   927.36    -0.212%
          2   10    6.0 s   923.28    -0.651%     <- 6 units of chi2, larger than the
                                                     Delta chi2 = 1 interval: NOT safe
  So screen at mult 5 / ppw 20 and refit and report at 20/50: 6.7x for 0.03%. Do not go to
  mult 2 -- the shift there is comparable to the width uncertainty itself. Re-measure this
  per project rather than assuming it; the commit established convergence out to a
  multiplier of 200, so 20 is deliberately conservative.
  Also note for the archive: a project that does not set the field now gets 20 instead of
  5, so every target-effect fit became ~3x slower with a4095a6. That is the price of
  correct narrow-resonance integration, but it should be a known price.
- 2026-09-13 -- A FREE ENERGY SHIFT IN THE .azr IS NOT A FREE PARAMETER IN A CAMPAIGN.
  `BoundsPolicy.shifts_free` defaults to False (local.py:39) because shifts are
  finite-differenced, and `default_bounds` then marks every shift column `fixed` (via the
  `Bounds.fixed` mask; lb/ub stay +-inf, so do not judge freedom by the bounds) regardless
  of the file's vary flag, and `polish()` optimizes only the columns not marked fixed. `run_polish` and `run_dataset_off`
  both use the campaign policy, so a "test of freeing the BandH shift" ran with the shift
  clamped at 0.0 for 201 iterations and I reported the shift-fixed result as a negative
  finding (13C+a/9-13-26_bandh_shift, void; see its readme). Corollary: segments 91, 94-96
  of 13C+a carried vary_shift=1 through every campaign and were pinned throughout.
  RULES: (1) to free shifts, put `"shifts_free": true` in the campaign policy at init --
  the CLI's `polish --shifts-free` frees them for that step only, and the dataset-off and
  structure paths ignore it; (2) only segments with a NON-ZERO shift systematic open, so
  set the systematic deliberately and zero it on any segment you want to keep pinned for an
  isolated test; (3) PRE-FLIGHT before submitting: build `default_bounds` on a single
  session with the campaign's policy and print `Bounds.fixed` for every `kind == "shift"`
  key (`9-13-26_bandh_shift_v2/preflight_shift_bounds.py`). Two minutes; it prints
  OPTIMIZED or excluded per shift and would have caught this before a node was spent.
- 2026-09-13 (later) -- the same shift test was void a SECOND time: `run_polish` had
  `shifts_free=False` as its own default and replaced the campaign policy's value with it,
  so `polish` without `--shifts-free` excluded the shift even though campaign.json said
  free. Fixed: `run_polish(shifts_free=None)` defers to the policy and logs
  `polish policy: shifts_free=...`; the CLI flag only ever forces True. Lessons that
  generalize: (1) a pre-flight must build its check through the SAME code path the job
  will run (`Campaign.policy` + the method's own argument handling), not from a hand-made
  object that merely looks equivalent -- mine passed twice while the job failed twice;
  (2) a fitted value of exactly 0.0000000e+00 for a continuous free parameter is never a
  result, it is a parameter that was not optimized; treat it as a failed test until the
  log line naming the policy in force says otherwise; (3) look for the round's
  "dead columns dropped" line -- its absence means the column had a live derivative, so
  an unmoved parameter was excluded upstream of the optimizer, not by the Jacobian.
- 2026-09-13 -- OPEN ISSUE, single-session path only: `SingleSession.residuals(x, jac=True)`
  on the full 13C+a model (34,951 points) with ONE energy shift free ran > 47 min at 100%
  CPU and was killed (`9-13-26_bandh_shift_v2/shift_effect_test.py`); the plan's own
  measurement was 158 s for the whole Jacobian with four free shifts. The SHARDED polish
  with the same shift free ran at its normal 105 s per 10 iterations (job 1440605), so
  per-shard finite differences of a shift are fine and campaigns are not affected. Until
  the single-session case is diagnosed: do not call a single-session Jacobian with free
  shifts (verify() only needs objective/score, which are unaffected); if you must, put
  it in the background with a cap and a kill, as here. The objective itself with the
  shift moved by hand is fast and correct (the +-10/+50 keV sensitivity table in
  9-13-26_bandh_shift_v3/readme came from that path).
- 2026-09-13 -- A POSITIONAL .sav IS INVALID AFTER A LEVEL ADD/REMOVE. AZURE2's param.sav
  names level parameters energy_N / width_N_c by position, so removing (or adding) a
  <levels> line renumbers every later level and `init --sav` puts the fitted values on the
  wrong levels: on 13C+a a seed with three removals applied gave a baseline of 7,199,392
  instead of ~378,020, with "theta2>2, E out of window x25" as the tell. The exported .azr's
  <levels> block already carries the fitted R-matrix values (round-trip 1e-9), so after a
  structure change seed ONLY the norms/shifts, from a .sav reduced to its segment_* lines
  (`grep '^\s*segment_' export.sav > norms.sav`). Also: csh `eval` of a step string drops
  everything after "#" in a label like 9/2+#2 and splits quoted --note values -- write job
  steps as explicit quoted lines and dry-run them into Python's argv before qsub.
  GUARD (same day): `Campaign._seed_vector` now compares the .sav's highest energy_N index
  with the model's; on a mismatch it logs "WARNING ... a positional .sav is invalid after a
  level add/remove" and seeds norms/shifts only. Verified end to end: the positional .sav
  fed to the 92-level model warned ("indexes 95 levels but the model has 92"), seeded 95
  segment parameters, and reproduced the correct baseline 378,052.431 exactly, where the
  unguarded seeding had given 7,199,392.
- 2026-09-13 -- EVERY STRUCTURE MOVE IS JUDGED AGAINST A SAME-BUDGET CONTROL POLISH, now
  including `run_relocate` (it had none: a relocation whose level walked straight back to
  its start still showed a "gain" of 309, which was 100 evaluations of warm-start drift).
  `run_remove`, `run_level_add`, the sign round and `run_relocate` all now polish the
  unmodified incumbent with the same budget first and judge the move against that. Rule:
  a move's gain is `control - moved`, never `incumbent - moved`; if you see the latter in
  any new move, it is wrong. The 13C+a 9.864 MeV 9/2+ relocation is the worked example
  (`9-13-26_structure2/readme`).
- 2026-09-13 -- THE REDESIGNED LEVEL SEARCH IS WIRED INTO `rmfit levelscan` with two new
  flags. `--differential-window KEV` (default 250): after the GN step the candidates are
  re-ranked on the DIFFERENTIAL segments within that half-width of the candidate energy
  (`levelscan.scoring_segments`; Ex ranges from one cached single session), because the
  total objective barely distinguishes J^pi; 0 restores ranking on the total. `--grid M,P`
  (default 5,20): the screen VARIANT's <targetInt> lines get resonanceWidthMultiplier M and
  pointsPerWidth P (`Campaign.set_grid_setting`), measured 6.7x cheaper than the converged
  20/50 for 0.03% -- the base model is untouched, so verify/export/report stay converged;
  'none' keeps the model's own setting. Also new defaults: `--scales` 0.3,0.6,1.2,2.5,5,10
  (x2 steps, reach matters more than density) and 2 GN slots reserved per J^pi group.
  Prior-partition templates come from the group's own existing levels automatically.
- 2026-09-13 -- REDESIGNED SCREEN, FIRST END-TO-END RESULT (13N, the removed 3/2- 3.5032 MeV
  level, `rmfit levelscan --ex 3.5032 --grid 5,20 --differential-window 250`): grid 5/20
  applied to both <targetInt> lines of the screen variant; 1008 patterns frozen-evaluated
  in 193 s (the same model took a 13N worker 69 CPU-min two days ago at the converged
  grid with the old pattern set); differential re-rank ran on 5 segments; the removed 3/2-
  ranked FIRST, and it was found through the PRIOR-PARTITION family (the partition of the
  group's own existing 3/2- at 20 MeV), gain +1,784,301 global / +353,977 local. The
  relocate regression also passes with its new control polish. So both mechanisms the
  13C+a gate failure asked for -- a grid that reaches the level, and a ranking that sees
  J^pi -- work on the fixture. First real use: 13C+a at 9.862 MeV (job 1441224).
- 2026-09-14 -- REDESIGNED SCREEN, FIRST REAL USE (13C+a at Ex 9.862 MeV, job 1441224,
  `13C+a/9-13-26_levelscan_986/`): NULL. 16,668 patterns frozen-evaluated in 4.0 h at grid
  5/20 (0.87 s per pattern; the cost is now the pattern COUNT, not the evaluation), 2,416
  below the incumbent; best after the GN step +30.8 on the global objective (1/2+ prior
  partition) and +3.0 on the 46 differential segments within +-250 keV (7/2+ single channel
  @9.902) -- both noise against tau ~378. The 7/2+ refit that followed sat at 377,782.9 after
  30 released iterations (incumbent 377,786.07) and threw Brune "Denominator less than
  zero" warnings for the added level throughout: a level the data do not want.  Read it
  as a limit of the FROZEN screen, not as "no level here": BandH's 9.857 MeV peak is real
  (seen by all ND 2021 angles and 28 Heil angles) and the model's 9/2+ 9.864 sits 5 keV
  away, so any new level here is the non-perturbative two-level case that the frozen
  single-level template cannot represent (the 9/2+ cannot readjust inside the screen).
  The move that CAN see it is `structure add` with a pinned energy and the released polish
  (9-14-26_add_52m_986, ENSDF 5/2- 9861.74 keV), judged against its control polish.
  Cost lesson for the ceiling of <= 30 min per energy: 16,668 patterns is 11 J^pi x
  6 widths x (single-channel + equal + dominant + prior families) x signs x 3 offsets;
  the offsets and the full sign enumeration of the equal split are where to cut.
- 2026-09-14 -- THE TEMPLATE SCREEN IS TWO-STAGE (`template_screen(stage_scales=(1.2, 5.0),
  refine_top=3)`, CLI `levelscan --stage-scales 1.2,5 --refine-top 3`; `--stage-scales none`
  restores the exhaustive stage).  Stage 1 evaluates every partition x sign SHAPE at the two
  stage widths and the candidate energy only; stage 2 takes the best `refine_top` shapes per
  J^pi to every width in `scales` and every energy offset.  The shape ranking is what carries
  the information; width and a 40 keV offset only scale it.  Counted on the 13C+a 9.862 MeV
  variant (11 J^pi, 12-24 channels each, 1-14 priors): 19,836 -> 2,852 evaluations, ~41 min
  at 0.87 s instead of 4 h; `--max-signs 8` takes it to ~34 min (the equal-split family's
  32 random sign patterns are the least informative shapes for n >= 12).  13N regression
  (`tests/test_levelscan_13n.py`, now with the 5/20 screening grid applied to its variant
  as the campaign does -- at the converged grid the 13N variant is ~7 s/evaluation, at 5/20
  ~0.2 s): the removed 3/2- still ranks first, gain +1,771,289 vs +1,784,301 exhaustive;
  30-70 evaluations per group instead of 180; whole test ~4 min on the login node.
  The per-J^pi line the screen now logs ("<jpi>: n channels, S shapes; N1 stage-1 + N2
  refinement evaluations in T s; best frozen ...") is the cost accounting the 4 h run lacked.
- 2026-09-14 -- FIRST LEVEL ACCEPTED BY `structure add` ON 13C+a AFTER THE FROZEN SCREEN SAID
  NULL at the same energy (job 1441491, `13C+a/9-14-26_add_52m_986/`): candidate energy
  9.8617 MeV pinned, `jpi: null` so all 12 J^pi were scanned on the +-0.5 MeV window model
  with a local polish each (1/2+ -2,487, 9/2+ -1,104, ..., the compilation's 5/2- only -27);
  the 1/2+ went to the full model: seeded 376,078 -> pinned 375,975 -> released 375,519
  (150 evaluations) vs a same-budget control 377,786 -> 377,786: delta +2,266 against tau
  756, rho 0.42, theta^2 1.9e-4 -> accept.  Fitted: Gamma_n(L=0) 2.92 keV, Gamma_alpha(L=1)
  39 eV.  The gain sits in Cierjacks n-total (-1,239) and Heil (a,a) (-1,176), not in the
  BandH 9.857/9.867 pair that motivated the energy (unchanged, -76/+34/+35 sigma).  Two
  lessons.  (1) The window scan's LOCAL POLISH sees what the frozen screen plus one damped
  GN step cannot: the same energy screened null (best +30.8) four hours earlier.  A frozen
  template is a lower bound on a level's worth, never an upper bound; a null screen is not
  a reason to skip the add when independent evidence (a peak the data resolve) points at
  the energy.  (2) `jpi: null` with the pinned energy is the right way to let the angular
  data assign J^pi: it disagreed with the compilation's tentative (5/2-) and won by a
  factor 90 on the window.  Record such disagreements as evaluation decisions for the human.
  Current best: 375,519.344, 93 levels.  Post-add sign round over 1/2+ = job 1441540.
- 2026-09-14 -- THE TEMPLATE SCREEN'S AMPLITUDES WERE WRONG BY ~100x ON 13C+a; BOTH "NULL"
  SCREENS AND THE FIRST GATE FAILURE WERE THIS BUG, NOT THE GRID.  Diagnosis chain, each step
  a measurement: (1) the accepted 1/2+ 9.862 dropped into the incumbent at its fitted widths
  with nothing else moved is worth +1,723 frozen (of the +2,266 released), so a frozen screen
  should have seen it; (2) the screen's own variant with the dummy 1/2+ at those widths through
  the screen's amplitude conversion gains +0.4, and flipping the sign changes nothing -- the
  amplitude column was inert; (3) the engine's x0 for the dummy's 1 eV reference width was
  5.5e-6 where the group's existing levels imply ~0.013 per sqrt(eV): the .azr -> rwa step
  goes through a Brune transform that fails for a 1 eV level among keV-wide neighbours
  ("**WARNING: Denominator less than zero while transforming", exactly one per dummy: 12 in
  the 9.862 run, 338 lines in its log, 576 in the first gate's); (4) the dummy set directly to
  the fitted rwa (0.0255 / 0.00364) gains +1,722.7 in the variant -- identical to (1).
  FIX (`levelscan.py`, backup `.bak-2026-09-14`): gamma = gamma_W sqrt(Gamma_c / Gamma_W) with
  `ev.wigner_amplitudes()` and `ev.wigner_widths(x_inc)` (eps-probe, no transform of the
  dummy); closed channels (no Gamma_W: an ANC) are excluded from the partitions instead of
  being given "widths"; the screen logs per group the x0-based/Wigner-based ratio and flags
  it outside 0.3..3.  The 13N fixture never showed this (its dummies transform fine, ratio
  ~1), which is why every test passed while the real screens were blind.  RULE: whenever a
  screen returns null at an energy the data visibly resolve, drop the hypothesised level in
  at plausible widths through the SAME code path and check the frozen gain is nonzero before
  trusting the null -- and grep the log for "Denominator less than zero" first.
  Consequences: the 2026-09-12 "width grid reach" finding stands as geometry but was NOT the
  binding defect of the first gate run; the 9.862 null (16,668 patterns, 4 h) is void; the
  gate (a) re-run and the 9.862 screen must be repeated with the fix.
- 2026-09-14 -- FIX VALIDATED ON THE REAL CASE (login node, 6 shards, grid 5/20, 1/2+ only at
  9.862 MeV on the 92-level incumbent): 1/2+ now has 3 OPEN particle channels (the 9 closed
  ones that used to be given "widths" are left out, so 42 shapes instead of 84); best frozen
  template = the prior partition of the 1/2+ @6.372 at 3 keV, +1,621; after the GN step
  +1,902 (the released add gave +2,266); 132 frozen evaluations in 285 s, GN of 6 in 441 s.
  Note the differential re-rank put that same candidate at local +5 / -217: the angular
  segments within +-250 keV do NOT carry this level (its gain is Cierjacks n-total and the
  Heil angles as a whole), so the LOCAL statistic is for deciding J^pi among candidates that
  the GLOBAL GN gain already flags, never for deciding whether a level exists.  The GN gain
  predicted the released gain to 16%.
  Tests after the fix: `test_levelscan_13n.py` (its exact-parameter check now converts through
  the screen's Wigner route and asserts the dummy's open channel has a Gamma_W; logs the
  x0/Wigner ratio, which is 1.0 on 13N for every group but 11/2+ at 0.12 and 9/2- at 1.6) and
  `test_levelscan_grid_13n.py` both pass; the 3/2- gain is unchanged (+1,771,289).
- 2026-09-14 -- GATE C' (`13C+a/9-14-26_levelscan_986_fixed`, job 1441597): the repaired
  screen at 9.862 MeV on the 92-level model ranks the accepted 1/2+ FIRST on the global GN
  gain (+1,902; released +2,266) -- 1,504 patterns, 26 min all-in (budget 30).  The LOCAL
  differential re-rank does not put it in the top 5, and instead ranks a 9/2+ 2 keV from the
  existing 9/2+ 9.864 first (+1,693 local / +1,814 global).  Two rules from this.  (1) The
  GLOBAL GN gain is the existence statistic; the LOCAL one discriminates J^pi among
  candidates the global gain already flags (a narrow n s-wave level is paid for by n-total
  and by every angle of (a,a), not by the +-250 keV window).  `run_levelscan --refit N`
  refits by the local order today; make it refit the top-N by global gain and report the
  local rank alongside.  (2) A candidate ON TOP of an existing level of the same J^pi is not
  noise: it says the existing level's strength or partition is wrong there and is tested as
  a pinned add with the released polish (the relocate move is the wrong tool -- it keeps
  the partition).  Wigner ratios on 13C+a: 0.012 (1/2+) down to 3.7e-6 (11/2+).
  Done the same morning: `Campaign.run_levelscan` records `best_global` (by GN gain, with the
  local gain alongside) next to the local-order `best`, notes "best per J^pi by GLOBAL gain
  ... | local order ...", and `--refit N` takes the top N by global gain.  13N smoke test
  through the CLI (`levelscan --ex 3.5032 --grid 5,20 --stage-scales 1 --refine-top 2`):
  304 patterns in 54 s, 3/2- first on both statistics (+1,784,313 global / +354,070 local).
- 2026-09-14 -- THE SCREEN'S LOCAL TOP CANDIDATE, A 9/2+ 2 keV FROM AN EXISTING 9/2+, WAS
  ACCEPTED BY THE CRITERIA (job 1441606, `13C+a/9-14-26_add_92p_986`): delta +3,898 vs tau
  751, rho 0.31, theta^2 1.6e-3; 375,519 -> 371,621 (94 levels).  The fitted level is 330 eV
  wide, 1.2 keV below the 3.8 keV 9/2+; BandH's +34 sigma point at 9.8613 went to -0.5.
  LESSON: the acceptance criteria (delta, tau, rho, theta^2, holdout) test whether a level
  HELPS, not whether it is a level.  A candidate that lands within a width of an existing
  level of the same J^pi is a LINESHAPE hypothesis, and the method must test the
  alternatives before adopting it: energy resolution / target convolution on the segments
  that drive the gain (BandH has none here; only segments 52 and 92 carry <targetInt>), a
  single point, or the existing level's own partition.  To build: `structure add` should
  flag "within 2 Gamma of a same-J^pi level" in its verdict and the campaign should stage
  a resolution test (`<targetInt>` beam spread grid, doublet removed, control polish) as
  the competing move.  Until then such acceptances are recorded as PROVISIONAL in the
  readme and left to the evaluator.
  Built the same morning: `run_level_add(near_same_jpi_kev=20)` checks the base model for a
  same-J^pi level within 20 keV of the candidate, logs "PROVISIONAL: ... lineshape
  hypothesis" in the verdict and the readme note, and records `near_same_jpi` in the
  ledger round (13N smoke test through `structure add`: accept, no flag, exit 0).  The
  competing resolution test is still to build.
- 2026-09-14 -- FIRST LEVEL FOUND FROM RESIDUALS ALONE (`13C+a/9-14-26_levelscan_resid`, job
  1441631, `levelscan --top 3 --refit 1`, 3.6 h): the three residual-run candidates
  (10.366, 10.452, 10.413 MeV) all screened to ~10.41 MeV; best global GN gain 7/2+ +3,706;
  refit accepted at delta +5,053 vs tau 1,669, rho 0.28: 371,621 -> 366,568 (95 levels).
  It is an n2 p-wave level (n + 16O*(6.13 MeV 3-), 12.6 + 5.6 keV) 140 keV above the n2
  threshold; the gain is the 13C(a,n2 gamma) data (48.2 -> 17.8 chi2 per point).  Per-energy
  cost at 10.4 MeV: 2,488 patterns / 37 min + GN 7 min (more open channels than at 9.86).
  Two things to fix from this run: (1) the accepted level carries formal theta^2 ~0.6 in
  channels with ~0 physical width (high-L / barely open) -- the theta^2 sanity uses formal
  amplitudes, so a level can pass with meaningless amplitudes in dead channels; zero and
  fix channels whose Wigner width is below ~1 eV at the level energy before the released
  polish; (2) `run_levelscan` screens three residual energies that turn out to be one
  feature -- cluster candidates whose +-40 keV screens converge on the same energy before
  spending 40 min on each.
  Built the same afternoon (campaign.py, backup `.bak-2026-09-14b`): (1) `run_level_add(
  dead_channel_ev=1.0)` zeroes and FIXES the new level's particle channels whose Gamma_W
  at the level energy is below 1 eV (closed or barely open; photon channels exempt), logs
  them, and leaves them out of tau's k -- the 7/2+ 10.413 sign round otherwise spent its
  polishes flipping four such channels at an unchanged objective; (2) `run_levelscan(
  merge_kev=60)` / `levelscan --merge-kev` drops residual candidates within 60 keV of a
  stronger one before screening (13N: 3.572 dropped next to 3.541; on 13C+a it would have
  saved the third 40-min screen at 10.413).  13N smoke through the CLI: levelscan with
  clustering exit 0; structure add ran through the window scan and pinned polish before
  the login-node time cap.
- 2026-09-14 (DeBoer) -- BROAD STRUCTURE FIRST.  "Don't get stuck investigating the
  inconsistencies between alignments of very narrow resonances ... focus on the broader
  resonances, widths greater than approximately 10 keV, to get a fit to the wider energy
  range first."  In practice: `levelscan --gamma 20 --scales 0.5,1,2,4` (10-80 keV trial
  widths, so sub-10-keV levels are simply not in the template grid), candidates ranked by
  broad residual runs, and no lineshape/resolution/doublet tests of narrow levels until the
  end stage -- record them PROVISIONAL and move on.  The 9.86 MeV morning (a 3 keV 1/2+, a
  330 eV 9/2+ partner, per-point BandH residuals) is the example of what NOT to spend a
  node-day on at this stage.
- 2026-09-14 (evening) -- STATE AND HANDOVER.  Campaign best: 366,565, 95 levels
  (`13C+a/9-14-26_signround_72p`).  Today's chain: 377,786 -> 375,519 (1/2+ 9.862) ->
  371,621 (9/2+ 9.862, PROVISIONAL) -> 366,568 (7/2+ 10.413, n2 level) -> 366,565 (7/2+
  sign round, empty).  The second residual search found nothing above tau at 10.37-10.49
  after the GN step and was terminated from outside during its GN step; a parallel session
  (DeBoer's) now runs the BROAD-structure search on the same seed
  (`9-14-26_levelscan_broad`: --gamma 20 --scales 0.5,1,2,4) under his instruction to fit
  Gamma >~ 10 keV resonances over the whole range before narrow alignment -- follow that
  ordering in future campaigns: broad screen first, narrow candidates near same-J^pi
  levels provisional.  Two nodes were briefly in use today (the other session's high-Ex
  window job and mine): check `qstat` for OTHER rmfit17O jobs before every submission,
  not just one's own.
- 2026-09-14 (evening) -- READ THE DRIVER'S campaign.log, NOT THE SGE STDOUT FILE.  The
  `-o rmfit_*.log` file is the job's stdout on NFS and can lag the driver by an hour (job
  1443295: stdout mtime 17:47 while campaign.log showed the refit verdict at 18:44).  I read
  the stale stdout, saw "nothing for 53 min", sampled CPU (one core busy -- the control
  polish's Jacobian, as it turned out), reproduced a suspect call standalone (fine), and
  qdel'd the job a few minutes after it had already recorded its verdict, during the
  report/export of an unchanged incumbent.  Nothing was lost, but the rule for every watch
  and every "is it hung?" question is: `campaign.log` (written directly by the driver) and
  the ledger are the truth; the SGE log is a convenience.  Monitors should tail
  campaign.log.
- 2026-09-14 (night) -- THE BROAD-ONLY SEARCH PAID OFF AT ONCE (`13C+a/9-14-26_levelscan_broad2`,
  job 1443649, `levelscan --ex 10.846 10.753 9.160 --gamma 20 --scales 0.5,1,2,4 --refit 1`):
  a 5/2- at 10.886 MeV, Gamma ~360 keV (n2 s-wave 220 keV, n0 f-wave 96 keV), accepted at
  delta +13,859 vs tau 1,097, rho 0.18: 366,565 -> 352,703 (96 levels).  Screen +3,999 after
  GN predicted a released +13,859 -- for a BROAD level the GN step underestimates by 3.5x
  (the whole neighbourhood re-shapes: the 5/2- 10.940 next to it changed by 100s of keV),
  so for broad candidates a screen gain of ~+2,000 already deserves the refit.  The gain is
  spread over every broad dataset (MANA (a,a) -6,231, ND 2021 -2,486, Heil -2,165, Cierjacks
  -1,754): that is what "broad structure first" buys.  Cost: 3 energies x ~40 min + refit
  ~70 min = 3 h per accepted level.  The dead-channel rule fired for the first time here
  (5 channels with Gamma_W < 1 eV zeroed and fixed).  Also the first time DeBoer's
  restriction changed the ranking: at 10.453 the narrow-capable screen had preferred a
  5/2- n2-only ~20 keV candidate (refit marginal); the broad grid at 10.846 found the
  360 keV one.
- 2026-09-14 (night) -- THE BROAD-ONLY SEARCH FOUND THE CAMPAIGN'S LARGEST GAIN
  (`13C+a/9-14-26_levelscan_broad2`, job 1443649, `levelscan --ex 10.846 10.753 9.160
  --gamma 20 --scales 0.5,1,2,4 --refit 1`): a 5/2- at 10.886 MeV, ~360 keV wide (n0
  f-wave 96 keV, n2 s-wave 220 keV, n2 d-wave 43 keV), screen +3,999 after GN, refit
  delta +13,859 vs tau 1,097, rho 0.18: 366,565 -> 352,703 (96 levels).  Paid for by MANA
  (a,a) -6,231, ND 2021 (a,n0) angles -2,486, Heil -2,165, Cierjacks n-total -1,754.  The
  earlier narrow-capable screens (Gamma 1.5-50 keV) at 10.366/10.452 never reached this
  energy; the 10-80 keV grid at the next residual energy found it in one screen.  Notes:
  (1) the dead-channel rule fired on first use (5 channels zeroed) and the level came out
  clean; (2) the existing broad 5/2- 54 keV away re-shaped with it (n2 s-wave 0 -> 152
  keV) -- two broad same-J^pi poles sharing strength is what a broad feature looks like in
  R-matrix, not a doublet to be "resolved"; theta^2 0.95 on the n2 s-wave widths says the
  strength is at the Wigner limit, an evaluator's question (background pole?) for later;
  (3) 9.160 MeV (15.6k chi2 in runs) was a clean null at broad widths: not every large
  residual run is a broad level.  Per-energy cost 34 + 6 min at 10.8 MeV.
- 2026-09-14 (late) -- SIGN ROUNDS NO LONGER PROPOSE FLIPS OF ~ZERO AMPLITUDES (search.py,
  backup `.bak-2026-09-14`): `sign_round` feeds the proposal generators a vector in which
  |rwa| < sign_eps * gamma_W (open channels a released polish left at ~0 -- typical of a
  just-added level; not caught by the dead-channel rule, which only fixes Gamma_W < 1 eV) is
  zeroed, and logs the excluded keys.  Flipping ~0 is a no-op that "improves" by noise and
  then costs a staged polish: the 7/2+ and 5/2- post-add rounds each spent ~3 h on such
  flips (0 accepted).  Screening/polishing still start from the true incumbent.  13N
  `tests/test_search_13n.py`: all checks ok, round unchanged (70,972.321 -> same).
  Diagnosed jointly with the parallel session, which also proposed fixing those columns
  after the released polish -- not done, because a later round may legitimately want them.
- 2026-09-15 -- A DETACHED (nohup) SUBMITTER CANNOT qsub ON THIS CLUSTER: "job rejected: job
  does not provide an AFS token".  The 9-15-26_levelscan_broad3 job staged at 01:41 was
  never submitted and nobody noticed for 15 h because the script logged "submitted" after
  the failed qsub.  Submit from a live session's shell or Monitor (both carry the token),
  check the qsub output for "Your job N", and never print "submitted" unconditionally.
- 2026-09-15 -- SECOND BROAD LEVEL FROM THE SAME SEARCH FAMILY (`13C+a/9-15-26_levelscan_broad4`,
  job 1447357, explicit energies 10.543/10.909/10.172/9.774 at Gamma 10-80 keV): a 9/2+ at
  10.909 MeV accepted, delta +1,976 vs tau 842, rho 0.34: 351,676 -> 349,630 (97 levels);
  gain again from MANA (a,a) -1,086 and Heil -398.  Three observations worth carrying:
  (1) the WINDOW scan read -19,832 while the full-model gain was +1,976 -- a truncated
  window model exaggerates a level's worth by an order of magnitude; never quote the window
  delta as the level's value, only use it to rank J^pi.  (2) tau fell from ~1,700 to 842
  because the dead-channel rule zeroed 6 of the level's channels: the threshold correctly
  follows the level's effective parameter count.  (3) theta^2 = 0.98 on this level and 0.95
  on the 5/2- 10.886/10.940 pair -- three broad levels at the Wigner limit in one 200 keV
  band says the model is short of poles there (background-pole question for the evaluator),
  not that three resonances are well determined.  Also: 10.3-10.5 MeV is now exhausted at
  broad widths (broad3: 7/2+ @10.449 delta +676 vs tau 737, marginal), and 9.774 and 10.172
  are clean nulls.
- 2026-09-16 -- A SIGN FLIP IN A CHANNEL THE DATA DO NOT SEE IS A NO-OP, AND SIZE DOES NOT
  TELL YOU WHICH ONES THOSE ARE.  Three post-add sign rounds in a row (7/2+ 10.413, 5/2-
  10.886, 9/2+ 10.909) spent 2-3.5 h polishing flips that left the objective unchanged.
  The 2026-09-14 tiny-amplitude exclusion (|rwa| < sign_eps*gamma_W) caught some but not
  these: measured on the 9/2+ round, `g[9/2+#6 p4 L4 S0.5]` has theta^2 = 0.98 -- at the
  Wigner limit -- and flipping it moves the objective by +0.24.  Reason: on 13C+a only
  pairs 1, 2, 3 and 5 appear in any active `<segmentsData>` line; pairs 4, 6, 7, 8, 9 are
  carried by every level but observed by nothing, so their signs reach the observables only
  through second-order multiple-scattering terms.  FIX (`search.py`, `SearchConfig.
  flip_floor = 0.01`): a proposal must move the FROZEN objective by at least
  `flip_floor * tau` (~5 on 13C+a) to earn a staged polish; the screen already computes
  that number, so the filter is free, and it logs what it dropped.  Regression test
  `tests/test_flip_floor.py` (no engine, milliseconds) pins the behaviour on the six
  measured no-op responses.  General rule: judge a proposal by the response it produces,
  not by the size of the parameter it changes.
- 2026-09-16 -- WHEN SINGLE-LEVEL ADDS STOP PAYING, MAP THE CHI^2 AND AUDIT THE CAPS BEFORE
  SCREENING ANOTHER ENERGY.  Three refits at 10.449 MeV on 13C+a came back marginal (+676,
  +617 against tau ~735) and two more energies screened null, while the residual-run ranking
  kept nominating the same band -- because those runs are large in POINT COUNT, not because
  one level is missing.  Two cheap diagnostics (one single session each, ~5 min) said what
  the searches could not.  (a) chi^2 in 200 keV Ex bands with the dominant datasets per band:
  the 13C+a misfit is SPREAD at 12-20 per point over Ex 9.0-11.2 (2,000-point bands, tens of
  thousands of chi^2) -- a single level returning +400 cannot matter there -- plus one sharp
  outlier, Ex 5.6-5.8 = 7,529 chi^2 over 60 points (125/pt), ALL in one n-total set and
  ~5,500 of it in six points (+74 sigma at 5.7321): a narrow n-channel feature, correctly
  deferred under "broad first".  (b) A theta^2 audit: 5 physical levels at Ex 10.54-10.94
  sit AT the cap (theta^2 0.83-1.00) and 5 background poles at Ex 12-27 sit AT theirs
  (1.7-3.0 of cap 3) -- the R-matrix has no room left in those J^pi groups, which is why
  capped adds buy nothing.  The physical reading is that the region is short of POLES, the
  archive's 16O+p lesson (an added l=0 background pole was what fixed that fit).  Running
  `13C+a/9-16-26_cap_diagnostic` (job 1447768): the same fit re-polished with the caps
  raised (physical 3, background 10) to measure what the caps cost -- a diagnostic, never
  adopted.  RULE: a band map plus a cap audit costs ten minutes and should precede any
  further level search once two consecutive searches come back marginal or null.
- 2026-09-16 -- THE MODEL'S BACKGROUND POLES WERE INERT, AND THAT IS WHY THE FIT WAS PINNED
  AT THE WIGNER LIMIT.  Every J^pi group of 13C+a carries a level at Ex = 50 MeV whose
  <levels> lines have levelFix = 1 and gamma = 0 in EVERY channel -- fixed and zero, so the
  engine creates no free parameter for them (`ev.level_energy` lists the level, but no key
  of any kind belongs to it: check with `[k for k in ev.keys if k.level == lab]`, not with
  the file's channelFix field, which is 0 and looks free).  Ten poles doing nothing, in a
  campaign that spent days looking for missing strength between 10.5 and 11 MeV.  Two
  separate traps here: (1) a level whose widths are all zero is invisible to a gradient
  polish even when free -- the cross section is quadratic in the amplitude, so the
  derivative vanishes at gamma = 0 and the polish leaves it at zero forever (the same
  reason the template screen needs a finite trial width); (2) `levelFix` (field 4) fixes
  the whole level, while `channelFix` (field 11) is per channel -- reading only the latter
  says "free" when the level is frozen.  DIAGNOSTIC that finds this in one session: count
  the free keys per level and list the levels with none.  The fix is a file edit plus a
  seed: `set_fixed(False)` on the pole and a physical width (200 keV used here) in the
  channels the data observe, then a normal polish + sign round (`13C+a/9-16-26_bgpoles`,
  job 1447775).  Generally: before concluding "the region needs another level", check that
  the poles already in the model are actually being fitted.
  Follow-up the same morning: seeding those poles took two attempts.  A `.sav` is POSITIONAL
  (`width_N_c`) and N is the ENGINE's level index, NOT the level's position in the file --
  the 3/2+ pole is level 47 in file order and `width_48_*` to the engine.  The first job
  filtered the .sav on the file index, so the real entries survived and zeroed every seed,
  and the init came back at exactly the incumbent's objective (349,619.415).  Nothing
  crashed; only the identical number gave it away.  RULES: (1) when seeding a parameter the
  .sav also carries, filter the .sav by the names a SESSION reports (`snapshot`'s
  `parameter.name` for the columns you seeded), never by an index you computed from the
  file; (2) after staging any seed, evaluate it once and check the objective MOVED before
  spending a node on it -- `13C+a/9-16-26_bgpoles/filter_sav.py` and the check in that
  readme are the pattern (seed verified at 349,549.427 vs 349,619.415 before resubmitting
  as job 1447790).
- 2026-09-16 -- level tests accepted broad absorbers on the 13C+a high-Ex window
  (`13C+a/9-14-26_highE_structure`): a 9/2- 12.738 (+222k) and a 7/2- 11.818 (+320k)
  with three or four channels pinned EXACTLY at the theta^2 cap (1.0, then 0.5) and the
  rest collapsed to ~1e-35; their Brune transformation fails ("Denominator less than
  zero") and `export`/`bake` then dies with "different free-parameter set". `classify`'s
  blowup test only fires above the cap; tightening the cap does not help (the polish
  drives the amplitudes to whatever cap exists). `classify` now takes `at_bound` (the
  new level's width keys the released polish left on a bound, from `FitResult.bounds_hit`)
  and returns "unphysical" -- not adopted -- for any at-bound width or theta^2 within
  0.1% of the cap. Test: `rmfit/tests/test_classify_gate.py`. Structures 4-7 of that
  campaign are the evidence; the search was rerun from structure 3.
- 2026-09-16 -- TWO SAME-J^pi LEVELS WITHIN ~1 keV MAKE A FIT UNWRITABLE.  The 13C+a fit at
  336,944.52 verified fine in-session but `bake` refused it: "the baked model has a different
  free-parameter set than its source".  The polish had drifted a free 9/2+ level to 0.66 keV
  from a pinned 9/2+ (the PROVISIONAL doublet partner added two days earlier), and the engine
  reading the WRITTEN file merges levels that close -- the group loses a level, its labels
  shift (9/2+#6 -> #5) and eight keys vanish.  In-session evaluation never sees this because
  the model is never re-read.  Three lessons.  (1) A same-J^pi pair at the same energy is a
  REDUNDANCY, not a doublet: measured here, separating them costs +1,700 at 1.1 keV, +5,400
  at 2 keV, +80,000 at 5 keV -- the fit had tuned them as one structure, so the cure is to
  remove one, not to pull them apart.  (2) `bake`'s key-set guard must delete its half-written
  files (patched): they were left on disk and evaluated at 403,867 while looking exactly like
  a successful export -- I quoted that file as the result before checking.  (3) After any
  campaign that adds a level near an existing one of the same J^pi, check the minimum
  same-J^pi energy gap before trusting an export; `structure remove` on the incumbent is the
  fix and needs no export to run.
  Corollary measured the same afternoon: the FROZEN removal screen cannot see a degenerate
  pair.  On the 13C+a fit the two 9/2+ levels 0.66 keV apart ranked NOWHERE in the ten
  cheapest removals (the cheapest were 1/2-#1 at +0.0, the freshly activated 9/2+ background
  pole at +836, 7/2-#4 at +1,153), because zeroing one level's widths without letting its
  partner absorb the strength is catastrophic by construction.  Frozen screens rank levels
  by what the model loses when a level is switched OFF; a redundancy only shows up when the
  partner is refit.  So for a suspected same-J^pi degeneracy, test the specific level with a
  targeted remove-and-refit, not with the ranking.

## 2026-09-17 — audit the pair table before fitting (13C+α)

Pair 8 of the 13C+α model (α + ¹³C* 3.854 MeV) had `z2 = 8` on every level line since
2016 — an oxygen Coulomb barrier on a carbon channel — through every fit in the archive.
It cost nothing in the fitted window only because the channel was also inert (all
amplitudes exactly zero, so never moved by any polish); above Ex 12 MeV it is a factor
3–12 in penetrability.  Nothing in AZURE2, pyazr or rmfit checks a pair's charges and
masses against each other.  Rule: at campaign init, print one row per pair —
`ir m1 m2 z1 z2 sepE e3 j3 pi3 radius` — from the `<levels>` block and compare with the
other pairs of the same particle (same target nucleus ⇒ same z2, same m2, same sepE
minus e3) before anything is fitted; also list channels with γ = 0 and channelFix = 0
(nominally free, actually dead — see the background-pole entry of 09-16).  A `rmfit
audit` command that does both is the right home; until it exists, the awk one-liner is
in `13C+a/8-10-26_more_data_claude/readme` (2026-09-17 19:40 entry).

## 2026-09-18 — seed dead channels in units of γ_W, never in keV (13C+α)

Seeding zero amplitudes with a physical width (`gamma = 200e3` eV in the `.azr`) is only
safe in well-penetrable channels.  Applied to the n₁/α₁ pairs of thirty background levels
it put the model at +166,747 (and +136,525 at 20 keV): for a high-L channel with
P ~ 10⁻⁶, 20 keV is a reduced width of thousands of MeV.  Do it in two steps instead:
a 1 eV physical seed so the channels exist as engine keys (zero amplitudes are not keys),
then set each key to `frac · γ_W` (θ² = frac², 0.04 is plenty) in the engine and `bake`
the vector — the baked `.azr`/`.sav` pair is consistent and round-trip checked, so no
`.sav` filtering is needed.  Always evaluate the seed before submitting: a seed should move
the objective by less than the polish's first step, not by 50%.  Script:
`13C+a/9-18-26_nbg/seed_nbg2.py`.

## 2026-09-18 — leave-one-dataset-out: totals are not enough (13C+α)

The 14-row table (`13C+a/9-18-26_loo`) separated cleanly into distorting sets (extra
recovery 2,400–37,000) and consistent ones (190–1,000 = the drift of a 100-evaluation
refit with no control polish — polish the reference first and quote that floor).  But
EVERY distorting set sat in the same energy region, so the table said "region", not
"dataset", and could not say which measurements disagree with which.  `run_dataset_off`
now also logs and stores (`by_dataset` in the round stats) the per-dataset χ² of the rest,
frozen vs refitted — the pairwise tension map.  Two reading rules: extra/own > 1 (a set
that returns more than its own χ² when removed) means it pulls against a specific partner,
look for the partner in the breakdown; and rank by extra recovery PER POINT, not total —
the leverage is in small precise sets (here 145 points returned 35 per point against 2–4
for the 10,000-point sets).  Test: `rmfit/tests/test_dataset_off_13n.py` (slow, ~15 min:
13N recomputes its capture integrals for the reduced variant).

## 2026-09-18 — a per-dataset table must reproduce a number you already know

`diagnose.dataset_table(ev, score)` double-counted any segment that a `Sharded` evaluator
splits into pieces (one row range per piece, χ² per base segment): Heil on 13C+α came out
at exactly 2.000× its known χ².  Caught only because the first use of a new breakdown was
compared, row by row, with the incumbent's REPORT table before being interpreted.  Fixed;
regression test `tests/test_dataset_table_split.py` uses a stub evaluator (runs in a
second) and was run against the OLD code first to prove it fails there.  Scope: the
"worst datasets" line at `init` and `compare()` on a sharded evaluator; accept decisions
were unaffected (their ρ is per segment from verified single-session scores).  Rule: the
first output of any new diagnostic gets checked against an independently known total
before anyone reads physics into it.

## 2026-09-18 — after a session restart, qsub may lose its Kerberos cache

"job rejected: job does not provide an AFS token" from an INTERACTIVE shell: the restarted
session inherited a KRB5CCNAME pointing at a cache file that no longer exists (`klist`:
"No credentials cache found") although the AFS token itself was fine.  Fix: point
KRB5CCNAME at the session's own still-valid cache (`ls /tmp/krb5cc_<uid>_*`, check each
with `KRB5CCNAME=FILE:<path> klist`; the one whose start time matches `tokens`' expiry is
this login's) for the qsub call, or ask the user to `kinit`.  Also: long login-node
computations launched with run_in_background die with the session -- launch them with
nohup and a done-file, and watch the done-file.  And `qsub -hold_jid <running job>` queues
the next step behind the current one without needing the session alive at hand-over.

## 2026-09-18 — new data in a channel the model never used (13C+α, (α,α₁γ))

Adding the first data to observe a particle pair (here 157 points of 13C(α,α₁γ), exit pair 4)
has three traps, all hit today:
1. **The data cannot act until the channel is seeded**: every amplitude in that pair was
   exactly zero (frozen cost 59 χ²/pt; a polish alone would never move).  Seed only the ZERO
   channels of the physical levels in the data's range.
2. **Do not create the engine keys with a placeholder physical width.**  1 eV in an L = 4–6
   α channel just above threshold is an enormous reduced width; the Brune transform of the
   whole level failed and its *fitted* amplitudes arrived corrupted (+40k whatever the seed
   size).  Write each new channel as θ²·Γ_W *of that channel* (Γ_W from the engine; it does
   not depend on the amplitude) straight into the `.azr`, use a norms-and-shifts-only
   `.sav` (a new channel shifts the positional `width_N_c` names), evaluate two seed sizes,
   keep the lower.  Result here: seeds alone −1,276; polish 321,168 → 314,231 against 322,444
   with no strength in the pair.
3. **Any structural edit of a level nudges its other amplitudes through the Brune
   transform** (~1 % here): amplitudes that sat exactly at a θ² cap land marginally over it
   and `init` reports a sanity flag.  Harmless (the bounded polish clips them), but check
   what the flag lists before believing or dismissing it.
Also: before trusting a secondary-γ dataset's absolute scale, test its detectors against
each other — a γ ray from a J = 1/2 state is isotropic, so the angles must agree (here
χ²/dof 3.1, 10–15 % scale differences).  And `levelscan.scoring_segments` now counts
`differential-cm` (isDiff 4) data as differential (`tests/test_scoring_segments.py`).

## 2026-09-18 — campaign directories go directly under the reaction directory

`rmfit export` writes `<reaction>.azr` only when the campaign directory's PARENT is the
reaction directory (`campaign.reaction_name()` looks for a "+" in the parent's name).  A
campaign nested one level deeper exports to the BASE FILE'S OWN STEM — i.e. it overwrites
its seed (`13C+a_seed.azr`) — and anything chained on `<reaction>.azr` finds nothing.  Found
when a five-step continuation aborted at its first hand-over (its guard worked; the refit
was intact in the overwritten seed).  Rule: one dated directory per fit, directly under the
reaction directory — which is the archive's convention anyway — and a chained job must
check for the previous step's export and stop if it is missing.

## 2026-09-19 — one bad point can hold a converged fit hostage (13C+α, MANA)

A per-ANGLE χ² table of a 14-detector dataset showed one segment at 4× its neighbours; it
was a single point, a factor 6.6 below both energy neighbours in the same detector and
below the adjacent detectors at the same beam energy (pull −61, χ² 3,746).  Dropping it
(evaluator's decision, new data file, old file untouched because a running job read it)
removed its 3,746 — and the refit then fell ANOTHER 3,743, in datasets that do not contain
the point (other channels' angular distributions, the neutron total), although the model
with the point in was converged.  Rules: (1) routinely print χ² per segment AND the largest
single-point χ² per segment; a segment whose maximum is most of its total is a data
question, not a physics one; (2) never edit a data file a running job reads — new name, new
directory; (3) predict the frozen effect of a data change before spending a node on it (here
the prediction matched to the digit), so that everything beyond it is attributable to the
refit.

## 2026-09-19 — levels added by the search had FIXED energies (method gap)

`structure.LevelCandidate.energy_fixed` defaults to True and `levelscan` passes True, so the
added level is written with `levelFix = 1`.  The "released polish" of `run_level_add`
therefore releases the rest of the model but never the new level's energy, and the accepted
structure keeps the flag, so no later polish moves it either.  Found only because an accepted
level read out at exactly 10.9 MeV, and then every level added since 09-14 turned out to sit
at a round number.  The plan had specified a released energy (±max(3Γ, 50 keV)).  Until the
code path is repaired (pinned stage by BOUNDS, released stage with the energy free in a
narrow window; 13N test), follow every accepted add with a job that frees the new level's
energy (flip token 4 of its channel lines) and polishes.  Read-out rule: after any accepted
structure change, print the new level's parameters and ask whether each one could have moved.

## 2026-09-19 — run_level_add now fits the new level's energy; and restarts are not the limit

**Repair.** `run_level_add` writes the full-model variant with a FREE energy for the new
level (whatever `cand.energy_fixed` says; set `cand.release_energy = False` to opt out), keeps
that energy out of the pinned stage's free mask, and in the released stage bounds it to
±max(50 keV, 3Γ) (Γ from θ²·Γ_W of the level's own channels), capped at `policy.e_window`.
The verdict line and the round stats report `energy_fit` / `e_moved`.  `classify(...,
e_window=)` calls a level "relocated" when it ends within 2 % of the EDGE of that window (the
old fixed 50 keV would have rejected a legitimate 60 keV move of a 300 keV-wide level; the
default keeps the old behaviour).  Test: `tests/test_level_add_energy_13n.py` (slow, ~1 h:
the 13N 3/2⁻ 3.5032 offered 30 keV high comes home to 3.5036).  Levels accepted BEFORE this
date have `levelFix = 1`; free them by flipping token 4 of their channel lines and polish
(on 13C+α four such levels were worth 2,453, shifts 4–23 keV).

**Restart test.** A polish that stops at max_nfev "still descending" IS continued adequately
by a new polish: one 280-evaluation call gained 13.6 where a restarted 200 had gained 16.
The large per-iteration gain seen at the end of a polish right after a structural release
is the tail of that release, not a lost descent.  No change of practice needed.

## 2026-09-19 — AZURE2's 1 keV level-merge tolerance is permanent; rmfit now guards against it

`JGroup::IsLevel` merges same-J^pi levels closer than 1 keV when a model is READ.  A change to
1e-6 (09-16) was reverted at deBoer's direction (09-17, azure2-eval SKILL.md) and pyazr was
rebuilt from the reverted source on 09-19.  A 13C+α lineage that still carried a 9/2⁺ pair
0.67 keV apart then read back as a different model (301,445 → 327,514; 583 keys instead of
586) and its `export` failed with "free parameters of the base model are not in keys".
What to take from it:
* A fit with two same-J^pi levels within 1 keV CANNOT be saved and re-read.  `local.sanity`
  now fails on such a pair (`merge_pairs`) and warns below 1.5 keV (`close_pairs`), using the
  free energies from the vector, not the file (`tests/test_sanity_same_jpi.py`).
* The bake guard's message "N free parameters of the base model are not in keys", with all N
  under one level label, means the level NUMBERING shifted: a level vanished on read.  Look
  for a sub-keV same-J^pi pair before anything else.
* Shared binaries change under you.  When a result stops reproducing, check the mtime of
  `pyazr/_azure2*.so` and `~/bin/AZURE2` and read the other skill's log before debugging your
  own code — and when two sessions work one tool, a deferred decision in one (“keep both
  levels, sort it out later”) can be overtaken by a decision in the other.
* A one-level template cannot test for a second state on top of an existing same-J^pi level
  (the Brune transform is singular; the window scan returns +inf).  Test the EXISTING level's
  width and channel partition instead.

- 2026-09-20 -- `add_level_variant` NEVER SEEDED OMITTED CHANNELS, whatever `LevelCandidate`'s
  docstring said (`seed_fraction` was an unused argument): only the channels the scan pattern
  names reach pyazr `add_level`, which writes every other channel gamma = 0, FIXED.  A level
  found through "channel 7 only" had one free channel for ever, could not carry inelastic
  strength, and was judged with those channels frozen (13C+a high-Ex window: 5 of 7 added
  levels had (a,a1) fixed at zero, two had a single free channel; DeBoer asked why narrow
  resonances were not fitting the (a,a1 gamma) data).  New, OPT-IN: `seed_omitted=True`, or
  `setenv RMFIT_SEED_OMITTED 1` in the job (no caller changes; default = old behaviour).
  `omitted_open_channels` adds the group's other PARTICLE channels that are open with
  Gamma_W = 2 P_L 3hbar^2/(2 mu a^2) >= 10 x the 1 eV dead rule -- decided offline with
  pyazr.transform.penetrability, so none of them is later zeroed through bounds into an
  exactly-zero free key (bake lost two such keys on 09-18) -- seeded at min(1 eV, 1% Gamma_W),
  i.e. theta^2 <= 0.01 per channel (a fixed 1 eV on a barely-open high-L alpha channel is a
  huge reduced width and breaks the level's Brune transform: the low-energy session's 09-18
  entry).  COST: the extra channels count in k_new, so tau rises (~1k for a one-channel
  level, ~18k with ~22 channels on 13C+a).  Test: `rmfit/tests/test_seed_omitted.py`.
  Companion traps from the same week, both in memory and the 13C+a readmes: an amplitude at
  exactly 0 in a channel no level of its J^pi feeds is a dead Jacobian column and never moves;
  and a channel that starts with zero strength must not enter with a FREE norm (the norm
  falls to its floor: right shape, 1/5 scale, for a week).

## 2026-09-20 — after a login-node REBOOT the Kerberos caches are gone

The cache I had pinned for qsub (`/tmp/krb5cc_<uid>_...`, see the 09-18 entry) vanished when the
login node rebooted (uptime 4 h; new kernel): "job does not provide an AFS token" again.  Do not
hard-code a cache path across days.  Before each qsub: `ls -t /tmp/krb5cc_<uid>_*`, take the
newest whose `klist` shows a krbtgt and whose time matches `tokens`' expiry, export KRB5CCNAME,
then submit.  If none is valid, ask the user to log in / kinit -- nothing else will work.
Also check the job id actually came back before writing it into a readme (a failed qsub prints
only "Exiting." and an awk over qstat then returns an empty string).

- 2026-09-20 (later) -- "EFFECTIVE PARAMETERS" IN tau WHEN CHANNELS ARE SEEDED (DeBoer's decision).
  With `RMFIT_SEED_OMITTED=1` a new level carries ~22 extra free channels; counting them all in
  tau = max(s^2 k ln N, 0.2% obj) put tau at ~18k on 13C+a -- above EVERY level the low-energy
  campaign accepted that month (+13.9k, +5.1k, +2.8k, +2.0k, the last a 1/2+ feeding a1 with a
  6 keV s-wave: exactly the case the seeding exists for) -- although most seeded channels end
  where they started.  `structure.effective_k`: the scan pattern's channels always count; a
  seeded channel counts only if the fit USED it -- partial width (theta^2 x Gamma_W) >= 1% of the
  level's summed width, or theta^2 >= 0.05.  The width share is the primary test on purpose: the
  seed is min(1 eV, 1% Gamma_W), so on a channel with Gamma_W = 6 MeV a final 5 keV is theta^2
  ~ 1e-3 yet dominates a narrow level.  `run_level_add` computes tau from k_effective, logs
  "k for tau: pattern / effective / all" with the three taus when anything was seeded, and
  stores k_pattern, k_effective, k_all, tau_pattern, tau_all in the round's stats, so a verdict
  can be re-read under another convention.  Nothing seeded => identical to the old k and tau.
  Tests: `rmfit/tests/test_effective_k.py`; end-to-end `test_campaign_13n.py` passed with the
  switch off and on (job 1460168) -- but 13N's groups have ONE particle channel, so that run
  seeded nothing; the non-empty branch was exercised offline only.  WATCH the first real use.

## 2026-09-20 — a template screen can manufacture a level (13C+α, 7/2⁺ at 10.42)

The residual-driven screen reported a 7/2⁺ worth +1,693 (after the GN step) on top of an
existing 7/2⁺.  Its winning template was an "equal split" of 10 keV over all 15 open channels —
667 eV each — although five of them have Wigner limits of 0.4 eV … 1e-15 eV (high L, barely
above threshold): θ² of 1e3 … 1e17.  The screen evaluates templates FROZEN AND UNCAPPED, so the
"gain" came from impossible reduced widths.  The honest add with only the channels that can
carry width (n₀, α₀, n₂ L=1) gained 28 with the energy pinned; released, the candidate shrank to a
~35 eV level 18 keV away (+678, "collapse", provisional).  Warning signs that were there from
the start and should stop a refit: the candidate was absent from the LOCAL (differential) top
five while leading globally, and it sat within 0.1 keV of a same-J^pi level.  To fix in
`levelscan`: build templates only from channels above `run_level_add`'s dead-channel rule (or
cap θ² ≤ 1).  Separately: the `collapse` gate (max θ² < 1e-4) mistakes a NARROW level for a
vanished one when Γ_W is 0.3–1 MeV; gate on the physical width instead.  And before testing
"can the existing level absorb it?", print Γ_W per channel: here the partition question was
answered in ten minutes on the login node (the other channels cannot hold 50 eV between them).

## 2026-09-20 (evening) — both fixes from the entry above are installed

1. **Templates the refit can reach** (`levelscan.py`).  `template_screen(dead_channel_ev=1.0)`:
   a particle channel enters the templates only if its Wigner limit Γ_W ≥ `dead_channel_ev`
   (the same rule `run_level_add` uses before its released polish); the ones left out are logged
   per J^pi ("N channel(s) with Gamma_W < 1 eV left out of the templates").  Every template
   amplitude is clipped at the policy's θ² cap (`theta2_cap` below `bg_energy`,
   `bg_theta2_cap` above), so a frozen template can no longer hold a reduced width the bounded
   refit could not.  Pure helpers: `carrying_channels`, `template_amplitude`, `align_priors`.
   By-product: the partition priors taken from a group's existing levels were truncated
   (`fractions[:n]`) rather than aligned, so removing a channel silently moved a prior's strength
   onto the wrong channel; they are now matched by (pair, L, S).
2. **"collapse" means vanished, not narrow** (`structure.classify`).  New arguments
   `level_gamma_ev`, `collapse_width_ev=1.0`, `narrow_width_ev`: collapse needs max θ² < 1e-4
   AND a summed physical width below 1 eV.  New outcome **`narrow`**: the add passes every test
   but its width is below the campaign's floor — recorded with its gain, not adopted.  The floor
   is `BoundsPolicy.min_level_width_ev` (default 0 = off); a broad-first campaign sets e.g.
   `"min_level_width_ev": 10000` in its `--policy`.  The verdict line now prints `Gamma … eV`
   and the round's stats carry `level_gamma_ev`.  With no width passed the old verdicts are
   reproduced exactly.
Tests: `test_template_channels.py`, `test_classify_width.py` (both FAIL on the old code — run
that way on purpose), the seven quick regressions, `test_levelscan_13n.py` (the dead-channel
rule fires on 13N's 9/2-, 11/2+, 11/2- groups, Γ_W 0.004–0.5 eV; the removed 3/2- still ranks
first) and `test_level_add_energy_13n.py`.  Backups `*.bak-2026-09-20b`.  NOT yet exercised on
13C+α: the first residual-driven scan after this date is the real test — compare its top
global candidate with the local (differential) top five before spending a refit on it.

## 2026-09-21 — first real use of the capped templates; a level add can be a SPLIT

13C+α, `9-20-26_levelscan_resid2`.  With dead channels out and θ² capped, the 10.37/10.46 MeV
residual energies screen NULL (+73, +75 vs τ 610) and global and local rankings agree on every
leader — the 09-19 +1,693 was the artefact.  Cost unchanged (~1,990 patterns, ~53 min per energy).
New failure shape the criteria do not catch: a 1/2⁻ added at 9.8637 passed everything (+1,230,
ρ 0.39, Γ 11 keV > the 10 keV floor) but the existing 1/2⁻ 5.6 keV away gave up exactly that
much n₀ width (41.0 → 32.6 keV; sum conserved) and 96% of the gain sat in ONE dataset with 1 keV
steps.  That is one level's lineshape being described by two — check, before adopting any add:
(i) nearest same-J^pi level and its width before/after, (ii) the share of the gain in the top
dataset.  To build: a "splits" outcome in `run_level_add`.  Also: `tail -F` on an NFS log misses
lines — watches must poll (`grep -c` in a sleep loop).

## 2026-09-21 (morning) — "splits" is now a verdict; the screen's "nan" explained

* **`splits`** (`structure.split_evidence`, `classify(splits=…)`, wired in `run_level_add`): returned
  only where `accept`/`narrow` would be, when a same-J^pi neighbour with a free energy is
  unresolved from the new level (separation ≤ max of the two widths) AND lost ≥ 50% of the new
  level's width AND ≥ 5% of its own in the released polish.  Recorded, NOT adopted — the evaluator
  decides.  The old "PROVISIONAL: within 20 keV …" line stays, but it only warned: the campaign
  adopted the level anyway, and a watch that greps for verdicts never shows it.  Also logged now:
  "N% of the summed dataset gains come from <dataset>" (stats `top_dataset_share`); a share near
  100% in a finely stepped dataset is the lineshape signature.
* **Replay before installing a rule.**  The first version (no 5% test) was run on the real case
  with one live `SingleSession` (variant + the two ledger vectors, ~3 min): it named a 29 MeV-wide
  background level 2.2 MeV away as the partner (it shed 16 keV, 0.05% of itself).  Unit tests built
  from my own idea of the case would never have shown that.  Overlap alone is also wrong: the
  accepted 5/2⁻ 10.886 overlaps the 5/2⁻ 10.940, which GAINED width.
* **"nan" in "best after GN"** was not a failure: with 12 J^pi groups × `gn_per_jpi` 2 = `gn_top`
  24, every GN slot is reserved and a good frozen row of an already-served group is never
  stepped; it sorted by its frozen value and printed the placeholder.  Now printed as
  "(frozen, not GN-stepped: +gain)".  If such a row ranks in the top five, raise `--gn-top`.
* Housekeeping: `pkill -f <pattern>` from a tool shell matches the shell's own command line and
  kills it — kill by pid.  Install order when a peer job may start a new `rmfit.cli` step at any
  moment: callee first (structure → levelscan → campaign), so every mix of old and new is valid.
Tests: `test_split_evidence.py` (fails on the old code), classify tests, 13N levelscan and
level-add regressions.  Backups `*.bak-2026-09-21`.

## 2026-09-21 (night) — J^pi reassignment study; two-channel templates; when a flag is a grid gap

**The study** (`13C+α/9-21-26_jpi_study`): for each broad level (Γ ≥ 10 keV, 36 of them) the
incumbent with that level REMOVED is a sub-campaign (fitted norms/shifts kept via a
norms-and-shifts-only `.sav`); `rmfit levelscan --ex <E> --gamma <Γ> --scales 0.5,1,2
--ex-offsets 0 --refit 0` then ranks all J^pi at the level's own energy and width, ~30 min per
level on one node; `summary.py` reads the campaign logs into one table.  Read it as: the level's
own J^pi must win; the GLOBAL ranking is meaningful only near the optimum, so a large removal
(hundreds of k) makes the global numbers a poor discriminator and the LOCAL (differential)
ranking is the one that carries J^pi — EXCEPT for a level with no α width (3/2⁻ 9.359: 271 keV
n₀, 44 eV α₀), which the (α,n₀) angular data cannot see at all: there the local gains are noise
(≤ 1.3k of a 276k removal) and the global ranking decides.  Rule now in `summary.py`: the local
ranking counts only when its best gain ≥ max(1,200, 10% of the own global recovery).
**Before believing a flag, run the drop-in.**  The 9/2⁺ 11.070 (199 keV: 44% n₀ g-wave, 56%
n₂ p-wave, nothing else) was flagged: its own J^pi recovered +2,724 of a 33,979 removal.  On the
screen's own templates variant, with the campaign's incumbent vector mapped on and all dummies at
zero (the file-derived rwa of a variant with dummies is garbage — the dummies' failed Brune
transforms corrupt their whole J-groups; start from `ev.x0` and overwrite the mapped keys, keep
the dummies' energies, zero their widths), the level's own partition dropped into the dummy
recovers 99%; "n₂ only" gives +102 (= the screen's number), "n₀ only" −132k.  The level lives in
the PRODUCT of two channels and the grid had no two-channel shape; `grid_coverage` said 0.107
because the nearest prior puts 23% into α₀ L=5, a channel the level does not carry and the
elastic data reject — distance understates a gap when a prior leans on a low-Γ_W channel.
**Fix**: `_patterns(..., weights=Γ_W of the carrying channels, pair_top=4, pair_splits=(0.5,
0.7))` adds a "pair" family: pairs among the top-4 channels by Γ_W plus the priors' channels,
50/50, 70/30, 30/70, both relative signs (≤ 60 shapes more; 74 → 134 on 13C+α); wired in
`template_screen`.  `tests/test_pair_family.py` (distance 0.107 → 0.007 for the 9/2⁺ case);
13N regression ALL OK (inert there: one particle channel per group).  Installed after the study
job ended (a running job re-imports `levelscan.py` at every step — never edit under it);
rescreen of the two affected levels queued behind it, gated on an install marker.
  First real result of the pair family (09-22 03:56): the 9/2⁺ 11.07 group's best template became
  "pair n₀ L4 + n₂ L1 50/50" at 94% of the removal, frozen (was 0.3%); with the grid able to
  represent the level, a 7/2⁻ (n₀ f + n₂ d, opposite signs) turned out 716 better — a real near-tie
  the first screen could not have shown.  Two reading rules that came out of the 36-level table:
  (i) a "runner-up within 1.5× globally" mark is meaningless for a huge removal (a 317 keV level's
  fixed-shape template stays 130k above the incumbent however right its J^pi is) — compare
  RECOVERIES, not distances; (ii) `summary.py` must read a level's LAST screen (a rescreen appends
  a round to the same campaign.log).  Stage 2 = `rmfit structure add --candidate i --jpi <J^pi>`
  on the removed-level campaign with a `candidates.json` holding the runner-up (widths mapped
  onto that J^pi's allowed L, same pair/magnitude/sign) AND the own J^pi as a same-budget control;
  each campaign judges against its own (removed) baseline, so the study's verdict is by hand.
  Stage-2 false start (09-22 04:49): `candidates.json` channels written as `gamma_eV` (what the
  `LevelCandidate` docstring said) while `add_level_variant` read `spec["gamma"]` — every add
  returned `scan_failed` in 4 min and the driver, testing only the exit code, marked them done.
  Fixed both: `structure.py` accepts `gamma` / `gamma_eV` / `gamma_ev` (docstring now says
  `gamma`), and a driver marks an add done only when the printed JSON `outcome` is not
  `scan_failed`/`crashed`.  Rules: `rmfit structure add` exits 0 on a failed scan — test the
  outcome, never the exit code; and never edit live rmfit code in the same breath as a job
  submission (the resubmitted job imported `structure.py` three seconds after a non-atomic
  rewrite of it — it happened to import cleanly; write to a temp file and `mv`).
  Second stage-2 lesson (09-22 06:53): `run_level_add` ADOPTS an accepted structure as the
  campaign's base model, so two adds in one campaign are sequential, not alternatives — the
  "own J^pi control" would have been added on top of the accepted runner-up.  Alternatives that
  must start from the same model each need their own campaign directory (init from the same seed).

- 2026-09-22 -- "BROAD" VERDICT (DeBoer).  Report figures of the 13C+a high-Ex window showed its
  broad levels fitting the AVERAGE (a,n0) cross section, not resonance shapes: of 17 free levels
  in Ex 11.2-13.3 only one was under 300 keV, the two original seed levels had grown to 7 and
  26 MeV, the four 12.5 MeV background poles to 0.7-3.2 MeV.  theta2_cap alone does not prevent
  it -- Gamma_W in open neutron channels is several MeV -- and nothing bounded a level's TOTAL
  width.  Mirror of the min-width rule: `BoundsPolicy.max_level_width_ev` (default 0 = off);
  `classify(..., broad_width_ev)` returns "broad" (not adopted) where "accept" would be, after
  "splits" and "narrow"; `run_level_add` passes it from the policy.  Used with 500000 eV in
  `13C+a/9-22-26_highE_narrow` after removing every window level wider than 500 keV and the
  12.5 MeV poles (chi2 jumps several-fold by design; the scan rebuilds with narrow templates).
  Note it only governs LEVEL TESTS; polishes can still widen existing levels -- watch survivors.
  Test: `rmfit/tests/test_classify_broad.py`.  Also from the same week: confine scans to an Ex
  range with a two-stage job (levelscan --refit 0, filter the readme's "levelscan at Ex" lines,
  then levelscan --ex <energies> --refit N) -- there is no range option in levelscan.
- 2026-09-23 -- THRESHOLD GUARD on free level energies.  In `13C+a/9-22-26_ntot_fix` the polish
  converged (1,403,526 -> 311,680) but the export failed in `bake -> ev.transform` with
  "z is not finite in log_Gamma": a 3/2- level had walked from 10.278 to 10.2730 MeV, 0.4 keV
  above the n2 threshold (10.2726), and the Brune transform diverges there (log-Gamma of a channel
  at zero energy).  The fit itself never saw it -- the rwa-space objective is finite on the
  threshold, only the physical-parameter transform is not -- so the failure shows up at export,
  after the node time is spent.  Bisection recipe that found it (`halve.py`): swap halves of the
  vector (widths / energies) between the seed and the fit until the transform fails, then single
  keys.  Fix: `BoundsPolicy.threshold_margin_mev` (default 0.003); `default_bounds` reads the
  particle-pair thresholds from the evaluator's .azr (`particle_thresholds`: sepE + residual
  excitation, photons excluded, cached) and `keep_off_thresholds` cuts every free energy window so
  it stays on the seed's side of each threshold by >= the margin; a level already inside the
  margin is moved to it by the polish's `bounds.clip`.  Cost on this model: 217 (311,680 ->
  311,897) for moving the level 3 keV up.  Tests: `tests/test_threshold_guard.py` (pure) and the
  three new lines in `tests/test_bounds_13n.py`.  Backups `local.py.bak-2026-09-23`.

## 2026-09-24 — the level search never estimated a width; `--gamma auto` now does (13C+α high-Ex)

DeBoer, on the 500 keV `max_level_width_ev` cap of the narrow campaign: "there likely are
some broader resonances in this region... do any of your diagnostics try to estimate the
width of the resonance type structures in the data?"  They did not.  The template screen
sampled `--gamma` × `--scales` (5 keV × 0.3..10 = 1.5–50 keV on that campaign) and the best
scale was only a tag in the log: it sat at the TOP of the grid 95 times out of 120 screens.
The released fits in the same window (9-22-26_highE_refresh, before the cap existed) came out
at 0.63, 0.63, 0.70, 1.2 and 2.1 MeV; one 7/2⁺ at 17 keV.  The cap set against absorbers was
also excluding the widths the region shows — the absorber signature is channels on the θ²
cap and a failing Brune transform (the `unphysical` gate), never the width alone.

Installed (`rmfit/levelscan.py`, `campaign.py`, `cli.py`; `tests/test_cluster_width.py`):
- `levelscan.cluster_width(runs)` → `width` (χ²-weighted MEDIAN extent of the cluster's
  residual runs in DIFFERENTIAL segments, ≈ FWHM of a peak or one lobe of a dispersion
  shape), `width_total` (the same over total/integrated segments) and `span` (full extent),
  MeV.  Three rules learned on the first tries, each from a real run: a run in an angular
  distribution at one energy has zero extent (`ne` < 3 distinct energies: skipped); the
  Cierjacks n-total's single 791-point run over 560 keV owned a χ²-weighted MEAN (hence
  median, and totals reported apart); a Ge (α,α₁γ) run covering 88% of its segment is a
  normalization deficit (`cover` > 0.5: skipped).  `ExCandidate` carries all three; the
  `residual-run candidates` log line and the readme `levelscan at Ex` note print them.
  Measured on the narrow incumbent (`9-11-26_add_MANA_an0/diag_cluster_width/`): differential
  widths 75 keV at 11.07–11.26, 114 keV at 12.48–12.86, 275–345 keV at 11.60–11.77, ~590 keV
  at 12.18 and 13.08; the n-total shows 560–580 keV-wide deficits at 10.9–11.46 and
  12.24–12.82.  So the region holds 100–600 keV structures, and the 500 keV cap was cutting
  into them.
- `levelscan --gamma auto [--gamma-bounds 10,1000]`: the template Γ of each screen is that
  cluster's run width clamped into the bounds, so `--scales` × Γ is centred on the data.
  `--ex` candidates carry no runs and fall back to 10 keV.  The numeric `--gamma` is unchanged.
- Policy for a window that shows this: `max_level_width_ev` ~2e6, keep `theta2_cap` and the
  unphysical gate as the absorber guard.  A 6 MeV "level" (the 3/2⁺ at 11.644 → 11.944,
  +486k, half from n-total) is still a background pole, not a resonance.

## 2026-09-24 (later) — stripping broad levels to "rebuild narrow" fails; rank residuals by density

13C+α high-Ex: `9-22-26_highE_narrow` removed the 13 broad window levels (> 500 keV) and let the
search rebuild with narrow ones.  4/4 level tests rejected (relocated at 6 MeV, broad 3.3 MeV, two
unphysical with channels on the θ² cap); export 1.51M vs 390k for the model that kept them.  With
the broad levels gone the largest residual IS the missing broad strength, and every candidate grows
into it regardless of template width or J^π.  **Build the broad component first, add narrow
structure on top** (DeBoer's call, `9-24-26_highE_addnarrow`).

Second cause, independent of the base: clusters were ranked by raw χ² in the residual runs, so a
12–125-point secondary-γ run never outranks a 791-point n-total run and the narrow γ features were
never screened.  Installed (`levelscan.py`, `campaign.py`, `cli.py`; `tests/test_cluster_rank.py`):
- `levelscan --rank density`: clusters ordered by χ² per MeV of structure,
  `cluster_density = score / max(width or span, 20 keV)`; `ExCandidate.density` carries it.
- `levelscan --seg-match REGEX`: only segments whose data-file name matches make candidates
  (`NoFeed` = the secondary-γ sets on 13C+α).
- `levelscan --ex-range lo,hi`: candidates outside are dropped BEFORE the `--top` cut — one pass
  with `--refit` replaces the stage-1 scan / `pick_window_energies.py` / stage-2 re-screen, which
  cost 2.5 h per energy screened twice.
On the addnarrow seed the γ-only density order gave 11.818 (20 keV wide), 12.334, 12.516, 12.438,
12.195, 12.797 MeV; the old order of all segments was led by 1.9 MeV-span n-total clusters.
- 2026-09-24 -- FIGURES OF UPOS DATA NEED DATA-MODE CURVES.  A <segmentsTest> extrapolation can never be UPOS
  (ESegment(ExtrapLine) forces isUPOS_ = false), so an extrapolated curve at a gamma detector's angle is the primary
  particle's angular distribution, not the secondary-gamma one the fit uses -- up to x2 off on 13C(a,n2 gamma), and
  it misled a whole day of reading on 09-22.  Evaluate UPOS curves in data mode on a pseudo-data grid instead
  (13C+a/9-22-26_ntot_fix/figures/upos_curves.py), and always check a figure's curve against the AZUREOut_*.out fit
  column at the data points before sending it.  Draw data at their fitted energy shifts, as the GUI does.

## 2026-09-25 — data-definition fixes outrank structure; measure them frozen first

On 13C+α, three data-definition corrections in three days outweighed every
structural move of the previous week:
- **(n,total) sum**: a neutron total re-summed over all open channels changed what
  the fit must reproduce above Ex 10.2 MeV.
- **(α,n) total yield**: a 4π yield redefined as the sum over n₀–n₄ gained
  −33,506 (309,355 → 275,849). The gain spread over the (α,n₀) angular
  distributions, not just the redefined set.
- **Outlier**: one point was removed after comparing two total-yield sets.

**Method.** Before any structure round, check that each total-type segment sums
every channel open in its range, and evaluate each fix frozen first (CLI mode 1
with the `.sav`). A frozen gain is enough to justify the refit. Afterwards, check
whether normalizations elsewhere relax toward 1: the ND 2020 angle-dependent
factors went from 1.7–1.9 to 1.0–1.2 here, which confirms the fix was real
rather than absorbed. See azure2-eval "Checking what the data and the figures
mean".

## 2026-09-25 — a level pinned at the threshold guard: test both sides, same budget

A free level energy that ends every polish on the guard's edge (3 keV from a
particle threshold) is telling you it wants the other side. It cannot sit *on*
the threshold, because the Brune transform diverges there. So run two variants
from the same seed with the same polish budget: the level where it is, and the
level moved to the other side by the same margin. Judge by verified objectives
and by the data the threshold feeds. Case: 13C+α 3/2⁻ #8, 10.276 against
10.270 MeV around the n₂ threshold 10.273, `13C+a/9-25-26_32m_n2thr`. **Result:** with 200 evaluations
each, the level held above the threshold stayed pinned (273,992.0). Moved below, it walked freely to
10.2625 MeV and won by 969 with the same parameter count (273,023.2). That matches, within 0.5 keV, the
(3/2⁻) level of an independent (α,n₁)-only fit. The guard had held it on the wrong side for three days.
**Also:** `rmfit export` names its file after the reaction directory (`13C+a.azr`) only when the campaign
directory sits directly in it (`Campaign.reaction_name`). A variant one level deeper (`<test>/control/`) falls
back to the base file's stem, so the export OVERWRITES THE SEED, whatever the seed is called. Copy the seed
aside (`*_ORIGINAL.*`) before the job, or pass `export --path`. Corrected 2026-09-26; the 09-25 note blamed
the file name.

## 2026-09-25 — zero amplitudes never move: audit closed channels below thresholds

An amplitude that is exactly 0 has zero slope, so no polish ever frees it,
whatever its "free" flag says (azure2-eval, zero-Jacobian trap). Level
searches place levels with widths only in their open channels, so a level just
below a threshold enters the fit with no sub-threshold (ANC) amplitude in that
channel. On 13C+α all 14 n₁ entries between 9.36 MeV and the n₁ threshold were 0,
while the (α,n₁) data are missed right at threshold.

**Method.** List zero entries of channels whose threshold lies within a few
hundred keV above the level. Seed the nearest ones at +0.05 γ_W, only channels
that are zero, with the 09-18 rules and `13C+a/9-25-26_n1_anc/seed_anc.py`. Run a
same-budget comparison with a control. **Result:** −583 for four amplitudes (passes τ narrowly). Three became real; two
changed sign through zero, and the f-wave one returned to ~0. The (α,n₁) data moved for the first time
(677 → 566). But the steep rise at threshold is still missed: sub-threshold tails on nearby levels do
not replace a level at the right energy. Always look at the figure, not only the χ².

## 2026-09-25 — two copies of one conversation: name one owner

A conversation resumed in a second terminal runs as two sessions with the same
history and name. On 09-24 both adopted the same fit, regenerated the same
figures and submitted the same continuation job, which was caught only by
cross-session messaging. **Rule:** when `ListAgents` shows a peer with this
session's name, agree immediately which copy owns the node jobs and the
Overleaf pushes. The other copy only reads, and the owner is told before any
job. When the user picks one, the other stands down explicitly: no watches, no
jobs, no edits.

## 2026-09-25 — width caps by energy window (`BoundsPolicy.theta2_cap_windows`)

`theta2_cap_windows = [[e_min, e_max, cap], ...]` (MeV, by each level's seed energy, first match wins) overrides the
θ² cap for the levels in a window only; everything else keeps `theta2_cap` / `bg_theta2_cap`. All cap lookups
(bounds, sanity, level adds, template clipping) go through `BoundsPolicy.cap_for(e)`. DeBoer asked for this on 13C+α:
"don't relax the reduced width limit for all of the levels in the fit, just those near the n1 threshold" — a
global relaxation would also free every background pole and every saturated channel at once, which tests
something else. Pass it in the `--policy` JSON, e.g. `"theta2_cap_windows": [[9.95, 10.60, 10.0]]`. Test:
`tests/test_cap_windows.py`. Backups `*.bak-2026-09-25`.

## 2026-09-25 — copying a campaign directory: repoint the ledger's structure paths

Parallel screens on copies of a campaign (13C+α `9-24-26_highE_addnarrow/fanout/`): all eight jobs
died in seconds with "campaign not initialized: run `rmfit init`".  `Campaign.structure_id` matches
`campaign.json`'s base model against the ledger's `structures.azr`, which is an ABSOLUTE path into the
original directory.  After copying, run
`UPDATE structures SET azr = replace(azr, '<orig>/', '<copy>/')` on the copy's ledger.sqlite and check
`Campaign('.').incumbent()` resolves before submitting.  Candidate vectors (`candidates.npz`) may stay
pointing at the original: they are only read.  A refit job held with `-hold_jid` on the screens then
starts at once when they all fail: stop it, and check the main ledger for stray rounds.

## 2026-09-27 — per-level theta^2 cap for added levels; the narrow loop

- `BoundsPolicy.new_level_theta2_cap` (0 = off): `Campaign.run_level_add` now wraps `_run_level_add`
  and gives the NEW level its own cap through a ±0.5 keV `theta2_cap_windows` entry at its seed energy
  (`cap_for` keys on the seed energy in the model file, which is unique per added level).  In force for
  the whole test; kept in campaign.json only on accept (adopt_structure writes the policy), removed
  on reject or exception.  Test `tests/test_new_level_cap.py`.  DeBoer set 3 for the narrow levels of
  13C+α (theta^2 = 3 is 3x Wigner: watch for absorbers; the unphysical gate tests against the new cap).
- Narrow loop (13C+α `9-24-26_highE_addnarrow/refit_narrow_loop.py`): a level added for a narrow
  feature that is accepted but ends > 2.5x the feature's measured width has "run off" to fill broad
  strength -- keep it, re-screen the same feature on the new model, add again (<= 3 attempts).
- Re-screens must not put any offset within ~1.5 keV of an existing level (AZURE2 merges < 1 keV):
  check against FITTED energies from the ledger, not only the model file's seeds.
- `pkill -f <script>` from a shell whose own command line contains <script> kills that shell too.

## 2026-09-28 — peak-weighted chi² and a per-channel width bound for narrow-level tests

DeBoer's alternative to the run-off loop: fit the strong narrow resonances first with the chi² weighted
to their peaks, then fill in the underlying structure.  Two additions, both off by default:
- `ObjectiveSpec.point_weights = [[ex_lo, ex_hi, w], ...]`: chi² × w on the data rows whose excitation
  energy lies in a window (`Evaluator.row_excitation()`, single and sharded); penalties unweighted,
  per-segment tables plain; `spec.plain` is False so the row path is used.  Verdict deltas under
  weights are NOT comparable to τ — judge on a plain polish afterwards.
- `BoundsPolicy.width_bound_windows = [[e_lo, e_hi, gamma_max_ev], ...]`: each free partial width of a
  level seeded in a window is bounded, θ_c² ≤ Γ_max/Γ_W,c, on top of the θ² cap.
Comparison protocol (`13C+a/9-24-26_highE_addnarrow/compare_peakweight.py`): campaign copy → adopt the
pre-run-off structure → same screen candidate under the variant → release (plain, unbounded, 150 evals)
→ plain control → plain delta and the level's width after release.  Copies via `fanout/make_copy.sh`.

## 2026-09-28 (later) — three traps from the run-off comparison

- **`ledger.add_structure` returns the EXISTING structure when the name repeats.**  A re-test of a
  candidate tested before (same J^π, same seed energy → same `add_<jpi>_E<..>.azr` name) inherited the
  old structure's candidates: `incumbent()` then returned the OLD run-off fit and every "release"
  polished that.  Fixed: `_run_level_add` gives a repeat a `_r<n>` suffix.  Check `stats["structure"]`
  is a NEW id whenever a candidate is re-tested.
- **A width on a `width_bound_windows` bound is not an absorber.**  The unphysical gate counted every
  `bounds_hit`; with a per-channel width bound the level is SUPPOSED to sit there.  `local.width_bound_pins`
  separates width-bound pins (logged) from cap pins (unphysical).
- **Any script that opens a SingleSession/Sharded needs `if __name__ == "__main__":`** (spawn-mode
  multiprocessing re-imports the main module) — bitten twice now (09-24 diagnostic, 09-28 doublet removal).
- Result of the comparison's valid part: a per-channel bound of 2× the feature width keeps ~95% of the
  gain at 12.68 (264 keV vs a 1.12 MeV run-off) and ~60% at 13.23 (136 keV vs 854 keV).  DeBoer's method:
  narrow levels fitted at their widths (bound permanent on accept), underlying structure afterwards.

## 2026-09-29 — csh job scripts: never put JSON in a shell variable

`set POL = \`python3 -c "...json..."\`` then `--policy "$POL"` died with "Missing '}'." — csh
brace-expands `{"a": 1, "b": 2}`.  Literal single-quoted JSON on the command line has always worked;
anything that passes through a variable or backquotes does not.  `rmfit init --policy @file.json`
(also `--search`, `--flags`, `--objective`) now reads the JSON from a file — use it for long policies
(cap windows).  Third false start of one campaign in a night, after two multiprocessing-guard slips:
dry-run every new helper script on the login node for a few seconds (the guard error appears at once)
before putting it in a job.

## 2026-09-30 — four lessons from the 12C+α campaign (`12C+a_onefile/9-8-26_rmfit_campaign/`)

- **`bg_theta2_cap` 3.0 was far too tight for α+¹²C elastic; ~100 was right.**  With 3.0, every
  `polish --tighten` stage 2 and every `structure add` pinned polish started at 1.7e8-2.4e8 and ended
  3.5-4.5× WORSE than the baseline (worst dataset each time: the high-statistics yield-ratio set
  ND_yields_2013_ratios).  Six code-level hypotheses were tested and refuted first (shard assembly,
  worker state, score vs residuals(jac=True), seed size, ...); the cause was policy: this channel has
  always needed large background amplitudes, and a small-error, many-point ratio set makes forcing them
  into the box catastrophic.  Re-`init` with `bg_theta2_cap: 100` (physical `theta2_cap` unchanged):
  stage-2 start fell 50-250×, tighten converged, and a structure-add delta went from -383,539 to -1,130
  (a trustworthy verdict).  Tell: the SAME background-pole keys pinned at the cap in every failed
  attempt, and the worst dataset is the one most sensitive to them.  Ask the evaluator how large a
  background contribution the reaction has historically needed before blaming the code.
- **`export` can fail `bake`'s 1e-4 verify after a tighten polish.**  Structure 8 (built on a tightened
  base, background amplitudes near the 100 cap) failed at rel 1.15e-4 / levels-only 1.20e-4, EXACTLY
  reproducible (systematic, not noise).  A direct `bake(..., verify_rtol=5e-4)` passed, and a CLI mode-1
  run with the exported `.sav` reproduced the source objective to 1.7e-6.  `export` has no rtol argument;
  call `bake` directly, and confirm with a mode-1 run using the `.sav`.
- **An in-memory policy override still auto-adopts on "accept".**  A diagnostic `run_level_add` with
  `c.policy = dataclasses.replace(c.policy, theta2_cap=10)` (campaign.json untouched) returned "accept"
  and `adopt_structure` switched `base_azr` to the unvetted structure; twice.  The level's θ² sat exactly
  at the override (9.99999998, then 99.99999972 at cap 100) -- railing at any cap is the tell that the
  level wants unbounded width (background-like), not a real acceptance.  After a diagnostic override,
  restore `base_azr` in campaign.json by hand and confirm `structure_id`/`incumbent()`.
- **`structure add` for a candidate above the data range crashes.**  The window scan keeps only segments
  within ±`half_width` (0.5 MeV, not exposed on the CLI); with zero segments `Sharded._assemble` fails
  (`IndexError` at `self.n_rmatrix = len(rblocks[0])`).  Find the true data ceiling by intersecting each
  segment's declared energy window with its data file (a segment's nominal `maxE` can be a leftover
  broad cut), and for tail-effect candidates call `Campaign.run_level_add(cand, half_width=2.5)` directly
  (`12C+a_onefile/9-8-26_rmfit_campaign/run_level_add_wide.py`).
