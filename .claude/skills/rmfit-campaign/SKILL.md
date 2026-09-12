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
fixed); `e_window` is the energy freedom of a level per polish (MeV); `norm_range`
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
  no "initial step" problem, but the same flat surface). MINUIT (CLI mode 2) cannot fit
  such an energy at all (10% initial step). Remedies to build (Phase 7 with the level
  search): per-level energy windows of a few times max(Gamma, resolution) for narrow
  levels in `default_bounds`, and a relocation move to the nearest residual peak;
  meanwhile map narrow levels by direct (E, Gamma) scans (`9-9-26_seg92_9halfplus/grid2d.py`).
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
- 2026-09-11 -- perturbation size: 8 flips on an already sign-converged 13C+a fit landed only
  1.0% above it and its rounds/crossover found nothing (a wasted 12 h); 12 flips on the same
  model landed 22% above it, a genuinely different basin. For a donor basin or a real
  recovery test use enough flips to move the objective by >= 10%; check the `perturb` line's
  verified value before committing the follow-up rounds.
