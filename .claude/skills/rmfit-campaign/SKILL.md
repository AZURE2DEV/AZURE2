---
name: rmfit-campaign
description: Drive, resume or judge an rmfit campaign -- the resumable global R-matrix fitting procedure on top of pyazr (sharded evaluation, bound-tightening polish, interference-sign rounds, crossover, level add/remove/relocate, levelscan, dataset-off tests, ledger, STATUS.md, export for review) and the idle-time audits run while its jobs run. Use when the task names rmfit, a campaign directory, STATUS.md, a sign round, a tightening polish or levelscan. For hand-driven fits use azure2-fit; for the search logic for missing levels use azure2-resonance-search.
---

# The rmfit process

`rmfit` does what an evaluator does by hand -- polish, flip interference signs, add, remove or move levels, switch
suspect data sets off -- as resumable steps that a script runs on one cluster node and a fresh session picks up from
`STATUS.md`. Every candidate parameter vector and its score are kept in a ledger, so the state of a fit is read, never
recomputed. The package is not installed; it is imported from its checkout:

```bash
export RMFIT_HOME=<directory that contains the rmfit/ package>   # its README is rmfit/README.md
export AZURE2_ROOT=<AZURE2 checkout or frozen build whose pyazr rmfit must import>
export PYTHONPATH=$AZURE2_ROOT:$RMFIT_HOME
python3 -m rmfit.cli <campaign_dir> <command> [options]
```

`$CAMPAIGN_ROOT` below is the reaction directory (e.g. one per compound system) that holds the fit directories.
Read this file top to bottom before the first command on a new reaction; to resume a campaign, sections 1 and 3
suffice. Detail lives in three references:

- `references/level-tests.md` -- levelscan, `structure add`, verdicts, and how to read screens and level tests;
- `references/policy-and-objective.md` -- `BoundsPolicy`, `ObjectiveSpec`, covariances, sharding, Park mode;
- `references/cluster-jobs.md` -- job scripts, frozen builds, interpreters, tokens, processes, quotas.

Siblings: hand-driven fits, MIGRAD, starting parameters before a first global fit, `.sav` handling: `azure2-fit`.
Where to look for missing levels and the width-floor recipe for a feature no search finds: `azure2-resonance-search`.
Whether a result can be trusted (pulls, Wigner limits, radius): `azure2-check`. Data covariances: `azure2-uncertainty`.
Figures of a result for a person: `azure2-plot`.

## 1. Read the state, never recompute it

A campaign directory carries everything a fresh session needs:

- `STATUS.md` -- best verified objective per structure, recent rounds and log lines. Regenerate with `status`;
  `status --json` gives a compact dump.
- `readme` -- the lab notebook. rmfit appends one paragraph per round; add your own with `log "..."`. Read it before
  touching anything.
- `campaign.log` -- every line the driver printed. It and the ledger are the truth for "is it hung?" and for every
  watch. The scheduler's stdout file (`rmfit_*.log`) is written through NFS and can lag the driver by an hour.
  Watches poll (`grep -c` in a loop); `tail -F` on an NFS file misses lines.
- `ledger.sqlite` + `candidates/*.npz` -- every candidate vector and score. The incumbent is the best candidate with a
  `verified` score (a fresh single-session evaluation of the full model). Sharded numbers rank; verified numbers are
  reported.
- `campaign.json` -- base model, engine flags, shard count, bounds policy, search and objective configuration.
- `<reaction>.azr` + `.sav` -- the reviewable export of the incumbent.
- `variants/` -- shard files, window models, with/without-level files (disposable).

Check the queue first (`qstat -u $USER`, and other users' jobs on the same campaign). One campaign command at a time
per directory: commands share the ledger and `variants/`. The login node is for `status`, `report`, `log`, `export`
and short tests (at most 6 workers, at most 30 min). Everything else runs as a job. When two sessions work on one
campaign, name one owner of the jobs; the other only reads.

## 2. Starting a campaign on a new reaction

Prerequisites, done in the reaction's fit directory, not in the campaign:

1. A `.azr` that runs with the CLI, with `data/`, an `output/` directory and a `readme`. Know the CLI flags the
   reaction uses (`--use-brune`, `--ignore-externals`, `--no-long-wavelength`, ...): they are not stored in the
   `.azr`, they differ per reaction (see the reaction directory's notes), and rmfit needs them as `--flags` JSON.
2. Reproduce the last finished fit before building on it: load its `param.sav` by name in pyazr, confirm that the
   objective matches `chiSquared.out` (χ² + normalization penalty), bake it into a clean file (`rmfit.engine.bake`),
   and check the baked file with the CLI in mode 1 with a blank external parameter file. A bake checked with the
   `.sav` as external file proves nothing (the `.sav` overrides `<levels>`). Never take parameters from a running or
   crashed fit, only from the last completed one. Purge inactive level lines.
3. Merge new data sets onto the baked file (`Project.add_dataset`), capped at the energy range the level scheme
   covers (`energy_window`), normalizations free with a systematic error taken from the paper. Look at every
   new file: one negative or spiked cross section can dominate the objective.
4. Audit the input file and the data (section 3b, items 1 and 2) before the first `init`.
5. Without a converged fit to start from, build starting parameters level by level first (`azure2-fit`, "First fits");
   then `init` a campaign on the result.
6. Record all of it in the reaction directory's `readme`.

Rules for `init`:

- **One dated directory per fit, directly under `$CAMPAIGN_ROOT`.** `export` names its file after the parent directory
  (it looks for a `+` in the name). A campaign nested one level deeper exports to the base file's own stem and
  overwrites its seed. A chained job checks that the previous step's export exists and stops if not.
- **Long JSON goes in files**: `init --policy @policy.json` (also `--search`, `--flags`, `--objective`).
- **Seeds from a `.sav` are matched by name and seed norms and shifts only after a structure change.** A `.sav` names
  level parameters by the engine's level index (`energy_N`, `width_N_c`), so a level added or removed renumbers every
  later level. `init --sav` compares the highest `energy_N` with the model and, on a mismatch, warns and seeds
  normalizations and shifts only. The exported `.azr` already carries the fitted R-matrix values (round trip 1e-9);
  reduce the `.sav` to its `segment_*` lines for such a seed. Never read a `.sav` by row position: newer files start
  with a `parametrization` row.
- **Freeing a normalization, or adding data, is a structure change**: edit a copy of the base file
  (`Project.set(<data set>, norm="1 fit prior 1 <percent>%")`, once per data set of that file),
  copy the `.sav` alongside, and `init` a fresh ledger. Crossover donors from a campaign without those norms still
  work: missing norm and shift keys are filled from the incumbent.
- **Energy shifts are free in a campaign only if the policy says so** (`"shifts_free": true` at `init`). Shifts are
  finite-differenced, so by default every shift column is marked fixed whatever the `.azr`'s vary flag says, and
  only segments with a nonzero shift systematic open. `polish --shifts-free` frees them for that one step only;
  `dataset-off` and `structure` follow the policy. Pre-flight: build `default_bounds` through `Campaign.policy` and the
  method's own argument handling, and print `Bounds.fixed` for every `kind == "shift"` key.
- **Pick the background cap with the evaluator.** `bg_theta2_cap` must allow the background the reaction has always
  needed (α+¹²C elastic needs about 100; 3 made every tightened start 3-4× worse). The tell of a cap that is too tight:
  the same background keys pinned at the cap in every failed attempt, and the worst data set is the one most sensitive
  to them.
- **A relvector covariance and a free normalization on the same rows are refused** (`init` names the file and the
  segments). Use the t0 treatment in `references/policy-and-objective.md`.
- **Evaluate the seed before any job**: sanity, bounds and the clipped objective under the policy the job will use. A
  start-up clip that costs more than a fraction of a per cent means a window or a cap is wrong. A staged seed must move
  the objective; an objective identical to the incumbent's means the seed did not reach the model.

## 3. The loop

Every step is a job (`qsub -v RMFIT_CMD="..." <job script>`, or `-hold_jid <id>` to queue it behind the previous
one). Every step ends by rewriting `STATUS.md` and appending to the readme; read both before the next step.

| step | command | purpose | done when |
|---|---|---|---|
| 1 | `polish --nfev 150 --tighten` | walk amplitudes beyond the θ² caps into the box (×10 stages), then converge | all θ² inside the caps; gain < 0.05 % per stage |
| 2 | `round --rounds 2` (repeat) | interference-sign rounds: screen every canonical flip (singles, exposure pairs, same-channel pairs, optional channel patterns), polish the top-k sign-locked then released, greedy union, verify | improving-on-screen and accepted counts fall; stop after two rounds with nothing accepted |
| 3 | `polish --nfev 150 --shifts-free` | warm-start continuation with energy shifts free | gain < 0.05 % |
| 4 | `structure add --candidate i`, `levelscan` | level tests (`references/level-tests.md`) | all candidates classified |
| 5 | `structure dataset-off --file <substr>` for the worst χ²/N sets | deactivate and refit; reports the extra recovery of the other sets beyond the set's own χ² | a table for the evaluator |
| 6 | `round --rounds 2` after any structure change | signs change with the structure | as in 2 |
| 6b | `crossover --donor <cid> [--donor-dir <campaign>]` | apply another good fit's sign pattern per channel, J^π group or all, and polish every proposal | every donor pattern polished; union verified |
| 7 | `export`, then review | bake the incumbent to `<reaction>.azr` + `.sav`, verified in a fresh process | CLI mode 1 with the `.sav` as external file reproduces the totals |

Other commands: `relocate --level 9/2+#1 --ex 8.4654 [--window 5]` kicks an existing level's energy to `--ex` (MeV),
re-polishes with that energy confined to ±window (keV), and adopts if it beats the control by τ -- for a narrow level
that drifted off its peak. `touch <dir>/STOP` ends a `round --rounds N` loop after the current round. `verify <cid>`
re-scores a candidate on a fresh session. `log "..."` appends to the readme. `report` writes `REPORT.md` (per-data-set
χ²/N, θ² table, structure decisions). Every level test and relocation writes a review figure to `figures/`.

Rules of the loop:

- **Every move is judged against a same-budget control polish of the unmodified incumbent.** The gain of a move is
  `control − moved`, never `incumbent − moved`: the incumbent keeps improving by hundreds per polish on its own.
  A control whose start objective differs from the incumbent's is broken; check it before reading a delta.
- **Signs.** Sign rounds skip amplitudes with |γ| < `sign_eps`·γ_W and proposals whose frozen response is below
  `flip_floor`·τ (default 0.01): a flip in a channel no data set observes moves the objective by ~0.2 whatever the
  amplitude's size. Judge a proposal by the response it produces, not by the size of the parameter it changes.
  Frozen screens favour small amplitudes; the Gauss-Newton (GN) screen is the ranking that matters, and multi-level
  patterns need several backtracking GN iterations (`search.gn_iters_multi`) and their own polish quota
  (`search.multi_top`). Some patterns rank badly on every screen and win only after a polish.
- **Count before submitting a sign round.** With `enumerate_channels` a 1,170-parameter model produced 7,407 frozen
  pre-screen patterns (11 h before the first log line); cap `max_patterns` (32 gave 3,643).
- **A second basin is the best proposal source.** Sign-converge, run the `perturb` benchmark, then
  `crossover --donor <old incumbent>`, then polish. For a donor basin perturb 2-4 flips inside one J^π group; flips
  scattered across groups produce a dead basin (8 flips: 1 % above, nothing found; 12 flips: 22 % above and
  unrecoverable). If the sign distance to the target falls by one flip or less in a round, stop the arm.
- **Data changes are measured frozen first** (CLI mode 1 with the `.sav`). A frozen gain justifies the refit; anything
  beyond it is attributable to the refit. Data-definition fixes (a total re-summed over every open channel, a 4π yield
  redefined as a sum, one outlier removed) have outweighed a week of structure moves.
- **Never edit a data file that a running job reads**: new name, new directory.
- **Bake-time guards.** `bake` refuses a model whose re-read key set differs from its source (a level vanished on
  read: look first for a same-J^π pair closer than 1 keV, which the engine merges), nudges free widths at exactly 0 to
  1e-6 (the engine writes no key for a zero width), and separates broad same-J^π pairs closer than 1.5 keV to 2.5 keV
  when the move costs less than 1 % of the broader level's width (`separate_max_cost`, χ²). Narrow pairs are refused;
  remove one member with `structure remove`. Half-written files of a refused bake are deleted; an export that exists
  is not proof of success -- read the log line.
- **Exports carry nominal normalizations in the `.azr`**; the fitted ones are in the `.sav`. Anything that re-evaluates
  an export seeds norms from the `.sav` by name.

## 3b. While jobs run: audit and read

A polish or a sign round takes hours. While a job runs and nothing is actionable, run the three reviews below. They
are read-only and cost no node. Most large gains of a long campaign come from this work, not from more minimization:
a free unpenalized normalization, a data file with statistical-only errors, a width window that swallowed a
neighbour, a wrong elastic/inelastic partition that only independent data exposed.

**When.** At the first idle window of every campaign, and again when data are added, a structure change is adopted or
a fit stalls. Mandatory when a fitted normalization leaves 0.7-1.4; a frozen or pinned part of the model carries a
large share of the objective; two fits of the same data disagree in absolute scale; a level sits at a bound or at a
threshold guard; a new data set is about to enter the fit.

**How.** One background subagent per review, on a cheaper model than the session's, at most three at a time, each with
a self-contained prompt (paths, file-format facts, what is known and must not be re-reported), each writing one report
into `$CAMPAIGN_ROOT/<date>_audit/` next to a `readme` that says why it was run. Reviews never edit a fit directory and
never submit jobs. Their findings are candidates: verify each against the files (number, line, segment key) before
telling the evaluator, and change nothing in a fit on a report alone (data treatment is the evaluator's decision,
section 5).

**1. Input-file audit** (the `.azr` of the current seed):
- particle pairs: masses, charges, spins, parities, excitation and separation energies, radius identical on every
  level line and right against the nuclear data. Print one row per pair (`ir m1 m2 z1 z2 sepE e3 j3 pi3 radius`) and
  compare pairs of the same target nucleus (same z2, same m2, same sepE − e3). Thresholds reproduce to ~1 keV.
  Nothing in AZURE2, pyazr or rmfit checks charges against masses (a carbon channel carried z2 = 8 for years);
- channel sets: every (pair, ℓ, s) allowed by angular momentum and parity, none missing at low ℓ, the same set on
  every level of a J^π group;
- levels: same-J^π pairs within 1.5 keV; free amplitudes exactly 0 (they never move); levels within 3 keV of a
  threshold; fixed/free flags against the campaign's rules. Count the free keys per level and list levels with none:
  a level with `levelFix` (field 4) set is frozen whatever `channelFix` (field 11) says, and a level with every width 0
  is invisible to a gradient polish even when free;
- segments: file exists; entrance/exit keys match the file; energy window against the data's range and the exit
  threshold; `isDiff` against the angle column and its frame; sum tails complete for every channel open in the window;
- segment tails: list every segment with tokens after the file name other than `0 0` (field meanings in `azure2-setup`)
  and record which are composite SUM/RATIO and which are UPOS (secondary-γ) segments. No extrapolation reproduces a
  UPOS observable: curves, scans and post-hoc tests on it run in data mode on a pseudo-data grid;
- normalizations: every `varyNorm 1` with `normError 0` is a free unpenalized norm -- list them with fitted values and
  ask whether each is intended (one at 1.72 hid a 72 % over-prediction of σ_tot); list fitted norms outside 0.7-1.4;
  check every norm and shift setting against the evaluator's recorded decisions;
- `<targetInt>` keys bind the intended segments; sharded convolution rules in `references/policy-and-objective.md`.

**2. Data audit** (every file of an active segment):
- integrity: non-numeric rows, zero or negative values or errors, duplicated rows;
- units and scale: ratios between overlapping sets of the same reaction, band by band; an absolute anchor (a total
  cross section, a unitarity constraint) compared with every independent measurement of the same quantity;
- frames: lab energy of the right projectile (test against thresholds), lab vs c.m. angles and cross sections against
  the segment type; a converted file is checked by evaluating it frozen in both frames with identical errors (the χ²
  must agree exactly);
- uncertainties: statistical-only where the source gives a systematic one; an advertised error floor that does not
  hold; implausibly small or constant errors;
- double counting: the same measurement under two names or in two overlapping windows. Data handed over between
  sessions are compared by value at common energies, never merged by file name;
- energy scale: for a set with narrow structure, scan the shift on the frozen model curve before any refit, and seed a
  freed shift at the scan optimum (`azure2-fit`, energy scale);
- per-segment χ² and the largest single-point χ² per segment: a segment whose maximum is most of its total is a data
  question (one point at pull −61 held 3,746 of χ² and another 3,743 in other data sets after the refit);
- against the source: re-fetch the EXFOR entry (`azure2-setup`) and compare values, units and error columns row by row.

**3. Literature review**: earlier evaluations of the same compound system (data, channels, radii, energy limits,
treatment of normalizations and backgrounds), measurements of channels the fit does not yet constrain, method papers on
the difficulty the campaign is stuck on. Every claim carries a citation; a paper that could not be read is listed as
such; nothing is quoted from memory. The output is a short list of actionable ideas, each with a cheap test.

**Cheap in-session checks** (minutes, no node):
- decompose a summed observable into its partial cross sections and compare each with independent data -- a total can
  be right with the wrong partition;
- the clip check before every submission (section 2);
- a coverage plot of what EXFOR holds against what the fit uses;
- when two consecutive level searches come back marginal or null: χ² in 200 keV E_x bands with the dominant data sets
  per band, and a θ² audit (levels and background poles at their caps). A misfit spread at 12-20 per point over MeV is
  not one missing level; groups whose levels sit at their caps are short of poles, not of resonances.

## 3c. A structure no search finds

The width-floor recipe (localise, floored level add with local release over every J^π, cross-checks, full free refit,
propagate) is in `azure2-resonance-search`.

## 4. Judging results

- **Thresholds scale with s² = objective/N**, because error bars are underestimated (s² ~ 15-20 on a large α-induced
  data set). A flip needs a verified gain ≥ max(0.1 % of the objective, 5 s² ln N) beyond the control polish. A level
  needs ≥ τ = max(s² k ln N, 0.2 % of the objective), θ² inside the caps, an energy that stays in its window, and a
  redistribution ratio ρ ≤ 0.5 (ρ = χ² lost by some data sets / gained by others). k counts the level's effective
  parameters (`references/level-tests.md`).
- **The verified global delta decides.** ρ and the top-data-set share are reported so the evaluator sees the split;
  read the absolute losses too -- ρ > 0.5 built from losses of ~200 is not a redistribution problem. A candidate aimed at
  one data set that improves another counts.
- **Outcomes**: accept / marginal / worse / reverts (the sign walked back in the released polish; never retried on that
  structure) / drift (walked back, but better than the control by ≥ τ: kept, joins the union) / unphysical / collapse /
  narrow / broad / splits / redistributes (a flag, not a veto) / blowup / relocated / live-harmful / crashed. Level-test
  outcomes are defined in `references/level-tests.md`.
- **Convergence of the sign search**: improving-on-screen counts decay geometrically (16 → 4 → 3 → 1). Two rounds without
  an acceptance = converged for the present structure.
- **`sanity`** flags θ² over the cap, |E| > 60 MeV, norms outside (0, 100), a non-finite Brune transform, parameters at a
  bound, same-J^π pairs under 1 keV (`merge_pairs`, fail; `separable_pairs` if bake can separate them) and under
  1.5 keV (`close_pairs`, warn). Any |value| > 1e6 or |E| > 1000 MeV means catastrophic cancellation in a
  total-cross-section χ², i.e. a fake improvement. Read what a flag lists before believing or dismissing it: a
  structural edit nudges a level's other amplitudes through the Brune transform by ~1 %, and amplitudes at a cap land
  marginally over it (harmless).
- **A fitted value of exactly 0.0 for a continuous free parameter is never a result**: it was not optimized. Find the
  log line naming the policy in force, and the round's "dead columns dropped" line (its absence means the column had a
  live derivative, so the parameter was excluded upstream of the optimizer).
- **After any accepted structure change, print the new level's parameters** and ask whether each one could have moved.
- **Read `REPORT.md`'s per-data-set table after every structure step.** A set whose χ²/N is far above the rest, or
  whose removal recovers more than its own χ² in the others, is an evaluation question (normalization, energy
  calibration, resolution), not a fitting question. Dataset-off tables: extra/own > 1 means the set pulls against a
  specific partner (find it in the per-data-set breakdown, `by_dataset`); rank by extra recovery per point, not total
  (small precise sets carry the leverage); polish the reference first and quote the drift of a no-control refit
  (190-1,000 on a 35k-point model) as the floor. When every distorting set sits in one energy region, the table says
  "region", not "data set".
- **Check every fitted norm against the set's fit/data ratio** (at a sensible minimum fit/data ≈ norm). A norm that
  moved the other way points to a data-scaled term in the set's covariance acting as a hidden normalization with
  Peelle's bias (`azure2-uncertainty`). Distrust norms on the `norm_range` bounds while χ²/N ≫ 1: they absorb shape
  misfit.
- **Robust objectives** stand in for uncertainties nobody has examined. Once the evaluator has corrected the error bars,
  use plain χ² (`references/policy-and-objective.md`).
- **The first output of any new diagnostic is checked against an independently known total** before anyone reads
  physics into it (a per-data-set table once double-counted split segments, exactly 2.000×).

## 5. What the human decides

rmfit hands these to the evaluator:

- normalization scale and systematic error of a new data set, and any rows dropped;
- an angle- or energy-dependent trend in fitted normalizations (the norms absorb it silently);
- whether a data set that distorts the rest stays in the fit;
- the energy range of the model (data above the level scheme's range are capped, not fitted);
- a level whose gain is clear but whose J^π or energy contradicts the literature, and every PROVISIONAL, splits,
  narrow or broad verdict;
- error-bar corrections (record each in the readme with the formula, statistics before/after and file checksums).

Review: `$AZURE2_BINARY calculate --parameter-file=<reaction>.sav <flags> <reaction>.azr` fills
`output/` with the fit (the `.sav` carries the fitted norms); then make the figures from it (`azure2-plot`).

## 6. Diagnostics and benchmarks

- `keys.canonical_signs(x, keys, fixed)`: the sign pattern modulo the exact symmetries (all amplitudes of one level; all
  amplitudes of one particle pair in every level; fixed nonzero amplitudes tie a level's gauge to a pair's). Valid for
  particle channels only.
- `keys.sign_distance(x1, x2, keys, fixed)`: the minimal number of flips between two patterns modulo those symmetries.
- `perturb --out <dir> --flips 8 --jitter 0 --ekick 0 --polish-nfev 150`: the recovery benchmark. Flips high-exposure
  signs of the incumbent, polishes to a genuine worse minimum in a sibling campaign, records the target; later
  verifications there log sign distance and gap. Success = back to the target within 1e-4. Run it once per reaction
  before trusting the method there.
- An oracle script that applies a known target's sign pattern per channel, per group or all at once and polishes
  separates "finding the pattern" from "reaching the basin from the pattern".
- `diagnose.dataset_table`, `theta2_table`, `residual_runs`, `compare`, `window_delta` feed `report` and the level
  tests.
- Timing a slow evaluator: time each shard's evaluation and Jacobian separately; the slowest shard sets the wall time
  and its contents say why.

## 7. Known limits

- Validated mainly on α-induced reactions with particle channels only (`ignore_externals`). With photon channels the
  pair-flip symmetry and the shard machinery (`intEC.*` per shard output directory) are less tested: check sharded ==
  single-session at three vectors before trusting a round (`init` checks the baseline).
- Recovery from a start with jittered magnitudes and kicked energies is untested; run `perturb` on every new reaction.
- Energy shifts are finite-differenced: keep them fixed during exploration and free them in the continuation polish.
  Held shift columns are skipped in the Jacobian.
- Norms of split segments are shared across shards; segments with target effects are never split, except as described
  for sharded convolution.
- `rmfit init` scores its baseline on a single session: a check via `init` is not a sharded check.
- Cost on a 35k-point model with 330 free parameters on 24 shards: objective 0.6-1.5 s, Jacobian ~20 s, a 60-iteration
  polish ~9 min, a sign round 2-4 h, a level test ~1.5 h, a levelscan energy ~30-50 min. The TRF SVD (~10 s per
  iteration, single thread) is the largest remaining cost after shard balancing.

## 8. Gotchas

- **One pyazr session per process**: two in one process desynchronise the parameter tables. rmfit runs one process per
  session and `bake()` verifies in a fresh process. Do not use `save_fit` for campaign files.
- **`bake()` writes the fixed parameters' current physical values.** A fixed parameter is fixed in reduced-width space,
  so its observed width moves with the level's free amplitudes (Brune's 1 + Σγ²dS/dE). To keep a "frozen" level's
  observed widths, fix its closed-channel amplitudes too. rmfit builds before 2026-10-24 (expires: 2027-02-01) kept the
  seed's observed values for fixed widths; such exports are a different model: re-export from the ledger vector.
- **New levels are appended at the end of the file and removed by deactivation**, so parameter keys stay valid across
  structures. Exception: a test that re-adds a candidate under the same name reuses the old structure unless it gets a
  `_r<n>` suffix; check `stats["structure"]` is a new id.
- **A released polish lets amplitudes cross zero**: sign patterns are not preserved, and a basin is defined only up to
  the amplitudes near zero.
- **Carry a fit between models by key** (`rmfit.keys.take` + `Campaign.verify_candidate`, then `polish --candidate`),
  never by writing it into `<levels>`: AZURE2's read of large seeded widths and the Brune transform of amplitudes at the
  Wigner limit change the model.
- **A diagnostic `run_level_add` with an in-memory policy override still adopts on accept.** Restore `base_azr` in
  `campaign.json` afterwards and confirm `structure_id`/`incumbent()`. θ² sitting exactly at whatever cap is set means
  the level wants unbounded width, not an acceptance.
- **`total 1`** (optical-theorem neutron total) is in rmfit's keyword table; an engine without it computes elastic
  only, silently. Check the build passes `tests/total_optical`.
- **Delete `output/intEC.*` whenever segments change** (`azure2-calculate`).
- **Multi-nucleus models** label levels `<J^π>#<n>@<nucleus>` (e.g. `g[1/2-#2@21Ne p1 L1 S0.5]`); single-nucleus labels
  are unchanged. Ledgers, checkpoints and `param_priors` written for multi-nucleus models by
  rmfit builds before 2026-10-30 (expires: 2027-02-01) no longer match by name. Sign-round groups, `free_mask(jpi=)`, level
  add and levelscan still work per J^π only.
- **`export` fails `bake`'s 1e-4 verify** occasionally after a tighten polish with amplitudes near the cap (rel 1.2e-4,
  reproducible). Call `bake(..., verify_rtol=5e-4)` directly and confirm with a CLI mode-1 run using the `.sav`.
- **Park mode** (`"use_park": true`) has its own rules: `references/policy-and-objective.md`.
- **Process, quota and job-script traps** (`pkill` with `+`, quotas, interpreters, tokens, `__main__` guards):
  `references/cluster-jobs.md`.
- **Job mail**: `-m abe` (start, end, abort), so a silent crash is noticed.
- **Changing rmfit itself**: every edit gets a test in `rmfit/tests` that fails on the old code, passes `test_keys.py`
  and a sharded fixture before it stays, and is installed callee-first through a temporary file and `mv` (running
  drivers re-import modules at every step). After an install, read the first new job's log past its evaluator start-up.
