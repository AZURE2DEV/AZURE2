---
name: rmfit-campaign
description: The rmfit process for global R-matrix fitting with AZURE2 -- start a campaign on any reaction, drive or resume one, use idle time while jobs run for input-file / data audits and literature review (sharded pyazr evaluation, bound-tightening polish, interference-sign rounds, level add/remove rounds, dataset-off diagnostics, ledger + STATUS.md), judge results, export for GUI review, and refine the method. Use whenever the task mentions rmfit, a campaign directory, STATUS.md, sign rounds, a tightening polish, a level-add test, fitting a new reaction "the rmfit way", or continuing a fit that a previous session left running on the cluster.
---

# The rmfit process

`rmfit` (this repository; the module table is in README.md) does what an
evaluator does by hand -- polish, flip interference signs, add or remove levels, switch
suspect data sets off -- as resumable steps that a script runs on one cluster node and a
fresh session can pick up from `STATUS.md`. It is not installed as a package:

```bash
export PYTHONPATH=<directory containing rmfit>:$AZURE2_ROOT
python3 -m rmfit.cli <campaign_dir> <command> [options]
```

Read this file top to bottom before the first command on a new reaction; for resuming a
campaign, "Read the state" and "The loop" are enough. The design and its open items are in DESIGN.md.

## 0. Two backends: pyazr 2 and pyazr 3

rmfit runs on the pyazr that `AZURE2_ROOT` provides (`rmfit/backend.py` reads its version): pyazr 2 (AZURE2 `dev`,
`.azr` text projects) or pyazr 3 (AZURE2 `feature/upgrade`, project format 3 only).  A campaign is pinned to one engine
build, so it stays on the backend it started with.  What pyazr 3 changes for an operator:

- **Convert first.** A format-1/2 `.azr` is converted once, `AZURE2 --convert OLD NEW` or `pyazr.convert_azr(old,
  new)`; read the converter's report.  A parameter file of the old project is renamed with
  `AZURE2 --convert-params OLD.par NEW.par --project OLD.azr`.
- **Data sets are named records**, numbered in file order for the outputs; rmfit's segment keys are those numbers.
  Edit a project through `pyazr.Project` (or `rmfit.azr3.Model` for level and data-set moves), never by position in
  a line.
- **Normalizations, shifts and their priors are the engine's** (a record's `norm 1 fit prior 1 5%`, a shared
  `normalization` or `energy_scale` record); rmfit takes every prior from the session's parameter rows.
- **The engine folds** target, beam and detector effects: such a data set is never split across shards, and there
  are no convolution pieces.  An extrapolation is folded only if its own record names a target, a beam or a detector
  (the audit reports `extrapolation-folded`; the old `<targetInt>` trap is gone).
- **Model data** (a measured scattering length, a width, a branching ratio), a prior covariance (single session
  only), abundance sums and decay priors are rows of rmfit's objective, taken from the engine.
- **Not on pyazr 3**: `rmfit.mixture` (several compound nuclei in one project replace it).  `rmfit uq` samples with
  rmfit's own NUTS (`rmfit/sampling.py`); a `uq.bands` grid is the fields of an `extrapolation` record plus
  `quantity_band` and `label`.
- **Audit codes on pyazr 3**: `engine-validate` (the engine's validator: masses, separation energies, data files),
  `model-setup`, `norm-prior-tight` (pyazr 2: `norm-error-units`), `extrapolation-folded` (pyazr 2: `targetint-test`).

## 1. Read the state, never recompute it

A campaign directory carries everything a fresh session needs:

- `STATUS.md` -- best verified objective per structure, recent rounds and log lines.
  Regenerate with `status`; `status --json` for a compact dump.
- `readme` -- the lab notebook (rmfit appends one paragraph per round; add your own
  with `log "..."`). Read it before touching anything.
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
   `--no-long-wavelength`, ...): they are not stored in the `.azr`, they vary per reaction
   (record them in the readme), and rmfit needs them as `--flags` JSON.
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

### First fits: the build phase of `rmfit drive`

The build phase of `rmfit drive` (section 3a) encodes how an evaluator starts a new model.  One line per step, with the
state that enforces it:

1. **Fix the scale with one representative data set** -> `SCALE_ANCHOR` (spec `data_sets[].role`: `anchor` = norm fixed at
   the nominal, `free` with `prior_sigma`; `drive.Driver._s_scale_anchor`).  Why: with every norm free, a first fit traded
   missing strength for norms of 0.4-0.7.
2. **Seed the lowest-L channel of every pair on every level** at `model.seed_frac` (0.01) γ_W, placeholders for free widths
   at exactly zero, photon channels at `model.photon_seed_ev` -> `SEED_CHANNELS` (`_s_seed_channels`).  Exercised on a
   from-scratch model (20 background channels seeded).
3. **Every channel that opens inside the fit range is in the model** -> `PAIR_COMPLETENESS` lists the ENSDF states of the
   residual nuclei below the window top that have no pair (`_s_pair_completeness`).  It only reports: the pair is still
   added by hand.  Never cap a data window at a threshold to avoid adding the pair.
4. **Level by level, from low energy up** -> `LADDER` (`model.ladder_Ex`; skipped without it, and not yet exercised).
5. **Then add levels where a resonance has no known level** -> `CANDIDATES` (ENSDF, broad levels first, every entry with a
   `hypothesis`) and, inside the loop, `LEVELSCAN` / `STRUCTURE_ADD` (section 3a).
6. **Only then a global fit** -> `INIT` (`Campaign.init --spec`, baseline verified, engine pinned).

Steps 1, 2, 5 and 6 have run end to end on from-scratch models; 3 and 4 only report.  A campaign initialised before the
drive (`init` first, as the bench does) skips steps 1-2, so data-set roles and seeds belong in the model before `init`.

Then, in a new dated campaign directory (`<reaction>/<M>-<D>-<YY>_rmfit_campaign/`):

```bash
mkdir <dir> && cp <base>.azr <dir>/ && cp -r data <dir>/data && mkdir <dir>/output
python3 -m rmfit.cli <dir> init --base <dir>/<base>.azr [--sav <base>.sav] --nshards 24 \
    --flags '{"use_long_wavelength": false, "ignore_externals": true, "use_brune": true}' \
    --policy '{"theta2_cap": 1.0, "bg_theta2_cap": 3.0, "bg_energy": <Ex above the data>, "e_window": 0.3, "norm_range": [0.2, 5.0]}' \
    --search '{"top_k": 6, "pairs_per_group": 4, "stage1_nfev": 30, "stage2_nfev": 60, "union_nfev": 100}'
```

`--flags` is passed straight to `pyazr.azure2()`, so the parametrization is chosen there:
`"use_brune": true` (Brune, the archive default), `"use_park": true` (Park: the fit vector
holds the observed reduced width amplitudes, `Gamma = 2 P gamma^2`; see the AZURE2 skill's
"Park's parametrization"), `"use_brune": false` (standard Lane-Thomas). One basis per
campaign -- the ledger's vectors, `.npz` snapshots and `.sav` companions are all in it.

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
(a polish raised the plain chi2 by 70 %, a data set was abandoned and the norms snapped
to nominal). Both methods are "statistically incorrect but useful when error bars are
underestimated". The objective
is a property of a campaign: candidates verified under different objectives are not
comparable, so a change of objective means a new campaign directory (seed it with the
previous export). `REPORT.md` and `segment_scores` keep the plain per-dataset chi-squared
whatever the objective. Test: `tests/test_objective_13n.py`.

`init` verifies the baseline (one full-model evaluation), registers the structure, and writes `STATUS.md` and the job
script `run_crc_rmfit`.  The scheduler is a site setting (`rmfit/site.py`, `~/.rmfit/site.json` or `$RMFIT_SITE`):
`"scheduler": "sge"` (the default; a csh script for 24 slots, with the queue from `sge_queue` -- no `-q` line
if unset -- and the parallel environment from `sge_pe`, default `smp`), `"slurm"` (a bash script with
`#SBATCH` lines from `slurm_partition`, `slurm_account`, `slurm_time`, `slurm_extra` and the campaign's shard count as
`--cpus-per-task`) or `"shell"` (a bash script with no scheduler, for a workstation or a node you already hold: start it
with `nohup ./run_crc_rmfit >/dev/null 2>&1 &`, `RMFIT_CMD="..."` in front for a command; it runs the campaign's shard
count of processes where it is started, writes `rmfit_shell.log`, refuses a second run on the same campaign through a
lock directory, and the drive script runs its legs one after another instead of resubmitting).  All of them run the
site's `setup` lines before Python -- modules, a Python with numpy and scipy on PATH -- written in the script's shell:
csh for SGE, bash for Slurm and shell.  Without them the job uses the environment it starts with.  Job mail goes to `RMFIT_MAIL`
or `"mail"`; with neither, the scheduler mails the submitting user.  The environment variables `RMFIT_SGE_QUEUE` and
`RMFIT_SGE_PE` override the file.  Example for an SGE site:

    {"mail": "you@example.org", "sge_queue": "long", "sge_pe": "smp", "setup": ["module load python/3.12 gsl"]}

for a workstation:

    {"scheduler": "shell", "setup": ["source ~/venvs/rmfit/bin/activate"]}

and for a Slurm site:

    {"scheduler": "slurm", "mail": "you@example.org", "slurm_partition": "compute", "slurm_time": "2-00:00:00",
     "slurm_extra": ["--mem=0"], "setup": ["module load python/3.12 gsl"]}

Re-running `init` keeps an existing configuration unless new values are given.

Policy meanings:

- `theta2_cap` is the Wigner-limit box for physical levels; `bg_theta2_cap` is the box for background poles (levels with
  Ex >= `bg_energy`, whose energies stay fixed).  3.0 is a starting guess, not a physics default: check where the
  baseline's background amplitudes already sit (`init`'s cap audit) before trusting it.  One alpha elastic fit needed ~100.
- Photon channels have no Wigner limit, so the theta2 caps never touch them.  `photon_width_max_ev` (1e4) and
  `bg_photon_width_max_ev` (1e5) cap |Gamma_gamma| of physical levels and background poles instead.  They are generous
  stops against runaways (a from-scratch background pole once polished to 1e12 eV), not Weisskopf limits: set them per
  reaction (0 = unbounded).
- `e_window` is the energy freedom of a level per polish (MeV); `norm_range` multiplies the nominal normalization.
- `narrow_window` (optional, e.g. `{"gamma_kev": 5, "factor": 3, "min_kev": 1}`) keeps a level whose free partial widths
  sum to less than `gamma_kev` within ±max(`min_kev`, `factor`·Γ) of its current energy instead of ±`e_window`, so a
  polish cannot carry a narrow level away from its peak (one narrow level drifted 180 keV in a single tightening polish).
- Shard count = cores.  Sharding is exact (chi-squared is additive over segments), so more shards only buy speed.

Cost check on a new model (`campaign.log` after `init`): one sharded objective should
be a few seconds and a Jacobian tens of seconds. If a segment with target-effect
convolution dominates, rmfit folds it itself (`model.ConvSpec`);
composite "sum" segments are allowed as convolution pieces.

### 2b. Spec, audit, decisions

Before `init`, write the campaign's **spec** (`spec.json`, template `rmfit/templates/SPEC_TEMPLATE.json`, schema in
`rmfit/spec.py`): every data set with its role (`anchor` = norm fixed at the nominal; `free` with `prior_sigma`; `excluded`
with the reason; `screen`), the fit window, the frozen block, the pairs that must exist, bounds / objective / search, what the
driver may do on its own (`pre_authorized`: `adopt_level` true / false / rules, `flip_adopt`, `dataset_off_adopt`, ...),
the stop rules and the **standing decisions** with their owner and date.  Then:

    python3 -m rmfit.audit <base>.azr --out AUDIT.md       # pre-flight, no engine session: zero free widths, 1 keV merges,
                                                           # nominal norms != 1, normError as a FRACTION, FLAG columns as
                                                           # errors, tiny errors, duplicates, spikes, thresholds, <targetInt>
                                                           # on test segments, legacy .sav names, output path, engine pin
    python3 -m rmfit.cli <dir> init --base ... --spec spec.json   # flags/policy/search/objective/nshards from the spec unless
                                                           # given; HIGH disagreements with the .azr refuse the init
    python3 -m rmfit.cli <dir> spec check                   # any time: SPEC_CHECK.md; spec render -> SPEC.md; audit -> AUDIT.md

A HIGH audit finding is fixed (or recorded as accepted in the spec's standing decisions) before the baseline.  During the
campaign, anything the spec does not answer goes to **`DECISIONS.md`** (`rmfit decisions` lists, `rmfit decide D-003
"adopt" --by evaluator` answers); a level add whose `accept` the spec does not pre-authorize is recorded as `accept-pending`
with a queued decision and is NOT adopted (automatic adoption applies only to campaigns without a spec).

The prompts for the three idle-time reviews of section 3b are generated:

    python3 -m rmfit.cli <dir> audit --kind all            # audits/{input,data,literature}/PROMPT.md + empty FINDINGS.md
    # hand each PROMPT.md to a background subagent (Agent tool, a cheaper model, at most three at once); it writes FINDINGS.md
    python3 -m rmfit.cli <dir> audit --findings            # unverified (candidate) rows of every FINDINGS.md
    python3 -m rmfit.cli <dir> audit --status data 1 confirmed --note "line 4007 checked"   # after YOU verified it against the files

The prompt carries the segment table, the active data files, the machine AUDIT.md findings, the spec's data-set roles and
standing decisions, the answered DECISIONS and the readme head, so the subagent does not re-report what is known; its table
has a `Status` column that only the driving session changes (candidate -> confirmed | rejected | decision).  Rules of 3b
still apply: findings are candidates until verified, nothing in a fit changes on the strength of a report, data treatment
is the evaluator's.  Without a campaign: `python3 -m rmfit.audit <model>.azr --kind data` writes `audits/` next to the model.

## 3. The loop

Every step is a job (`qsub -v RMFIT_CMD="..." run_crc_rmfit`, or `-hold_jid <id>` to queue
behind the previous one; on Slurm `sbatch --export=ALL,RMFIT_CMD="..." run_crc_rmfit` and
`--dependency=afterany:<id>` -- `--export` splits on commas, so JSON options go in a file as `@file`).  The examples
below use qsub; read sbatch on a Slurm site. Every step ends by rewriting `STATUS.md` and appending to the
readme; read both before the next step.

| step | command | purpose | when done |
|---|---|---|---|
| 1 | `polish --nfev 150 --tighten` | walk amplitudes beyond the theta^2 caps into the box (x10 stages), then converge | all theta^2 inside the caps; gain < 0.05% per stage |
| 2 | `round --rounds 2` (repeat) | interference-sign rounds: screen every canonical flip (singles, exposure pairs, same-channel pairs, optional channel patterns), polish the top-k sign-locked then released, greedy union, verify | improving-on-screen and accepted counts fall round over round; stop after two rounds with nothing accepted |
| 3 | `polish --nfev 150 --shifts-free` | warm-start continuation with energy shifts free | gain < 0.05% |
| 4 | `structure add --candidate i` (per entry of `candidates.json`) | forward addition: window J^pi scan at pinned energy, pinned then released full-model polish, control polish without the level, classify, adopt if accepted | all candidates classified |
| 5 | `structure dataset-off --file <substr>` for the worst chi2/N sets | deactivate + refit; reports the extra recovery of the *other* sets beyond the set's own chi2 (the data-set-off or "Botek" test) | a table for the evaluator |
| 6 | `round --rounds 2` after any structure change | signs may change with the structure | as in 2 |
| 6b | `crossover --donor <cid> [--donor-dir <other campaign>]` | basin crossover: apply another good fit's sign pattern per channel / J^pi group / all and polish every proposal; use it whenever the ledger (or a sibling campaign, e.g. a recovery benchmark) holds a second verified basin | every donor pattern polished; verified union |
| 7 | `review [--grid-check]`, `export`, then GUI review | bake the incumbent to `<reaction>.azr` + `.sav`, verified in a fresh process | CLI mode 1 with the `.sav` as external file reproduces the totals |
| 8 | `uq --uq uq.json --nprocs 4 --threads 6` (queue; section 3d) | uncertainty quantification of the incumbent: export, saturated poles held, Laplace covariance, NUTS chains, R-hat / ESS / divergence gates, linear vs sampled bands -> `UQ.md` | gates PASS; the band table says where the linear band is adequate |

`candidates.json` (write it before step 4): `{"candidates": [{"energy": Ex_MeV, "jpi":
"5/2+" or null to scan, "energy_fixed": true, "channels": [{"pair": n, "L": l, "S": s,
"gamma": width_eV_with_sign}, ...], "source": "...", "note": "..."}]}`. Sources: ENSDF/TUNL
(`pyazr.nds.fetch_levels`), the reaction's literature tables, and residual runs
(`diagnose.residual_runs`) in the fit's worst segments. Omitted channels are seeded at a
small fraction of the Wigner limit.

**Default J^pi of `levelscan` / `structure add`.**  With no `--jpi` (or `"jpi": null`) these use the model's own
active J^pi groups (`structure.default_jpis`), and a requested J^pi that no particle pair can couple to raises an error
(`structure.resolve_jpis`; e.g. 1/2+ on an integer-J compound nucleus).  The log names the list and warns about any NEW
group: its channels are partial waves the model lacks, and their hard-sphere terms change the elastic channels even at
zero width.  In every case, check the template screen's log line `template screen at Ex ...: incumbent <obj>`: it must
equal the campaign's verified incumbent.

Residual-driven level search (`rmfit levelscan [--ex E ...] [--top K] [--gamma keV]
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
the single-level signature is gone (the compensation limit) and only the refit
of the top candidates can recover it. Template amplitudes are the engine's own conversion
of a tiny reference width scaled by sqrt(Gamma); candidate energies are rounded to 1 eV
before writing (see the AZURE2 skill's add_level gotcha).

Every level test and relocation writes a review figure to `figures/` (data with the
incumbent and the candidate calculation around the energy: excitation functions per
segment, angular distributions for single-energy segments) -- the evaluator's visual check.

Relocation of an existing level (`rmfit relocate --level 9/2+#1 --ex 8.4654 [--window 5]`):
energy kick to `--ex`, then a full re-polish with that level's energy confined to
+-`--window` keV, verified and adopted if it beats the incumbent by tau. Use it for a
narrow level that has drifted from its peak (see `narrow_window` above for keeping it
there afterwards). 13N test: a level shifted by 60 keV comes back and recovers the full
model's objective.

Other commands: `touch <dir>/STOP` ends a `round --rounds N` loop after the current
round; `verify <cid>` re-scores a candidate; `log "..."` appends to the readme;
`report` writes `REPORT.md` (per-dataset chi2/N, theta^2 table, structure decisions).

## 3a. The closed loop: `rmfit drive`

`rmfit <dir> drive` runs sections 2-4 as a resumable state machine instead of a session choosing every step.  It needs
`<dir>/spec.json` (section 2b): the spec's `budget`, `stop_rules` and `pre_authorized` are the thresholds, and anything the
spec does not answer becomes a `DECISIONS.md` entry.

    cp <base>.azr <base>.sav <dir>/; cp -r data <dir>/data; mkdir <dir>/output; write <dir>/spec.json
    python3 -m rmfit.cli <dir> drive --base <dir>/<base>.azr --sav <dir>/<base>.sav --dry-run   # records the base, prints the plan,
                                                                                               # writes run_crc_drive
    qsub <dir>/run_crc_drive                     # one leg per job: states run until DONE / WAITING / STOP or DRIVE_MAX_MIN (600),
                                                 # then the job resubmits itself (CONTINUE)
    python3 -m rmfit.cli <dir> drive --dry-run   # any time: where it is, what comes next, what it waits for (also DRIVE.md)
    python3 -m rmfit.cli <dir> decide D-003 "adopt" --by evaluator; qsub <dir>/run_crc_drive     # after a WAITING leg
    touch <dir>/STOP_DRIVE                       # stop cleanly after the current state

On Slurm the leg asks for `--time` = DRIVE_MAX_MIN + 30 min (the last state and the resubmission); when the site's
`slurm_time` is shorter, DRIVE_MAX_MIN is shortened to fit and the script says so.  The resubmission runs `sbatch` on the
compute node; where a site forbids that, resubmit from a login node.

**Build phase** (section 2's prerequisites and "First fits"):
`PIN_CHECK AUDIT REPRODUCE SPEC_CHECK SCALE_ANCHOR SEED_CHANNELS PAIR_COMPLETENESS LADDER CANDIDATES INIT RECALL`

- `AUDIT` runs the machine checks.  On a new campaign it then writes the three subagent audit prompts (section 3b) and waits
  for their verified `FINDINGS.md` (decision `audit-subagents`: 'done' or 'skip'; `pre_authorized.skip_subagent_audits`
  for fixtures).  Each row left as 'decision' becomes one wait decision.  Why: a from-scratch model once fitted a
  lab-frame error for an hour before the data audit found it.
- `REPRODUCE`: the archived `chiSquared.out` must reproduce.
- The spec's data-set roles are applied to a working copy `<stem>_drive.azr`; free widths at exactly zero are seeded at
  `model.seed_frac` (0.01) gamma_W; ENSDF states of the residual nuclei that open below the window top without a pair are
  listed; `candidates.json` is built from ENSDF, broad levels first, when it is absent.  Every candidate needs a
  `hypothesis`.
- `RECALL` (`structure.recall_levels`) removes the broadest known levels one at a time and asks the level-add procedure to
  find them again with no J^pi given.  Below `recall_min_pass`, every add of the run is queued as a decision instead of
  adopted.

**Loop** (section 3 steps 1-7 with section 4's stopping rules):
`TIGHTEN SIGN POPULATION SHIFTS_FREE LEVELSCAN STRUCTURE_ADD SIGN_POST DATASET_OFF REVIEW EXPORT CYCLE`

- **Polish** rounds run until the gain is below `polish_gain_rel_min`.  Size `budget.polish_nfev` to the model: 150 is far
  too few for 200 free parameters.
- **Sign rounds** start only once the polish has converged (`sign_after_converged_polish`) and run until
  `sign_rounds_without_accept`.  On a moving incumbent they are deferred to the next cycle (in one benchmark 18 rounds on an
  unconverged incumbent accepted nothing).  They are deferred at most once and never in the last cycle; after that they run
  on the unconverged incumbent and the readme says the polish budget is too small (one from-scratch model with 162 free
  parameters never converged in 6 x 300 evaluations).
- **POPULATION** is the basin-crossing stage (`Campaign.run_population`, CLI `population`).  It makes
  `budget.restarts_per_cycle` sign-perturbed restarts of the incumbent (6-12 random high-exposure flips, at most 3 per
  J^pi, 20 % jitter) and polishes each sign-locked, then released, to convergence in rounds of `restart_nfev`
  (0 = `polish_nfev`), together with a control.  Then `crossover_passes` passes of `crossover` moves run from the best
  member with every other member as donor; each pattern is polished and the winner polished to convergence.  Members are
  judged on polished objectives only.  An improvement is a new basin, and SIGN_POST finishes its last signs.  A round
  gaining less than `population_gain_rel_min` switches the stage off until the next structure change
  (`population_rounds_max` caps it).
  Why: in a development benchmark the sign rounds stopped in a basin 13.5 % above the reference, at sign distance 17.  The
  reference's pattern, applied to the incumbent, started at three times its objective and polished down to the reference,
  while no single, pair, channel or group flip gained more than 2.6 even with a 600-evaluation polish.  Local proposals
  do not cross basins; a population of polished basins and crossover between them does (six restarts and one crossover
  pass recovered the reference).
- **Structure search** follows the spec's `structure` section.  Candidates go broad first: a known width below
  `broad_first_kev` (10) waits until the broad queue and the residual scan are exhausted, then the narrow-alignment phase
  (`narrow_phase`) takes them.  When the broad queue is empty, LEVELSCAN (`levelscan_per_cycle`, `levelscan_top`,
  `levelscan_gamma_kev`; null = the run's width) screens the residual-run clusters with one level of every J^pi and appends
  the best (J^pi, pattern) to `candidates.json`, with the screen as hypothesis and the next two J^pi as alternatives.
  Every add polishes the neighbourhood first (`local_release_mev`, `local_nfev`).  Level adds continue until
  `level_adds_marginal_in_a_row`; adoption follows `pre_authorized.adopt_level`.
- **DATASET_OFF** runs the data-set-off test (the "Botek test" in the logs) on the worst chi2/N sets.  The result is
  recorded; a data set is never excluded by the drive.
- **REVIEW**: HIGH flags wait for the evaluator.
- **CYCLE**: another pass while a cycle gains `cycle_gain_rel_min` and `cycles_max` allows; then UQ if `spec.uq.enabled`
  (section 3d).

`DRIVE.md` is rewritten after every state; `drive_state.json` is the checkpoint (`--state S` jumps, `--reset` starts over
and keeps the old file).  Engine states refuse the login node (`--local` is for fixtures).

What the supervising session does: read `DRIVE.md`, `STATUS.md`, `DECISIONS.md` and `REVIEW.md`; relay the pending
decisions to the evaluator with the evidence files; `decide`; resubmit.  It does not pick a step the spec answers, and it
does not answer a data-treatment decision itself (section 5).

Limits: no channel lines or pairs are added automatically (PAIR_COMPLETENESS and SEED_CHANNELS report what is missing);
LADDER needs `model.ladder_Ex`; multi-compound-nucleus models are outside the sign rounds and levelscan.  The drive has run
on the 13N fixture (`tests/test_drive_13n.py`, ~70 min on 3 slots; rerun it after every change to the state machine), on
the bench cases of section 7 and on from-scratch models, but not yet on a full evaluation.

## 3b. While jobs run: audit and read (idle time is work time)

Minimization is only part of the work.  A polish or a sign round takes hours and needs no attention between milestones.
Whenever a job is running and nothing is actionable, do NOT just wait for the next log line: run the three reviews below.
They are read-only and cost no node, and on a large evaluation every large gain of the merged fit came from this kind of
work rather than from more minimization (a free unpenalized normalization on the neutron total data, a data file with
statistical-only errors, a width window that swallowed a neighbouring level, a wrong elastic / inelastic partition that
only independent data exposed).

**When.**  At the first idle window of every campaign, and again whenever data are added, a
structure change is adopted, or a fit stalls.  Mandatory, not optional, when: a fitted
normalization leaves 0.7-1.4; a frozen or pinned part of the model carries a large share of the
objective; two fits of the same data disagree in absolute scale; a level sits at a bound or at a
threshold guard; a new data set is about to enter the fit.

**How.**  One background subagent per review (Agent tool, `model` set to a cheaper model than the
session's so the reviews do not consume the main budget), at most three at a time, each with a
self-contained prompt (paths, file-format facts, what is already known and must not be
re-reported), each writing ONE report into `<reaction>/<M-D-YY>_audit/` next to a `readme` that
quotes why it was run.  The reviews never edit a fit directory, never submit jobs.  A report's
findings are CANDIDATES: verify each one against the files yourself (numbers, line, segment key)
before telling the evaluator, and change nothing in a fit on the strength of a report alone --
data treatment stays the evaluator's decision (section 5).  Start from the prompts that `rmfit audit --kind all`
writes (section 2b).

**1. Input-file audit** (the `.azr` of the current seed):
- particle pairs: mass, charge, spin, parity, excitation and separation energy, radius identical
  on every level line and right against the nuclear data; thresholds reproduce to ~1 keV;
- channel sets: every (pair, L, S) allowed by angular momentum and parity, none missing at low L,
  the same set on every level of a J^pi group;
- levels: same-J^pi pairs within 1.5 keV; free amplitudes that are exactly 0 (they never move);
  levels within 3 keV of a threshold; fixed / free flags against the campaign's rules;
- segments: file exists; entrance / exit keys match what the file holds; energy window against
  the data's range and the exit threshold; isDiff against the angle column and its frame; "sum"
  tails complete for every channel open in the window;
- segment tails: list every segment with tokens after the file name other than `0 0`
  (azure2-eval SKILL, ".azr file anatomy", has the field meanings and an awk one-liner) and
  record which are composite SUM / RATIO segments and which are UPOS (secondary-gamma)
  segments. A UPOS segment's observable is the secondary gamma, not the primary cross section,
  and no extrapolation reproduces it, so curves, scans and post-hoc tests on it must run in data
  mode on a pseudo-data grid (once 39 UPOS and 138 RATIO segments in one model);
- normalizations: every `varyNorm 1` with `normError 0` is a FREE UNPENALIZED norm -- list them
  with their fitted values and ask whether each is intended; list fitted norms outside 0.7-1.4;
  check every norm / shift setting against the evaluator's recorded decisions;
- `<targetInt>` keys bind to the intended segments; in sharded runs a convolved segment must be angle-integrated or at
  ONE fixed angle (minA == maxA); an angular distribution with a
  `<targetInt>` line still needs a single session.

**2. Data audit** (every file of an active segment):
- integrity: non-numeric rows, zero or negative values or errors, duplicated rows;
- units and scale: ratios between overlapping sets of the same reaction, band by band; an
  absolute anchor (a total cross section, a unitarity constraint) compared with EVERY independent
  measurement of the same quantity, not only the one in the fit;
- frames: lab energy of the right projectile (test against thresholds), lab vs c.m. angles and
  cross sections against the segment type; a converted file is cross-checked by evaluating it
  frozen in both frames with identical errors (the chi-squares must agree exactly);
- uncertainties: statistical-only where the source gives a systematic one; an advertised error
  floor that does not hold; implausibly small or constant errors;
- double counting: the same measurement under two names or in two overlapping windows;
- energy scale: for a set with narrow structure, run a post-hoc shift scan on the frozen model
  curve BEFORE any refit (none / constant / a+b*sqrt(E) / a+c*E, norm re-optimized; azure2-eval
  SKILL, "Fitting a data set's energy scale"); a shift freed afterwards must be SEEDED at the scan
  optimum -- MIGRAD from zero stalls in a local minimum next to zero (-193
  from zero vs -795 seeded in one case);
- against the source: re-fetch the EXFOR entry (`pyazr.nds`, AZURE2's `nds-explorer` skill) and compare values, units and error columns row by row.

**3. Literature review**: earlier evaluations of the same compound system (what data, channels,
radii, energy limits, how they treated normalizations and backgrounds), measurements of channels
the fit does not yet constrain, and method papers on the specific difficulty the campaign is
stuck on.  Every claim carries a citation; a paper that could not be read is listed as such;
nothing is quoted from memory.  The output is a short list of actionable ideas, each with a cheap
test, not a survey.

**Cheap in-session checks that belong to the same habit** (minutes each, no node):
- decompose a summed observable into its partial cross sections by extrapolation and compare each
  partial with independent data -- a total can be right with the wrong partition;
- the clip check before every qsub: sanity + bounds + clipped objective of the seed under the
  policy the job will use; a start-up clip that costs more than a fraction of a
  per cent means a window or a cap is wrong;
- a coverage plot of what EXFOR holds against what the fit uses.

## 3c. A structure you can see but no search finds: the width-floor recipe

The recipe comes from a broad resonance in the backward angles of one excitation function that every tool had declared
null: levelscan at three energies (all J^pi, 50-400 keV), sign rounds over the groups there, a seed / re-partition test,
and a first local-release level add in which every J^pi candidate collapsed to eV widths and came out "marginal".  The
level finally accepted (a width of about 200 keV in the exit channel, a large gain over its control) had been within reach
of the very first test.  It was lost because a released polish lets a broad seed shrink to nothing while the neighbours
re-absorb the misfit, and a collapsed candidate is then scored as "nothing there".  The recipe that found it:

1. **Localise and characterise the deficit before fitting anything.**  Map chi2 per 0.2 MeV of Ex by data set (from the
   incumbent's per-point residuals), then per 10 keV with the MEAN PULL per band and per angle group.  A broad deficit
   shows as the model below the data by 2-3 sigma over >= 50 keV with the same pull sign, at every angle or in one angle
   range.  A narrow spike (two 10 keV bands, the pull changing sign) is a lineshape or resolution question, not a level
   search.  Write the anatomy into the readme first.
2. **Width-floor level add, every J^pi, with the local release and a same-budget control.**  `structure add --local-release
   0.5 --local-nfev 80 --pinned-nfev 60 --nfev 150` with the policy's `total_width_floor_windows: [[E-0.0005, E+0.0005,
   100000]]` (a floor on the TOTAL width over the level's open channels; the older per-channel `width_floor_windows` forces
   the width into one named channel) and `setenv RMFIT_SEED_OMITTED 1`, so that every open channel of the group exists and
   can take width.  Seed the lowest-L channels broad (e.g. 150 keV in the exit channel, 30 keV in the entrance channel).
   Use one campaign directory per J^pi (an accepted level is adopted into its base, so groups must not share one), three
   groups per node job.  The floor is one-sided on the sign the amplitude starts with and never above 0.95 of the cap; it
   is skipped for an amplitude that is exactly 0 (the control's "level off" state).  A control whose start differs from
   the incumbent is broken (for example because it carries the floored level): check it before reading any delta.  The
   verdict is the delta against the control, with tau from k; an empty "at_bound" list at the end means the data, not the
   floor, chose the width.
3. **Cross-checks that make the verdict believable.**
   (a) Run the same recipe on the other plausible hypotheses, e.g. widen an existing neighbour of the same region (a floor
   on its width) instead of adding a level.
   (b) A data-only window multi-start on the data set that shows the resonance (±0.5 MeV window, norms fixed, random
   sign patterns, every J^pi with a level, 250 evaluations each, ~1 h on one node).  It is a SCREEN that must agree with the
   full-data result on J^pi and width, never a start vector: vectors transplanted from it into the full model came back far
   worse.  When the screen ties two J^pi, the full data break the tie.
   (c) Attribute the gain per data set and per angle: it should sit where the deficit was.
4. **Full free refit of the accepted fit** (every parameter, no floor, 400 evaluations): the level must hold on its own,
   with nothing else moving much.  That export is then the incumbent.
5. **Propagate.**  Any model that freezes this block of levels gets a NEW seed with the refreshed levels (copy energies and
   channel widths by (J^pi, energy order), add the new levels as fixed, keep the fixed flags, write a partial .sav, check
   the frozen objective, recentre the policy windows) and a refit.  A splice cost is normal, because the free levels above
   the cut were tuned to the old block.
6. **Then the next region**, chosen from the band map, with the SAME recipe.  Expect most J^pi to be rejected quickly (the
   high-spin ones by orders of magnitude: a 100 keV floor is above their Wigner limit).  A candidate that improves a data
   set other than the one you aimed at still counts: the goal is to fit all the data.

Cost: eleven J^pi = 4 node jobs of 6-8 h; the anatomy and the window screen are login-node work.

## 3d. Uncertainty quantification: `rmfit uq`

`rmfit <dir> uq [--uq uq.json|JSON] [--nprocs N] [--threads T]` runs `rmfit/uq.py`, as a job (`qsub -v RMFIT_CMD="uq --uq
uq.json --nprocs 4 --threads 6" run_crc_rmfit`) or as the drive's last state when `spec.uq.enabled` is set:

1. **EXPORT** the incumbent to `uq/<name>.azr` + `.sav`.
2. **SETUP**: Park's J of every level from a Park session on the same file.  Levels with J < `j_saturated` (0.01) are
   saturated -- their widths sit at the channel-radius ceiling, the data do not set them, and the Brune Hessian is
   singular along them -- and are HELD.  Laplace covariance of the sampled block.
3. **SAMPLE**: `pyazr.sampling` NUTS, `chains` x (`warmup` + `draws`), dense metric from the Laplace covariance,
   `target_accept` 0.9.
4. **GATES**: split R-hat <= 1.01, minimum bulk ESS >= 400, divergences <= 1 %, E-BFMI >= 0.3.  A failure re-runs once
   with 2x the draws (subject to the guards below); after that it is the decision `uq-gates`.
5. **BANDS** on the spec's `uq.bands` grids or the model's own `<segmentsTest>` energy scans: the LINEAR band from the
   Laplace covariance on the free R-matrix block, the same band times sqrt(chi2/nu) printed beside it and NOT adopted, and
   the SAMPLED band from `band_draws` posterior draws (percentiles 2.5/16/50/84/97.5).
6. **REPORT**: `UQ.md`.

Norms are sampled (the full posterior) unless `sample_norms` is false; `uq.priors` {engine name: [mean, sigma]} adds
Gaussian priors beyond the objective's own penalties.

**Bounds.**  The chains run inside the campaign's theta^2 box: every sampled particle amplitude stays at
|gamma| <= sqrt(cap) gamma_W with the policy's cap for its level (`bounds: "policy"`, the default).  NUTS runs on
x = c + h tanh(u) with the Jacobian; draws are stored in x, and the metric in u comes from the Hessian of -log p(u) at the
mode (the Laplace covariance divided by (dx/du)^2 froze the chains of a validation run).  The box is needed because the
data alone need not bound a background pole: in one validation without it, an amplitude left its Wigner limit in two
thirds of the draws and emcee's acceptance fell to 0.11.  `bounds: "none"` lifts the box; `theta2_cap` gives one cap to a
bare .azr + .sav that has no policy.  The box is the campaign's (the polish ran in it); a different prior on an amplitude
is the evaluator's.

**Re-run guards.**  There is no automatic re-run when split R-hat > `rerun_rhat_max` (1.5) -- chains in different regions
do not merge with longer runs -- or when the re-run would exceed `rerun_max_hours` (12, estimated as extend_factor x the
attempt's wall time).  The stage then reports FAIL with the reason and the drive asks `uq-gates`.  Each chain prints its
rate and the hours left after 100 draws.  (The guards came from a run with R-hat 2 whose second attempt would have taken
~60 h.)

**Reading the bands.**  A half-width / linear sigma near 1 and a small median offset mean the linear band is adequate
there.  A ratio well above 1, or an offset of a sigma, means the posterior is skewed along a direction the Jacobian misses
(a saturating width, a sign-ambiguous amplitude).  Everything is conditional on the model and the spec's data decisions.

Cost: a gradient is ~2 objective evaluations; a model with 195 free parameters needs ~25-50 gradients per draw.
Pitfalls: a script that spawns the chains must guard `if __name__ == "__main__"` (spawn re-imports it; without the guard it
recurses).  `pyazr.bands.extrapolation_bands` opens its session with pyazr's DEFAULT flags, so `uq` trims the model and
opens the extrapolation session with the campaign's flags itself.  Validation: section 7.

## 4. Judging results

- Acceptance thresholds scale with s^2 = objective / N because errors are
  underestimated (s^2 ~ 15-20 on a large model). A flip needs a verified gain
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
  (typically 16 -> 4 -> 3 -> 1). Two rounds without an acceptance = converged for
  the present structure.
- `sanity` flags theta^2 over the cap, |E| > 60 MeV, norms outside (0, 100), a
  non-finite Brune transform, and parameters at a bound (listed in the readme entry).
  Any |value| > 1e6 or |E| > 1000 MeV means catastrophic cancellation in
  total-cross-section chi2, i.e. a fake improvement.
- Read `REPORT.md`'s per-dataset table after every structure step: a data set whose
  chi2/N is far above the rest, or whose removal recovers more than its own chi2 in
  the others, is an evaluation question (normalization, energy calibration, resolution),
  not a fitting question.
- **Run `rmfit <dir> review` before judging a fit by eye**: `REVIEW.md` + `figures/review_*.png`
  hold the checks the evaluator makes in the GUI -- per-point summary per data set (fraction within 1 / 2 sigma, median |r|, how
  many points carry half of the chi2 and where), the top chi2 carriers, every norm against its segment's own optimum n*
  (a free norm on the other side of 1 from n* was pulled by a group / covariance / prior), n* trends with angle or energy
  inside one experiment (a shape problem dressed as a scale), residual runs with the level inside or missing,
  residuals that vanish when blurred over 50 keV (the fit follows the average), MeV-wide levels spanning visible
  unfitted peaks, a zigzag of the calculated curve (bumpiness), levels narrower than 2 % of the data spacing under a
  target integration (`--grid-check` re-scores on a converged grid), levels seen by one data set only.  Engine-free on
  any finished CLI run: `python3 -m rmfit.review <file.azr> [--output DIR]`.  Flags are candidates, not verdicts: look
  at the figure before writing a readme line, and leave data treatment to the evaluator.
- **Check every fitted norm against the set's fit/data ratio** (the norm multiplies the data and its errors, so at a
  sensible minimum fit/data ≈ norm). A norm that moved the other way (in one case a norm of 1.32
  with the fit at 0.64 of those data) usually points to a data-scaled `vector` in that set's covariance acting as a
  hidden normalization with Peelle's bias, not to an rmfit bug. Scan that one norm with the plain and the covariance
  χ²; the fix is to write analysis variants as `relvector` terms (README, "Data covariances and normalization groups"). Also distrust norms
  sitting on the `norm_range` bounds while χ²/N ≫ 1: the priors are negligible there and the norms absorb shape misfit
  (see "First fits" in section 2).

## 5. What the human decides

rmfit hands these to the evaluator rather than deciding them:

- normalization scale and systematic error of a new data set (and any rows dropped);
- an angle- or energy-dependent trend in fitted normalizations (the norms absorb it
  silently);
- whether a data set that distorts the rest stays in the fit;
- the energy range of the model (data above the level scheme's range are capped, not
  fitted);
- accepting a level whose gain is statistically clear but whose J^pi or energy
  contradicts the literature;
- error inflation (the Birge factor sqrt(chi2/nu) that `UQ.md` prints beside the plain band), priors on R-matrix
  parameters, and what to do with chains that fail their gates after the automatic rerun (`uq-gates`).

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
- An oracle test applies a known target's sign pattern per channel / per group / all at
  once and polishes: it separates "finding the pattern" from "reaching the basin from the pattern".
- `diagnose.dataset_table`, `theta2_table`, `residual_runs`, `compare` feed `report`.
- **Bench suite**: `python3 -m rmfit.bench list|score|run` on `rmfit/bench/cases/<name>.json`.
  `score` applies the evaluator's recovery yardstick from the ledgers alone (data chi2 <= reference, objective within 0.5 %, zero wait
  decisions, no HIGH/MEDIUM review flag the reference lacks; 'exact' 1e-4 as a secondary column) and appends to `bench/results.tsv`;
  `run` writes the SGE job (reference init, `perturb`, spec, `drive`).  Supersedes "success = within 1e-4" above for drive runs: a
  benchmark post-mortem showed 1e-4 selects one of two near-degenerate background configurations.

## 7. What the benchmarks show, and the known limits

The recovery yardstick (section 6): data chi2 at or below the reference's, objective within 0.5 %, no review flag the
reference lacks, zero wait decisions.  The tier-0 case shipped with rmfit (`bench/cases/13N_tier0.json`, AZURE2's
`tests/13N` fixture: four high-exposure sign flips) passes exactly in about 25 minutes on 3 slots.  On larger references
(hours to a day on 8 slots) the development benchmarks showed:

- **Sign rounds alone do not cross basins; the population stage does.**  Single, pair, channel and group flips can leave a
  run in a basin a few per cent above the reference; sign-perturbed restarts polished to convergence plus crossover between
  them recovered the reference's interference pattern, including branches an evaluator found by hand.
- **The yardstick, not a 1e-4 match**: near-degenerate background configurations can differ by a fraction of a per cent with
  the data chi-squared on the right side; "exact" is a secondary column.
- **A from-scratch model needs every stop an inherited model has**: photon widths without a cap run away; seeds for free
  zeros, the anchor, the caps and the frozen block are spec items, not defaults to trust.
- **UQ can fail its gates on a posterior that is genuinely hard** (a rarely visited tail, an unbounded background pole);
  the stage then reports FAIL and the evaluator decides, which is the intended behaviour.

Still true (limits):
- No channel lines or pairs are added automatically; LADDER is unexercised; multi-compound-nucleus models are outside the sign
  rounds / levelscan; density / Doppler / thickness / screening kinds are single-session only (no shards).
- The structure search has been validated on recall only; a level no reference model had has not been found by it yet.
- Energy shifts are finite-differenced.  Data covariances, normgroups, shiftgroups: `covariance.py`, README "Data
  covariances and normalization groups" (relvector terms are scaled by one reference model, `cov_reference.npz`).
- Every campaign is pinned to one engine (`campaign.json["engine"]`, `EnginePinMismatch` unless
  `RMFIT_ALLOW_ENGINE_CHANGE=1`); pin a frozen copy of the engine for long campaigns, not a checkout that changes.

## 8. Gotchas

- Never open two pyazr sessions in one process (the parameter tables desynchronize); rmfit runs one process per session,
  and `bake()` verifies in a fresh process.  Do not use `save_fit` for campaign files.  `bake()` skips inactive levels.
- `bake()` also writes the FIXED parameters' current physical values: a fixed width is fixed in rwa space, so its observed
  value moves with the level's free amplitudes.  To freeze a level's observed widths, fix its closed-channel amplitudes
  too.  Test: `tests/test_bake_fixed_13n.py`.
- `total 1` (the optical-theorem neutron total, AZURE2 commit d37d4e6) is in rmfit's keyword table and in bake's keyword
  guard; AZURE2 builds before d37d4e6 ignore it silently and compute the elastic cross section only.
- Delete `output/intEC.*` whenever the segments change (the stale-cache rule of the AZURE2 skill).
- New levels are appended at the end of the file and levels are removed by deactivation, so every existing parameter key
  stays valid across structures.
- A released polish lets amplitudes cross zero: sign patterns are not preserved by polishing, and a "basin" is only
  defined up to the amplitudes that are near zero.
- Frozen screens (flip, no refit) are cheap but favour flips of small amplitudes; the GN screen (one damped Gauss-Newton
  step on the group's columns) is the ranking that matters.  Coordinated multi-level flips need several backtracking GN
  iterations (`search.gn_iters_multi`) and their own polish quota (`search.multi_top`).
- Job mail: `-m abe` (start, end, abort), so that a silent crash is noticed.
- **`pkill -f` / `pgrep -f` patterns are regexes, and every reaction name has a `+`.**  Find the pid with a fixed-string
  match (`ps -eo pid,cmd | grep -F "<name>.azr"`), kill BY PID, confirm with `kill -0 <pid>` failing, and only then remove
  the log.  Launch background programs with `< /dev/null` and a size cap.  **When a batch of jobs dies at once, check the
  group's disk quota (`df -h <group filesystem>`) first.**
- **`init --sav` and FIXED parameters.**  AZURE2's parameter-file reader sets every parameter it names, fixed level
  energies and widths included; rmfit seeds only the free ones and evaluates the .azr's fixed values.  Both readings
  occur in practice: a CLI fit's `param.sav` holds the fixed values the fit ran with, while a script-built base means the
  .azr's values and may come with a stale .sav.  So `init` logs every disagreement with the objective both ways and keeps
  the .azr's values by default; `--sav-fixed apply` (for a CLI fit) writes `<stem>_sav.azr` (+ .sav) with the file applied
  and makes it the base (`Campaign.apply_sav_fixed`, test `tests/test_init_sav_fixed_13n.py`).  Keep the archived
  `chiSquared.out` in the base's output directory so that init's reproduce check runs (bench: `chisq_out`).
- **REVIEW's HIGH flags and the drive.**  The init baseline is reviewed once (`REVIEW_baseline.md`); a HIGH flag the base
  model already carried is reported, not waited on, and a flag answered 'continue' is not asked again (`drive_state`
  `review_continued`).  A flag that the evaluator's own model carries (a broad pole the evaluator kept, with "no extra
  levels" a standing decision) still queues `review-high` when the start did not have it: the spec cannot tell the drive
  that a standing decision answers it, so answer it with the evidence.
- **Close same-J^pi pairs: a report, not a constraint.**  `review` flags two levels of one J^pi whose gap is below a
  fraction of the NARROWER total width (MEDIUM `level-close-width`; default 10 %, `review.CLOSE_GAP_FRACTION`;
  `rmfit <dir> review --close-fraction f`, `python -m rmfit.review <azr> --close-fraction f`, spec budget
  `"close_fraction"`; 0 = off).  Nothing in rmfit enforces a minimum spacing: a centroid-plus-bounded-gap
  reparametrization would be a model constraint whose scale is the evaluator's to choose, and a MEDIUM flag does not make
  the drive wait.  The flag is a question, never a verdict: it has marked a same-J^pi pair 1 keV apart, with widths of 240
  and 280 keV, that merge tests found NOT redundant (one data set needs both).  The test of a flagged pair has two halves,
  both from the same start: a merged refit, AND a refit of the pair started from the merged solution.  A merged model that
  only wins against a two-level fit stuck in a worse minimum proves nothing.
- **Park-mode campaigns (`"use_park": true`).**  (1) The J > 0 wall is carried as residual rows (`engine.ParkWall`, test
  `tests/test_park_wall_13n.py`).  (2) `--sav` seeds are applied by NAME with no basis conversion: convert a Brune `.sav`
  through a CLI mode-1 run in Park mode first.  (3) `theta2_cap` bounds the observed theta^2 in Park mode and the formal
  one in Brune mode.
