# rmfit — agent-operable global R-matrix fitting on AZURE2

`rmfit` drives the in-process AZURE2 engine (`pyazr`) to do what an evaluator does by
hand, on a schedule a script can follow and an LLM agent can resume:

- **sharded evaluation** — the data segments are split across N single-thread engine
  sessions (chi-squared is additive over segments), so one objective or Jacobian of a
  35,000-point model costs seconds instead of tens of seconds;
- **bounded local polishes** in reduced-width-amplitude space with the analytic
  Jacobian (Wigner-limit bounds, energy windows, sign locks);
- **sign rounds** — canonical single/pair interference-sign flips per J-group, screened
  with one Gauss-Newton step, polished with the flipped sign locked, then merged one at a
  time and re-verified on a fresh full-model session;
- **structure rounds** — ATARI-style backward elimination and forward addition from a
  candidate bank with J^pi scans at pinned energies, holdout checks and outcome
  classification (accept / collapse / blowup / live-harmful / relocated / reverts);
- **energy-ladder continuation** for from-scratch fits;
- a **ledger** (SQLite + npz) of every candidate with provenance, and a `STATUS.md`
  a fresh session can read to continue a campaign.

Process guide (the Claude skill `rmfit-campaign`): `.claude/skills/rmfit-campaign/SKILL.md` in this repository,
linked as `~/.claude/skills/rmfit-campaign`.  It lived in the AZURE2 repository until 2026-11-07.

**This copy runs on pyazr 3.0** (`Session`, `Objective`, `Evaluation`, `Parameters`, `Project`; project format 3
only).  It is rmfit of AZURE2 `dev` (527f63de, 2026-11-15) re-written where it touched the engine or the project
file; the search, the ledger, the driver's logic and the criteria are unchanged.  What differs:

- every project is format 3 (`AZURE2 --convert OLD NEW` first); rmfit reads it through `azr3.py` (the data sets
  from the engine's description and the records' fields, levels and pairs from `pyazr.Project`) and never parses
  a line of it;
- the rows are the engine's: data rows from `Objective.residuals`, the priors of normalizations, shifts and the
  other nuisances from the session's parameter rows, and the engine's model rows (model data such as a measured
  scattering length, a prior covariance, Park's wall with them, decay priors, abundance sums) as one block
  (`engine.ModelRows`); a single session's objective is `Session.objective()` exactly;
- shards are record edits (`azr3.write_variant`); a data set the engine folds (it names a target, a beam or a
  detector) or one of a shared energy scale is never split, and there are no convolution pieces;
- `bake` writes its files with `Session.save_fit` and keeps rmfit's own fresh-session check;
- `audit` takes masses, separation energies and file problems from the engine's validator;
- `uq` samples with `sampling.py` (dev's NUTS, kept here because pyazr 3.0 has no sampling module);
- `mixture.py` is gone (superseded by several compound nuclei in one project).

The notes below this table are the history of rmfit on the old project format; where they name
`<segmentsData>`, `<targetInt>`, `AzrModel` or a segment-line keyword, the format-3 counterpart is the record
field of the same meaning (`norm`, `shift`, `covariance`, `target`, `beam`, `detector`, a `normalization` or
`energy_scale` record).

Layout (design rationale and roadmap: `DESIGN.md` in this directory; the plan file the first campaign readme cites is gone):

| module | role |
|---|---|
| `keys.py` | structured parameter keys (energy, width, norm, shift, shift_sqrt, shift_fn, abundance, density, doppler, thickness, screening) and index maps between model variants; canonical sign patterns; `sign_distance` (minimal flips between two patterns modulo the gauge symmetries) |
| `azr3.py` | the project as rmfit reads and edits it on pyazr 3.0: `SegmentTable` (data sets), `write_variant`, `Model` (level and data-set edits through `pyazr.Project`) |
| `model.py` | model variants written as new `.azr` files (shards, single-session copies) |
| `sampling.py` | NUTS, R-hat and ESS (dev's `pyazr/sampling.py`), `AzurePosterior` on `Session.objective()` |
| `engine.py` | `SingleSession` (authoritative) and `Sharded` evaluators; `ObjectiveSpec` (plain chi2, reduced-chi2 per data set, Sivia loss) |
| `local.py` | bounded trust-region polishes, staged fits, sanity checks |
| `ledger.py` | candidate store, rounds, status |
| `moves.py` | physics-informed proposals |
| `search.py` | sign rounds, basin-crossover rounds, the shared screen/polish/union machinery |
| `structure.py` | level add / remove / dataset-off hypotheses |
| `windows.py` | energy-ladder driver |
| `levelscan.py` | residual-driven level search: candidate energies from residual runs, template screen of every J^pi x width partition x sign pattern at each (`rmfit levelscan`); relocation of an existing level with a narrow energy window (`rmfit relocate`) |
| `diagnose.py` | per-dataset tables, residual scans, runaway detection, figures |
| `cli.py` | `rmfit ...` subcommands |
| `site.py` | site settings, separate from any campaign: environment variable > `$RMFIT_SITE` or `~/.rmfit/site.json` > default; `mail` (env `RMFIT_MAIL`) is the job scripts' mail address, none if unset; `scheduler` (env `RMFIT_SCHEDULER`) `sge` (default: `templates/run_crc_*.csh`), `slurm` (`templates/run_slurm_*.sh`, bash, with `slurm_partition`, `slurm_account`, `slurm_time`, `slurm_extra`) or `shell` (`templates/run_shell_*.sh`, bash with no scheduler: run in place, log file, lock directory, drive legs looped; the scripts set `RMFIT_SHELL_JOB`, which counts as a job); `sge_queue` / `sge_pe` (env `RMFIT_SGE_QUEUE` / `RMFIT_SGE_PE`): the SGE `-q` line (none if unset) and the `-pe` PE (default `smp`); `setup`: the lines every job script runs before Python (modules, Python), in the template's shell (csh for SGE, bash for Slurm and shell); job detection sees SGE `JOB_ID`, Slurm `SLURM_JOB_ID` and `RMFIT_SHELL_JOB` (2026-11-14) |
| `enginepin.py` | the engine a campaign was started on (AZURE2_ROOT, `_azure2` .so md5, commit); refusal to run on another (2026-11-07) |
| `spec.py` | the fit specification `spec.json` / `SPEC.md`: data-set roles and norms, window, frozen block, bounds, pre-authorized actions, standing decisions; `rmfit spec check\|render\|derive`, `init --spec` (2026-11-07) |
| `decisions.py` | `DECISIONS.md`, the queue of questions outside the spec (`rmfit decisions`, `rmfit decide`); a level add the spec does not pre-authorize ends `accept-pending` with an entry here |
| `review.py` | `rmfit review` / `python3 -m rmfit.review <file.azr>`: the post-fit "GUI eye" as checks on the per-point table (from `AZUREOut_*.out` + `chiSquared.out` + `parameters.out`, engine-free on any archived fit, or from a campaign's incumbent through a single session) -> `REVIEW.md` + `figures/review_*.png`: per-point summary in the archive's wording, top chi2 carriers, norm vs the segment's own optimum n* (opposite direction, fixed norms the data would move), n* trends with angle / energy inside one experiment, residual runs, average-vs-structure (blurred residuals, MeV-wide levels over visible unfitted peaks), zigzag of the calculated curve, levels below the data spacing (target-integration grid), one-data-set levels; `--grid-check` re-scores on a converged grid |
| `drive.py` | `rmfit drive`: the closed loop (DESIGN.md Phase 2.1, 2026-11-08) -- a resumable state machine, one checkpoint per state in `drive_state.json`: build phase PIN_CHECK, AUDIT, REPRODUCE (base + .sav vs the archived chiSquared.out), SPEC_CHECK, SCALE_ANCHOR (spec roles -> `<stem>_drive.azr`), SEED_CHANNELS (free widths at zero -> seed_frac x gamma_W, baked), PAIR_COMPLETENESS (ENSDF states open below the window top without a pair), LADDER (spec `model.ladder_Ex` only), CANDIDATES (ENSDF, broad first, hypothesis required), INIT; loop phase TIGHTEN, SIGN, SHIFTS_FREE, STRUCTURE_ADD, SIGN_POST, DATASET_OFF, REVIEW, EXPORT, CYCLE with the thresholds from `spec.json` `budget` / `stop_rules` / `pre_authorized`; anything outside the spec -> `DECISIONS.md` (a `wait` default ends the leg WAITING); `DRIVE.md` for the supervisor; `run_crc_drive` (templates/run_crc_drive.csh) runs one leg per job and resubmits on CONTINUE; engine states refuse the login node |
| `audit_subagent.py` | `rmfit audit --kind input\|data\|literature`: the three idle-time reviews of skill section 3b as self-contained subagent prompts (`audits/<kind>/PROMPT.md`, from `templates/audit_*.md`: segment table, data files, machine findings, spec roles and standing decisions, readme head); parser of their `FINDINGS.md` (`--findings`, `--status KIND N confirmed`) |
| `audit.py` | `rmfit audit`: engine-free pre-flight checks of the model and data -> `AUDIT.md` (zero free widths, 1 keV level merges, nominal norms != 1, normError as a fraction, FLAG columns as errors, tiny errors, duplicates, spikes, cross-file duplicates, thresholds, `<targetInt>` on test segments, legacy .sav names, output path, engine pin) |
| `covariance.py` | data covariance matrices (`covariance <file>`) and normalization groups (`normgroup <label>`): the engine's whitening redone on the evaluator's rows, the normgroup leader map (2026-10-22) |


Running: `rmfit` is not installed as a package; put its parent on the path:

    PYTHONPATH=<directory containing rmfit>:$AZURE2_ROOT python3 -m rmfit.cli <campaign_dir> status

Tests (13N model): `python3 rmfit/tests/test_spine_13n.py`, `python3 rmfit/tests/test_search_13n.py`.  They
were seconds when written (2026-09-05); since the mid-September target-integration tolerance change a 13N
session opens in ~5 min and evaluates in ~8 s (about 100 sub-points per point, 370k external-capture
amplitudes), so the spine test's polish section now runs 1-2 h -- give it a long timeout or a quiet node.

Several compound nuclei in one project (a natural or compound target) are the engine's own: `keys.py`
knows a parameter kind `abundance` (`ParamKey.nucleus = (Z, A)`), its prior is one of the nuisance rows and
the engine's abundance-sum term comes with its model rows (`engine.ModelRows`).  The earlier joint fit of
several one-nucleus projects (`rmfit.mixture`, 2026-10-12; superseded on 2026-10-13 and since then runnable
only against one old build) is not part of rmfit on pyazr 3.0.  Two 2026-10-13 regressions, both fixed the same day: `from_segments` was first called with
`self.azr` from `Sharded`, which only has `base_azr` (killed five sharded jobs),
and `sign_distance` enumerated pair nodes without the stage-2 member tag (`KeyError` on every call;
`keys.py.bak-2026-10-13c`).  Any edit here must pass `tests/test_keys.py` (seconds) and a `Sharded` fixture
(`tests/test_levelscan_13n.py`, queue job template `tests/_queue_sharded/run_crc_sharded`) before it stays.

Level labels with several compound nuclei in one project (2026-10-30).  The engine
has one J-group per (J^pi, nucleus), each numbered from 1, so in a natural-target project with two odd-A nuclei
(21Ne and 23Ne in an archived campaign) the plain label `"1/2-#2"` named two levels: before
this fix `level_energy` kept ONE energy per label (`default_bounds` then took width caps and floors, the background
test and the total-width floors from the other nucleus's level) and two free energies shared the key `E[1/2-#1]`.
Now every level label of a multi-nucleus model reads `"<J^pi>#<n>@<nucleus>"` (`"1/2-#2@21Ne"`; `keys.level_label`,
the nucleus from the pairs of the J-group's widths, `jgroup_nucleus_tags`); a single-nucleus model keeps the plain
label byte for byte, so ledgers and checkpoints are unchanged.  `ParamKey.ordinal` / `.jpi` / `.nucleus_tag` and
`keys.level_group` (`"1/2-@21Ne"`) parse it; the same-J^pi spacing checks (`local.py`), the backward-elimination
group count, `Campaign._level_energy` (active levels of that J^pi AND nucleus) and
`structure.deactivate_level_variant(nucleus=...)` use the nucleus; `build_keys` refuses duplicate energy/width
keys.  Still per J^pi only (conservative, frees or selects both nuclei's groups): the sign-round groups
(`moves.py`, `search.py`) and `free_mask(jpi=...)`.  Level ADD / levelscan variants are not multi-nucleus aware
(they pick pairs 1 and 2).  Test `tests/test_level_labels_multicn.py` (fixture `tests/fixtures/multicn_labels`,
13N + 15N; Single and Sharded, ~10 s); backups `*.bak-2026-10-30-cnlabel`.

Data covariances, normalization and shift groups (2026-10-22, `covariance.py`).  AZURE2's segment keywords
`covariance <file>` (dev 4b71177), `normgroup <label>` (dev 3e13fb2) and `shiftgroup <label>` work in both evaluators:

- Workers never see `covariance`: `SingleSession` runs a sibling `<model>.nocov.azr`, shard files drop it,
  and the evaluator whitens the assembled data rows with the engine's own algorithm (same rows in file-row
  order, same D = sqrt(C_ii), same B components, V and T; `W = (I + V T V^T) L_B^-1`), so the rows equal
  the engine's `residual_jacobian` rows to rounding and a group split across shards keeps its correlations.
  `SingleSession` compares its chi2 with an engine session on the model as written at its FIRST
  evaluation and stops on a mismatch (RMFIT_COV_CHECK=0 skips it) -- which is also what happens with a
  pyazr built before 4b71177, whose engine ignores the keyword.
- `relvector` terms: the model is frozen at the evaluator's first evaluation (as the engine freezes at a
  run's first evaluation); `ev.refreeze_covariance()` re-freezes at the next one.  Each `rmfit` CLI step is
  a new process, so every polish starts from a model frozen at its own starting point (Chiba-Smith).
- Row weights (`dataset_weight`, `point_weights`) act on the plain rows BEFORE whitening (inflated
  errors, correlations kept); the Sivia loss is refused on covariance rows.  `rows_chi2(r)` is the plain
  chi2; `point_chi2(r)` gives each row's share `r_i (M^-1 r)_i` (adds up to chi2; levelscan/diagnose use
  it), `point_residuals(r)` the plain normalized residuals before whitening.  Per-segment chi2 are the
  engine's partition shares (can be negative).
- `normgroup`: single sessions leave it to the engine (the worker's one-evaluation residuals now use the
  leader's norm for a follower); shard files drop it, a follower's pieces run at their own fixed norm `n0`
  and `Sharded` rescales them to the leader's `n`: `r = (n0/n)(r0 + q) - q`, `q = d/sigma`.
- Keep covariance files under `data/`.  Shards need none (the evaluator reads it from the base model's
  directory) and `single_variant` links one that sits next to the model, but the structure, window and
  levelscan variants link only `data/`.
- `shiftgroup`: every member follows the constant shift of the group's first active segment, so `build_shards`
  puts all pieces of a group into ONE shard, which keeps the keyword (exact there; a finite-difference step of
  the group's shift moves every member).  A group with a free shift is one shift-solo unit.
- `shiftfngroup` (AZURE2 dev 62392fc; rmfit 2026-10-24): every member uses the shift function of the group's first
  active member, which declares it, so a group's pieces go to ONE shard (as a shiftgroup) and no member is ever split
  (a split owner would declare the function twice).  Shard files keep file order, so the owner stays first.  Test
  tests/test_shiftfngroup_sharded.py (fixture tests/fixtures/shift_fn_groups; needs AZURE2_ROOT at a >= 62392fc build).
- Free `sqrtshift` coefficients and `shiftfn` parameters (2026-10-22; before, `build_keys` raised "unknown parameter
  kind" and rmfit could not open such a model): keys `('shift_sqrt', seg)` and `('shift_fn', seg, term=j)`, the
  engine's penalties as rows (`Penalties.from_segments`, equal to pyazr's `penalties()`), bounds and `free_mask`
  like constant shifts (`BoundsPolicy.shifts_free`, `shift_nsigma`); a segment with any free shift term is never
  split and counts as shifted for the shift-solo shards.
- Held finite-difference columns are not computed (2026-10-30, the evaluator "go ahead with the column mask"): the engine
  finite-differences every free shift, sqrt(E)-shift, shift-function and density column (two evaluations each), also
  the ones polish's bounds hold.  `Evaluator.residuals(x, jac=True, cols=...)` now names the columns the caller reads;
  the workers skip the others (pyazr `residual_jacobian(skip=)`, AZURE2 dev 753c6f0 or later), which come back zero
  on the data rows.  `polish` asks for its free columns, `gn_step` for its own.  Older pyazr builds compute every
  column as before.  Before this, a multi-nucleus polish with 8 held shift terms spent ~16 of ~18 evaluations per
  Jacobian on columns it discarded (7.5 h without a finished Jacobian).  Test tests/test_jac_skip_columns.py.
- Fitted target densities (`densityfit` / `densitygroup` on a `<targetInt>` line; AZURE2 branch
  dev f272a35, 2026-10-29): key `density[t<n>]` (`ParamKey.effect` = the
  line's effect number, a group's first line), its prior as a row read from the base .azr
  (`AzrModel.target_densities`), bounds [0.01 N0, max] in `default_bounds` (`SingleSession.density_info`).
  SingleSession only: `build_keys` raises for a shard (whether a shard keeps every `<targetInt>` line is untested).
  The density is in atoms/cm^2 (~1e17): `polish`'s `x_scale="jac"` handles the scale.
- Fitted Doppler temperatures and sample thicknesses (AZURE2 2026-11-05: `doppler T vary sigma [Tmax]` on a
  `<targetInt>` line; `transmission|yield n vary sigma` and `thicknessgroup` on a segment line): keys `doppler[t<n>]`
  (`ParamKey.effect`, as a density) and `thickness[seg<k>]` (`ParamKey.seg`; a thicknessgroup follower has no free
  parameter, so only the leader's key appears). Priors as rows read from the base .azr (`AzrModel.dopplers`,
  `AzrModel.sample_thicknesses`); bounds in `default_bounds`: T in [T0/4, the layout temperature], n in
  [1e-6 n0, inf) (`SingleSession.doppler_info` / `thickness_info`). SingleSession only: `build_keys` raises for a shard
  (a split transmission segment and a thicknessgroup across shards are not designed). The fitted values travel in
  the .sav (pyazr `save_fit`), the .azr keeps the prior centres. Needs a pyazr that knows kinds 9 and 10.
  Test tests/test_doppler_thickness_kinds.py.
- Electron-screening potentials (`screening Ue vary sigma` on a `<targetInt>` line; AZURE2 dev since 2026-10): key
  `screening[t<n>]` (`ParamKey.effect`), its prior as a row from `AzrModel.screenings`, bounds [0, inf) as the
  engine's. SingleSession only. Until 2026-11-05 `build_keys` stopped at this kind. Test tests/test_screening_kind.py.
- levelscan equal-split sign templates (`_patterns`): until 2026-11-06 all 2^n sign vectors were listed before
  `--max-signs` of them were sampled, so any J^pi with >= ~30 open channels was OOM-killed (300-650 GB at n = 30-31). Now the same indices are drawn from 2^n and decoded:
  identical templates, bounded memory. Test tests/test_patterns_signs_memory.py.
- Conv pieces are made only for segments without any energy shift: until 2026-10-22 a fixed nonzero shift,
  a `sqrtshift`, `shiftfn` or a shift group's shift was silently dropped on the conv grid line.

Tests: an archived campaign (`test_covariance_rmfit.py`: rows and Jacobian = the engine's,
Single = Sharded, gradients, refreeze, weights, bake; the old/new bit-identity dumps; the ND2021r scale run).

## Licence

Copyright (C) 2026 the rmfit authors.

rmfit is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.  It is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License in `LICENSE` for details.  The same
licence as AZURE2 (https://github.com/AZURE2DEV/AZURE2), whose `pyazr` engine rmfit drives.
