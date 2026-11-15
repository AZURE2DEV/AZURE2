# rmfit design

rmfit runs an evaluator's global R-matrix fit with AZURE2 as resumable, checked steps that a script executes on a cluster
node and an LLM agent (or a person) can supervise and resume from plain-text state files.  This document describes what it
is for, how it is built and what is still open.  How to use it in a session is in the `rmfit-campaign` skill
(`.claude/skills/rmfit-campaign/SKILL.md`); the module table is in README.md.

## Principles

1. **Human input at the front and at the end.**  The evaluator declares the objective before the run (`spec.json`: data
   sets and their roles, normalizations and priors, fit window, model scope, what may be adopted without asking) and
   reviews a candidate after it.  The middle is a closed loop that does not stop to chat.
2. **Data treatment stays with the evaluator.**  Normalization scales and systematics of a data set, dropped rows, shape-only
   treatment, outliers: these are recorded in the spec with a date and an owner, never decided inside a fitting step.  A
   question the spec does not answer becomes an entry in a decision queue (`DECISIONS.md`) and the run waits.
3. **Judge per point, not by total chi-squared.**  A large chi2/N is often carried by a few isolated points; the review
   reports the fraction of points within 1 and 2 sigma, the median |r| and where the chi-squared sits.
4. **Audit inputs before fitting.**  Input errors (units, frames, error columns, normalization-error units, swapped columns)
   are found by audits, not by the fit.
5. **Autonomy is measured.**  A benchmark of degraded starts against known references scores every change to the method.
6. **One engine per campaign.**  A campaign records the engine build it was started with and refuses another unless told,
   because objectives are not comparable across builds.

## Components

| module | role |
|---|---|
| `engine.py`, `model.py`, `_worker.py` | exact evaluation sharded over single-thread pyazr sessions (chi-squared is additive over segments); Jacobian-cost shard balancing |
| `local.py` | bounded trust-region polishes in reduced-width-amplitude space with the analytic Jacobian; Wigner-limit caps, energy windows, sign locks; sanity checks |
| `search.py`, `moves.py`, `mixture.py` | interference-sign rounds (single, pair, channel, group flips screened by one Gauss-Newton step) and the population stage (sign-perturbed restarts polished to convergence plus crossover between them: the step that crosses basins) |
| `structure.py`, `levelscan.py`, `windows.py` | level add / remove / relocate with a required hypothesis per candidate; residual-driven level search; energy-window variants |
| `keys.py` | sign-pattern keys and distances between solutions |
| `covariance.py` | data covariance matrices, normalization groups, shift groups, model-scaled (relative) terms with a fixed reference model |
| `spec.py`, `decisions.py` | the declared objective and the decision queue |
| `audit.py`, `audit_subagent.py` | machine pre-flight checks (pair masses, dead channels, normalization conventions, data integrity, target energies, legacy parameter names, output paths) and templated prompts for three agent audits (input file, data, literature) |
| `review.py` | the automated "GUI eye": per-point summary, top chi-squared carriers, normalization direction and trends, residual runs, broad levels averaging narrow peaks, levels below the resolution, grid exploitation, close same-J^pi pairs |
| `drive.py` | the closed loop (state machine, below) |
| `uq.py` | posterior sampling and uncertainty bands as the last stage |
| `campaign.py`, `ledger.py`, `cli.py` | campaign directory, SQLite ledger of every candidate, STATUS/REPORT files, command line |
| `enginepin.py` | engine fingerprint and pin |
| `bench/` | recovery benchmark: cases, scorer, results |

## The drive

`rmfit <dir> drive` is a resumable state machine; one call is one leg (one scheduler job), and a leg stops at DONE, at
WAITING (a decision is queued), at a wall-time cap (the job resubmits itself) or on request.

Build phase (a new model): engine check, audit, reproduce the starting objective, spec check, scale anchor (one
representative data set's normalization fixed), seed channels (lowest L per pair), pair completeness, energy ladder,
candidate bank.  Loop: initialize, tightening polish, sign rounds, population stage, free energy shifts if
pre-authorized, level search and structure adds (broad levels first, narrow alignment later; a recall test on known
levels before level adds are trusted on a new reaction), sign rounds after an adoption, data-set-off test on the worst
sets, review, export; another cycle while it gains, then UQ if enabled.  Every state has a budget; stop rules end a
stage on small gains.  The decision kinds that stop the drive are the evaluator's: adopting a level (unless
pre-authorized), switching a data set off, a new high-severity review flag, failed UQ gates.

## Uncertainty quantification

`rmfit uq` (also the drive's last state): export the incumbent; hold saturated poles (Park's J below a threshold: their
widths sit at the channel-radius ceiling and the data do not set them); Laplace covariance at the mode; NUTS chains with a
dense metric; automatic gates (split R-hat, effective sample size, divergences, E-BFMI); linear and sampled bands on
requested energy grids; a report.

- **Bounds.**  Particle amplitudes are sampled inside the campaign's theta^2 box (the box the polish ran in), by a tanh map
  with its Jacobian; the metric in the mapped coordinates comes from the Hessian at the mode.  Data alone need not bound a
  background pole; an unbounded posterior can run away along it.  A different prior on an amplitude is the evaluator's.
- **Gates.**  A failed gate re-runs with longer chains, except when R-hat is far off (chains in different regions do not
  merge with longer runs) or when the re-run would exceed a wall-time limit; then the stage reports FAIL with the reason and
  the drive asks the evaluator.  Error inflation (the Birge factor) is reported, never applied.
- **Basins.**  Two fits of equal quality can differ in an observable by more than either fit's band (for example through
  one background pole that is saturated in one fit and free in the other).  UQ on the incumbent alone does not see this; run
  it on the equal-quality alternatives the population stage finds.

## Benchmarks

A case names a reference fit and degraded starts (perturbed signs, jitter, removed levels, or the start before a key human
decision).  The score: data chi-squared at or below the reference's, objective within 0.5 %, no review flag the reference
lacks, zero wait decisions; node-hours and human turns are recorded.  Reference objectives are re-measured on the pinned
engine.  The shipped tier-0 case uses AZURE2's own `tests/13N` fixture.

## Wide energy ranges: hypothesis under test

For models with many hundreds of parameters, the working hypothesis is:

1. structure search in overlapping energy windows (a level or two of overlap), levels below a window fixed with the broad
   ones free;
2. quantities the windows share declared in the spec before either runs (the normalization of a data set that crosses the
   boundary fixed or tied to one window; the J^pi of an overlap level fixed);
3. frozen snapshots, never synchronization: nothing is reconciled before the merge;
4. the merge as a drive run of its own (concatenated levels, sign rounds and the population stage across the junction,
   genuine conflicts queued as decisions);
5. the merged full-range fit as the stopping test, and UQ on the merged model and its equal-quality alternatives.

Below roughly 200 free parameters a full-range fit is the default.  The hypothesis is scored at each merge with the
benchmark yardstick, the decisions queued, the cost, and what the evaluator changed afterwards.

## Open items

- A scheduler layer: the site file (`site.py`: `~/.rmfit/site.json` or `$RMFIT_SITE`) chooses SGE, Slurm or plain-shell templates and
  holds the mail address, the SGE queue and parallel environment, Slurm partition / account / time / extra options and
  the setup lines both kinds of template run.  Still open: the Slurm templates have not yet run under a real Slurm
  scheduler (checked with bash and a simulated job environment); SGE resource requests other than the slot count
  (memory, run time) go in no setting yet.  The plain-shell kind (`"shell"`) has no queue: nothing stops a run from
  taking more cores than are free.
- Multi-compound-nucleus models in sign rounds and the level search; density, Doppler, thickness and screening kinds in
  shards.
- Channel lines and particle pairs are not added automatically; the energy ladder is unexercised.
- The level search has been validated on recall; a level no reference model had has not yet been found by it.
- Competing-resolution test for near-duplicate levels; differential-convolution shards.
- An engine build identifier compiled into pyazr, so a pin can name a commit directly.
