# Running a CLI fit: convergence, gradients, priors, long jobs

The CLI fit (mode 2; mode 4 fits and then extrapolates) runs MIGRAD on the full objective (data χ² plus normalization,
shift and nuisance penalties) and writes `param.sav`, `chiSquared.out` and `fitStatus.out` when it ends. This page says
how to tell a converged fit from a stopped one, how to use the analytic gradient and its diagnostics, how nuisance
priors behave, and how to run fits of days on a batch cluster. Behaviour reference:
`docs/source/user_guide/fitting.rst`, `docs/source/reference/command_line.rst`,
`docs/source/reference/output_files.rst`. Scripted prompts and the parameter-file ordering of the CLI:
`azure2-calculate`.

## Convergence

- MIGRAD has no call limit unless `--max-calls=N` sets one. Every fit ends by printing
  `MIGRAD converged: valid minimum, EDM ...`, or `WARNING: MIGRAD did NOT converge: stopped at the call limit ...` /
  `... invalid minimum ...`, and writes `output/fitStatus.out` (`minimizer`, `fitter`, `status
  converged|call_limit|iteration_limit|stopped|invalid|failed`, `fitter_status`, `valid`, `edm`, `edm_goal` (`none`
  except MIGRAD), `calls`, `call_limit`); every fitter writes it, LM too. Read it before calling a fit
  converged; `chiSquared.out` and `param.sav` are written either way, at the lowest objective the fit reached (at the
  call limit: the lowest MIGRAD evaluated). LM ignores `--max-calls` (`call_limit N (--max-calls) ignored`
  and a `Note:` line).
- builds before 2026-10-31 (expires: 2027-02-01) stopped MIGRAD and NLopt at a fixed 50,000 function calls without a
  message; one 145-parameter refit stopped there in most of its J^π fits and looked converged. On such a build a log
  ending near 50k Minuit calls means "not converged": restart from `param.sav` until χ² stops dropping.
- The `Iteration: N` counter in the log is the engine's own count of objective evaluations (`AZURECalc.cpp`, printed
  every 10th call); it can exceed Minuit's count, and MIGRAD checks its budget only between its own steps, so runs
  overshoot a limit (fits with 300–650 free parameters logged 56k–136k evaluations before stopping).
- On a large model (hundreds of free parameters) a flat stretch is not a minimum: long plateaus broken by discrete drops
  recur throughout a run, including near its end. Do not infer convergence from a plateau of a few hundred evaluations;
  trust the explicit message, or report the result as final for this run while noting that a warm-started refit may
  improve it.
- Results are written only when the fit ends; a killed or crashed run loses everything after its last `param.fit`
  snapshot (below).

## The analytic gradient (`--use-gradient`)

- The engine supplies a Gauss–Newton G2 from the residual Jacobian to seed MIGRAD (log line
  `Seeding MIGRAD with the analytic Gauss-Newton G2.`) and takes gradient and G2 from one Jacobian per step (about
  1.3 s per step for 79 free parameters on 11k points at 8 threads). `Analytic G2 unavailable ...` means the old
  numerical seed is in use.
- **Compare a test and its control with the same minimizer**, preferably `--use-gradient`. On one model the analytic
  fit converged about 490 below a numerical-gradient fit that MIGRAD had declared converged early on its noisy
  gradient, so a numerical control overstated a level-add gain (−902 against the numerical control, −545 against the
  analytic one).
- A log whose `Iteration` counter stalls during the numerical-gradient seeding phase is not proof of a hang: Minuit's
  parallel numerical-gradient threads do not advance it.
- builds before 2026-10-10 (expires: 2027-02-01): `--use-gradient` was no faster than the numerical fit, single χ²
  evaluations ran on one core, and penalized norms (the GUI's inert nuisance flag on them), target-integrated A_y data
  or many free shifts made the fit fall back to finite differences; the finite-difference step for norms was
  1e-6·(|x|+1), wrong by up to 14× for a thick-target norm of 2 × 10⁻⁵.
- On models with target effects the cost is dominated by re-shifting: every sub-point's Coulomb functions are
  recomputed when a segment shift changes. Check the grid density of the target integration (`azure2-targets`); one
  model used 200–300 points per width where 10 was converged (`azure2-targets`).

Diagnostics (environment variables; standard error):

| variable | effect |
|---|---|
| `AZURE_GRAD_CHECK=1` (or `exit`) | compare sampled analytic derivatives, and the Jacobian path, with central finite differences; `exit` stops after the check |
| `AZURE_GRAD_CHECK_N=12`, `_IDX`, `_SEGS`, `_H` | number of sampled parameters, which parameters or segments, the step |
| `AZURE_GRAD_TIMING=1` | per-call timing breakdown |
| `AZURE_GRAD_ADJOINT=1` | force the adjoint gradient instead of the Jacobian path |
| `AZURE_TE_WEIGHTS_CHECK=1`, `AZURE_TE_WEIGHTS_FD=1` | check the target-effect combiner weights |
| `AZURE_POOL_CHECK=1` | compare the pooled, parallel χ² with a fresh serial one |
| `AZURE_SERIAL_CHI2=1` | evaluate χ² serially |
| `AZURE_TRACE_FCN=1` (or `fresh`) | print every χ² with its parameter vector, bit-exact, to compare runs; `fresh` also checks each against a fresh model copy |

The finite-difference reference of `AZURE_GRAD_CHECK` uses the single-threaded χ², about 10× slower than a fit
evaluation. Differences of 10⁻⁴–10⁻³ near narrow levels are finite-difference truncation (they shrink as h²); larger
ones are a defect to report.

## Nuisance priors on level parameters

- A `<parameterSettings>` entry with nuisance 1 on a level energy or width adds ((p − p₀)/σ)² to χ², with p₀ and σ given
  as physical values (MeV; eV for widths; fm^-1/2 for ANCs) and converted to reduced widths
  (`docs/source/user_guide/fitting.rst`, `docs/source/user_guide/chi_squared.rst`). The `--use-lm`/GSL penalty rows
  apply them too; pyazr and the MCMC never apply `<parameterSettings>` priors (the mode-6 `<mcmc>` priors act on the
  physical values directly).
- builds before 2026-10-27 (expires: 2027-02-01) evaluated the penalty on the wrong parameter whenever any fixed level
  parameter preceded the flagged one (a fixed ground-state energy is enough): a parameter at its own nominal value
  paid a large penalty and MIGRAD pulled another parameter toward that nominal (+12,460 instead of +4 in one case; a γ
  width at its nominal paid +273,466). On such builds do not use nuisance entries on level parameters. Test on any
  build: `calculate` with the external parameter file, and compare the log's `Total Chi-Squared` (with penalties) with
  `chiSquared.out` (data only); `tests/nuisance_level_index` guards it.
- Limits and priors written by the Qt GUI, and their level numbering: `references/parameters-and-sav.md`.

## Long fits on a batch cluster

- **Memory.** A mode-2 fit holds max(4, `OMP_NUM_THREADS`) + 1 copies of the model and data (CompoundSystem +
  WorkingCopy), plus
  one for the `param.fit` snapshots. Each copy is the whole data set including every target-integration and convolution
  sub-point, so heavy target-effect models are large: 4,600 points with many target-integrated data sets took about 3.4
  GB per copy, 80–100 GB at 24 threads, steady for days. Measure the peak (`maxvmem` in the scheduler's accounting) and
  reserve it per slot.
- builds before 2026-10-05 (expires: 2027-02-01) leaked memory in glibc's per-thread arenas (6–17 MB per evaluation
  at 24 threads, until the node killed the fit); if a long fit's memory climbs, check the build date first. Builds of
  that period also occasionally segfaulted (exit 139) on target-effect models.
- builds before 2026-11-06 (expires: 2027-02-01) that already shared point tables between segments could crash
  (use-after-free) in threaded fits whose segments share energies and shifts (shift groups, shift functions). A run
  that finished is fine; rerun a crashed one rather than suspecting the model. Single-thread runs were not affected.
- **When a fit is killed the job script keeps going**: only the AZURE2 process dies ("Killed"/exit 137 for out of
  memory; "Bus error" when its binary was replaced on a network file system; exit 139 for a segfault). The job can
  restart from the last snapshot, losing at most 100 evaluations. A restart loop: after each pass treat
  `<outdir>/chiSquared.out` as "finished" (it is written only at the end of a fit, never by the snapshots); otherwise
  copy `param.fit` to `restart_<n>.sav` and run mode 2 again from it. Guard it: at most 15 restarts, and stop when two
  passes die within 10 minutes of starting (the node is out of memory; resubmit elsewhere).
- A restart resets Minuit's covariance, so the first few hundred evaluations rebuild the gradient (brief χ²
  excursions are expected). If the last snapshot fell on a bad probe, verify it with a `calculate` run before
  trusting the restart (`param.fit`: `docs/source/reference/output_files.rst`).
- **The autosaved point can be a line-search or gradient probe** (one crashed run's last autosave printed
  χ² = 8.5 × 10⁴⁰). Check the log's χ² at the last multiple of 100 evaluations and restart from `param.fit` only if
  that value is on the plateau. Otherwise restart from a finished fit's `param.sav`, or from the interim
  `parameters.out` (physical values, written into `<levels>` and verified with a `calculate` run). Taking only the
  `segment_*` lines of `param.fit` is always safe (a norms-only first stage handing over to the second).
  Builds before 2026-09-30 (expires: 2027-02-01) wrote the unfitted input energies and widths into `param.fit`, so a
  restart began far off (290,000 against 26,646).
- Run fits on batch nodes, never on a login node.
