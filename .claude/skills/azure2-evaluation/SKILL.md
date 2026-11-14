---
name: azure2-evaluation
description: The procedure of a whole AZURE2 R-matrix evaluation, from the data audit to the final tables - stages, the gate each must pass, which skill owns each step, what is recorded (provenance, decisions, held-out data), when to stop, and what the final report contains. Use to plan, run, resume or review an evaluation of a reaction (S factor, rate, resonance parameters) with AZURE2 and pyazr. A single fit is azure2-fit, errors azure2-uncertainty, a level search azure2-resonance-search, an rmfit campaign rmfit-campaign.
---

# An AZURE2 evaluation, from the data to the tables

An evaluation is a claim: these data, with this model, give these numbers (S(E), S(0), N_A⟨σv⟩(T₉), resonance
parameters) with these uncertainties, and the rest of the data agree. This skill is the order of the work, the gate
each stage must pass before the next starts, and what must be written down so that another person can rerun and
judge it. Each step is done with the skill named in it; this skill holds no recipe of its own.

## The evaluation directory

One directory per evaluation, under version control, with:

- `README.md`: the question, the reaction, the energy range, the quantities to deliver, the decisions (one line
  each, dated by stage, with the evidence file), and the status of each stage.
- `data/`: the data files as fitted, each with its provenance (EXFOR entry or paper, table, digitisation, what was
  changed: frame, units, a split of the errors) and its covariance file (`azure2-uncertainty`, its
  data-covariance reference).
- `held_out/`: the data sets kept out of every fit until the last stage (chosen in stage 1, never looked at in a
  fit; see stage 6).
- `model/`: the project (`.azr`, format 3), one file per adopted stage (`stage3.azr`, ...), each a fitted project
  written by `Session.save_fit` (it carries its values: no loose `param.sav` is the fit).
- `runs/`: the outputs of each script (JSON, NPZ), named by stage; `results/`: the final tables and figures.

## Stages and gates

| Stage | What | Owner skill | Gate (all must hold before the next stage) |
|---|---|---|---|
| 1. Data audit | inventory, frames, units, error budgets, overlaps, normalisations; covariance per set; held-out choice | `azure2-setup`, `azure2-uncertainty` (covariance), user-level `data-auditor` | every set has provenance, a split of its errors (point-to-point / correlated / scale / energy scale) and a covariance file or a written reason; `write_covariance.py --check` passes on every file |
| 2. Model | pairs, thresholds, levels from ENSDF, channels, radii, external capture, apparatus (targets, resolution), THM | `azure2-setup`, `azure2-targets`, `azure2-thm`, `azure2-calculate` | `AZURE2 --validate` clean; the model at the start reproduces a published curve or the previous evaluation's χ² where one exists; grid convergence of the target integration (`azure2-check`) |
| 3. First fit | staged fit, nuisances with priors, multiple starts | `azure2-fit` (`lsq_fit.py`, `engine_fit.py`) | the fit converged (status, EDM); the lowest of several starts (`lsq_fit.py --starts`) and `migrad-gradient` started there agree to 1e-5 relative; the snapshot reproduces the objective when reopened |
| 4. Trust | χ² reproduction, pulls, residual structure, χ²/ν per set, Wigner limits, Park norms, norms against priors | `azure2-check` | no pull pattern unexplained; no width above its Wigner limit without a stated reason; every norm within 3σ of its prior or explained |
| 5. Structure | missing or surplus levels, background poles, interference signs, channel radius | `azure2-resonance-search` (`level_test.py`), `azure2-fit` | each change accepted on the significance rules below, on all data; the channel-radius scan (`azure2-check`) recorded: its minimum and its χ² profile |
| 6. Held-out check | refit nothing; predict the held-out sets | `azure2-check`; `level_test.py --hold-out` of `azure2-resonance-search` | held-out χ²/N within the range the fitted sets show; otherwise back to stage 4 with the finding written down |
| 7. Uncertainties | Fisher information, covariance (HESSE, LM), profiles, posterior (NUTS, emcee), bands, systematic axes, model average | `azure2-uncertainty` (`fisher.py`, `profile_scan.py`, `sample_nuts.py`, `draw_bands.py`) | covariance and posterior agree for the quantities quoted, or the posterior is quoted; R-hat < 1.01, ESS ≥ 400 for each quoted quantity; every systematic axis has a variant |
| 8. Deliverables | S(E), S(0), rates, resonance parameters, tables, figures, GNDS/ENDF | `azure2-calculate`, `azure2-plot`, `azure2-exchange` | every number in the tables traces to a run file; figures checked by eye |

Large campaigns (hundreds of parameters, many sign and level rounds on a cluster) run stages 3 and 5 with rmfit:
`rmfit-campaign`.

## Significance rules for structure changes (stage 5)

A change of the level scheme is accepted when, refitted on all the data against the incumbent refitted with the
same budget:

1. the objective falls by more than the F-test's 99 % critical value for the added parameters (`level_test.py`
   prints p; accept below 0.01), **and** BIC falls;
2. the held-out sets do not get worse (`--hold-out`);
3. the new level is physical: θ² below 1 in its channels, a Park overlap J ≥ 0.5 (`azure2-check`), and an energy
   inside or near the data.

Removal is the same test read the other way. A change that only one data set wants is recorded as provisional,
with that set named. Two levels of one J^π that exchange strength (a shared Δχ²) are one decision, not two.

## What is recorded

- **Provenance** of every data set (stage 1) and of every model choice (radius, background-pole energy, level
  J^π from ENSDF or from the fit), in `README.md`.
- **Decisions**: one line each, with the stage, what was decided, why, and the run file that shows it (for example
  "stage 5: 2⁻ level at 19.0006 MeV accepted, Δ = 89 359 for 3 parameters, F-test p below double precision, BIC
  lower by 89 340; `runs/stage5/scan.json`", the numbers of `level_test.py` on `tests/upgrade/T3`).
- **Held-out data**: which sets, chosen when, and the stage-6 result. A held-out set that was used to choose
  anything is no longer held out: record it.
- **The fit of each stage**: the project written by `save_fit`, the objective, the data χ², N, the free
  parameters, the fitter and its status, the build (`pyazr describe`).
- **The options** of the calculation (basis, external capture, long-wavelength) as `<config>` keys of the project,
  not as command-line switches only.

## Stop conditions

Stop and write down why, instead of going on, when:

- a stage's gate fails twice after its fix;
- the objective's minimum depends on the start (multiple starts end more than Δ = 1 apart) and no physical reason
  picks one: report both;
- a quantity to deliver depends on a model choice by more than its statistical error and no data decide the choice
  (radius, background-pole energy): report it as a systematic axis (`azure2-uncertainty`), do not pick the
  favourable variant;
- the held-out data disagree with the fit: the model or the data are wrong; this is a finding, not a tuning target;
- the posterior is not closed (a non-positive Hessian eigenvalue, chains drifting along one direction): a prior
  is needed, and choosing it is a decision for the evaluator.

## The final report

1. The question and the delivered quantities with their uncertainties (statistical, each systematic axis, total).
2. The data table: set, reaction, quantity, range, points, scale error, covariance, fitted norm, χ²/N, held out or
   not.
3. The model: pairs, radii, levels (J^π, E_x, widths or ANCs, θ²), background poles, options.
4. The fit: objective, data χ², N, free parameters, fitter, convergence; the decisions of stage 5 with their
   statistics.
5. The uncertainties: method (covariance, profile, posterior with its R-hat and ESS), systematic axes and their
   combination, what is not included.
6. The held-out check.
7. The figures (data with the fit and pulls, the S factor with its band and decomposition, the rate ratio to a
   previous evaluation).
8. Reproduction: the commands, in order, from the data to the tables.

## Pitfalls

- An evaluation that changes the data after looking at the fit (a "bad" point removed because it pulls) must say so
  and show the fit with and without it.
- Comparing two fits with different data sets, objectives, bases or builds compares nothing (`azure2-fit`).
- A quantity quoted from the covariance where the posterior is not Gaussian (a width near 0, an interference
  minimum) is wrong in the tails; compare with the posterior first (`azure2-uncertainty`).
- The previous evaluation is a model, not data: fit it only as a prior (prior covariance), never as points.
