# Phase 2: one objective — design

ROADMAP Phase 2: an `inference/` target with one `ParameterVector` and one `Objective` made of a list of terms (value,
gradient, residuals, Jacobian); every penalty of `EData`, `CNuc`, `ESegment` and pyazr's `penalties()` becomes a term,
and the CLI fit, the gradients, the API and the output writer call it. Gate: bit-identical χ², gradients equal to
finite differences, no penalty computed in Python.

## 1. What exists (inventory of 2026-11-09, `feature/upgrade` 33538a5c)

Abbreviations: `AC` = `src/AZURECalc.cpp`, `ED` = `src/EData.cpp`, `MC` = `src/AZURECalcMCMC.cpp`,
`API` = `api/src/AZUREAPI.cpp`, `PY` = `pyazr/azure2.py`.

| # | term | value | where |
|---|---|---|---|
| T1 | data χ² per segment | Σ((f − n d)/(n e))² | `AC` 631–649; residuals `AC` 2105 |
| T1a | THM profiled norm | χ² = S_dd − S_md²/S_mm, norm not a parameter | `ESegment.cpp` 312, 323 |
| T1b | THM experiment group | shared norm and background profiled at the group's last segment | `EDataThm.cpp` 687 |
| T2 | data covariance | rᵀM⁻¹r per group, as a correction to Σr² | `ED` 4743 |
| T3 | normalization prior | ((n − n₀)/(n₀ σ%/100))², every segment with σ ≠ 0 (fixed too), not THM | `AC` 650–655 |
| T4 | energy-shift priors | ((ΔE − ΔE₀)/σ)² + (b/σ_b)² + Σ(a_j/σ_j)², varied terms only | `ESegment.cpp` 945 |
| T5 | nuisance priors | ((p − nom)/err)², `useAsNuisance` levels and cbkg | `AC` 2339 |
| T6 | Park wall | Σ_{J<0} (J/10⁻³)² | `CNuc.cpp` 2867 |
| T7 | abundances | per-nucleus priors + per-sample sum-to-one | `CNuc.cpp` 570, 590 |
| T8 | screening, density, Doppler, thickness priors | ((x − x₀)/σ)² on the stored (clamped) value, fixed ones too | `TargetEffect.h`, `ESegment.h` 372 |
| T9 | measured scattering lengths | Σ_{Re, Im} ((a − a_m)/σ)² | `ScatteringLength.cpp` 209 |
| T10 | prior covariance | (p − μ)ᵀC⁻¹(p − μ) | `PriorCovariance.cpp` 273 |
| T11 | gamma_decay priors and the remainder's wall | Σ ((q − c)/w)² over the record numbers with a prior and the remainder 1 − Σb; ((Σb − 1)/10⁻³)² for Σb > 1 | `GammaDecayMap.cpp` (`gamma_parametrisation.md` §2.6); rows in `AppendModelRows`, gradient direct (`GradAccum::direct`) |

The paths that evaluate them, and which they include:
1. the Minuit FCN `operator()` (`AC` 561–742);
2. `Chi2On`/`Chi2Value` (finite-difference columns, LM cost; `AC` 815–882);
3. the adjoint `Gradient()` (`AC` 920–1240);
4. `GradientFromJacobian` / G2 (`AC` 1385–1485);
5. the GSL residuals (`AC` 1820, 2169–2200);
6. the MCMC (`MC` 52–250, 703);
7. the API (`API` 990–1168);
8. pyazr `objective()` (`PY` 2107);
9. the `chiSquared.out` writer (`ED` 2725–2830).

Paths 1–3 include all ten terms. In the other paths:

| term | 4 | 5 | 6 MCMC | 7 API | 8 pyazr | 9 writer |
|---|---|---|---|---|---|---|
| T1a/b THM profiling | Y | Y | **no** | Y | Y | Y |
| T3 norm | Y | row | prior | **no** | Python | all segments, THM too |
| T4 shift | Y | row | prior | **no** | Python | Y |
| T5 nuisance | **no** | Y | **no** | **no** | **no** | **no** |
| T6 Park | gradient only | **no** | wall (−∞) | Y | engine | cached |
| T7 abundance | Y | priors only (sum dropped) | split | Y | engine | cached |
| T8 target priors | Y | Y | prior | **no** | Python (fixed skipped) | Y |
| T10 prior covariance | rows | rows | replaces the diagonal priors | stored only | rows | cached |

The same list of terms is written out about eight times: `operator()`, `Chi2On`, `Gradient`,
`GradientFromJacobian`, `G2`, `PrepareFreeParams` (by name), `IndicesWithOwnPrior`/`HasNuisanceParameters` (by name),
MCMC `BuildAutoPriors` (by name), the API with pyazr, and the writer. Penalty totals are stored into whichever
`EData` was evaluated (a pooled clone during a fit), and the writer re-reads or recomputes them. pyazr's
`chi2_and_grad` has no gradient for any of the terms it adds in Python.

## 2. The design

**`inference/` target** (library `azure2_inference`, depends on physics/ and on the engine through adapters only,
like physics/'s `CNucSystem`):

- `ParameterVector`: the full vector in today's order (`BuildParamIndexMap`, `AZUREGrad.h` `ParamKind`), with kind,
  label, value, fixed flag, limits and the packed↔full map. It is built once per configuration, and every path reads it
  instead of matching names (`PrepareFreeParams`, `IndicesWithOwnPrior`, `BuildAutoPriors`).
- `Term`, an interface: `value(p)`, `gradient(p, g)` (adds into g), `residuals(p, r)` and `jacobian(p, J)` when the term
  is a sum of squares (`num_residuals()`, 0 if not), and `name()` for the writer. Its parameters are packed indices.
- `Objective`: the ordered list of terms, `value = Σ term.value`, the same for gradients, residuals and Jacobians. The
  terms are T1–T10, one class each. T1 (with T1a/b and T2) wraps today's data evaluation.
- The output writer prints one line per term from the same objects, with no cached totals in `EData`.

**Order of the work, each step bit-identical on the paths that already agree:**
1. `ParameterVector` from `BuildParamIndexMap`; the name-based lookups replaced by kinds. Check: identical outputs and
   identical free-parameter lists in every path.
2. T3–T10 as terms, used by paths 1–5: `operator()` and `Chi2On` sum the objective; `Gradient`, `GradientFromJacobian`
   and the GSL rows take each term's gradient or rows. The penalty lists are deleted from `AC`. Check: χ² bit-identical
   in paths 1–3 on the 179 projects × 3 bases; gradient against finite differences on T1–T3; outputs identical.
3. The writer reads the terms (no cached totals). Check: `chiSquared.out` identical, except where decision D4 changes
   it.
4. The API returns the objective's terms (value and gradient), and pyazr's `penalties()` becomes a call. Check:
   `tests/pyazr/output_penalties_test.py` and pyazr's χ² identical, except under decisions D3 and D5.
5. MCMC reads the same terms as log-priors (−χ²/2 plus its normalization constants). This needs decisions D1 and D2.
   **PI, 2026-11-09: the MCMC shares the fit's configuration** ("The MCMC should share the same configuration with the
   frequentist"). The posterior is exp(−χ²/2) of the fit's objective, term for term (nuisance priors and the prior
   covariance included), from the same `.azr`. The MCMC differs only where stated: J ≤ 0 lies outside the support
   (D2), and its Gaussian normalization constants. Its name-based `BuildAutoPriors` is replaced by the terms'
   enumeration.

## 3. Decisions (the paths disagree today; unifying them changes χ² somewhere)

**PI, 2026-11-09: every proposal D1–D7 below accepted** ("I agree with all your proposals. The fitting was getting
messy").

- **D1, THM profiling in MCMC.** The MCMC uses a plain Σr² with the THM norm unprofiled (it keeps whatever the clone
  holds), so it disagrees with the fit on THM projects. Proposal: the MCMC uses the profiled T1a/b, as the fit does.
- **D2, Park.** The fit uses a quadratic wall, the MCMC −∞, and the GSL residuals nothing. Proposal: one term, a
  quadratic wall in fits and GSL alike. The MCMC keeps −∞, a support restriction stated as such.
- **D3, prior covariance and diagonal priors.** The fit, API and pyazr apply both on the same parameters, which counts
  them twice. The MCMC replaces the diagonal priors. Proposal: as the MCMC does, the covariance replaces the diagonal
  priors of its parameters everywhere (χ² of projects with `priorCovariance=` changes).
- **D4, the norm prior of profiled THM segments in the writer.** The writer counts a norm term the fit never
  minimized. Proposal: the writer follows the fit.
- **D5, fixed target-effect parameters.** C++ penalizes a fixed screening, density, Doppler or thickness value away
  from its centre, while pyazr skips it; a fixed norm pays, a fixed shift does not. Proposal: a fixed parameter pays
  no prior (it is a constant). This changes χ² only by a constant, but the reported values change.
- **D6, nuisance priors (T5) outside the CLI.** They are missing in the API, pyazr, MCMC and the writer. Proposal:
  every path includes them (a term like the others).
- **D7, abundance sum-to-one in the GSL/LM step.** It is dropped there with a warning. Proposal: it becomes a residual
  row like the rest (it is a sum of squares).

Steps 1–2 change no χ² and need none of these decisions. Steps 3–5 apply them, each decision in its own commit with
its own comparison, so that every changed output is attributable to one decision.

## 4. As built (branch `feature/upgrade-2-objective`, 2026-11-09)

| commit | what | outputs |
|---|---|---|
| c99a8e07 (2a) | `Objective`, one ordered list of terms; `operator()` and `Chi2On` sum it | identical |
| 802eb6fb (2b) | `ForEachDiagonalPrior`, one enumeration for `GradientFromJacobian`, `G2`, the least-squares rows (no name matching) | identical |
| 615c71f8 (2c) | nuisance and own-prior bookkeeping by `ParamIndexMap` kind | identical |
| 82ac146a (D5) | a fixed parameter pays no prior | constants only |
| 64aba582 (D4) | no norm prior for a profiled THM norm in the writer and pyazr | THM `chiSquared.out` |
| 5834693d (D3) | the prior covariance replaces its parameters' diagonal priors (`ReplaceDiagonalPriors`) | `prior_covariance` |
| 3cd21e02 (D1) | the MCMC likelihood is `Objective::DataChi2` (THM norms profiled) | THM MCMC |
| 5cf9cad9 (D2, D7) | Park's wall and the abundance sum-to-one as least-squares rows (`AppendModelRows`) | Park and abundance LM/GSL fits |
| 32f067df (step 4, D6) | the API's objective, gradient and breakdown; pyazr computes no prior; nuisance priors in the API and the writer (pyazr 2.11.0) | pyazr objective |
| ab9e13e1, d35f6d6b (step 5, D6) | the MCMC's priors by kind from the fit's enumeration; nuisance priors in the MCMC | MCMC with nuisance priors |

The test changes are in commits of their own: 1c11fe1e (`prior_covariance` (c), D3), bf1c1f45
(`pyazr_azr_parse_cache`, step 4), 58be9a73 and e8468cdc (`pyazr_objective_terms`; e8468cdc restores
`objective_test.py`, which 86ac2970 had overwritten), c6caf712 (`nuisance_prior` (d), (e)).

**Evidence.**
- **2a + 2b:** 179 × 3 identical. The fit tests of 14 projects run through an archiving wrapper (every engine call of
  their `check.sh`, MIGRAD numerical and analytic, LM and GSL; 215 runs) give identical outputs against the build
  before. ctest 119/119, `run_tests.sh` 133/133, threads identical.
- **2c, D3–D5:** 179 × 3 identical. The one test that moved, `prior_covariance` (c), is re-baselined for D3: the
  difference with and without the covariance is 2C⁻¹(p − μ) minus the replaced shift prior's gradient, 40.
- **D1, D2, D7:** 179 × 3 identical. Of 179 fit-test runs, five changed, all as decided:
  - the three THM MCMC runs (D1). `mcmc_classification`'s THM chain had run to a level energy of −3.2 × 10⁷ MeV with
    the unprofiled norm; it now sits at 6.29 MeV;
  - `natural_target`'s LM and GSL fits (D7): the full objective at their optimum goes from 32.8674 to 32.8673, now that
    the sum-to-one term is minimized.
- **Step 4:** 179 × 3 identical; ctest 130/130 (with the slow tests) after `pyazr_azr_parse_cache` was updated: the
  penalties are now the engine's, on the model it loaded, so an edit of the `.azr` on disk no longer moves them within
  a session, as it never moved the χ²; `run_tests.sh` 133/133; threads identical. `pyazr_objective_terms`:
  objective = the sum of its terms, its gradient = finite differences on all 19 parameters of T2 (2.8 × 10⁻⁸; norm,
  shift, density, abundance and covariance priors non-zero there).
- **Step 5:** of 148 fit-test runs (MCMC projects included) none changed; ctest 130/130, `run_tests.sh` 133/133.
  `nuisance_prior` (e) fails on the step-4 build and passes on step 5.

**Not done (open).** The writer still reads totals stored by the final evaluation (`Set...Penalty`) rather than the
terms themselves (step 3 as designed). The `inference/` target is not a separate library yet: the terms use the
engine's `EData` and `CNuc` directly, and become a library when Phase 4 makes the data plain values. The MCMC's
level-energy spread and its user priors keep their code.

## 5. T9, the model data (Phase 4 step S15, decision M10)

T9 is one term, `ModelDataTerm` in `src/Objective.cpp` (`ModelOnly() = true`, at T9's place in the list, so the sum
keeps its order), over the rows of `ModelDataSet` (`include/ModelData.h`): first the lines of `<scatteringLengths>`
(Re, then Im, each when its σ > 0), then the `datum` records of `<modelData>` in file order (one row, two for a
scattering length with an imaginary part). Each row is (q_k(θ) − d_k)/e_k; the term is Σ r_k² in row order. With
`<scatteringLengths>` alone the rows, their order, their χ² and every output are those of the former
`ScatteringLengthTerm` (the 179 × 3 comparison of the branch).

Every path asks the objective: the value (`Objective::Chi2`, `ModelChi2`; the MCMC likelihood `ModelDataChi2`), the
least-squares rows (`ModelDataResiduals`, the GSL path; `AppendModelDataRows`, LM and the API, after the data rows and
before the prior covariance, where the scattering-length rows were), and the gradient (`ModelDataGradient`, the
Minuit adjoint path and the API). The six hand-written scattering-length sites of `AZURECalc`, `AZUREAPI` and
`AZURECalcMCMC` call these.

- **Scattering lengths**: `ScatteringLength::ConstraintValue`; gradient and rows exact, through the level-matrix
  adjoint at the threshold (`ScatteringLength::ConstraintGradient`).
- **Level quantities** (`resonance_strength`, `resonance_energy`, `total_width`, `partial_width`, `branching_ratio`,
  `anc`): from the observed energy and widths of the current fit values, `CNuc::ObservedOfFitValues` (TransformOut's
  two steps without messages: formal → Brune by `brune_from_formal` in the formal basis, then `observed_levels`; the
  N_f^−½ that TransformToBrune sets is put back, so nothing a calculation reads changes). The formulas are
  `physics/ResonanceQuantities.h`. Gradient and rows: central differences of the observed quantities over the level
  parameters of the data's J^π groups (all groups with external capture), each parameter set and set back by
  `CNuc::SetFitLevelParameter` (bit for bit), steps 10⁻⁵ max(1, |E|) MeV and 10⁻⁵ max(|γ|, 10⁻³) MeV^½; Park's shift
  functions are recomputed after. An exact chain rule through the formal → Brune iteration was not written (no
  adjoint exists for it).

Totals: the writer gets the `<scatteringLengths>` part (`Total-Scattering-Length-Chi-Squared`) and the `<modelData>`
part (`Total-Model-Data-Chi-Squared`, written only with the section); pyazr's `penalties()["scattering_length"]` is
the whole term ("model data"); `pyazr.model_data` gives the rows of the section as records.

## 6. The ensemble sampler samples this objective (2026-11-13)

**Measured before the change** (branch start ac2488a0; a scratch program built on the engine, `work/sampler/measure/`,
outside the repository; −½ χ² of `AZURECalc::operator()` with the fit's limits manager against the sampler's
`CalculateLogPrior` + `LogLikelihood` / `LogLikelihoodPhysical` after `UpdateParameterVectors`, `BuildAutoPriors` and
the prior covariance's set-up, at the start, at a displaced point and at a point outside a limit). Line numbers are
those of `src/AZURECalcMCMC.cpp` (`MC`) and `src/Objective.cpp` (`OB`) at that commit.

| term or limit | fit (`OB`) | sampler, fit basis (RWA) | sampler, observed basis | lines (`MC`) |
|---|---|---|---|---|
| T1 data χ² (THM profiled) | 23–80 | yes, `DataChi2` | yes, the same | 97, 172 |
| T2 data covariance | 84–90 | yes, `DataChi2` | yes | 97, 172 |
| T3 norm prior | 59–68 | yes, Gaussian from `ForEachDiagonalPrior`, normalised | yes, the same | 345–357, 552–579 |
| T4 shift priors | 71 | yes, as T3 | yes | 345–357 |
| T5 nuisance priors | 92–101 | yes, on the amplitude, normalised | yes, on the width as written | 431–472 |
| T6 Park wall | 104–114, quadratic | −∞ for J < 0 (D2) | −∞ for J < 0 | 72, 160 |
| T7 abundances | 141–151 | per-nucleus priors normalised; sum-to-one in the likelihood | the same | 100, 175, 345–357 |
| T8 screening, density, Doppler, thickness, apparatus priors | 117–138 | yes, as T3 | yes | 345–357 |
| T9 model data | 155–169 | yes, `ModelDataChi2` | yes | 101, 176 |
| T10 prior covariance | 184–194 | yes, normalised, replaces diagonal priors | refused on a level parameter | 510–558 |
| T11 `gamma_decay` priors and wall | 173–181 | **no** | **no** | — |
| MINUIT limits of any parameter | MIGRAD stays inside | **no** (never read) | **no** | — |

At the three points, sampler − (−½ χ²) (a constant where every term agrees: the Gaussian priors' normalisation):

| project (basis) | start | displaced | outside a limit |
|---|---|---|---|
| `gamma_decay` record, no prior (Brune, Park) | 0 | 0 | b(6.7931) = −0.1: 0, finite (−41554.0) |
| the record, Σ b > 1 (b(6.7931) = 0.6) | | | 20161.77: the wall missing |
| the record, priors on Γ_γ, b(6.1763), remainder | 3.0617 | 0.8956 | Γ_γ = −0.01: 81.117, finite |
| T2 (Brune; many priors, prior covariance) | 0.881047 | 0.881047 | abundance −0.001: 0.881047, finite |
| T2, density 10¹⁵ (its lower limit 2×10¹⁵) | | | 0.384797, finite |
| 13N (Brune, Park, formal) | 0 | 0 | (no limited parameter free) |
| `mcmc_physical_norm` 13N_fast (Brune) | 2.076794 | 2.076794 | |
| `energy_shift_sqrt`, `thm_coherent_background/toy`, `identical_pp_res` (Park) | 0 | 0 | |
| `nuisance_prior`'s model (b) (Brune) | 3.686232 | 3.686232 | |

The same numbers in both bases. So: every term of the fit is in the sampler, with the same value up to the Gaussian
normalisation, except T11 (the record's priors and wall, a difference that moves with the record's numbers) and the
limits, which the sampler never applied to any parameter. Outside a limit the fit and the sampler also differ where a
term reads a clamped value: the density prior is evaluated by the fit on the value held inside its limits, by the
sampler on the sampled value (T2 above, 0.496 of −½ χ²).

**The rule.** The sampler's log-probability is −½ × the fit's objective, taken from the same object and the same
enumeration of terms, plus a documented constant, and −∞ outside the support:

    ln P(x) = −½ Σ_terms χ²_t(p(x)) + C,          C = Σ_priors −ln(σ √(2π)) − ½ ln det(2π C_prior)

Each term of `Objective::Terms()` says how the sampler takes it (`ObjectiveTerm::AddForSampler`,
`ObjectiveTerm::ExcludesForSampler`):
- by value, in the evaluation (the default of the interface, so that a new term reaches the sampler with no edit of
  `AZURECalcMCMC.cpp`): T1 and T2 (the data, THM norms profiled), the abundances' sum-to-one, T9, T11;
- as normalised Gaussian priors from their enumeration, in the log-prior (`ForEachDiagonalPrior`, the nuisance
  priors, the prior covariance), because the posterior should be a normalised density in those parameters and the
  log-prior is analytic (evaluated before the model): T3, T4, T5, T7's per-nucleus priors, T8, T10;
- as a support: T6, Park's J > 0 (D2: −∞ in the sampler, a quadratic wall in the fit).

The limits: the sampler gives −∞ to a point outside the MINUIT limits of a parameter whose kind's row says so
(`ParamKindRow::samplerSupport`, `include/ParameterKind.h`). Today that is the `gamma_decay` record's numbers (Γ_γ ≥ 0,
b_f ∈ [0, 1]; δ has none): there the limits define the parametrisation (gamma_parametrisation.md section 2.6). For the
other kinds the flag is off, so the posterior of a project without a record does not change (below).

The columns of `samples.mcmc`: `logprior` is the normalised Gaussian priors (the enumerated ones and the user's);
`loglikelihood` is `logprob − logprior`, −½ × every term taken by value, so it holds T11 for a project with a record.

**What remains the sampler's own**, and why none is a term of the fit:
- the user's priors of the `<mcmc>` block (`SetPriors`): a configuration of the sampler, which the fit does not read;
- the normalisation constants: they do not depend on the parameters (the fit minimises, the sampler needs a density
  only up to a constant, but the `logprior` column is a normalised log density, as documented in the user guide);
- the profiled norms: not parameters in either (T1a/b; the same function);
- the observed-parameter basis: the sampler samples the observed widths (or the record's numbers) and evaluates the
  objective at the amplitudes they give; there is no Jacobian. The prior is flat in the observed widths, which is the
  choice of that basis, not a transformation of the fit's flat prior in the amplitudes;
- a nuisance prior on a width in the observed basis is a Gaussian on the width as written, in the fit a Gaussian on the
  converted amplitude (on purpose: each basis puts the prior on its own coordinate);
- the walkers' initial spread and the level-energy spread.

**Differences left for the director** (the posterior of a project without a record would change; not done here):
1. The limits of every other kind (abundances ≥ 0, thickness ≥ 0, densities and Doppler temperatures inside their
   layout range, apparatus parameters, a reduced width's limit, `<parameterSettings>` limits): the sampler accepts
   points outside them (abundance −0.001 in T2: finite). Proposal: `samplerSupport = true` on those rows (a one-column
   change), after a decision on the observed basis, where a level parameter's limit is set on the amplitude and the
   sampler moves the width.
2. Outside a limit, the fit's terms that read a clamped value (the density prior, T8) and the sampler's prior on the
   sampled value differ (0.496 in the example); this disappears with item 1.
