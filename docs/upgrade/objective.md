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
