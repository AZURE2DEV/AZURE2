---
name: azure2-uncertainty
description: Uncertainties of an AZURE2 evaluation - parameter covariance from the Hessian or the least-squares Jacobian, MINOS and profile likelihood, MCMC sampling on the objective, bands on cross sections, S factors and reaction rates (pyazr.bands, --covariance-band), prior parameter covariance (sequential update), data covariance matrices for every data set (covariance keyword, relvector, Peelle's puzzle), systematic axes (radius, background, model choices) and their combination. Use when a result needs an error. Fitting is azure2-fit; whether the fit can be trusted at all is azure2-check.
---

# Uncertainties of an AZURE2 evaluation

An evaluation reports observables y (σ in b, S in MeV b, N_A⟨σv⟩ in cm³ mol⁻¹ s⁻¹) computed from fitted parameters
p (level energies E_λ in MeV, reduced-width amplitudes γ_λc in MeV^½, normalizations, energy shifts). Their
uncertainty has three parts: the statistical part carried by the data covariance into the parameter covariance C_p;
the data systematics (scales, energy scales) carried by nuisance parameters; and the model part (channel radius a_c,
background poles, level scheme, kernel choices) that the data do not determine. This skill gives the procedure for
each and for their combination. The objective itself is defined in `docs/source/user_guide/chi_squared.rst`.

## When to use, and when a sibling owns the task

- Use when a number or a curve needs an error band, when data are prepared for a fit (covariance matrices), or when
  a previous evaluation enters as a prior.
- First run `azure2-check`: an uncertainty on an untrustworthy fit is meaningless. Fitting and nuisance parameters:
  `azure2-fit`. Bands drawn in a figure: `azure2-plot`. THM model averaging: `azure2-thm`. A review of a finished
  error budget: the user-level `uncertainty-audit` skill; a referee-style read: `physics-referee`.

## 1. Data covariance: one matrix per data set, before the fit

No data set enters a fit without a covariance file, or a written reason in the fit's README why not (a set used only
for a plot; an alternative view of points already in the fit). What C_data holds:

- In C: everything that moves the points of a set relative to each other in a known pattern: statistics
  (`diagonal`, or `block` per spectrum fit), per-detector solid angle or efficiency (`relvector` per detector),
  efficiency shape (`relvector` per shape parameter), analysis variants (`relvector` per variant), extra scatter
  (second `diagonal`).
- Out of C: the overall scale of a set (a free normalization with a prior, or a `normgroup`), and energy-scale
  uncertainties (an energy shift with a prior). Putting the scale in C as well counts it twice.
- A multiplicative term is a model-scaled `relvector`, never a data-scaled `vector`: a data-scaled term pulls the fit
  low (Peelle's puzzle).

Engine side: `docs/source/user_guide/chi_squared.rst` ("Correlated data: covariance matrices"), file format in
`docs/source/reference/data_formats.rst`, outputs (`covarianceModel.out`) in `docs/source/reference/output_files.rst`.
The full recipe for own and literature data, the checks and the traps: `references/data-covariance.md`.

## 2. Parameter covariance at the minimum

- **CLI.** A fit with the band on (`--covariance-band`, or "y" at the band prompt) writes `output/covariance.dat`:
  the covariance of the free level energies and reduced-width amplitudes only (normalizations, shifts and other
  nuisances left out), in `param.sav` order and space, with a `#columns` header. From MIGRAD it is HESSE's matrix;
  from `--use-lm`/`--use-gsl-lm` it is (JᵀJ + priors)⁻¹. `--scale-covariance` multiplies it by χ²/ν when that exceeds 1
  (Birge ratio; never deflates). MINOS writes asymmetric errors to `param.errors`. Details:
  `docs/source/reference/output_files.rst`, `docs/source/reference/command_line.rst`.
- **pyazr.** `pyazr.bands.load_covariance(path)` reads it. With your own least-squares fit, C_p = (JᵀJ)⁻¹ from
  `m.residual_jacobian(x)` (whitened residuals; priors and penalties included when you fit `m.objective`).
- **Validity.** C_p is the linear (Gaussian) approximation at p̂. It fails for parameters near a bound, for saturated
  background poles (`azure2-check`), and for widths whose sign or branch is not determined. Then profile or sample.

## 3. Profile likelihood

For a quantity q (one parameter, or a derived number such as S(0) or a resonance strength): fix q at a grid of
values, refit all other parameters at each (warm-start from the neighbour), and read the interval where
Δχ² = χ²(q) − χ²_min ≤ 1 (68 %, one parameter of interest). With χ²/ν > 1, scale Δχ² by it or state that you did
not. Use the profile when the χ² surface is visibly non-parabolic or the interval touches a physical bound
(θ² = 1, Γ = 0).

## 4. Sampling on the objective

- Sample exp(−χ²/2) with χ² = `m.objective(x)` (data, normalization and shift penalties, priors), not
  `calculate_chi2_rwa` (data only). Worked samplers: `pyazr/examples/fit_emcee.py`, `fit_zeus.py`; the built-in
  sampler: `docs/source/user_guide/mcmc.rst`. One engine per worker process; keep workers × `OMP_NUM_THREADS` at most
  the cores you hold.
- Discard burn-in and check the integrated autocorrelation time (a chain shorter than about 50 autocorrelation times
  is not quoted) before using a posterior (`docs/source/user_guide/mcmc.rst`, "Analysing the chain").
- **Memory in long chains.** A sampler that calls `transform_rwa` per step (e.g. with a fixed total γ width) or that
  samples a_c (a new radius per sample) must stay flat in memory. Diagnose with `m.cache_stats()`:
  `input_files_read/served/held` show whether input files are re-read per rebuild (the engine keeps no Coulomb
  memo, so its counters are 0).
  Builds before 2026-10-28 (expires: 2027-02-01): `transform_rwa` kept every Coulomb function at the moving level
  energies, about 5,000 entries per call per worker (4 GB per worker after 1,400 steps of a 12C(α,γ) chain).
  Builds before 2026-11-04 (expires: 2027-02-01): every new radius added about 100 kB to the memo and every rebuild
  re-read the `.azr` and all data files. On such builds plan the memory reservation, recycle workers
  (`multiprocessing.Pool(..., maxtasksperchild=N)`), or write the radius into the `.azr` per sample. Regression test:
  `tests/pyazr/coulomb_memo_transform_test.py`.

## 5. Bands on observables

Linear propagation: σ_y² = gᵀ C_p g, with g_j = ∂y/∂p_j the sensitivity of each point.

```python
from pyazr import azure2, uncertainty_bands, extrapolation_bands, sensitivities
m = azure2("model.azr")
bands = uncertainty_bands(m, keys=[3, 7])        # {segment key: Band}; params default to output/param.sav
lo, hi = bands[3].interval(1.0)                   # also .relative
ext = extrapolation_bands("model.azr", keys=[1])  # trimmed private model: cost scales with what you ask for
```

- Evaluate at `param.sav`, not at the `.azr`: the `.azr` need not carry the fitted parameters.
- `method="analytic"` (one reverse-mode adjoint per point, about two forward evaluations for the whole matrix) is the
  default; `"central"` finite differences check it (`pyazr/examples/sensitivities.py --check`).
- `quantity="sfactor"` gives S-factor bands. Bands on the CLI's own grid: `--covariance-band` writes `.band` files,
  which `uncertainty_bands` reproduces to six decimals (`pyazr/examples/uncertainty_band.py`).
- Elasticities ∂ln y/∂ln p (`pyazr/examples/sensitivities.py`) say which parameters a data set constrains.
- **Reaction rates** have no band output. Propagate by recomputing N_A⟨σv⟩(T_9) for parameter draws (posterior
  samples, or draws from N(p̂, C_p)) and quote percentiles per temperature. `pyazr/examples/reaction_rate.py` computes
  the rate integral in Python for that purpose; its grid must resolve the Gamow window at the lowest T_9 and every
  narrow resonance (`azure2-calculate` owns rates).
- A band from C_p is statistical only (plus nuisances that were free in the fit). Systematic axes come next.

## 6. Prior parameter covariance (sequential update)

`priorCovariance=<file>` in `<config>` adds (p − μ)ᵀ C⁻¹ (p − μ) for the parameters the file lists (SAMMY's PUP). The
workflow: fit A with `--covariance-band`; `pyazr.prior_from_fit("output/param.sav", "output/covariance.dat",
"prior.dat")`; fit B with `priorCovariance=prior.dat`. The file's optional `#parametrization` tag must match the run;
fixed parameters are marginalized; a malformed file is an error. A per-parameter prior on the same parameter is
added in a fit and replaced in MCMC (a note is printed). `m.objective()` includes the term, `calculate_chi2_rwa()`
does not. Builds that predate the keyword ignore it without a message: check the log for the prior's note.
Reference: `docs/source/user_guide/chi_squared.rst`, "Prior parameter covariance".

## 7. Systematic axes and their combination

1. List the axes the data do not decide: a_c per partition, background-pole energy and number, level-scheme choices
   (J^π of a doubtful level, an extra level), kernel choices (target integration against Gaussian, stopping table),
   data-set inclusion (leave one set out), and for THM the vertex model.
2. For each axis refit a variant (complete fit, own directory, same data). Report each variant's χ², number of free
   parameters k, and the quantities of interest with their statistical errors.
3. Combine with `pyazr.model_average(variants, weights="aic" | "bic" | "chi2" | "flat", rescale=None | "best")`:
   mean m = Σ w_i m_i, statistical variance S² = Σ w_i s_i², model variance M² = Σ w_i (m_i − m)², total T² = S² + M².
   Weights compare models of the same data; with χ²/ν > 1 use `rescale="best"`, or Akaike weights put all the weight
   on one variant. Quote S and M separately. Module docstring: `pyazr/modelavg.py`.
4. A variant that changes the result by more than its statistical error on an axis with no physical preference is
   a systematic to report, not a variant to discard.

## Pointers

- `docs/source/user_guide/chi_squared.rst`, `docs/source/user_guide/mcmc.rst`, `docs/source/user_guide/pyazr.rst`,
  `docs/source/reference/output_files.rst` (`covariance.dat`, `param.errors`, `covarianceModel.out`).
- `references/data-covariance.md`.
- Examples: `pyazr/examples/uncertainty_band.py`, `sensitivities.py`, `fit_emcee.py`, `fit_zeus.py`,
  `reaction_rate.py`.
