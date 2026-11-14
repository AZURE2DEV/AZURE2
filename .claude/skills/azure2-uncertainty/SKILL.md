---
name: azure2-uncertainty
description: Uncertainties of an AZURE2 evaluation - Fisher information (J^T J, sloppy directions, information per data set), parameter covariance (HESSE, LM), profile likelihood and MINOS, posterior sampling with NUTS on the gradient and emcee, convergence (R-hat, ESS), bands from draws and from the covariance, data covariance files (writer, checks, the two chi2 splits), prior covariance, systematic axes and model averaging. Use when a result needs an error. Fitting: azure2-fit; whether to trust the fit: azure2-check; a review of a finished budget: uncertainty-audit.
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

**The tool.** `scripts/write_covariance.py OUT --file data/a.dat --file data/b.dat [terms]` writes the file from
components: the error column (`diagonal`, on by default), `--extra-diag F`, `--common F` (one term over every
row), `--common-per-file F`, `--shape PATH` (N fractions), `--matrix PATH` (a full N × N matrix, absolute units);
the fraction terms are model-scaled `relvector`s unless `--absolute`. It then checks the file as
`--check COV` does: row counts against the data files, the block part B and C positive definite (B must be: the
program whitens with it), the eigenvalue range, √C_ii against the error column (the program warns above 1 %).
Attach it with `Project.set(name, covariance="cov.txt")` on every data set it covers.

**Two χ² splits, one meaning.** With a covariance joining data sets, `Evaluation.data[i].chi2` (the program's
split, r_i (M⁻¹r)_i summed over the set) and `sum(r[d.rows]**2)` of `Objective.residuals` (the whitened rows) differ
set by set (48.8874 against 48.8867 on the self-test's two joined T3 sets) and agree in sum; only the sum over the
joined sets has a meaning. A per-set χ²/N of joined sets is not a test of either set.

## 2. Fisher information and the parameter covariance

With r the objective's rows (`r, J = f.residuals_and_jacobian(x)`, whitened data rows, then one row per prior;
sum r² = objective) the Fisher information is F = JᵀJ; the Gauss-Newton Hessian of the objective is 2F, so F⁻¹ is
the covariance at Δχ² = 1. `scripts/fisher.py PROJECT --x-json fit.json` prints and saves:

- the spectrum of the scaled F (D^−½ F D^−½, D = diag F): its condition number and the **sloppy directions**, the
  eigenvectors of the smallest eigenvalues by their largest components (on T3: energy_1 and width_1_1 of the 2⁺
  level, eigenvalue 0.20);
- per parameter its error √(F⁻¹)_jj and the **share of its information from each data set**,
  F_i,jj / F_jj with F_i = J_iᵀ J_i over the rows of set i (and from the priors): which data set constrains which
  parameter; with `--json` also the errors with one set left out;
- **what a better measurement adds**: `--improve SET=k` gives the errors with the errors of SET divided by k
  (F + (k² − 1) F_SET); a new measurement at other energies needs its points in the project (a data set of
  expected values) and the same script;
- `--compare-engine lm` (or `migrad-gradient`): the engine's errors over √(F⁻¹)_jj. At the same point F⁻¹ is the
  engine's LM covariance (JᵀJ + priors)⁻¹ to 6e-15 (self-test); MIGRAD's HESSE includes the second-order term and
  differs where the model is non-linear over the errors.

The covariance is the linear (Gaussian) approximation at x̂. It fails near a limit (a width at 0, θ² at 1), for
saturated background poles (`azure2-check`), and for widths whose sign is not determined: profile or sample.

**CLI.** A fit with the band on (`--covariance-band`) writes `output/covariance.dat` (free level energies and
amplitudes, `param.sav` order); `--scale-covariance` multiplies it by χ²/ν when that exceeds 1 (Birge ratio);
MINOS writes asymmetric errors to `param.errors` (`docs/source/reference/output_files.rst`).
`pyazr.bands.load_covariance(path)` reads `covariance.dat`; `s.fit_with_engine(method).covariance` is the same
from Python.

## 3. Profile likelihood and MINOS

For one parameter: hold it at a grid of values, refit all others, and read the interval where
Δ objective ≤ 1 (68 %). `scripts/profile_scan.py PROJECT --parameter NAME --points 9 --span 3 --minos` does it
(least squares on the rows, held column in `skip`, warm starts outward from the fit, crossing interpolated in
√Δ) and compares with MINOS (`fit_with_engine("migrad", minos=True)`): on the self-test's energy_4 the profile gives
−0.000569 / +0.000545 MeV against MINOS's −0.000568 / +0.000546. For a derived quantity (S(0), a strength) profile
by adding it as a datum with a tight error at each grid value (model data, `azure2-setup`) or sample. With
χ²/ν > 1, scale Δ by it or say that you did not. Use the profile when the surface is not parabolic or the interval
touches a limit.

## 4. Sampling the posterior

The target is `f.log_prob(x)` = −f(x)/2, −∞ outside the fit's limits and where Park's J < 0 (the program's
sampler's rule), with the priors of the project; never `s.evaluate(x).chi2` (data only). Procedure and lessons:
`references/sampling.md`.

1. Start at the best fit (MAP): `lsq_fit.py` / `engine_fit.py` of `azure2-fit`.
2. `scripts/sample_nuts.py PROJECT --x-json fit.json --chains 4 --warmup 300 --draws 500 --out nuts.npz`: exact
   Hessian at the MAP, whitening, NUTS on the analytic gradient (`scripts/nuts.py`); stops on a non-positive
   eigenvalue (an open direction: a prior is needed).
3. `scripts/sample_emcee.py ... --out emcee.npz`: the ensemble sampler as a second opinion, with its 50 τ rule.
4. `scripts/mcmc_stats.py --npz nuts.npz`: split R-hat < 1.01, ESS ≥ 400 per quoted quantity, after folding the
   level-sign symmetry; no unexplained divergences.
5. Bands and derived quantities from the draws (section 5).

One engine per worker process; keep workers × `OMP_NUM_THREADS` within the cores you hold. The program's own
ensemble sampler (MCMC tab): `docs/source/user_guide/mcmc.rst`.

**Memory in long chains.** A sampler that changes the radius per draw or transforms widths per step must stay flat
in memory; `s.cache_stats()` shows whether input files are re-read per rebuild (regression test:
`tests/pyazr/coulomb_memo_transform_test.py`). Recycle workers (`multiprocessing.Pool(..., maxtasksperchild=N)`)
on builds that grow; builds before 2026-11-04 (expires: 2027-02-01) grew by 100 kB per new radius.

## 5. Bands on observables

Linear propagation: σ_y² = gᵀ C_p g, with g_j = ∂y/∂p_j the sensitivity of each point. From draws (a posterior,
or N(x̂, C)): `scripts/draw_bands.py PROJECT --x-json fit.json --draws nuts.npz` (or `--laplace 400`) gives the
16/50/84 % quantiles per point beside the linear band; on the self-test they agree to 4 %.

```python
from pyazr import Session, uncertainty_bands, extrapolation_bands
s = Session("model.azr")
bands = uncertainty_bands(s, keys=[3, 7])         # {data set number: Band}; params default to output/param.sav
lo, hi = bands[3].interval(1.0)                   # also .relative
bands = uncertainty_bands(s, keys=[5], params=x, covariance=C)   # a point and a matrix of your own (R-matrix
                                                  # parameters only, e.g. F^-1 of fisher.py)
ext = extrapolation_bands("model.azr", keys=[1])  # trimmed private model: cost scales with what you ask for
```

- Evaluate at the fit: `param.sav` by default, or `params=` (a project written by `save_fit` carries its fit).
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
added in a fit and replaced in MCMC (a note is printed). `s.objective()` includes the term, `s.evaluate(x).chi2`
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
- `references/data-covariance.md`, `references/sampling.md`.
- Scripts: `fisher.py`, `write_covariance.py`, `profile_scan.py`, `nuts.py`, `sample_nuts.py`, `sample_emcee.py`,
  `mcmc_stats.py`, `draw_bands.py`; `scripts/selftest.py` runs them on T3 (52 s at one thread).
- Examples: `pyazr/examples/uncertainty_band.py`, `sensitivities.py`, `fit_emcee.py`, `fit_zeus.py`,
  `reaction_rate.py`.
