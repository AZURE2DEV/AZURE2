# Sampling the posterior of an AZURE2 fit: target, whitening, NUTS, emcee, convergence

The posterior of a fit is exp(−objective/2) over the free vector x (`Objective.log_prob`), the objective being the
program's: data χ² plus every prior of the project. A multichannel R-matrix posterior is hard to sample for two
separate reasons. The raw parameters are badly scaled (energies in MeV to 10⁻⁴, amplitudes, norms): a coordinate
problem, removed by whitening with the Hessian at the MAP. Some directions are not closed by the data at all
(background poles, levels at or above the data edge): a modelling problem, which only priors solve. The method
below is ported from the PI's user-level skill `rmatrix-posterior`, where it was worked out on a 202-parameter ⁷Be
evaluation; its numbers quoted here are from that study.

## The target

- `f.log_prob(x)` = −f(x)/2, and −∞ outside the fit's limits of every parameter kind (`Objective.bounds`: the
  limits a parameter is made with, then a `<parameterSettings>` block's) and where Park's overlap J of a level is
  negative. That is the program's own sampler's rule (`docs/source/user_guide/mcmc.rst`, "Limits").
- The gradient is −∇f/2 (`f.value_and_grad`). A point where the engine fails (`pyazr.EngineError`) must return
  −∞: in NUTS it ends the trajectory as a divergence; it must never raise or return a large finite penalty.
- The priors are part of the model. A project with only normalization priors can leave widths of poorly constrained
  levels free to wander: on ⁷Be, Newton from the adopted fit lowered the objective by 13 and stalled on a ridge
  with a negative Hessian eigenvalue, where the engine began to fail. Physicality priors (θ = γ/γ_W ~ N(0, 1) on
  physical levels, a broad scale such as 10 on background poles, a Park norm J ≥ 0.5) close those directions; they
  are written into the project as nuisance priors (`<parameterSettings>`) or a prior covariance, and stated in the
  paper.
- Temperature: T = 1 is the posterior of the quoted errors; T = χ²/ν is the Bayesian analogue of scaling the
  covariance by the Birge ratio. Quote which.

## Whitening at the MAP

x = x̂ + L z with L = V diag(λ^−½), where H/2 = V diag(λ) Vᵀ is the Hessian of −log p at the MAP x̂. In z the posterior
is close to N(0, I) wherever the Laplace approximation holds.

- **Exact Hessian, not Gauss-Newton, as the metric.** `sample_nuts.py --hessian exact` takes central differences of
  the analytic gradient (2n gradients). 2JᵀJ drops Σ rᵢ ∇²rᵢ, which is large where widths of poorly constrained
  levels enter non-linearly: on ⁷Be a 1σ step of the JᵀJ-whitened coordinates cost more than 1 in log p along 14 of
  202 axes (34 000 along the worst); with the exact Hessian 200 of 202 stayed within 0.25.
- **A non-positive eigenvalue means an open direction.** `sample_nuts.py` stops with exit 3 and names the
  eigenvalue: add or tighten a prior, or fix the parameter; do not floor it away.
- **The Laplace approximation is the start and the metric, not the answer** in many dimensions: importance
  sampling from it gave a Kish ESS of 1 out of 1500 draws on ⁷Be although every axis deviated little.

## NUTS (`scripts/nuts.py`, `scripts/sample_nuts.py`)

Multinomial NUTS (Betancourt 2017) with Stan's extra U-turn checks, dual-averaging step size (target acceptance
0.85) and Stan's windowed warm-up; in whitened coordinates the identity is already a good metric, and the
warm-up's metric update is `auto`: dense only when a window holds ten draws per dimension, diagonal otherwise. A
dense estimate from short windows is rank-deficient: in 202 dimensions it needed 914 gradients per draw against 15
with the diagonal update.

On the two-parameter T3 project of the self-test: exact Hessian, 2 chains × (200 + 400), 3.9 gradients per draw,
no divergence, R-hat 0.999, posterior sds 1.01 and 1.04 times the Laplace sds (ESS about 500). On the eight free
parameters of T3 (2 × (150 + 200)) about 13 gradients per draw, 63 s at one thread.

- Run chains as separate one-thread processes (`--chains 1 --seed k`), more than one; start them at draws of the
  Laplace approximation (the default) and some wide.
- **Fold the exact sign symmetry of each level before R-hat.** Changing the sign of every width of one level leaves
  every observable unchanged, so the posterior has mirror copies a chain cannot cross: on ⁷Be one chain sat in the
  copy of a background pole, giving R-hat up to 1.66 in 14 parameters; folding each draw to the sign of the level's
  largest width at the MAP brought it to 1.09. Fold by the level's widths (`s.parameters.table()` names each width's
  level) before `mcmc_stats.py`.
- `--max-depth 8` caps a trajectory at 255 gradients; hitting the cap often means the metric is wrong.

## emcee (`scripts/sample_emcee.py`)

The ensemble sampler on `f.log_prob`, walkers started at Laplace draws inside the limits. Report its
autocorrelation time with the 50 τ rule; in many dimensions do not judge it by its own τ: on ⁷Be emcee (416
walkers) reported τ = 10, 15, 29 stored steps after 0.95, 1.3 and 2.6 thousand steps, an estimate that grows with the
run, and after 4.4 h its widths were still 22 % (narrowest) to 94 % (median) of the NUTS widths. Compare it with a
NUTS run, parameter by parameter, before quoting it. It agrees on well-constrained quantities first (S(0) did on
⁷Be) and fails in the prior-dominated directions.

## Convergence (`scripts/mcmc_stats.py`)

- Split R-hat < 1.01 over at least four chains (after the sign fold); bulk ESS ≥ 400 for each quoted quantity.
- Divergences: none unexplained. A few at a hard limit (a width at 0, a Park wall) are expected; many in the bulk
  mean the step size or the metric is wrong.
- Gradients per draw: in whitened coordinates a few to a few tens; hundreds mean a bad metric.
- Effective samples per core-hour, NUTS against emcee: gradient calls × t_grad against log-p calls × t_logp,
  measured at fresh points on one thread.

## Bands and derived quantities from draws

Evaluate the quantity on thinned draws (`scripts/draw_bands.py --draws`), one evaluation per draw, and quote
percentiles. Where the posterior is close to Gaussian the band from draws equals the linear band
sqrt(gᵀ C g) (0.96 on the self-test's data set); where they differ, the draws are the band. Rates and
extrapolations are computed per draw on a session of the extrapolations (`azure2-calculate`).
