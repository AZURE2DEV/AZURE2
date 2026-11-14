# Fitting from Python with pyazr

pyazr gives the program's objective as one object, `f = s.objective()`, and the program's own fitters as one call,
`s.fit_with_engine(method)`. This page says what a hand-written fit on `f` must do to match the program's fit, and
the traps that stop a Python fit or send it wrong without an error. Behaviour reference:
`docs/source/user_guide/pyazr.rst` (Objective, Evaluation, the fit) and `docs/source/user_guide/chi_squared.rst`.
The tools are `scripts/lsq_fit.py` (least squares on the rows and Jacobian, staged, multiple starts) and
`scripts/engine_fit.py` (the engine's fitters, MINOS).

## Inspecting the model

```python
s = pyazr.Session("model.azr")
print(s.describe())               # pairs, levels, data sets, parameters, options in effect
print(s.parameters.table())       # one row per parameter: name, kind, unit, start, free, bounds, owner
f = s.objective()
f.names, f.x0, f.bounds           # the free vector: names, start, the fit's limits (-inf/inf where none)
print(s.evaluate(f.x0).table())   # objective, terms, one line per data set
```

## The objective and its rows

```python
ev = s.evaluate(x)                # every number at x, from one engine evaluation
ev.chi2                           # the DATA chi-squared (Total-Chi-Squared of chiSquared.out)
ev.objective                      # what the program's fit minimises: data + every prior and constraint
ev.terms                          # {name: value}, one entry per term of the objective
ev.data[i].chi2, .norm, .norm_chi2, .shift_chi2, .rows   # per data set
v, g = f.value_and_grad(x)        # analytic gradient
r, J = f.residuals_and_jacobian(x)   # rows with sum(r**2) == f(x), and dr/dx
```

- **Fit `f`, not `ev.chi2`.** The data χ² alone lets the normalizations absorb every discrepancy: a lower number
  the program's fit never reaches.
- The rows are the data rows (whitened where a data covariance joins points; data set i has `ev.data[i].rows`),
  then the model data, the prior covariance, Park's wall, the `gamma_decay` priors, the abundance sums and one row
  per nuisance prior. A least-squares fit on them sees the priors with their Jacobian. The data part of the
  gradient alone is `2 * J[:n].T @ r[:n]` with `n = ev.data[-1].rows.stop`.
- With a data covariance joining data sets, `ev.data[i].chi2` (the program's split) and `sum(r[d.rows]**2)` (the
  whitened rows' split) differ set by set; only their sum over the joined sets has a meaning.
- `f.residuals_and_jacobian` raises for a data set the engine has no analytic Jacobian for. Energy-shift,
  √E-shift, shift-function and density columns are finite differences (two passes of the segment, or of the whole
  model for a composite or mapped one); `skip=[cols]` leaves held columns undifferenced and zero in the data rows
  (a prior row keeps its derivative).
- `f.log_prob(x)` is −f(x)/2, −∞ outside the fit's limits and where Park's J < 0: the program's sampler's rule.
- **Only the fit vector is an input.** Every call takes the free vector of the basis the session fits in; physical
  values for a report come from `s.transform(x)` or the saved project.
- **Threads.** Evaluations spread over `OMP_NUM_THREADS`; with several worker processes keep processes × threads
  within the cores you hold. A session cannot be pickled: one session per worker process
  (`pyazr/examples/fit_emcee.py`).

## Traps that stop a Python fit silently

**A zero Jacobian column freezes `least_squares`.** A parameter no active data set sees (a γ width without capture
data, a channel nothing populates) has an identically zero column; with `x_scale="jac"` its scale is zero and TRF
reports convergence after one evaluation. On one model this froze every normalization at exactly 1.0 and left
χ²/N = 5.6; dropping the dead columns gave 1.01 from the same seeds. `lsq_fit.py` holds such columns and names
them. A norm penalty of exactly 0.0 beside free norms is the tell: `sum(d.norm_chi2 for d in ev.data)`.

**An amplitude that is exactly zero never moves.** χ² depends on γ through γ² near zero, so its derivative there
vanishes. The common case is a closed channel just above a level (all 14 n₁ entries of the levels below the n₁
threshold were 0 in one model while the (α,n₁) data were missed at threshold). Seed in units of the Wigner
amplitude (about 0.05 γ_W), touch only channels that are zero, and judge against a control with the same budget.

**An engine failure at a trial point aborts the fit.** A step can put a width where the Coulomb functions overflow
(`z is not finite in log_Gamma`) or the level matrix is singular; `pyazr.EngineError` is raised. Return large rows
and a zero Jacobian instead (`lsq_fit.py` does), so the trust region shrinks.

**One start is one branch.** Interfering amplitudes give several minima: four on `tests/upgrade/T3/T3_start.azr`
(755.22, 591.83, 580.33, 509.35), and the engine's LM and MIGRAD end on different ones from the same start. Run
`lsq_fit.py --starts 8 --jitter 3` (starts spread by three parabolic errors) and compare with
`engine_fit.py --method migrad-gradient` (`references/interference.md`).

**The Wigner limits are not bounds unless the project sets them.** `f.bounds` are the limits of the parameters
and of `<parameterSettings>`; the Qt GUI's "Populate Wigner Limits" wrote them there. Without them a fit can take
a width above the single-particle limit (`azure2-check` judges it).

## A model built from hand-chosen seeds

Free the parameters in the order the data constrain them, warm-starting each stage, and run the sequence two or
three times (`lsq_fit.py --free REGEX`, one call per stage, each from the previous stage's `--save-fit`):

1. the particle widths of the pair that carries the dominant reaction;
2. the other particle pairs;
3. the γ widths, at fixed particle widths;
4. the background poles;
5. everything, including level energies and normalizations (`--free 'norm$'` selects the data sets' norms).

Save each stage with `--save-fit` in the project's directory and start the next from it.
