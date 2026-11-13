# Fitting from Python with pyazr

pyazr exposes the model, the data χ², the full objective and an analytic residual Jacobian; it does not minimize. The
minimizer is the caller's. This page lists the calls, what a hand-written fit must do to match the engine's own fit,
and the traps that make a Python fit stop or go wrong without an error. Behaviour reference:
`docs/source/user_guide/pyazr.rst` and `docs/source/user_guide/chi_squared.rst`. Running observables and the two
parameter conventions (`*_rwa` against plain calls) are in `azure2-calculate`.

Why no built-in minimizer: `m.objective().residuals_and_jacobian` returns the whole Jacobian for about two forward evaluations
(reverse-mode adjoint), which beats Minuit's numerical gradients on this problem; binding `MnMigrad`/`MnMinos` would
mean mirroring the engine's objective, parameter limits and Wigner bounds inside the API, two fits that could disagree.
The CLI (modes 2 and 4) remains the reference implementation.

## Inspecting the model

```python
print(m.level_scheme)          # pairs → J-groups → levels → channels (read-only)
print(m.parameters.table())    # idx, name, kind, free, J^pi, E, L, S, pair, rad, value, wigner
print(m.datasets.table())      # per data segment: file, reaction, observable, E range, norm err
print(m.extrapolations.table())
m.pairs                        # PairSet: masses, charges, spins, sep/excitation E, channel radius
m.parameters.free / .widths / .energies / .norms / .shifts
m.parameters.by_physical_level()   # {LevelKey: ParameterSet}, one entry per resonance or pole
m.physical_levels()                # the LevelKeys, ordered (jgroup, level)
m.datasets.sys_errors(vary_only=True)   # per-segment normalization systematics (fractional)
```

Name pitfalls (`energy_<N>` against `param.sav`, the engine pair number, the introspection `RuntimeError`):
`references/parameters-and-sav.md`.

## χ², residuals and the objective

```python
ev = m.evaluate(x)          # one evaluation: every number at x
ev.chi2                     # total DATA chi-squared (Total-Chi-Squared of chiSquared.out)
ev.objective                # what the engine's own fit minimizes: data chi2 + every prior and constraint
ev.terms                    # {name: value}: each term of the objective, by the engine's names
ev.data[i].chi2             # per data set (its row of chiSquared.out); .n, .norm, .norm_chi2, .shift_chi2
f = m.objective()           # the objective as one callable for a minimizer or sampler
val, g = f.value_and_grad(x)            # analytic gradient (energies, widths, norms)
r, J = f.residuals_and_jacobian(x)      # its least-squares rows: sum(r**2) == f(x)
r[ev.data[i].rows]                      # the rows of data set i (one per point)
```

- **Fit against `m.objective()`, not `ev.chi2`.** The latter is the data χ² only. The engine also adds, per
  data set, ((n − n₀)/(n₀ · normErr/100))² and, for a free energy shift, (ΔE/σ_ΔE)²
  (`references/normalizations-and-shifts.md`), and every other prior. Minimizing the bare χ² lets the normalizations
  absorb every discrepancy, a lower number the engine would never reach (−480 on one model).
- `f(x)` equals the sum of the `Total-*` terms of `chiSquared.out`; `tests/pyazr/objective_test.py` checks it against
  the engine.
- The rows of `f.residuals(x)` are the data rows (whitened where a data covariance joins points), the model data, the
  prior covariance, Park's wall, the gamma_decay priors, the abundance sums and one row per prior of a parameter: an LM
  or TRF fit on them sees the priors with their Jacobian.  The data part of the gradient alone is
  `2 * J[:n].T @ r[:n]` with `n = ev.data[-1].rows.stop`.
- `f.residuals_and_jacobian` raises for an analytically unsupported segment. Energy-shift, √E-shift, shift-function
  and density columns are forward differences: a self-contained segment's shift costs two evaluations of that
  segment's points; any other (composite, mapped, density) two full forward passes. `f.jacobian(x, skip=[cols])`
  leaves the listed free columns undifferenced and zero, for a fit that holds them. `AZURE_PROFILE=1` reports
  "shift cols (segment)" and "shift cols (full)".
- **Use `m.active_datasets`** (one entry per engine segment, in result order) to label per-segment results; `m.datasets`
  also lists inactive segments. Use `m.datasets.by_key(p.segment_key)` rather than `datasets[key - 1]`: a segment key
  counts inactive segments too.
- **Threads.** Every evaluation (`evaluate`, the objective, its gradient and rows) spreads its points over
  `OMP_NUM_THREADS` threads; results equal the serial path to a few ulp (`AZURE_SERIAL_CHI2=1` restores it). With
  several worker processes (rmfit shards, emcee pools) keep processes × threads within the slots requested: a worker
  inherits the job's `OMP_NUM_THREADS`.
- **Only the fit vector is an input.** pyazr's calls take the free vector `x` of the basis the session fits in; a
  physical vector is no input (`transform_rwa` turns `x` into physical values for a report).

## Free-vector helpers

`Parameter.free_index` is a parameter's position in the free vector; threading it by hand is the usual source of
silent misalignment.

```python
m.n_rmatrix                       # where the R-matrix block ends: x[:m.n_rmatrix]
w = m.parameters.widths
w.indices()                       # their free-vector positions
w.take(x)                         # their entries of x
w.put(x, 0.0)                     # a copy of x with every width zeroed
m.datasets.by_key(p.segment_key)  # the segment a norm or shift belongs to
```

## A subset fit with `scipy.optimize.least_squares`

Fitting a subset with everything else frozen is the usual evaluation step:

```python
from scipy.optimize import least_squares

fit = np.array([p.free_index for p in m.parameters.widths
                if not p.fixed and p.radiation_type in ("E", "M")], int)

f = m.objective()

def resid(z):
    v = np.array(best, float); v[fit] = z
    return f.residuals(v)          # every row: data, priors and constraints

def jac(z):
    v = np.array(best, float); v[fit] = z
    return f.jacobian(v)[:, fit]

sol = least_squares(resid, best[fit], jac=jac, method="lm",
                    xtol=1e-14, ftol=1e-14, max_nfev=2000)
x = np.array(best, float); x[fit] = sol.x
```

- Fit in rwa space (the Jacobian is exact there) and convert for reporting with `transform_rwa`. Sanity check:
  `transform_rwa(best[:nR])[p.free_index]` equals `p.value` for every width.
- Free the relevant normalizations too (`m.parameters.norms`) when a subset fit would otherwise be absorbed by them.
- **Respect the Wigner limits** if they matter: `m.wigner_widths()` gives the bound per channel. The engine's
  `ParameterLimitsManager` enforces them in its own fit; a Python fit does not unless they are given as bounds.
  `Parameter.wigner_limit` is γ²_W = 3ħ²/(2μa²) in MeV, the Teichmann–Wigner unit; `teichmann_wigner()` returns it
  unchanged and `lane_thomas_wigner()` gives the ħ²/(μa²) unit. builds before 2026-10-27 (expires: 2027-02-01):
  `teichmann_wigner()` multiplied by 1.5 twice, giving θ² too small by 1.5.
- For asymmetric (MINOS-style) errors, point `iminuit` at `f = m.objective()` and `f.grad`: `minos()` then runs
  against the same objective the CLI uses (`azure2-uncertainty`).
- Save the result with `m.save_fit("fitted.azr", x)` (`references/editing-levels.md`).
- MCMC: `pyazr/examples/fit_emcee.py` and `fit_zeus.py` show the pattern: one in-process engine per pool worker
  (constructed at module level), priors per `p.kind`, `log_prob = -0.5*(chi2 + offset) + log_prior` with
  `offset = Σ log(2π σ²)` (`azure2-uncertainty`).

## Traps that stop a Python fit silently

**A zero Jacobian column freezes `least_squares`.** `f.jacobian` returns an identically zero column for any
parameter no active data set sees: a γ width in a model with no capture data, a channel of a level nothing populates.
With `x_scale="jac"` that column's scale is zero, and TRF reports convergence after one function evaluation while the
gradient elsewhere is still in the hundreds. On one model this froze every normalization at exactly 1.0 and left
χ²/N = 5.6; dropping the dead columns first gave χ²/N = 1.01 from the same seeds.

```python
J = fitter.jacobian(x0)[:, cols]
cols = cols[np.max(np.abs(J), axis=0) > 0]      # before least_squares
```

A norm penalty of exactly 0.0 alongside free normalizations is the tell: `sum(d.norm_chi2 for d in m.evaluate(x).data)`.

**An amplitude that is exactly zero never moves.** χ² depends on a reduced-width amplitude γ through γ² near zero, so
its derivative there vanishes and no gradient method moves it; "free" in the `.azr` means nothing at 0. The common case
is a closed channel: levels just below a particle threshold carry that channel as a sub-threshold (ANC-type)
amplitude, and if the levels were given widths only in their open channels, that contribution is absent. Audit before
blaming the level scheme: list the zero entries of channels whose threshold lies just above the levels (in one model
all 14 n₁ entries of the levels between E_x = 9.36 MeV and the n₁ threshold were 0, while the (α,n₁) data were missed
right at threshold). To test them, seed in units of the Wigner amplitude (about 0.05 γ_W), never as a physical width;
touch only channels that are zero, and judge against a control with the same budget.

**`f.jacobian` raises, and the exception aborts the fit.** A trust-region step can put a reduced width where the
Coulomb functions overflow (`RuntimeError: z is not finite in log_Gamma`) or where the level matrix is singular. Catch
it and return a large residual with a zero Jacobian; the optimizer then shrinks its region, as it would after a finite
but large χ².

## A model built from hand-chosen seeds

A fresh `.azr` is far from any minimum, and one global fit over sixty parameters walks into a poor one. Free the
parameters in the order the data constrain them, warm-starting each stage, and run the sequence two or three times:

1. the particle widths of the pair that carries the dominant reaction;
2. the other particle pairs;
3. the γ widths, at fixed particle widths;
4. the background poles;
5. everything, including level energies and normalizations.

A staged driver should carry the penalty rows, the dead-column filter and the exception guard above. Interfering
amplitudes can still land on the wrong branch: `references/interference.md`.
