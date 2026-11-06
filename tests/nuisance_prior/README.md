# nuisance_prior — a nuisance prior lands on the parameter it names

A `<parameterSettings>` row with the nuisance flag adds
`((p - nominal)/error)^2` to chi-squared for the parameter the row names.
The chi-squared function, its gradient and the least-squares drivers receive
the **full** parameter vector (fixed parameters included, as
`CNuc::FillCompoundFromParams` reads it), while `ParameterLimitsManager`
numbers the parameters by their position among the non-fixed ones.
`AZURECalc::CalculateNuisanceChiSquared` and `AddNuisanceGradient` used to
index the full vector with that non-fixed index, so as soon as one fixed
parameter preceded the prior-carrying one the prior was evaluated on another
parameter.  The prior is now looked up by the full-vector index
(`ParameterLimitsManager::NuisancePrior`), in those two functions and in the
Levenberg-Marquardt / GSL setup (`PrepareFreeParams`).

`check.sh` runs `tests/7Li_p_a` in Calculate mode with a prior on
`Level 2 Energy (MeV)` = `energy_2` = 27.494 MeV, a free parameter preceded by
the fixed zero widths `width_1_6`, `width_1_7`:

| case | prior | expected total chi2 |
|---|---|---|
| none | — | chi2_0 |
| (a) | 27.494 +- 0.01 | chi2_0 (rel 1e-12) |
| (b) | 27.504 +- 0.01 | chi2_0 + 1 (abs 1e-6) |
| (c) | (b), with an external `param.par` that marks `energy_1` fixed | chi2 of that run without the prior + 1 |

Case (c) covers a parameter list whose fixed flags differ from a freshly
filled one: the limits manager is applied to the parameters as `param.par`
left them, so a non-fixed index recomputed from `CNuc::FillMnParams` points one
parameter off.

The chi-squared is read from the `Total Chi-Squared:` line of the run's
stdout (printed with 12 significant digits), which includes the priors;
`chiSquared.out` holds the data and normalization terms only.

Before the fix (HEAD 431477f's `AZURECalc.cpp`):

    chi2: no prior 2137.82733857   (a) 7561338.18734   (b) 7566837.98734

(a) is `((0 - 27.494)/0.01)^2 = 7.559e6` — the prior read `width_1_6 = 0`.
With only the direct-index fix, (a) and (b) pass and (c) gives 4.86e6, the
prior read against `width_1_5`.

    ./tests/nuisance_prior/check.sh path/to/AZURE2
