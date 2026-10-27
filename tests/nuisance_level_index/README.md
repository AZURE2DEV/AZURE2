# nuisance_level_index

A nuisance prior on a level parameter penalizes that parameter.

`<parameterSettings>` can mark a level energy or width as a nuisance parameter, and the fit adds
`((p - nominal)/sigma)^2` to chi-squared (`AZURECalc::CalculateNuisanceChiSquared`; `AddNuisanceGradient`
adds the matching gradient term).  The `ParameterLimitsManager` indexes the settings by the parameter's
position among the *non-fixed* parameters, while the vector those functions receive is Minuit's *full*
vector, fixed parameters included.  Up to dev 89fa22f both read `p[nonFixedIndex]`: a prior on any level
parameter preceded by a fixed level parameter was evaluated on a different parameter, and MIGRAD then
pulled that other parameter toward the nominal value (audit 2026-10-25, item 14).  The least-squares
drivers (`--use-lm`, GSL) build their penalty rows by name and were right; `AZURECalc::Chi2On`, which
they and the gradient check use for the cost value, shared the MIGRAD mistake.  Both functions now index
the full vector with the parameter's actual position.

`check.sh` uses `tests/13N/13N.azr` (energy_1 fixed, width_1_1 free, energy_3 free; its settings block names
every free level parameter with nuisance 0) and edits one settings line per case, all in mode 1 (the run
log's "Total Chi-Squared" includes the penalty; `chiSquared.out` does not):

- (a) Level 3 Energy as a nuisance parameter with nominal `E_3 + 0.02` and sigma 0.01: the log total is the
  13N total (`tests/13N/expected/chiSquared.out`) + 4;
- (b) Level 1 Channel 1 Width as a nuisance parameter at its own value (1.62 eV, sigma 0.162 eV): penalty 0
  (the first free parameter, preceded by the fixed energy_1);
- (c) as (b) for Level 3 Energy at its own value with `--no-brune` (the nominal is converted to the formal
  energy): the total equals the `--no-brune` calculation without the prior;
- (d) a fit of model (a) with `--use-gradient` and `AZURE_GRAD_CHECK=exit`: the analytic gradient, penalty
  term included, agrees with finite differences to 1e-5 of the largest component (with a nuisance prior the
  analytic Gauss-Newton seed is off and the check runs after MIGRAD's numerical seeding steps, so it tests the
  agreement, not the start value);
- (e) as (b) for the gamma width of level 3 (channel 6, -0.48 eV, sigma 0.048 eV): penalty 0, the physical
  nominal being converted to the reduced-width amplitude it corresponds to.

dev 89fa22f: (a) +12,460 instead of +4 (the prior was evaluated on width_1_1 = 1.2719 in reduced units),
(b) +56 (on energy_1 = 0), (c) +12,018 and (e) +273,466 by the same mechanism.  The 2026-09-13 observation in the
azure2-eval skill (nuisance entries on ANCs and gamma widths at their nominal values adding 316, 3,306 and
3.3e9 to the 12C+alpha total) was this index mistake, not the physical-to-reduced conversion.
