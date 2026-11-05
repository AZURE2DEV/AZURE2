# prior_covariance — a multivariate Gaussian prior on fit parameters

`priorCovariance=<file>` in `<config>` (path relative to the `.azr`'s directory) adds

    (p - mu)^T C^-1 (p - mu)

to chi-squared over the parameters the file lists, on their values as the fit's parameter vector holds them (the
`param.sav` space).  This is SAMMY's "PUP" (sequential Bayesian update): a fit starts from an earlier evaluation's
parameters and their correlations.  The file: optional `#parametrization <0|1|2>` (param.sav's tag; another basis
than the run's is refused), the parameter names, their means, the n rows of the symmetric covariance; `#` comments.
Listed parameters fixed in the run are marginalized (Note); the term is reported as
`Total-Prior-Covariance-Chi-Squared` in `chiSquared.out` (15 digits).  `pyazr.prior_from_fit` writes the file from a
fit's `param.sav` and `covariance.dat` (whose `#columns` header, new with this feature, names its columns).
Introduced 2026-11-05.

Model: the p+p model of `tests/norm_groups` (five levels; its data are the model's exact cross sections at 40 and 90
degrees times 0.8 / 0.9 with 1 % errors); every project is built by `check.sh` (no committed models).  A copy of the
90-degree file with its errors multiplied by 1e8 gives a normalization the data do not constrain.

`check.sh path/to/AZURE2` (about 10 s):

- (a) a diagonal prior on two norms = their normError priors (switched off): same data chi2 and the prior term = the
  normError terms (1e-10 against the exact value; the 6-digit Total-Norm-Chi-Squared to 1e-5), at moved norms and
  after an `--use-lm` fit;
- (b) a correlated 2x2 prior on a norm and a reduced width: the term = d^T C^-1 d computed by hand (1e-10);
- (c) `AZURE_GRAD_CHECK` with every level parameter free and a correlated 3x3 prior on a norm, a width and an energy
  shift: analytic = finite differences (1e-6 of the largest), Jacobian path = adjoint (1e-10; the shift column
  1e-8 of its own), and (with prior) - (without) = 2 C^-1 (p - mu) on the three parameters, 0 elsewhere;
- (d) the unconstrained norm: a fit returns the prior mean, and with a correlated prior (rho 0.95) the conditional
  mean mu_B + C_AB/C_AA (A - mu_A) (0.024 away from mu_B), to 1e-4, `--use-gradient` and `--use-lm`;
- (e) refusals with an ERROR, a non-zero exit and an empty `output/`: unknown name, non-symmetric, not positive
  definite, size mismatches (means, rows, row length), a `#parametrization` of another basis, a missing file; the
  run's own tag in param.sav's form is accepted; a listed parameter fixed in the run is marginalized (Note; the term
  is the 1-parameter marginal, 1e-10);
- (f) `prior_from_fit` round trip: fit A (`--use-lm --covariance-band`) -> prior file (names = covariance.dat's
  `#columns`, means = param.sav, the tag) -> fit B from the `.azr`'s start with the prior and the same data returns
  A's widths (1e-6 relative with `--use-lm`, 1e-2 of each prior sigma with MIGRAD), the same chi2, a prior term ~0;
- (g) MCMC with RWA parameters: `samples.mcmc`'s logprior column is the prior's Gaussian log density (1e-10), which
  replaces the norm's normError auto prior (noted); the chain's mean and spread are the prior's; a prior on a
  level parameter is refused in a physical-parameter chain.

`tests/pyazr/prior_covariance_test.py` covers the pyazr side (penalties, objective, residual_jacobian rows,
prior_from_fit, AzrModel.prior_covariance); `tests/gui/runtime_options_test` the GUI's config round trip.
