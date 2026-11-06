# fixed_param_par -- param.par "fixed" flags under every minimizer

`check.sh` fits a one-level p+7Li -> a+a project with an external param.par
that fixes a width the .azr leaves free, under MIGRAD, `--use-lm` and
`--use-gsl-lm`. The fixed width must not move, the free parameters must, the
three minima must agree, and an LM fit with `--covariance-band` must save a
covariance over the two free R-matrix parameters only. Before the fix LM and
GSL-LM built their free set from the .azr alone: they fitted the fixed width
(chi2 9e-11 against MIGRAD's 12.59) and saved a 3x3 covariance.
