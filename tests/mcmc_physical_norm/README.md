# mcmc_physical_norm

Physical-parameter MCMC (mode 6, "Use Reduced Width Amplitudes: no", the CLI and GUI default)
must apply the sampled normalizations.

Before 2026-10-26 `AZURECalcMCMC::CalculateLogLikelihoodPhysical` filled the compound and the
energy shifts from the sampled vector but never called `FillNormsFromParams`: the pooled data kept
its starting norms while the sampler moved the norm parameters under their priors.  The norm
posterior was the prior, the R-matrix posterior was taken at frozen norms, and the periodic
`param.mcmc` snapshot recorded the starting norm (its tail came from `FillMnParams` on a fresh
clone).  `AZUREAPI::CalculateChi2Physical` (pyazr `calculate_chi2`) had the same omission
(`tests/pyazr/physical_norm_test.py`).

## Model

`13N_fast.azr` (`make_model.py`): `tests/13N/13N.azr` without its capture segments (the external
capture integrals cost 0.7 s per evaluation), every `<targetInt>` line inactive, and the
normalization of the `meyer_84` elastic segment free with a 5 % systematic error.  314 elastic
points; the chain (4 walkers, 100 steps, 0.5 % spread) runs in two seconds.

## Checks

The chain is random (no seed control), so the checks are relations that hold for any chain:

- **(a)** `param.mcmc` carries the norm of the best walker of the last step of `samples.mcmc`.
- **(b)** the log-likelihood recorded for that walker is −χ²/2 of its parameters evaluated with
  the **sampled** norm (a mode-1 calculation from `param.mcmc` with that norm; `param.mcmc` holds
  7 digits, so the agreement is ~1e-5 of χ²), and not with the start norm.

dev dfc64e1 (before the fix): (a) `param.mcmc` norm 1.0000000 against a sampled 1.0038; (b) the
likelihood matched the start-norm χ² to 4e-6 and the sampled-norm χ² only to 2e-3.
Fixed build: (a) 3e-8, (b) 1e-5 against the sampled norm, 2e-3 against the start norm.

```
bash tests/mcmc_physical_norm/check.sh path/to/AZURE2
```
