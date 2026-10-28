# parameter_settings_nan

`nan` or `inf` in a `<parameterSettings>` field must not end the run (audit
2026-10-25, 06_parameters_outputs_api C11, checklist item 30).

`ParameterLimitsManager::ReadParameterSettings` parses every numeric field with
its own `custom_stod`, which threw `std::invalid_argument` on `nan` and `inf`.
A fit error written from a singular covariance is exactly that, and nothing in
`AZUREMain` caught the exception, so the command-line run terminated before any
output. The parser accepts them now (any case, signed); a limit that is not
finite is reported and ignored.

## Checks (`check.sh path/to/AZURE2`)

- (a) `tests/13N` as it is reproduces `expected/chiSquared.out`;
- (b) the same with the fit_error field of one `<parameterSettings>` line set
  to `nan`: the run completes with the same total;
- (c) with `inf`: the same.

dev 2363ada: (b) and (c) terminate (`terminate called after throwing an
instance of 'std::invalid_argument'`) with no output.
