# modified_sfactor

Data sets and extrapolations of `quantity modified_sfactor` (format 3; `docs/upgrade/data_kinds.md` section 4.1):
the modified S factor of heavy-ion fusion S*(E) = E σ(E) exp(2πη + gE), g the record field `shape` (1/MeV).

**Models** (`models/`, data in `models/data/`), all on the levels of `tests/15N_p_a` (15N(p,α)12C):
- `pa_extrap.azr`: four extrapolations on one c.m. grid, 0.1–1.6 MeV: σ, S, S* with g = 0.46 /MeV, S* with g = 0.
- `fit_sigma.azr`, `fit_sstar.azr`: a fit of `data/sigma.dat` (the model's σ times 1 + 0.03 sin(3.7 i), 5 % errors,
  from `truth.py data`) and of the same rows written as S* (`data/sstar.dat`, value and error times
  E exp(2πη + 0.46 E) in mpmath).

**Truth** (`truth.py`, mpmath, no AZURE2 code): η = Z1 Z2 α sqrt(μc²/2E) with the masses `<levels>` prints and the
constants of `physics/Constants.h`.

**Checks** (`check.sh path/to/AZURE2`): S* against the σ block times E exp(2πη + gE) to 1e-10 (found 4.7e-11, the
output's 11 digits), its σ column equal to the σ block, the g = 0 block equal to the S block character for character;
the two fits with the same χ² and parameters within 1e-6 (found 3.4e-10); `--describe`; 1 against 2 threads; the
refusals (no `shape`, a `shape` in MeV, a neutral entrance).
