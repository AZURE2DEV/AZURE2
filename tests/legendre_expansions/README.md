# legendre_expansions

The five normalisations of `quantity legendre` (field `expansion`: `a_over_a0`, `absolute`, `rsl`, `2l2`, `4pi`;
`docs/upgrade/data_kinds.md` section 4.4, `experiment/LegendreExpansion.h`), for data sets and extrapolations.

**Models** (`models/`, data in `models/data/`), the levels of `tests/angular_coefficients` (13C(α,n)16O[6.13], two
positive-parity levels) and a 5/2⁻ level at 11.05 MeV, so that orders 0–4 are non-zero:
- `extrap.azr`: the coefficients L = 0…6 in the five expansions at E_lab = 5.9, 6.1, 6.3 MeV, and dσ/dΩ at the 24
  Gauss–Legendre nodes in cos θ_cm.
- `data.azr`: the same rows in the five expansions (`data/<expansion>.dat`, `truth.py data`: the projection times
  1 + 0.02 sin(1.3 i), 3 % errors).

**Truth** (`truth.py`): A_L = (2L+1)/2 Σ_i w_i dσ/dΩ(x_i) P_L(x_i) of the engine's own dσ/dΩ (exact for the model's
degree-8 polynomial in cos θ), each expansion formed from it.

**Checks** (`check.sh path/to/AZURE2`): the CLI's `.acoeff` against the projection to 1e-8 (found 4.3e-10, the
output's 11 digits); pyazr's values to 1e-8 (4.3e-13) and the expansions against each other to 1e-12 (5.1e-16); five
equal χ² for the same rows; 1 against 2 threads; the analytic gradient against finite differences; `--describe`;
refusals (rsl in b/sr, a fitted norm on rsl).
