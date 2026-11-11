# ratio_to_reference

Cross sections measured relative to a tabulated standard: `quantity ratio_to_reference` with
`reference <file> [unit <u>]` (format 3; `docs/upgrade/data_kinds.md` section 4.9), on the 15N(p,α)12C levels of
`tests/15N_p_a`.

**Models** (`models/`, data in `models/data/`):
- `extrap.azr`: σ, the ratio to `ref_cm.dat` (E_cm in MeV, σ in mb: 50 exp(−1/E)(1 + 0.3 sin 3E), 37 rows) and the
  ratio to `ref_one.dat` (1 b) on one c.m. grid.
- `fit_sigma.azr`, `fit_ratio.azr`: `sigma.dat` (the pseudo-data of `tests/modified_sfactor`) and the same rows over
  the table (`ratio.dat`, `truth.py data`).
- `folded.azr`: `sigma.dat` as σ and as the ratio to 1 b, both through a 10 keV lab beam spread.

**Truth** (`truth.py`, mpmath): the table read and interpolated linearly in ln σ against E_cm, independently of the
engine (which interpolates in E_lab, a linear map of E_cm).

**Checks** (`check.sh path/to/AZURE2`): the ratio against σ/σ_ref to 1e-10 (found 4.2e-11) and the ratio to 1 b
equal to σ; the two fits with the same χ² and parameters within 1e-6 (found 3.6e-8); the folded ratio to 1 b equal
to the folded σ; the analytic gradient with the fold (1.7e-6); `--describe`; 1 against 2 threads; four refusals.
