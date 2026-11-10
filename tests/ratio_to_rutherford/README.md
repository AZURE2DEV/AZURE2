# ratio_to_rutherford

Data sets and extrapolations of `quantity ratio_to_rutherford` and `ratio_to_mott` (format 3; `docs/upgrade/
data_kinds.md` section 4.3, step S14 of `data_model.md`): the model is dσ/dΩ (c.m.) over the engine's own
Rutherford or Mott cross section at the same point, each folded by the data set's apparatus.

**Models** (`models/`, data in `models/data/`):
- `aC_hs.azr`: α + ¹²C, a = 5 fm, every partial wave l ≤ 4 pure hard-sphere scattering (the one level has no width):
  dσ/dΩ and the ratio on a c.m. grid, the ratio at lab angles.
- `aa_hs.azr`: α + α (identical spin-0 bosons), a = 4.5 fm, every allowed l ≤ 6 hard-sphere: dσ/dΩ and the ratio to
  Mott.
- `aC_data.azr`: ratio data (`data/aC_rth_cm.dat`, from `truth.py data`: the hard-sphere ratio times
  1 + 0.02 sin(7.3 i), 3 % errors) on the hard-sphere model; `tests/pyazr/ratio_to_rutherford_test.py` makes its
  variants (multiplied by σ_R as `dsigma_domega`, a 5 keV shift, a 0.2 MeV beam spread).
- `aC_story.azr`: the worked example of `project_format3.md` section 9.2 (lab frame, five angles in one file, 3 %
  common normalization; three ¹⁶O levels) with synthetic data (`data/aC_rth.dat`: the model times
  1 + 0.03 sin(3.1 i), divided by 1.02, 4 % errors).  User page: `docs/source/examples/ratio_to_rutherford.rst`.

**Truth** (`truth.py`, mpmath, no AZURE2 code): f = f_C + f_HS from Coulomb functions at the channel radius; for
identical bosons f_C(θ) + f_C(π − θ) + 2 f_HS with even l.  `expected/aC_hs.dat`, `aa_hs.dat` are its values on the
engine's grids (`truth.py extrap`).

**Checks** (`check.sh path/to/AZURE2`): the hard-sphere values against the truth to 1e-9 (found: 3.5e-10 with the
Coulomb table, 2.8e-10 without, the output's 11 digits); the example's `--describe` and `--print-predictions` as
committed (`expected/story_*.txt`); 1 vs 2 threads identical; the CLI gradient check with a beam spread; the
identical-particle refusal.
