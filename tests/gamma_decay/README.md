# gamma_decay — a level's radiative decay as Gamma_gamma, branchings and mixing ratios

The `gamma_decay` record of `<gammaDecays>` (docs/upgrade/gamma_parametrisation.md) against the same level written
with partial widths. `check.sh` says what each part checks; run it as `./check.sh path/to/AZURE2`.

* `models/partial.azr`: the 14N(p,g)15O model of `tests/thick_target_resonance` (external capture off), level 5 the
  259.7 keV resonance (1/2+, E_x = 7.5567 MeV) with its five radiative widths free, and five data sets: the capture
  cross sections to the ground state and to 6.7931, 6.1763, 5.183 and 6.8594 MeV at 18 lab energies across the
  resonance (`models/data/`, made by `models/make_data.py` from the model with 3 % errors and a fixed scatter).
* `models/record.azr`: the same project with level 5 written as a record, made by `models/make_record.py` from the
  note's definitions (Gamma_gamma = sum, b_f = Gamma_f / Gamma_gamma, delta of Krane and Steffen), not by the engine.
* `models/sign_partial.azr`: a 3/2+ level of 15O formed in one channel (l = 2, s = 1/2) decaying to 5.183 MeV by M1
  (1 eV) and E2 (-0.25 eV), the case a of the note's section 2.3: the record (delta = +0.5) gives the same angular
  distribution, a2/a0 = -0.99282 as Krane, Steffen and Wheeler's formula for delta = +0.5.

Errors this check was shown to catch (each injected once into `src/GammaDecayMap.cpp`, then removed): the sign of the
L+1 width in the map (the sign case and its a2/a0 fail), a wrong remainder (identity, parameters.out and the
gradient fail), no pull-back to the level energy or to a particle amplitude (the gradient rows fail), no wall on the
sum of the branchings (the sum ends at 1.2), the sign of the remainder prior's gradient (the gradient rows fail).
