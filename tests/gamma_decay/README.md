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

* `models/sign_ad.azr`: the sign case with its angular distribution as data (`models/data/sign_ad.dat`, nine angles,
  3 %): the record's mixing ratio fitted from 0 must end at +0.5, and its gradient at 0 must not be 0.

* External capture on the record's channels (`check.sh` (h), note section 2.7): the twin is `models/partial.azr` with
  capture switched on (column 31) into the ground state and 6.1763 MeV by E1 (level 5's g_ext is -0.0049 and
  +0.114 MeV^1/2 against internal amplitudes of 0.0077 and 0.109) and into 6.7931 MeV by M1 (its g_ext is 0 for
  an s-wave level); the record is made from it by `make_record.py`, both are written by the check. Compared:
  every output file in both bases, also after moving level 5's proton width, its energy and the ground state's ANC;
  the record's widths against `TransformOut`'s; the gradient in every free parameter (level 5's energy and proton
  width, the ANCs of the ground state and of 6.1763 MeV, the record's numbers) at three points, the third with a
  second record on level 6 (8.7502 MeV), where Im g_ext is 1.3 % of the amplitude; MIGRAD with the gradient and LM
  from displaced starts against the twin's minimum and branchings; 1 against 3 threads; a ground-state width of
  3.85e-12 eV, below the floor 2P (Im g_ext)^2 / n = 6.4e-10 eV, handled as `TransformIn` handles it.
  `tests/pyazr/gamma_decay_test.py` part 4 reads and sets the record's parameters of this project in a session.

* The vector of observed parameters (`check.sh` (i), note section 2.8): a narrow resonance, the ¹⁴N(p,γ) model of
  `tests/thick_target_resonance` with Γ_p = 0.01 eV and a layer yield added, as twin and record: the copy of the
  compound that places the rate's breakpoints (mode 5) and the meshes' anchors (a fit with free norms) has the
  twin's widths (the grid's resonance line of `AZR_DEBUG_GRID`; before the fix the record's gave 0.870 eV of photon
  width for Γ_γ = 0.0385 eV), the rate is the twin's; a limit (900, 950) eV on level 5's proton width ends at 950 eV
  as the twin's, a limit on Γ_γ at its bound; the observed-parameter MCMC (mode 6, `no`) with only the record's
  numbers free: ln L moves with them (it was constant), `param.mcmc` holds the best walker's Γ_γ, and its
  −2 ln L is the χ² of a calculation at its parameters. `tests/pyazr/gamma_decay_test.py` part 5: Γ_γ × 1.2 set in
  the physical vector gives the χ² and the `update_segments` cross sections of the project written with it, and the
  twin's with its radiative widths × 1.2 (before: the χ² of the unchanged project).
* A final state's free energy (`check.sh` (j), note section 2.7): the external-capture record with a second particle
  pair (α + ¹¹C, threshold put at 5.0 MeV, in a 9/2⁺ group of its own with width 0), so that the 6.1763 MeV final
  state lies above a threshold and its energy is free; the gradient in every free parameter against central
  differences, that energy included.

Errors this check was shown to catch (each injected once into `src/GammaDecayMap.cpp`, then removed): the sign of the
L+1 width in the map (the sign case and its a2/a0 fail), a wrong remainder (identity, parameters.out and the
gradient fail), no pull-back to the level energy or to a particle amplitude (the gradient rows fail), no wall on the
sum of the branchings (the sum ends at 1.2), the sign of the remainder prior's gradient (the gradient rows fail), the derivative in a mixing ratio taken through
p/(2 Gamma) (0 at delta = 0: the gradient check at 0 and the LM fit from 0 fail), Minuit2's call count in fitStatus.out
(the 1-against-3-thread comparison fails). With external capture (h): no term in the final states' ANCs in the
pull-back (the gradient rows at all points and the LM fit fail); no term in the level's particle amplitudes (the
gradient rows fail); d gamma_int / d g taken as 1 (the third gradient point fails, 2.3e-3); gamma^ext from the
final state's set-up amplitudes instead of the current ones (gradient rows and both fits fail); gamma^ext not
subtracted (every comparison fails, and pyazr part 4's chi^2 against the twin); the energy's central difference
without gamma^ext (the energy rows fail, 8e-3).
For the observed-parameter paths and the final state's energy (each injected once, then removed): the record's
numbers not turned into partial widths by `WithPartialWidths` ((i): the rate's and the fit's grid widths fail);
`FillMnParams` taking the file's numbers instead of the compound's input numbers ((i): the MCMC's ln L spread and
`param.mcmc` fail; pyazr part 5: all six checks fail); a limit converted at 1.1 times its value in
`ConvertPhysicalLimitToReduced` ((i): the proton-width limit fails, both end at 990 eV); a limit on a record's number
not applied ((i): Γ_γ ends at 0.0385 eV); no term in the final state's energy ((j): rel. 0.12). A limit on a record's
number converted as if it were a width passes: the conversion's round trip gives back the record's number (section
2.8), which is the rule.
