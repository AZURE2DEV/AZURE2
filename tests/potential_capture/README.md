# potential_capture

The interior direct capture of a potential background (`wave ... background interior ... capture interior`;
`docs/upgrade/potential_background.md` section 5).

**Model** (`models/capture_a40.azr.in`): 3He(α,γ)7Be as a pure two-body potential model in the engine. Pair 1 = 3He + 4He
(atomic masses 3.0160293 u and 4.0026032 u) with the Gaussian potential VD of Tursunov, Turakulov and Kadyrov, Phys.
Rev. C 97 (2018) 035802 (arXiv:1707.05948v3, Tables I and II): α = 0.15747 fm⁻², R_C = 3.095 fm, depths s1/2 67.5,
p1/2 81.815179, p3/2 83.589554, d3/2 66.0, d5/2 69.0 MeV, in both halves (r_m = 15 fm). The entrance s1/2, d3/2, d5/2
waves have the interior background (B_bg = 0) with interior capture; the p waves have their depths only. The final
states are the potential's own p3/2 and p1/2 bound states: the separation energy of pair 1 is the p3/2 binding energy
with these masses, 1.4548867432 MeV (so 7Be(3/2⁻) is at E_x = 0), the 1/2⁻ level at E_x = 0.4177109916 MeV, and their
ANCs are the potential's (4.0931604450 and 3.5000439504 fm⁻¹ᐟ²), all from `make_reference.py`. The entrance J^π groups
have one level each at 40 MeV with zero widths (so that the channels exist and nothing but the potential acts). E1
in the long-wavelength form (the engine's default). The VD p3/2 depth does not bind at the paper's 1.5866 MeV with
these masses (the paper uses ħ²/2u = 20.7343 MeV fm² and nuclear masses); `make_reference.py --paper` reproduces the
paper's B = 1.5866 MeV and C = 4.34 fm⁻¹ᐟ².

**Reference** (`data/reference.csv`, `make_reference.py`): the direct capture cross section of the same potential,
computed without the engine: scipy DOP853 for the scattering waves (0 to 15 fm, matched to mpmath's Coulomb functions)
and the bound states (shooting, mpmath's Whittaker function beyond 30 fm), the overlap integrals from 0 to infinity,
and the standard E1 formula (the paper's Eq. (8), checked term by term in the script's docstring). About 15 minutes.

**`models/fit_level.azr.in`**: the same with a 1/2⁺ level at E_x = 3 MeV on top of the s-wave background (its particle
and γ widths free) and the ground state's ANC free, fitted to made-up capture points (`data/gs_pseudo.dat`; not data,
only a χ² away from its minimum for the gradient check).

`check.sh` runs the identity at a_c = 3.5, 4.0, 4.87, 5.5 fm, the same without `capture interior`, the refusal of a
final level away from the potential's eigenvalue, the gradient check and the thread check.

The models are kept as `.azr.in` (copied by `check.sh` into its scratch directory) so that `tests/upgrade/
predictions_test.py`, which lists every `.azr` under `tests/`, and its committed table are not touched by this research
branch; renaming them to `.azr` needs `predictions_test.py --write` (lines added for these two projects only).
