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
and the standard E1 formula (the paper's Eq. (8), checked term by term in the script's docstring). About 15 minutes
to 3 MeV and 45 more to 9 MeV (`python3 -I make_reference.py 3.5 4.0 ... 9.0` appends the rows' values; the file
holds both runs).

**`models/fit_level.azr.in`**: the same with a 1/2⁺ level at E_x = 3 MeV on top of the s-wave background (its particle
and γ widths free) and the ground state's ANC free, fitted to made-up capture points (`data/gs_pseudo.dat`; not data,
only a χ² away from its minimum for the gradient check).

`check.sh` (also the ctest `potential_capture`) runs the identity at a_c = 3.5, 4.0, 4.87, 5.5 fm for E_cm = 0.02-9
MeV, through the poles of the interior R-function (s wave at 7.59 MeV for 4.87 fm; s, d5/2, d3/2 at 5.14, 7.02, 7.53
MeV for 5.5 fm), the same without `capture interior`, a nonzero B_bg (s −0.05, d3/2 −0.27, d5/2 −0.25) and the energy
table off (`--no-ec-table`), which must all agree (the table divided by the dressed √P̃ and e^{−iφ̃} before
2026-11-20 and was off by up to a factor 28 there: `potential_background.md` section 5.8), the refusal of a final
level away from the potential's eigenvalue, the gradient check and the thread check. About 4.5 minutes.

The models are kept as `.azr.in` (copied by `check.sh` into its scratch directory) so that `tests/upgrade/
predictions_test.py`, which lists every `.azr` under `tests/`, and its committed table are not touched by this research
branch; renaming them to `.azr` needs `predictions_test.py --write` (lines added for these two projects only).

**`models/bound_only.azr.in`**: the model above with the p waves `only bound_state` (their potential gives the
bound-state shape of the interior capture only; potential_background.md section 5.9) and two p-wave levels with
widths (3/2⁻ at E_x = 8 MeV, Γ = 1 MeV; 1/2⁻ at 9 MeV, Γ = 2 MeV) and p-wave phase-shift extrapolations: `check.sh`
section 7 compares it with the same levels without any potential (identical p waves), with the potential in the p
waves (capture equal where that potential is negligible beyond a_c), the binding rule, and the record's round trips.

**`models/subtract_level.azr.in`** (written by `make_subtract_level.py`): section 9's model at 4.87 fm with the s-wave
box state of 2 nodes subtracted and put back as a 1/2⁺ level whose particle and photon amplitudes are the box state's
(potential_background.md section 5.11); check.sh section 10 compares it with the unsubtracted model.

**`models/e2m1.azr.in`**, **`data/reference_e2m1.csv`** (`make_reference_e2m1.py`): E2 and M1 (orbital) capture into
VD's bound 0d5/2 state through the interior capture of the s1/2 and d3/2 waves, the final level's ANC 10⁻³ so that its
own pole is negligible; check.sh section 11 (section 5.12).

**`models/bound_level.azr.in`** (written by `make_bound_level.py`, no engine needed): the model at 4.87 fm with the
s½ depth 104 MeV, where the s-wave potential binds its 2-node state by 0.126 MeV; that state subtracted and put back
as a ½⁺ level given by the bound state's observed energy, ANC and photon widths, all from scipy. `check.sh` section
17 compares its capture with the unsubtracted model's (potential_background.md section 14.9).
