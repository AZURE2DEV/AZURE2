# potential_background — the interior R-function of a potential as the background of a channel

The project files are format 3. Design note: `docs/upgrade/potential_background.md` (a research option of the
side branch `feature/upgrade-9-potential-background`, off unless a `potential` record has a `wave ... background
interior` field).

System: ³He + ⁴He (pair 1 of `tests/scattering_length/he3he4.azr`, atomic masses 3.0160293 u and 4.0026032 u) with the
Gaussian potential VD of Tursunov, Turakulov and Kadyrov, Phys. Rev. C 97 (2018) 035802, Tables I and II (Dubovichenko's
potential): V₀ exp(−α r²), α = 0.15747 fm⁻², depths per partial wave s₁/₂ 67.5, p₁/₂ 81.815179, p₃/₂ 83.589554,
d₃/₂ 66.0, d₅/₂ 69.0, f₅/₂ 75.9, f₇/₂ 84.8 MeV, Coulomb of a uniform sphere R_C = 3.095 fm. Each partial wave has a
`wave` field with its depth (in both halves, interior and exterior) and, in the identity models, `background interior`.

| model | a_c (fm) | content |
|---|---|---|
| `identity_a35.azr` | 3.5 | background in all seven partial waves, B_bg = 0, one level without width per J^π (so that the channel exists) |
| `identity_a50.azr` | 5.0 | the same with B_bg = −1 |
| `exterior_only.azr` | 3.5 | the potential in the exterior only (no background) |
| `bound_level_a40.azr` | 4.0 | as identity_a35, plus the potential's own p₃/₂ bound state (one node) as a 3/2⁻ level at its energy and ANC (`data/phases.csv`: E_c.m. = −1.454887 MeV, C = 4.0932 fm⁻¹ᐟ²), the box state with one node subtracted from the p₃/₂ background (`subtract 1`) |
| `bound_level_a55.azr` | 5.5 | the same with B_bg = −1 in the p₃/₂ wave |
| `bound_level_nosub.azr` | 4.0 | the same level without the subtraction |
| `fit_level.azr` | 4.0 | s-wave background of the potential VDa (77.0 MeV, α = 0.180 fm⁻²; same paper), a free 1/2⁺ level on top, the VD s-wave phases as data (`models/data/swave.dat`, 0.5° errors) |

`data/phases.csv` is written by `make_phases.py` (scipy DOP853 and mpmath's Coulomb functions; not the engine).

Checks (`check.sh`, about 0.3 s):
1. With no levels the channel's U is the potential's (section 2.3 of the note): the phase shifts of both identity
   models equal the reference to 1e-4° (observed 6.8e-6°), and each other to 1e-4° (1.1e-6°), though a_c and B_bg
   differ. A wrong sign of R⁰ in the dressing, or the background taken from another depth than the exterior's, fails
   this by degrees.
2. The background acts: without it the s wave misses the reference by 7.8° (at 4 MeV).
3. A level on top (section 3 of the note, rule 1): the potential's bound state as a level, its pole subtracted from
   the background, restores the potential's phase shifts (9.3e-5° at 4.0 fm, 3.1e-5° at 5.5 fm), while its formal
   reduced width differs (1.374 and 0.566 MeV^½): the observed energy and ANC are radius-independent, the formal
   pole is not. Without the subtraction the level gets γ = 0: the dressed W̃(a_c) vanishes at a bound state the
   background already holds, so Brune's transformation gives the level nothing. This checks the dressed shift
   function below threshold, the ANC factor and the subtraction together.
   (The scipy reference binds this p₃/₂ state at 1.4549 MeV with the engine's atomic masses, not at the paper's
   1.5866 MeV; the engine's own physics/ gives the same root, 1.4549 MeV.)
4. The analytic gradient on the level on top of a background equals finite differences (8.6e-10): the background
   changes channel functions only, which the adjoint reads as constants.
5. 1 and 3 threads give identical files.

How the models were made: `work/u9-potential-background/scripts/he3he4.py` of the branch's research directory
(pyazr's Project editor on a copy of `he3he4.azr`); edit the records to change them.
