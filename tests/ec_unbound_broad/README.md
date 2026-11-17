# ec_unbound_broad — capture into a broad unbound final state (the line-shape integration)

A capture data set or extrapolation with the field `line_window W` is integrated over the energy E_f of its unbound
final state within E_R +- W (`src/LineWindow.cpp`): the photon energy is E - E_x(E_f), the exterior of external
capture is the real scattering wave of the final channel at E_f, and the weight is the one-level line shape
rho(E_f)/N_f of the final level. Without the field a capture keeps the narrow-resonance treatment
(`tests/ec_unbound_final`).

**Models** (`models/`, dev's, written in format 3 by `AZURE2 --convert` with 0 notes): `a_*`, `b_*`: 12C + p with
capture pair 3 into the 13N 1/2+ state at E_x = 2.3655 MeV (E_R = 422 keV, Gamma_p = 39 keV as in Rolfs and Azuma,
NPA 227, 291 (1974); `b_*`: 39 eV), a = 4.48 fm, every other width 0, so pure E1 direct capture from hard-sphere
p waves; `d_*`, `g_*`, `h_file`: the same with a 5/2+ level at E_R = 1.60 MeV (Gamma_p 60 keV) that decays by E2
(0.05 eV) into the 2.3655 state (`d_ec*`: that width 0; `g_*`: external capture into pair 3 off; `g_off`: the pair's
energy 2 keV off the level; `h_file`: level and pair 20 keV higher). `_none` has no window, `_w1`, `_w3`, `_w20` one
of 1, 3 and 20 Gamma.

**Truth.**

- (a), (b) An independent calculation of the windowed over the unwindowed direct capture (damped radial integrals
  with Richardson extrapolation, mpmath/ODE Coulomb functions, a quadrature over the window; in the narrow limit the
  window fraction of the line plus the non-resonant continuum term): the numbers are in `check.sh`; 5e-5 (3e-4 for
  the 3 Gamma window).
- (c), (f) Internal capture alone is int rho/N_f ((E - E_f)/(E - E_R))^5 dE_f times the unwindowed one (scipy
  quad; the numbers in `check.sh`): 2e-5, with external capture into the pair on (as a difference of two runs) and
  off (the final level then found by J^pi and energy).
- (d), (f) The formal parametrization gives the Brune result (1e-9).
- (g) The line follows the final level's current energy: a parameter file that moves the level 20 keV gives the
  cross sections of the project with the level there (1e-6), and they differ from the unmoved ones by more than 2 %.
- (e) What is refused: `line_window` without a photon exit (`E_FIELD_NOT_APPLICABLE`); a final state that is not a
  level of the project (the run's ERROR).

**Origin.** dev's `tests/ec_unbound_broad` (2026-11-15). This engine's cross sections equal dev's on every model
(0 to 7e-11 relative; compared in a scratch directory, not a test).

**Run.** `./tests/ec_unbound_broad/check.sh build/src/AZURE2` (about 3.5 minutes: 48 final-energy nodes per point,
each with its direct radial integrals).
