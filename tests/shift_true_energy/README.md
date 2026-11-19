# shift_true_energy

An energy shift is a calibration of the beam energy, so **the fold of a point measured at the nominal lab energy E_n
with the shift s is the fold of the same point without a shift at the beam energy E_n + s(E_n)** (the rule of
`include/ShiftedFold.h`; `docs/upgrade/energy_grid.md` section 2.5; the PI, 2026-11-12).  For the target kernel this is
not a rigid translation of the kernel: the back edge of a layer solves int_{E_out}^{E_0} dE/eps = N from the true E_0,
an infinitely thick target reaches from E_0 down to 1 keV, the straggling grows along that window, and the stopping is
read at the nodes' own energies.  Until 2026-11-13 the engine translated the kernel rigidly: the stopping, the target
thickness and the window stayed those of the nominal energy, an error of (d ln eps/d ln E) s/E on a yield on a
resonance (0.9 % on `tests/shift_sqrt_subpoints/models/narrow_sqrt.azr`).

## Models (`make_models.py BUILD`; the outputs are committed)

The physics of `tests/stopping_yield`: 15N(p,a)12C across one isolated 1- level of 16O at E_r = 0.9 MeV c.m.
(Gamma_p 60 eV, Gamma_a 40 eV), E_r/f = 0.9605 MeV lab; four nominal energies 0.9585, 0.9620, 0.9750, 0.9840 MeV
(below the resonance, its rising edge, the plateau, the back of the layer).

- `identity.azr`: 36 cases, every combination of
  - quantity and target: `thin` (quantity sigma on a layer: the target-averaged cross section), `layer` (quantity
    yield on a layer), `thick` (quantity yield on an infinite target);
  - stopping: `formula` (eps = 2.4e-20 x^-1.5 MeV cm^2 per active nucleus), `srim` (`stopping evaluation`, SRIM-2013
    for p in TiN), `table` (`stopping table`, the same written by `azure2-stopping` from 1 keV to 1.2 MeV);
  - straggling `none` or `bohr`; beam `none` or a 0.8 keV lab spread;

  each with a shift, the kinds in turn: an offset of 3 keV with a prior (a varied shift, so the planner's reach is not
  0), b sqrt(E) with b = 0.003 MeV^1/2, the function a0 + a1 x (1 keV + 0.002 E); and next to each its twin
  `<case>_twin` without a shift whose data file holds E_n + s(E_n) to 17 digits.  Layers of 1e18 active N atoms/cm^2
  (about 24 keV).
- `shifted.azr`, `twins.azr`: the same data sets split in two projects, so the two members of a pair are on two
  meshes and share no node.
- `fit.azr`: a scan of 39 energies (0.950 to 1.045 MeV lab) of a 72 keV layer (3e18 atoms/cm^2, the formula
  stopping) with a 1.5 keV beam spread, across the rising edge, the plateau and the back edge; the shift an offset fitted
  with the prior 0 +- 3 keV from 0; the data the engine's own yields at the shift 2 keV with 1 % errors.  At the back
  edge the window's width (the true energy's thickness) enters the yield's slope in the shift.

The kernels whose shape depends on the energy beyond a translation (`make_kernel_models.py BUILD`, the outputs
committed; the owner's answer of 2026-11-18, `docs/upgrade/OWNER_ANSWERS.md` section 8): the same pairs, each with an
offset with a prior, b sqrt(E) or a0 + a1 x in turn, against their twins in the same project.

- `kernel_gauss.azr`: 15N(p,a)12C as above with a Gaussian beam spread sigma(E_lab) = a0 x^8 (2 keV at 0.96 MeV), alone
  (`gauss`: its 5 sigma support is that of the true energy) and in a layer's kernel (`gausslayer`), and a free-gas
  Doppler broadening at 1e6 K alone (`doppler`: its window in velocity grows as sqrt(E)); five energies 0.948-0.984 MeV.
- `kernel_profile.azr`: 12C(p,g)13N across the 1/2+ level (`tests/13N`'s levels), a beam profile (one skewed Gaussian,
  xi 0.47 MeV, omega 15 keV, alpha -1.5) with a 3 keV detector resolution inside each point's +-4 keV window and the
  detailed-balance weight; the twin's window is the shifted point's window moved by s(E_n).
- `kernel_udr.azr`: 17O(n,a)14C (`tests/udr_low_energy`'s levels) with a UDR kernel (its 0.2 ns delay, a 125 ns burst, a
  channel of 7.55 bins per decade whose width in time is t0(E) ln 10 / 15.1) at 0.10-0.30 MeV.
- `fit_gauss.azr`, `fit_profile.azr`, `fit_udr.azr`: a scan of 16 energies of each kernel, the offset fitted with the
  prior 0 +- 3 keV from 0, the data the engine's values at the shift 2 keV with 1 % errors.

## Checks (`check.sh path/to/AZURE2`)

| check | tolerance | measured | what limits it |
|---|---|---|---|
| 1. one project: every shifted data set = its twin | 1e-9 | 0 in all 11 printed digits, 36 pairs | the twin's energy E_n + s written to 17 digits against the engine's sum; the same nodes |
| 2. two projects: the same pairs | 3e-6 | 1.04e-6 (beam spread, SRIM and table); 2.5e-8 (beam spread, formula); 6.5e-8 (thick target); 0 (layers without a spread) | the two meshes' nodes: the beam spread's Gaussian averaged on two different depth grids |
| 3. the gradient of the fitted offset (Jacobian path, a central difference of the fold on the fixed nodes) = a central difference of chi^2 (`AZURE_GRAD_CHECK`) | 1e-5 | 3.9e-7; fit_gauss 2.5e-7, fit_profile 5.5e-10, fit_udr 6.7e-9 | the two differences' steps and rounding |
| 4. one project per kernel model: every shifted data set = its twin | 1e-9 | 0 in all printed digits, 15 pairs | as check 1 |

What each check catches (each error injected once and removed, 2026-11-13):
- the engine before the rule (stopping and thickness at the nominal energy, the window translated) fails check 1 by
  2.3 to 27 % on the thin targets and layers (the most at the back-edge point) and by 0.18 to 0.47 % on the thick
  targets: (d ln eps/d ln E) s/E;
- the window's growth dropped (`Unit::SupportThicknessChange` held at 0, the thickness still at the true energy): check
  1 fails on six thin-target cases (8.4 % without straggling or spread, NaN where the moved support no longer holds the
  window); the layers pass, since their back edge is solved on the nodes' own stopping inside a 2 % margin;
- the shift column of the gradient reading the stopping translated from the nominal energy (the rule before) while
  chi^2 reads it at the true energy: check 3 fails, 5.3e-2.

- the engine before the owner's answer of 2026-11-18 (the supports of the kernels below translated from the nominal
  energy, the beam profile's detailed-balance norm and the UDR table at the nominal energy) fails check 4 on 13 of
  its 15 pairs: the Gaussian with sigma(E) by 9.5e-5 to 1.7e-4, the Doppler window 1.7e-8 (the function shift; the
  other two below 1e-9), the beam profile 2.8e-3 to 4.4e-3, the UDR kernel 3.9e-3 to 8.7e-3; the `gausslayer` pairs
  pass with it too (the target kernel's support reaches 25 sigma of the spread, beyond any tail the change moves).

A kernel author must not read the point's nominal energy (`KernelPoint::nominalLabEnergy`) where the true one is meant:
checks 1, 2 and 4 are what stop that.  The beam profile reads it only to move its reconstructed-energy window by the
point's shift (the window is stored unshifted).
