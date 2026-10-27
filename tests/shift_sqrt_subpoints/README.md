# shift_sqrt_subpoints

Sub-points move rigidly with their data point under an energy-dependent shift.

A segment's energy shift is a calibration of the measured beam energy.  `ESegment::UpdatePointEnergyWithShift`
moved each sub-point of a target-integration / resolution grid by the shift at the sub-point's OWN energy,
which for a `sqrtshift` or a `shiftfn` (never for a constant shift) had two effects (audit item 11, findings 3
and 4 of 05_experimental_effects_chi2.md, fixed 2026-10-26):

- the integration window of the no-spread target-integration branch stretched by `1 + ds/dE`, a yield factor
  (1.9 % for b = 0.02 at 0.28 MeV in a flat region) that a fit could trade against b, and inconsistent with the
  beam-spread / straggling branch, which clips the window to its physical width;
- the resonance-anchored fine grid, placed by `EData::ReplaceSubPointGrids` with the data point's shift s(E0),
  sat off a narrow resonance a depth d inside the target by s' x d (b = 0.02, d = 38 keV: 0.7 keV, more than the
  +-0.2 keV lattice of a 10 eV resonance).

Every sub-point now moves by the data point's own shift (c.m.), so the window keeps its width and the fine grid
is exactly on the shifted resonance.

## Models (`make_models.py`, from `tests/thick_target_resonance`: 14N(p,g) across the 259.7 keV resonance,
Gamma_p 0.975 keV, 3.2e19 at/cm2, four yield points at lab 272-320 keV, no external capture)

- `sqrt`: one four-point segment with `sqrtshift 0.02 0 0`;  `const`: four one-point segments, each with the
  constant shift 0.02 sqrt(E_k) of its own point -- the same calibration point by point;  `shiftfn`: the same
  shift as `shiftfn "a0*sqrt(x)"`;  `dense`: `sqrt` on a uniform 2000-point grid (`useAdaptiveGrid=0`);
  `noshift`: b = 0;  `narrow_sqrt` / `narrow_const`: `sqrt` / `const` with the resonance's proton width
  100x smaller (9.7 eV).

## Checks (`check.sh path/to/AZURE2`)

- (a) `sqrt` == `const` (1e-9): the sqrt shift is the point-by-point constant calibration;
- (b) `shiftfn` == `sqrt` (1e-12);
- (c) `sqrt` == `dense` to 1 % (the adaptive grid resolves the shifted resonance);
- (d) `sqrt` differs from `noshift` (the shift is active);
- (e) `narrow_sqrt` == `narrow_const` (1e-9): the 10 eV resonance is on the fine grid after the shift.

dev 81abfd4 (before the fix): (a) 7.9e-4 off (the window stretch, diluted by the resonance's dominance of the
yield), (e) 2.1 (213 %) off: the fine grid missed the shifted narrow resonance; (b), (c), (d) pass.  Fixed build:
(a) 5e-11, (c) 2.8e-5, (e) 5e-11.
