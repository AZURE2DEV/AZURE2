# shift_sqrt_subpoints

Sub-points move rigidly with their data point under an energy-dependent shift.

A segment's energy shift is a calibration of the measured beam energy.  `Segment::UpdatePointEnergyWithShift`
moved each sub-point of a target-integration / resolution grid by the shift at the sub-point's OWN energy,
which for a `sqrtshift` or a `shiftfn` (never for a constant shift) had two effects (audit item 11, findings 3
and 4 of 05_experimental_effects_chi2.md, fixed 2026-10-26):

- the integration window of the no-spread target-integration branch stretched by `1 + ds/dE`, a yield factor
  (1.9 % for b = 0.02 at 0.28 MeV in a flat region) that a fit could trade against b, and inconsistent with the
  beam-spread / straggling branch, which clips the window to its physical width;
- the resonance-anchored fine grid, placed by `WorkingCopy::ReplaceSubPointGrids` with the data point's shift s(E0),
  sat off a narrow resonance a depth d inside the target by s' x d (b = 0.02, d = 38 keV: 0.7 keV, more than the
  +-0.2 keV lattice of a 10 eV resonance).

Every sub-point now moves by the data point's own shift (c.m.), so the window keeps its width and the fine grid
is exactly on the shifted resonance.

Since 2026-11-13 the shift is the true beam energy (`include/ShiftedFold.h`, the rule; `tests/shift_true_energy`): the
fold of a shifted point is the fold of the same point without a shift at E + s, so the target's stopping, thickness
and window are those of the shifted energy, not translated from the nominal one.  On these models (a 107 keV layer,
the stopping formula falling with energy, shifts of 10-11 keV) every yield moved by -0.89 to -0.96 %: the window's
1/eps at the resonance is eps(E_r - s)/eps(E_r) larger before (the director's measurement, 2026-11-12: 1.00898 at
272 keV, 1.00974 at 320 keV).  Checks (a), (b), (e) compare two shifted models and (c), (f) a model with itself
refined: their numbers did not move (0, 0, 3.10e-5, 0, 1.78e-4).

## Models (`make_models.py` (deleted; see "How the models were made"), from `tests/thick_target_resonance`: 14N(p,g) across the 259.7 keV resonance,
Gamma_p 0.975 keV, 3.2e19 at/cm2, four yield points at lab 272-320 keV, no external capture)

- `narrow_sqrt_nores`: `narrow_sqrt` with the proton width of the 1/2+ level at 7.5567 MeV set to 0 (check (g)).
- `sqrt`: one four-point segment with `sqrtshift 0.02 0 0`;  `const`: four one-point segments, each with the
  constant shift 0.02 sqrt(E_k) of its own point -- the same calibration point by point;  `shiftfn`: the same
  shift as `shiftfn "a0*sqrt(x)"`;
  `noshift`: b = 0;  `narrow_sqrt` / `narrow_const`: `sqrt` / `const` with the resonance's proton width
  100x smaller (9.7 eV).

## Checks (`check.sh path/to/AZURE2`)

- (a) `sqrt` == `const` (1e-9): the sqrt shift is the point-by-point constant calibration;
- (b) `shiftfn` == `sqrt` (1e-12);
- (c) `sqrt` == `sqrt` under `AZURE_REFINE_QUADRATURE=8` to 9.3e-5, three times the 3.1e-5 measured (2026-11-12);
- (d) `sqrt` differs from `noshift` (the shift is active);
- (e) `narrow_sqrt` == `narrow_const` (1e-9): the 10 eV resonance is on the fine grid after the shift;
- (f) narrow resonances are covered by the mesh: `narrow_sqrt` == itself refined eightfold to 5.4e-4, three times the
  1.78e-4 measured (2026-11-12; eightfold against 64-fold 3.2e-6);
- (g) the analytic plateau (`plateau.py`): `narrow_sqrt` minus `narrow_sqrt_nores` (the same model with the 1/2+
  level's proton width at zero: the non-resonant part, 0.07-0.15 % of the yield) equals

  ```
  sigma_r = (lambda_r^2/2) omega Gamma_gamma (m+M)/(M eps(E_r,lab) N) x [arctan(2(E-E_r)/Gamma_p) - arctan(2(E_back-E_r)/Gamma_p)]/pi
  ```

  at each point's true c.m. energy E, to 5e-4 (measured 2.06e-4, 8.3e-5 with the quadrature refined eightfold; the
  engine before the rule: 9.96e-3).  It replaces, as a grid-free check, the uniform
  grid removed with `useAdaptiveGrid`.  The strength is omega Gamma_gamma Gamma_p / (the sum of the particle widths),
  here omega Gamma_gamma, and the Lorentzian's width Gamma_p alone: the radiative widths do not enter the level matrix
  (`docs/source/theory/conventions.rst`), so the engine's resonance is that of a level whose total width holds the
  particle widths only.  The textbook omega gamma = omega Gamma_p Gamma_gamma / (Gamma_p + Gamma_gamma) is 0.39 % lower
  here (Gamma_gamma/Gamma_p = 3.9e-3).  With a constant stopping and no shift, refined 64-fold, the director measured
  1.1e-4 (2026-11-12); before the rule above the shifted model was 1.3 % off (0.9 % of it the stopping at the nominal
  energy).

Until 2026-11-12, (c) compared `sqrt` with `dense`, the same model on a uniform 2000-point sub-point grid
(`useAdaptiveGrid=0`), as a grid-free check.  The key and its uniform grid were removed: on `narrow_sqrt` that grid
gave 2.06e-8 for the first yield against the mesh's 2.18022e-6 (refined 2.17983e-6), and 5.09e-7 at 64 times the
points.  `dense` went with it; the mesh against itself refined is the check now, and (f) guards the narrow
resonance (with the level anchors taken out of the mesh plan, (c) fails at 0.68 and (f) at 0.98).

dev 81abfd4 (before the fix): (a) 7.9e-4 off (the window stretch, diluted by the resonance's dominance of the
yield), (e) 2.1 (213 %) off: the fine grid missed the shifted narrow resonance; (b), (c), (d) pass.  Fixed build:
(a) 5e-11, (c) 2.8e-5, (e) 5e-11.

## How the models were made

The committed models are the reference.  They were written by `make_models.py` in the old project format and converted
to format 3 with `tools/format3_tree_convert.sh` (`AZURE2 --convert ... --to 3 --exact-grids`); the generator
was deleted with the old format (it is in the repository's history: `git show f56e02a9:tests/shift_sqrt_subpoints/make_models.py`).  To change a
model, edit its records.  The generator's description of what it wrote:

```text
Models for the energy-dependent-shift sub-point test (audit item 11, 2026-10-26).

All from tests/thick_target_resonance (14N(p,g) across the 0.97 keV-wide 259.7 keV resonance, 3.2e19 at/cm2
target, no external capture), the x20 target-integration line, four yield points at lab 272-320 keV:

  sqrt    : one segment, sqrtshift b = 0.02 MeV^1/2 (fixed), constant shift 0
  const   : four one-point segments, each with the CONSTANT shift b sqrt(E_k) of its own point -- the same
            calibration point by point, so the yields must be identical if sub-points move rigidly with
            their data point
  shiftfn : one segment, shiftfn "a0*sqrt(x)" with a0 = b (the same shift through the other code path)
  noshift : the sqrt model with b = 0 (to show the shift matters)
  narrow_sqrt, narrow_const : sqrt and const with the resonance's proton width 100x smaller (9.7 eV), where
            the fine grid's position matters: the old engine's grid missed the shifted resonance

Run from this directory: python3 make_models.py  (writes models/ and models/data/).
```
