# target_thickness_exact

The thin-target energy loss is the exact range integral, not the first-order product.

Target integration computes `Y(E0) = (1/N) int_{E0-Delta}^{E0} sigma(E)/epsilon(E) dE`. Up to
dev 5ef20e5 `TargetEffect::TargetThickness` returned `Delta = epsilon(E0) N`. A stopping cross
section that rises through the target (every ion above the Bragg peak) then puts the back edge too
shallow, and since `1/epsilon` sits in the integrand a flat cross section came out low by
`(1/N) int dE/epsilon < 1`: 2 % for a 5 % target, 8 % for a 20 % one, growing towards low energy,
a shape distortion of a fixed-target excitation function (audit 2026-10-25, tier 3 entry 2).
`TargetThickness` now solves `int_{E0-Delta}^{E0} dE/epsilon = N` (Newton from the first-order
estimate, Simpson on 16 intervals; exact for a constant epsilon, so every suite with a constant
stopping cross section is unchanged).

## Checks (`check.sh path/to/AZURE2`)

Model from `tests/thick_target_resonance` (p + 14N, external capture off) with two data segments of
the total-capture cross section: segment 1 at lab 320-500 keV under one `<targetInt>` line with
`epsilon = a0 / x^0.8` (18 eV/1e15 at 360 keV) and `N = 4e18` atoms/cm^2 (a 20 % target at 360 keV),
no beam spread; segment 2 the bare cross section on a 0.1 keV lab grid, no target effect. The truth
is segment 2 integrated by hand over the exact window, whose back edge has the closed form
`E_exit = (E0^1.8 - 1.8 a0 N)^(1/1.8)` for the power law.

- (a) every yield within 0.5 % of the truth (1 % for the 320 keV point, whose window holds the
  259.7 keV c.m. resonance): the new build agrees to 0.04 % everywhere;
- (b) the first-order window gives 4-65 % less at the four background points (65 % at 360 keV,
  where the exact window reaches into the resonance tail), so (a) has teeth.

dev 5ef20e5: (a) fails at the four background points, by exactly the first-order deficit.
