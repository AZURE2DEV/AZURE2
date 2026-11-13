# Checking a target integration against an independent integral

Use this when a folded curve looks wrong. It separates a modelling error (wrong kernel, wrong N or ε) from an
engine error, by computing the same integral outside the engine.

## Recipe

1. Copy the project. Edit only the `<targetInt>` line in the copy (isActive, density, the Gaussian flag and σ_b,
   the target-integration flag) for each variant.
2. Bare curve: give the segment a fine dummy data file (0.2 keV spacing in E_lab across the region) with the
   effect switched off, and read the calculated values:
   ```python
   from pyazr.bands import active_indices
   bare = m.evaluate(m.params_rwa).data[active_indices(m.datasets)[key]].model
   ```
3. Reference integral: Y(E_0) = (1/N) ∫_{E_0−Δ}^{E_0} σ(E′)/ε(E′) dE′ with numpy on the bare grid, Δ from
   ∫ dE′/ε = N, ε from the line's stopping-power formula.
4. Compare with the engine's folded values at every data point.

## What a correct build gives

On a ¹²C(p,p) model with a 7.5 μg cm⁻² foil:

- engine against the numpy integral: ≤ 1.3 × 10⁻³ relative at every data point;
- density × 10⁻³ (Δ ≈ 3 eV): the bare curve to 2 × 10⁻⁴;
- 10 and 100 sub-points: agreement to 5 digits.

## Older builds

builds before 2026-09-22 (expires: 2027-02-01): a broad background pole could be the only level covering a
thick-target window, and the sub-point grid then collapsed to three points (start, E_R, end), weighting σ(E_R) by
half the window; one ¹⁴N(p,γ) yield came out 27 times too large at 264.3 keV. A model with few narrow levels near the
point is exposed. The numpy check above detects it.
