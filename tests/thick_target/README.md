# thick_target

A thick target (target integration with a density of 0) is an explicit mode of the effect, not a
mutation of it.

A `<targetInt>` line with target integration and a density of 0 means an infinitely thick target:
the yield per incident particle, ∫₀^E₀ σ(E)/ε(E) dE, not divided by a density. Up to dev 5df7024
`EData::ReadTargetEffectsFile` implemented that by overwriting the **shared** effect's density with
1e24 at the first point whose thickness came out below 1e-10 (audit 2026-10-25, item 17). Two
things followed. That first point kept a target thickness of 0, and the Gaussian branch of
`EPoint::IntegrateTargetEffectImpl` integrates from E₀ − thickness to E₀, so with a beam spread or
straggling on the line its yield was exactly 0 (13N segment 1, first point: 0 instead of 2.85e-14).
And the branch also fired for a *thin* target whose stopping cross section evaluated to zero or
below at some point, silently turning the rest of the line into a thick target with the yields
divided by 1 instead of by N.

Now `TargetEffect::IsThickTarget()` is set from the density when the line is read, each
thick-target point integrates from its energy down to 1 keV (its thickness is E₀ − 1 keV, so the
spread and straggling branches cover the whole range), the yield is not divided by a density, and a
non-positive stopping cross section on a thin target is an error naming the segment and energy.

## Checks (`check.sh path/to/AZURE2`)

Models are made at run time from `tests/13N/13N.azr`, whose line `"1,2"` is a thick target without
spread:

- (a) the file as it is reproduces the 13N total;
- (b) a 2 keV beam spread on the thick line: every yield of segments 1 and 2 is positive and the
  first point of segment 1 is within 1 % of its value without the spread;
- (c) straggling (coefficient 0.04) on the thick line: the same;
- (d) a thin target (density 1e18) whose stopping cross section a₀x + a₁ = x − 0.5 is negative
  below 0.5 MeV stops with an ERROR and writes no `chiSquared.out`.

dev 5df7024: (b) and (c) give 0 for the first point, (d) runs and writes negative yields
(chi2 436,536).
