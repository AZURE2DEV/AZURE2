# shift_subpoint_angle — energy shifts must not move the sub-points of a target or resolution integration

Sub-points (the energy grid of a convolution, beam profile or target integration) are built at their parent
point's c.m. angle and carry no data of their own (`EData::ReadTargetEffectsFile`, `EPoint(angle, energy, parent)`
puts the angle in both the lab and the c.m. slot). `ESegment::UpdatePointEnergyWithShift`, which moves every
point of a segment with an energy shift, re-ran the parent's lab-to-c.m. angle and cross-section conversions on
each sub-point. That took the parent's c.m. angle for a lab angle and moved every sub-point by the full lab-c.m.
difference: 2-5 degrees for p + light nuclei, 18 degrees for alpha + 12C. It happened for any lab-frame
differential segment in a sub-point target effect with an active shift: nonzero, varied, `shiftgroup`,
`sqrtshift`, `shiftfn`. Shift derivatives (finite differences through the same path) were wrong with it.

`template.azr` is the 2026-10-25 audit model, p + 16O with the 17F ground-state and 0.495 MeV capture pairs,
and the Harris 16O(p,p) lab-frame differential data at 11 angles (`data/Harris_pp.dat`). `make_models.py` gives
every model a 0.1 keV Gaussian convolution (20 sub-points) on all segments:

| model | segments | shift |
|---|---|---|
| a0 | 16O(p,p), 11 angles | 0 |
| a1 | the same | 1e-7 MeV |
| av | the same | 0, varied |
| g0 | a0 plus 16O(p,g0) lab-frame differential at 45, 90, 135 deg (made-up data, `data/pg.dat`) | 0 |
| g1 | g0 | 1e-7 MeV on the capture segments |

`check.sh` runs each in calculate mode and checks:

| | check | tolerance |
|---|---|---|
| (a) | a1 equals a0 | 1e-4 (the shifted grid is rebuilt, so the discretisation moves chi2 by ~1e-5; the bug is 7e-1) |
| (b) | av equals a0 | 1e-6 |
| (c) | g1 equals g0 on the three capture segments (gamma branch of the same code) | 1e-5 |
| (d) | a0 against `expected/a0_chi2` (the unshifted value, which the fix does not change) | 1e-3 |

**What it caught (fixed 2026-10-26).** Before the fix, (a) gave 23,096.7 instead of 13,484.5 (the audit's
reproduction) and (c) moved the capture segments by 1.8e-4. The fix sets each sub-point's angle to the parent's
c.m. angle at the shifted energy and skips the conversions, which is what `ReadTargetEffectsFile` does when the
sub-points are first built. (b) passes on the old build too: a shift of exactly 0 does not enter the conversion.
