# shared_effect_pairs

A `<targetInt>` line shared by segments of different entrance pairs converts its lab energies with
each point's own pair.

Every segment listed on a `<targetInt>` line, and every component of a SUM/RATIO or natural-target
segment, points to the same `TargetEffect`. Up to dev 1599725 the Gaussian sigma (and a `sigmalab`
expression, the beam profile and, since 5ef20e5, the Bohr straggling width) was converted lab to
c.m. once, with the mass ratio `m2/(m1+m2)` of the first segment on the line (`ConvertSigmaToCM`,
`ConvertBeamProfileToCM`), so every other pair sharing the line got a kernel scaled by the wrong
ratio: 0.43 % on the width for n + 20Ne against n + 22Ne, tens of keV on a beam-profile centroid
(audit 2026-10-25, tier 3 entry 3, checklist t3-22). The line now keeps its lab values and each
point converts them with its own factor and charges (`EPoint::SetTargetKinematics`, set by
`EData::ReadTargetEffectsFile` for data points and component points alike).

## Checks (`check.sh path/to/AZURE2`)

Model from `tests/natural_neon` (n + 20Ne pair 2, n + 22Ne pair 4): two plain angle-integrated
elastic data segments at the same 41 lab energies, 0.85 to 1.25 MeV, and a 30 keV lab beam sigma.

- (a) one line `"1,2"` equals one line per segment on both segments to 1e-9;
- (b) a split model whose 22Ne line carries the sigma scaled by f(20Ne)/f(22Ne), which is what the
  old build computed for the shared line, differs from the correct split on the 22Ne segment by
  more than 1e-5 somewhere (4e-5 here), with the 20Ne segment untouched.

dev 1599725: (a) fails on the 22Ne segment by the (b) amount. `tests/natural_neon` (three isotopes'
components under one beam-spread line) moved accordingly and is re-baselined.
