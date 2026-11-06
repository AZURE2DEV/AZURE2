# straggling_bohr

The straggling coefficient is a centre-of-mass quantity, and a coefficient of 0 means Bohr straggling.

`TargetEffect::StragglingSigma` gives the width `c_s sqrt(dE [keV]) keV` with the energy loss taken in
the c.m. (surface minus depth, both c.m.) and the result used as a c.m. width, with no frame factor,
while the beam sigma on the same line is converted lab to c.m. and the stopping cross section is a
lab function. The docs never said which frame `c_s` belonged to (audit 2026-10-25, tier 3 entry 1).
It is now documented as c.m.: nothing in the engine changes for an explicit coefficient
(`tests/thick_target_resonance`, c = 0.8, pins that), and a lab coefficient is entered as
`c_lab sqrt(m2/(m1+m2))`.

The default coefficient was 0.04 keV^1/2, several times below Bohr's value (0.3-0.5 for light ions in
light targets), so a line that kept it had essentially no straggling. Straggling on with a coefficient
of 0 (the new default) now means Bohr straggling, `Omega_lab^2 = 4 pi Z1^2 e^4 Z2 N x`, with `N x` the
areal density of target atoms traversed (`int dE / epsilon` accumulated over the sub-points), the
charges those of the entrance pair (`TargetEffect::SetBohrCharges`), and the c.m. width
`m2/(m1+m2) Omega_lab` (`TargetEffect::BohrStragglingSigma`). The unused `Straggling` class
(a placeholder "ERYA methodology" never instantiated) is gone.

## Checks (`check.sh path/to/AZURE2`)

Models from `tests/thick_target_resonance` (p + 14N, the 259.7 keV c.m. resonance) with one
`<targetInt>` line for both segments: 0.1 keV beam sigma, a constant stopping cross section
1.8e-20 MeV cm^2 and 1.1e18 atoms/cm^2 (a 20 keV target), eight points on the rising edge of the
thick-target yield (lab 276-300 keV). With a constant epsilon Bohr's width is exactly
`c_B sqrt(dE_cm)` with `c_B = sqrt(1000 f 4 pi Z1^2 Z2 e^4 / epsilon)` = 0.3075 keV^1/2.

- (a) straggling off, (b) on with coefficient 0, (c) on with `c_B` explicit: (b) equals (c) at every
  point to 1e-6 and differs from (a) by more than 1e-3 somewhere on the edge;
- (d) the 13N thick target (density 0, p + 12C integrated down to 1 keV) with Bohr straggling: every
  yield finite and positive.

dev c6ed589: (b) has no straggling, equal to (a), and up to 90 % off (c) (the lab 300 keV point, whose target back edge sits on the resonance).
