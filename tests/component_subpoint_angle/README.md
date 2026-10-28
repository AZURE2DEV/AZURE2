# component_subpoint_angle

The sub-points of a SUM/RATIO component sit at the component's own c.m. angle
(audit 2026-10-25, 01_kinematics_frames A3, checklist item 32).

Component points are copies of the base segment's points.
`ReadTargetEffectsFile` laid out a target effect's sub-points from the copied
point's c.m. angle, before `InitializeComponentSegments` converted the
component point's angle with its own pair, and never revisited them: a
component with a different exit pair (or a fixed-angle RATIO denominator,
whose raw lab angle was used) was integrated at the base pair's angle. The
sub-points are now put at the converted angle before their Legendre
polynomials are computed.

## Checks (`check.sh path/to/AZURE2`)

Model `tests/natural_target/models/M.azr` (p + 12C and p + 13C in one
project), test segments at 90 deg lab, output in the lab frame
(`kinflag.dat` = 1):

- P: 1 -> 1 and 2 -> 2 alone; Q: the SUM of 1 -> 1 and the 2 -> 2 component,
  no target effect; Qc: the same SUM with a beam-resolution convolution of
  sigma = 0.1 keV (a sub-point effect that changes nothing but the sub-point
  bookkeeping);
- (a) Qc equals Q to 1e-6: the second term's sub-points are at the p + 13C
  angle (94.4 deg), not the p + 12C one (94.8 deg);
- (b) Q equals P1 + P2 to 1e-5 (2e-6 observed): the lab SUM is the sum of the terms' lab cross
  sections. (The default c.m. output of such a SUM is the lab sum taken to the
  base pair's c.m. frame, so there the second term carries J_base / J_term,
  1.0005 here; that is the design, not a defect.)

dev 510e6d1: (a) differs by 0.5 %.
