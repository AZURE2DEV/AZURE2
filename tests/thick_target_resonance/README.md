# thick_target_resonance — thick-target yield across a narrow resonance

The 14N(p,g)15O model of the paper's straggling example, reduced to four
Formicola et al. yield points (lab 272, 282, 284 and 320 keV) across the
0.97 keV-wide 259.7 keV resonance, with a 3.2e19 at/cm2 target, 0.1 keV beam
sigma, straggling on (c = 0.8 keV^1/2), 20 sub-points. External capture is
switched off (every ecMultMask is 0) so the run takes about a second; the grid
the regression lived in does not depend on it.

The same data are read twice, as two segments that differ only in the
targetInt resonanceWidthMultiplier:

* segment 1: 5, the example's own setting. The window of a point more than
  5 widths above the resonance is covered only by the 10 MeV background pole
  (the 30 MeV level). Between a4095a6 and the fix, GenerateGrid anchored to
  that pole's 200 keV lattice and crossed the 100 keV window in one step,
  picking up E_R on the way: the grid [start, E_R, end] put sigma(E_R) on
  half the window and the yield jumped x34 on the plateau (chi2 983975
  instead of 3.608). With only the anchoring fixed, the smooth step still
  jumped into the resonance's region and skipped its upper wing: the
  264.9 keV point was 36% high (chi2 68.47).
* segment 2: 20, the default a4095a6 found converged. The regression did not
  show here; it pins the converged answer.

Since 2026-10-22 the zone of nodes anchored to a narrow resonance goes on past
its lattice (resonanceWidthMultiplier widths) with graded steps, each 10 % of
the distance to the resonance, until the step reaches the coarse one
(AdaptiveIntegrationGrid, tests/reference/adaptive_grid_reference_test). Before,
the lattice edge met coarse steps placed from the window's top edge, so the tail
was integrated poorly, and differently at every resonance position. Both segments
now agree with the 1000-sub-point integral: segment 1 at multiplier 5 to 0.55 %
at the first point and 0.02 % elsewhere (it was 15 % high), segment 2 to 0.08 %
and 0.003 %. The expected chi2 moved accordingly, from 3.607 / 15.418 to
15.031 / 15.591. The old segment-1 value of 3.6 came from the under-resolved
first point lying closer to the data.

Yields (total capture, c.m. keV):

| E_cm     | seg 1 (x5) | seg 2 (x20) | 1000 sub-points | before 2026-10-22, seg 1 / seg 2 | pre-regression 7e533b7, seg 1 / seg 2 | broken dev, seg 1 |
|----------|------------|-------------|-----------------|----------------------------------|---------------------------------------|-------------------|
| 253.7846 | 3.8641e-08 | 3.8461e-08  | 3.843e-08       | 4.405e-08 / 3.851e-08            | 3.973e-08 / 3.848e-08                 | 4.347e-07         |
| 263.0573 | 2.0396e-06 | 2.0394e-06  | 2.0393e-06      | 2.0473e-06 / 2.0394e-06          | 2.0444e-06 / 2.0396e-06               | 2.588e-06         |
| 264.8950 | 2.0776e-06 | 2.0773e-06  | 2.0773e-06      | 2.0853e-06 / 2.0774e-06          | 2.0835e-06 / 2.0775e-06               | 7.346e-05         |
| 298.3287 | 2.1509e-06 | 2.1504e-06  | 2.1504e-06      | 2.1645e-06 / 2.1505e-06          | 2.1603e-06 / 2.1507e-06               | 7.391e-05         |

Before the graded tails, segment 2 agreed with the pre-regression build to 0.1 %,
and with a 1000-sub-point integral to 0.2 %, at every point. Segment 1 was not
identical to it: b273abc (the resonance-anchored lattice) and a4095a6 (physical
widths) both changed that grid on purpose, and at multiplier 5 the first point,
6 keV below the resonance, was 15 % above the 1000-sub-point value (the
pre-regression grid happened to be 3 % above it).

2026-10-27 (audit tier 3 entry 2, `tests/target_thickness_exact`): the thin-target energy
loss is now the exact range integral, not `epsilon(E0) N`. This target is 3.2e19 at/cm2 under
a stopping formula whose value falls towards lower energy in this range (3.35e-21 MeV cm2 at
272 keV), so the first-order loss of 107 keV (39 % of the beam energy) becomes 102 keV, and the
yields move by a few 0.1 %: expected chi2 15.0305 / 15.5912 -> 15.0799 / 15.6401 (total 30.6218
-> 30.7201). Re-baselined for that reason only.
