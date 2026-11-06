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

Segment 2 agrees with the pre-regression build to 0.1%, and with a
1000-sub-point integral to 0.2%, at every point. Segment 1 is not identical to it: b273abc (the resonance-anchored lattice) and
a4095a6 (physical widths) both changed that grid on purpose. At multiplier 5
the first point, 6 keV below the resonance, is still under-resolved -- 15%
above the 1000-sub-point value, where the pre-regression grid happened to be
3% above it. That is the multiplier, not the regression; use 20.

Re-pinned 2026-09-26: the shift-function fix b6cc41b removed ~1e-9 noise from S(E) that ShftFunc::EnergyDerivative amplified into dS/dE, which enters the Brune transformation of sub-threshold levels. 19.0295 -> 19.0249 (2.4e-4).

Re-pinned 2026-09-26: the adaptive grid now continues each resonance's lattice past its uniform core in geometric tails (tests/thm_narrow_fold). Segment 1 (multiplier 5) is therefore converged too. Its first point went from 4.405e-08, 15 % above the 1000-sub-point integral, to 3.856e-08 (+0.3 %); its other points agree to 1e-4. Segment 2 went to 3.847e-08 / 2.0394e-06 / 2.0773e-06 / 2.1506e-06, within 0.1 % of the 1000-sub-point column. The earlier segment-1 chi2 of 3.61 was the grid error happening to sit closer to the data; with both segments converged they carry the same chi2. 19.0249 -> 30.8444 (segment 1 3.607 -> 15.28, segment 2 15.418 -> 15.565).

2026-10-27 on dev (audit tier 3 entry 2, `tests/target_thickness_exact`, dev 3207ae1): the thin-target
energy loss is the exact range integral, not `epsilon(E0) N`. This target is 3.2e19 at/cm2 under a stopping
formula whose value falls towards lower energy in this range, so the first-order loss of 107 keV becomes
102 keV and the yields move by a few 0.1 %. Re-pinned at the merge of dev 2d97b1f (November 2026) for that
reason only: 30.8444 -> 30.9094 (segment 1 15.2792 -> 15.3329, segment 2 15.5652 -> 15.5765). The merged
build with the first-order loss restored gives 30.8444 and every output file of the thm build byte for byte,
so the whole move is that fix (dev's own pin, 30.7201, is on dev's grid; see MERGE_NOTES.md).
