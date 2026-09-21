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

Yields (total capture, c.m. keV):

| E_cm     | seg 1 (x5) | seg 2 (x20) | 1000 sub-points | pre-regression 7e533b7, seg 1 / seg 2 | broken dev, seg 1 |
|----------|------------|-------------|-----------------|---------------------------------------|-------------------|
| 253.7846 | 4.405e-08  | 3.851e-08   | 3.843e-08       | 3.973e-08 / 3.848e-08                 | 4.347e-07         |
| 263.0573 | 2.0473e-06 | 2.0394e-06  | 2.0393e-06      | 2.0444e-06 / 2.0396e-06               | 2.588e-06         |
| 264.8950 | 2.0853e-06 | 2.0774e-06  | 2.0773e-06      | 2.0835e-06 / 2.0775e-06               | 7.346e-05         |
| 298.3287 | 2.1645e-06 | 2.1505e-06  | 2.1504e-06      | 2.1603e-06 / 2.1507e-06               | 7.391e-05         |

Segment 2 agrees with the pre-regression build to 0.1%, and with a
1000-sub-point integral to 0.2%, at every point. Segment 1 is not identical to it: b273abc (the resonance-anchored lattice) and
a4095a6 (physical widths) both changed that grid on purpose. At multiplier 5
the first point, 6 keV below the resonance, is still under-resolved -- 15%
above the 1000-sub-point value, where the pre-regression grid happened to be
3% above it. That is the multiplier, not the regression; use 20.
