# thm_energy_shift — an energy shift of a THM segment

`check.sh` runs `tests/17O` (17O(n,a)14C THM: 23 points from E_lab = -54 keV,
four of them below the n+17O threshold; the 33 eV 5- level at E_cm = 74 keV;
21 keV Gaussian fold) with a segment energy shift d, twice:

- **A**: `segment_1_energy_shift = d` in an external parameter file (the
  `param.par` of a plain run with that one line changed);
- **B**: the data file's lab energies moved by d, no shift.

The two are the same measurement, so the model at every point (column 4 of
`AZUREOut_aa=1_R=2.out`; column 1 is the shifted E_cm in both) must agree to
1e-4. It does to ~2e-6, for d = +6 keV and d = -4 keV.

What it guards (fixed 2026-09-28):

- `ESegment::UpdatePointEnergiesWithShift` moved only points and folding
  sub-points with original energy > 0. The guard (with the floor of
  `ShiftedEnergy`) keeps an ordinary segment's energies positive; a THM
  segment has points and sub-points below the entrance threshold, and those
  stayed put: the four points below threshold were 27 % off, and the fold of
  every point within reach of threshold mixed shifted and unshifted
  sub-points. For THM segments every point and sub-point now shifts, with no
  floor; ordinary segments are unchanged.
- `EPoint::RefreshSubPointGrid` did not rebuild the grid of a shifted segment,
  so the lattice built around the narrow level at the unshifted energies sat
  a whole shift away from it (30-40 % off near 74 keV). A THM segment's grid
  is now built in the frame of its unshifted data around the anchors moved
  back by the shift, quantized there, so it follows the level as the shift
  moves (angle-integrated THM segments; other shifted segments keep their
  grid).

About 15 s.
