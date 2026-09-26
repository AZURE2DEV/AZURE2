# thm_narrow_fold — Gaussian fold across a resonance much narrower than the Gaussian

THM (HOES) spectra are folded with the experimental resolution, sigma of order
20 keV, and can contain resonances tens of eV wide. `check.sh` takes the
`tests/17O` model (21 keV lab Gaussian; the 5- level at E_cm = 74 keV is 46 eV
wide) and checks three things.

* **A, default grid.** The `<targetInt>` fold at 38 energies (E_cm 20–131 keV)
  against a fold of AZURE2's own *unfolded* curve. That curve is sampled at
  23,664 energies, 0.5 eV apart within 2 keV of the level and 20 eV apart
  elsewhere, and integrated with the trapezoid rule over the same ±5 sigma
  window. The check requires every point within 1e-3 and the area under the
  folded curve within 5e-4.
* **A, the "5 50" grid of `tests/17O`.** The same comparison with
  resonanceWidthMultiplier 5 and pointsPerWidth 50.
* **B, a moved level.** An external parameter file puts the level 1.5 keV
  (33 widths) above the `.azr` value. The sub-point grid is built from the
  `.azr` when the data are filled, as it is when a fit later moves a level.
  The result must equal a run whose `.azr` carries the moved energy.

| check | before this test | now |
|---|---|---|
| A, default (20 widths for a file that omits it, before) | 4.9e-5 | 1.2e-4 (2 widths), area +8e-5 |
| A, 5 widths, 50/width | 5.2e-3, area +1.9e-3 | 1.1e-4, area +7e-5 |
| B, level moved 1.5 keV | 31 % | 0 (identical grids) |

On the 33 eV level of the La Cognata/Guardo 17O(n,α) analysis, the "5 50" grid
was 16 % off at 150 sub-points. A level moved by 1–2 keV left the grid 13–70 %
off. A fit that freed that level's energy stopped 1 keV short of the true
value, at chi2 25 where the true minimum is 0.

What changed (`AdaptiveIntegrationGrid`, `EPoint::RefreshSubPointGrid`):

* The lattice of each narrow resonance now continues past its uniform core,
  ±resonanceWidthMultiplier·Γ, in geometric tails. The tail step is
  proportional to the distance from the resonance, and the tails run until
  that step reaches the smooth step. Before, the core ended abruptly and the
  smooth step bridged the 1/(E−E_R)² wings with chords that lie above them.
* Where lattices of several resonances overlap, only the finest is kept.
* The lattice is anchored at the width rounded to a power of 1.25 and the
  energy rounded to a quarter of that width. The grid is therefore a
  piecewise-constant function of the parameters. It is rebuilt from the
  current fit parameters whenever those anchors change (Brune formalism;
  not for components, beam profiles, energy-shifted segments or models with
  external-capture levels).
* A pure Gaussian convolution whose `<targetInt>` line leaves out
  resonanceWidthMultiplier now uses 2 widths, not 20. This needs about a fifth
  of the sub-points, and memory scales with the sub-points. A value on the line
  still wins, and target integration keeps 20.

Runtime is about 15 s.
