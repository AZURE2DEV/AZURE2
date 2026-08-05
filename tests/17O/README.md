# 17O — Trojan Horse (THM) regression

Guards the modified R-matrix / half-off-shell path, which nothing else in
`tests/` touches. Without it, a change to the ordinary on-shell machinery can
break THM silently, because the two share `ESegment`, `EPoint`, `AZURECalc` and
the gradient code but diverge at the observable.

## The reaction

<sup>17</sup>O(n,α)<sup>14</sup>C, extracted from the three-body
<sup>2</sup>H(<sup>17</sup>O,α<sup>14</sup>C)H quasi-free breakup.

M. L. Sergi, C. Spitaleri, M. La Cognata *et al.*,
*Improvement of the high-accuracy* <sup>17</sup>O(n,α)<sup>14</sup>C
*reaction-rate measurement via the Trojan Horse method for application to*
<sup>17</sup>O *nucleosynthesis*, Phys. Rev. C **91** (2015) 065803 —
data as distributed with the THM example, DOI 10.1103/PhysRevC.95.025807.

23 points, one segment.

## Why THM needs its own case

In a Trojan Horse measurement the reaction of interest is induced by a nucleon
bound inside a carrier nucleus, so it proceeds **below the Coulomb barrier
without penetrability suppression**, and the extracted cross section is
*half-off-shell*: the entrance channel is off the energy shell, the exit channel
on it. Two consequences make it structurally different from every other segment
type in the suite.

* **The observable is not the on-shell cross section.** It carries the THM form
  factor, so the R-matrix amplitude is folded differently. A segment declares
  this with an `isDiff` offset of **+10** on top of its ordinary observable
  code, which is why `isDiff = 10` here means "THM, angle-integrated".
* **The analytic gradient does not cover it.** The adjoint differentiates the
  T-matrix observable, not the half-off-shell one, so THM segments are excluded
  from the analytic accumulation and their contribution to the energy, width and
  normalization derivatives comes from central differences instead. That split
  lives in `AZURECalc::Gradient` and is exactly the kind of thing a regression
  test should pin.

## What is checked

`run_tests.sh` runs a plain calculation and compares the total and per-segment
chi-squared against `expected/chiSquared.out`:

```
1, 12.3193, N = 23,  norm = 2.63145e-05
```

The file is stored **at its fitted minimum**. The eight reduced width
amplitudes come from `output/param.sav` of a fit by JS; refitting from them
moves chi-squared by less than 0.02, so they sit in a genuine minimum.

Two things about this case are worth knowing before changing it.

**The widths are reduced width amplitudes, not partial widths.** Every channel
line carries the `gammaIsRWA` flag, so the `gamma` column holds an amplitude in
MeV^(1/2) rather than a width in eV or an ANC. AZURE2 keeps that convention on
output too — `CNuc::TransformOut` returns the amplitude unchanged for such a
channel — so a fit written back does not silently flip the convention. This
matters when editing the file by hand or through `pyazr`: converting rwa to
physical before writing would corrupt it.

**The normalization is small because the THM cross section is in arbitrary
units.** A Trojan Horse measurement determines the *shape*; its scale is fixed
by matching to direct data, so the fitted norm absorbs the unit conversion. It
is a fitted quantity, not a physical cross-section scale.

The case is deliberately left at these parameters rather than refitted on every
run, so the number is a fixed reference. A change to it means the THM path
moved.
