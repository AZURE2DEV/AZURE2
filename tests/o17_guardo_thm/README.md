# o17_guardo — Trojan Horse (THM) regression

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
1, 469.236, N = 23,  norm = 6.3529e-05
```

The normalization is small because the THM cross section is in arbitrary units —
a Trojan Horse measurement gives the *shape* and its scale is fixed by matching
to direct data, so the fitted norm absorbs the unit conversion. It is a fitted
quantity here, not a physical cross-section scale.

The case is deliberately kept at the parameters it ships with rather than
refitted, so the number is a fixed reference. A change to it means the THM path
moved.
