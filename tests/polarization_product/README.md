# polarization_product — P dsigma/dOmega as a fitted observable

Ten digitized points of the neutron polarization from
<sup>11</sup>B(alpha,n)<sup>14</sup>N at E_alpha = 2.049 MeV, Niecke *et al.*,
Nucl. Phys. **A289** (1977) 408, fig. 6(b), entered as `isDiff = 8` on the
<sup>11</sup>B+alpha entrance channel.

This pins the observable introduced alongside the analysing power: the
outgoing vector polarization multiplied by the differential cross section,
which is the form most polarization measurements are published in.

It is deliberately a **spin-0 projectile** case. An alpha has no vector
analysing power, so `isDiff = 7` is identically zero here and these data
cannot be fitted that way at all — the only alternative is the inverse
reaction <sup>14</sup>N(n,alpha)<sup>11</sup>B via the time-reversal
relation. `isDiff = 8` fits them on the channel where they were measured.

What the recorded chi-squared guards, beyond "the numbers did not move":

* the exit-index Pauli trace in `AmplitudeMatrix::OutgoingPolarizationPy`,
  including the Clebsch-Gordan decomposition of the exit channel spin into
  ejectile and residual, and its sign;
* the multiplication onto AZURE2's own differential cross section, so the
  result carries the cross section's units rather than the amplitude
  matrix's bare spin sum;
* `ESegment`/`EPoint` treating code 8 as a centre-of-mass differential
  observable and propagating it to target-effect sub-points.

Identities checked when the observable was written, which this case will
break if any of them stops holding:

* for elastic scattering P_y and A_y are equal — 178 of 178 points of
  `tests/13N` agree to every printed digit, zero mismatches;
* P_y vanishes at theta = 0 and 180 degrees, since the polarization lies
  along k_in x k_out;
* an `isDiff = 8` run equals P_y times an `isDiff = 5` run to 5e-8;
* `P_y * sigma / N` is exactly 1/nEntrance, N being the polarization
  numerator — so the published product is the *bare numerator*, bilinear
  in the amplitudes. That is why its analytic adjoint
  (`OutgoingPolarizationNumeratorBar`) carries no 1/D^2 and is better
  conditioned than the analysing-power adjoint it is built from.

**The chi-squared is large (228 per point) and that is not a failure.**
The <sup>11</sup>B+alpha model in this file is not fitted to these data —
it predicts a polarization of the opposite sign over the back hemisphere.
The reference is a fixed number for regression purposes, nothing more.

The suite runs mode 1 with a blank external parameter file, so everything
the number depends on is in the `<levels>` block of the `.azr`.
