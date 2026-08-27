# target_effect_ranges — a convolution restricted to part of a segment

The 7Li_p_ay analyzing-power model with a 30 keV Gaussian convolution
applied to both Deineko segments -- but only inside the lab-energy window
1.95-2.55 MeV, with a 120 keV smoothstep blend at each edge and an
automatic-application tolerance of 0.002.

This pins the optional trailing tokens of the targetInt line
("lo-hi,..." ranges, blend width, auto tolerance): points outside the
window are bit-identical to an unconvolved run, points well inside are
bit-identical to a fully convolved run, edge points blend between the
two, and the automatic decision skips the integration wherever it would
change the observable by less than the tolerance -- so any discontinuity
it introduces is bounded by that tolerance. Files that do not use the
tokens are read and written exactly as before, in both directions.

The chi2 (3861.43 at the recorded parameters) sits between the bare
(3955.3) and fully convolved (3800.8) values, as it must.
