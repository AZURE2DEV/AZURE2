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

The chi2 (3861.51 at the recorded parameters) sits between the bare
(3955.3) and fully convolved (3800.8) values, as it must.

Re-pinned 2026-09-26: the shift-function fix b6cc41b removed ~1e-9 noise from S(E) that ShftFunc::EnergyDerivative amplified into dS/dE, which enters the Brune transformation of sub-threshold levels. 3861.43 -> 3861.51 (2.1e-5, segments up to 3.4e-5). About 1.3e-5 of that predates b6cc41b.

Re-pinned 2026-09-28: the dS/dE fix (ChannelFunc::DerivativeStep) starts the central difference of the shift function from a quarter of the distance to threshold, at most 1 keV, instead of 1 eV, where it was round-off dominated (~1e-7..1e-5 relative, and different with and without FMA -- Linux and Windows CI differed from an x86-64-v3 build by up to 1.6e-3 in tests/thm_options). 3861.51 -> 3861.54 (7.8e-6).
