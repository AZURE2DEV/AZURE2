# level_at_threshold -- a level energy at a particle threshold

Every place that evaluates the shift function, penetrability or dS/dE at a
level energy (CNuc::CalcBoundaryConditions, CalcShiftFunctions, TransformIn,
TransformOut, CheckRadiativeWidths; the THM vertex=constant shift
ShiftAtLevelEnergy; AZUREGrad's dS/dE table; the adaptive grid's width
estimate) and the external-capture entrance channel (EPoint, ECIntegral) goes
through `ChannelFunc` (include/ChannelFunc.h), which handles threshold once:
S continued through E = 0, P = 0 at and just above a Coulomb threshold
(eta > 100), dS/dE differenced on the continuous S where its stencil crosses
threshold.

`check.sh` puts the first level of three one-level projects (12C+12C ->
a+20Ne and elastic; p+7Li -> a+a; the THM project tests/18O_p_a_thm) at
E_res = -1e-9, 0 and +1e-9 MeV with the entrance reduced width given as an
amplitude, and requires every run to finish within 60 s with finite output and
the model of the three runs to agree to 1e-6 (measured 1e-8, 1.3e-7, 2e-9).
The binary before the change aborted (core dump) at E_res = 0 and +1e-9 for
p+7Li and the THM project.

The unit test tests/reference/channel_threshold_test (ctest
channel_threshold) checks the helper itself on p+7Li, 6Li+d, 12C+12C and
17O+n, l = 0-3.
