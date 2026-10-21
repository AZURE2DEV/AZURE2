# channel_space — the channel-space level-matrix solution

A synthetic p + 12C model built to exercise the channel-space route of
`AMatrixFunc` (`include/ChannelSpace.h`): J = 1/2+ with eight levels and
3/2- with seven, each with four channels (p0, the two channel spins of p1 to a
synthetic 2+ state of the target at 0.8 MeV, and a photon channel), plus a
one-level 5/2+ group that keeps the level matrix.  The data are the model
times 1 + 0.05 sin(7 E + segment) with 4 % errors (`data/*.dat`, lab frame):
elastic scattering at 90 and 150 degrees, the inelastic and the capture cross
section.  The 90-degree set has one point at the energy of the 4.70 MeV 3/2-
level, a pole of the Lane-Thomas level matrix, where the channel route must
keep that level out of R and solve it with the channels (deflation).

`run_tests.sh` compares the chi-squared with `expected/chiSquared.out`;
`check.sh` compares the channel route with `AZURE_LEVEL_SPACE=1` (outputs,
bilinear forms, gradients) and the gradients with finite differences, for the
Brune, Lane-Thomas and Park level matrices.  Section (d) extends the model with a `<scatteringLengths>` block (the
1/2+ length is evaluated at the p + 12C threshold on the channel route), an
angular distribution of 29 angles at each of 10 energies, and electron
screening with U_e free on the elastic segments; in a Lane-Thomas variant a
1/2+ level lies 0.5 keV above the threshold.  It checks that
`scatteringLengths.out` is identical on both routes, that the gradient stays
analytic with the constraints, that it (and d chi^2/d U_e on the adjoint and
residual-Jacobian paths) agrees with finite differences, and that it agrees
with the level-space gradient.
