#ifndef AZURE2_PHYSICS_COULOMBSCATTERING_H
#define AZURE2_PHYSICS_COULOMBSCATTERING_H

#include "azure2/physics/ChannelFunctions.h"

namespace azure2 {
namespace physics {

/*!
 * The Coulomb reference cross sections of elastic scattering (docs/upgrade/data_kinds.md section 4.3), in the c.m.
 * frame, b/sr, at a c.m. kinetic energy E (MeV) and a c.m. angle theta (deg, 0 < theta < 180):
 *
 *   Rutherford   dsigma_R/dOmega = eta^2 / (4 k^2 sin^4(theta/2))
 *   Mott         dsigma_M/dOmega = eta^2 / (4 k^2) [ 1/sin^4(theta/2) + 1/cos^4(theta/2)
 *                                   + 2 w cos(eta ln tan^2(theta/2)) / (sin^2(theta/2) cos^2(theta/2)) ],
 *                w = (-1)^(2s) / (2s + 1) for two identical particles of spin s
 *
 * (1 fm^2 = 0.01 b).  The Mott form is |f_C(theta) + f_C(pi - theta)|^2 for spin 0, and for spin s the average over
 * the channel spins S = 0 .. 2s of |f_C(theta) + (-1)^S f_C(pi - theta)|^2 with weights (2S + 1)/(2s + 1)^2, which
 * leaves the interference term with w.  eta and k are the engine's (eta and rho of ChannelFunctions.h, with the
 * kinematics of the pair: non-relativistic, or relativistic with kin.relativistic), so that the reference is the
 * Coulomb part of the engine's elastic cross section, |C_alpha|^2 pi/k^2 (GenMatrixFunc), with the same masses and
 * constants.  Non-relativistically the Rutherford form is (Z1 Z2 e^2 / 4E)^2 / sin^4(theta/2), e^2 = alpha hbar c.
 *
 * d_dE is the derivative with respect to E at fixed c.m. angle (b/sr/MeV), from the derivatives of eta and k: for a
 * fitted energy shift, and as the first-order size of a shift's effect on the reference.
 */
struct CoulombCrossSection {
  double value = 0.0;  ///< b/sr
  double d_dE = 0.0;   ///< b/sr/MeV
};

/// Rutherford cross section of the pair `kin` at c.m. energy E_MeV > 0 and c.m. angle theta_deg in (0, 180).
CoulombCrossSection rutherford_cross_section(const ChannelKinematics &kin, double E_MeV, double theta_deg);

/// Mott cross section of two identical particles of spin two_s / 2 (two_s >= 0) with the kinematics `kin`.
CoulombCrossSection mott_cross_section(const ChannelKinematics &kin, int two_s, double E_MeV, double theta_deg);

}  // namespace physics
}  // namespace azure2

#endif
