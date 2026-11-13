#ifndef AZURE2_PHYSICS_WIDTHS_H
#define AZURE2_PHYSICS_WIDTHS_H

#include "azure2/physics/ChannelFunctions.h"

namespace azure2 {
namespace physics {

/*!
 * Width scales of a channel, as the engine uses them (docs/upgrade/pyazr3.md section 2.7, "Widths"):
 *
 *   gamma_W^2 = 3 hbar^2 / (2 mu a^2)          the Wigner (Teichmann-Wigner sum-rule) limit of a reduced width (MeV),
 *                                              mu the reduced mass, a the channel radius;
 *   Gamma_sp(E) = 2 P_l(E) gamma_W^2           the single-particle width at c.m. energy E > 0 (MeV), P_l the
 *                                              penetrability at rho = k a (ChannelFunctions.h);
 *   P_gamma = (|E_gamma| / hbar c)^(2L+1)      the photon "penetrability" of a multipole L (fm^-(2L+1)), the factor
 *                                              that turns a reduced photon amplitude into a width (Lane and Thomas,
 *                                              Rev. Mod. Phys. 30 (1958) 257, section XIII).
 *
 * The expressions are the engine's, operation for operation (AChannel::SetWignerLimit and the photon branch of the
 * level parameters call these), so the engine's numbers do not move.  Pure functions: no state.
 */

/// gamma_W^2 = 1.5 (hbar c)^2 / (mu c^2 a^2) (MeV) for a reduced mass mu_u (u) > 0 and a radius a_fm (fm) > 0; NaN
/// otherwise.  The convention of the engine's Parameter.wigner_limit (3 hbar^2 / (2 mu a^2), not hbar^2 / (mu a^2)).
double wigner_limit(double mu_u, double a_fm);

/// Gamma_sp = 2 P_l(E) gamma_W^2 (MeV) of channel (kin, l, a) at c.m. energy E > 0 (MeV), P_l from the Coulomb
/// functions of the given algorithm (coulsm, the engine's).  NaN for E <= 0, a <= 0, mu <= 0 or a failure of the
/// Coulomb functions.
double single_particle_width(const ChannelKinematics &kin, int l, double a_fm, double E_MeV,
                             coulomb::Method method = coulomb::Method::coulsm);

/// (|E_gamma| / hbar c)^(2L+1) (fm^-(2L+1)) for a photon of energy E_gamma (MeV) and multipolarity L >= 1.
double photon_penetrability(double E_gamma_MeV, int L);

}  // namespace physics
}  // namespace azure2

#endif
