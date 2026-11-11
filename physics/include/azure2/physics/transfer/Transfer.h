#ifndef AZURE2_PHYSICS_TRANSFER_TRANSFER_H
#define AZURE2_PHYSICS_TRANSFER_TRANSFER_H

#include <complex>
#include <vector>

#include "azure2/physics/transfer/BoundState.h"
#include "azure2/physics/transfer/DistortedWaves.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The DWBA amplitude and cross section of the transfer A(a, s)F of a cluster x, a = s + x, F = A + x, to a bound
 * state of F (docs/upgrade/transfer.md sections 5-6).  Post form, the interaction the binding potential V_sx(r_sx)
 * of x in a (Coulomb included, as FRESCO's KIND=7, IP1=0), no remnant term (FRESCO's IP2 = 0):
 *
 *   T = < chi^(-)_f(r_sF) I_F | V_sx(r_sx) | I_a chi^(+)_i(r_aA) >,                                         (6.1)
 *
 * I_a = A_a [phi_{l_a j_a}(r_sx) x J_s]_{J_a} and I_F = A_F [phi_{l_F j_F}(r_xA) x J_A]_{J_F} the overlap functions,
 * coupled as FRESCO couples them ((l, s_x) j, J_core; J_composite), A the spectroscopic amplitudes (signed; their
 * isospin factors and any sqrt(N) of identical nucleons are in them, as for FRESCO), phi = (u/r) Y_lm (u of
 * BoundState.h).  With r = r_xA and rho = r_sx, alpha = m_A/m_F, beta = m_s/m_a (masses in u):
 *
 *   r_aA = r + beta rho,   r_sF = alpha r + rho,   d^3r d^3rho = d^3r_aA d^3r_sF / (1 - alpha beta)^3.         (5.2)
 *
 * The distorted waves are those of DistortedWaves.h without spin-orbit term (the spins of a, s, A, F enter only
 * through the overlaps): entrance a + A at E_i, exit s + F at E_f = E_i + Q.  The cross section, summed over the
 * final and averaged over the initial projections,
 *
 *   dsigma/dOmega = (mu_i mu_f / (2 pi hbar^2)^2) (k_f/k_i) / ((2J_a + 1)(2J_A + 1)) sum_{M} |T_M|^2   (mb/sr),  (6.6)
 *
 * mu_i, mu_f the reduced masses of a + A and s + F; the spin sum is done explicitly over the Clebsch-Gordan
 * coefficients of the two overlaps (Eq. (6.5) of the note).  The orbital amplitude O_{m_F m_a}(theta) and its
 * partial-wave and kernel form are Eqs. (6.2)-(6.4) of the note.
 *
 * Zero range: V_sx(rho) phi_{00}(rho) = D0 delta^3(rho) (l_a = 0 only; D0 in MeV fm^3/2, BoundState.h's
 * ZeroRangeConstant), then r_aA = r, r_sF = alpha r (FRESCO's KIND=5 with P2 = 0: no finite-range correction).
 */

/// One overlap: the bound state of x relative to the core (its l, 2j), the spectroscopic amplitude and the spins.
struct Overlap {
  BoundState state;            ///< u(r) of x in the composite; its potential is the post-form interaction
  double amplitude = 1.0;      ///< A (signed)
  int two_J_core = 0;          ///< 2J of the core (s for the projectile, A for the residual)
  int two_J_composite = 0;     ///< 2J of the composite (a, F)
};

/// The post-form interaction: the whole binding potential V_sx of the projectile overlap (FRESCO's default,
/// KRPOT = KBPOT), or its nuclear part V_sx - V_C (FRESCO with KRPOT a potential without the Coulomb term).
enum class PostInteraction { BindingPotential, NuclearPart };

/// The reaction.
struct TransferReaction {
  double mass_a_u = 0.0, mass_A_u = 0.0, mass_s_u = 0.0, mass_F_u = 0.0;  ///< masses (u)
  int two_s_x = 0;     ///< 2 s_x, the intrinsic spin of x
  Overlap projectile;  ///< x in a (core s)
  Overlap residual;    ///< x in F (core A)
  PostInteraction interaction = PostInteraction::BindingPotential;
};

/// The numerical settings of the finite-range amplitude.
struct TransferGrid {
  double step_fm = 0.1;           ///< step of the r_aA, r_sF grid: a multiple of both waves' steps
  double r_max_fm = 20.0;         ///< largest r_aA, r_sF (<= the matching radius of both waves)
  double vertex_range_fm = 10.0;  ///< V_sx(rho) u_a(rho) taken as 0 beyond this r_sx
  int angle_nodes = 48;           ///< Gauss-Legendre nodes in cos(r_aA, r_sF)
};

/// The finite-range post-form dsigma/dOmega (mb/sr) at the c.m. angles theta (deg) of s relative to a.
std::vector<double> FiniteRangeCrossSection(const TransferReaction &reaction, const DistortedWaves &entrance,
                                            const DistortedWaves &exit, const TransferGrid &grid,
                                            const std::vector<double> &theta_deg);

/// The zero-range dsigma/dOmega (mb/sr) with the constant D0 (MeV fm^3/2); needs l_a = 0.  The radial integral runs
/// on the residual overlap's grid.
std::vector<double> ZeroRangeCrossSection(const TransferReaction &reaction, double d0_MeV_fm32,
                                          const DistortedWaves &entrance, const DistortedWaves &exit,
                                          const std::vector<double> &theta_deg);

/*!
 * The zero-range orbital integral of the note's Eq. (6.7),
 *
 *   O_m(theta) = int d^3r chi^(-)*_f(c_f r) (u(r)/r) Y*_lm(r^) chi^(+)_i(c_i r)   (fm^3/2 for u in fm^-1/2),
 *
 * for m = -l ... l (index m + l), k_i along z and k_f at theta (deg) in the xz plane; the post-form zero range is
 * c_i = 1, c_f = alpha; the engine's prior-form THM weight M(E) (src/ThmDistortion.h) is c_i = beta, c_f = 1 with
 * l = 0 and M = sqrt(4 pi) O_0 (its phi has no Y_00).
 */
std::vector<std::complex<double>> ZeroRangeOrbitalIntegral(const BoundState &state, const DistortedWaves &entrance,
                                                           const DistortedWaves &exit, double c_i, double c_f,
                                                           double theta_deg);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
