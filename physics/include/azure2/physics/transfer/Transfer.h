#ifndef AZURE2_PHYSICS_TRANSFER_TRANSFER_H
#define AZURE2_PHYSICS_TRANSFER_TRANSFER_H

#include <complex>
#include <vector>

#include "azure2/physics/transfer/DistortedWaves.h"
#include "azure2/physics/transfer/Interaction.h"
#include "azure2/physics/transfer/Overlap.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The DWBA amplitude and cross section of the transfer A(a, s)F of a cluster x, a = s + x, F = A + x, to a bound
 * state of F, a bin or a scattering state of x + A (docs/upgrade/transfer.md sections 5-6, 15-17).  Written here for
 * the post form with the binding potential V_sx(r_sx) of x in a as the interaction (Coulomb included, as FRESCO's
 * KIND=7, IP1=0) and no remnant term (FRESCO's IP2 = 0); the prior form and the remnants are other sums of the
 * terms of Interaction.h (Eq. (17.1)) in the same kernel:
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
 * The distorted waves are those of DistortedWaves.h: entrance a + A at E_i, exit s + F at E_f = E_i + Q.  Without a
 * spin-orbit term the spins of a, s, A, F enter only through the overlaps; with one, the waves of a + A depend on
 * (l_i, j_i), j_i = l_i + J_a, and those of s + F on (l_f, j_f), j_f = l_f + J_s (the waves' spin must then be J_a,
 * J_s of the projectile overlap; the target's spin stays a spectator), and the kernel integrals are computed once per
 * pair of waves (note section 19).  The cross section, summed over the final and averaged over the initial
 * projections,
 *
 *   dsigma/dOmega = (mu_i mu_f / (2 pi hbar^2)^2) (k_f/k_i) / ((2J_a + 1)(2J_A + 1)) sum_{M} |T_M|^2   (mb/sr),  (6.6)
 *
 * mu_i, mu_f the reduced masses of a + A and s + F; the spin sum is done explicitly over the Clebsch-Gordan
 * coefficients of the two overlaps and of the distorted waves (Eqs. (19.5)-(19.6) of the note; without spin-orbit they
 * reduce to Eq. (6.5)).  The orbital amplitude and its partial-wave and kernel form are Eqs. (6.2)-(6.4) of the note.
 *
 * Zero range: V_sx(rho) phi_{00}(rho) = D0 delta^3(rho) (l_a = 0 only; D0 in MeV fm^3/2, BoundState.h's
 * ZeroRangeConstant), then r_aA = r, r_sF = alpha r (FRESCO's KIND=5 with P2 = 0: no finite-range correction).
 */

/// The reaction.  The overlaps are those of Overlap.h (a bound state, a bin or a scattering state of F); the
/// interaction is a sum of the terms of Interaction.h (default: the post form without remnant, Eq. (6.1)).
struct TransferReaction {
  double mass_a_u = 0.0, mass_A_u = 0.0, mass_s_u = 0.0, mass_F_u = 0.0;  ///< masses (u)
  int two_s_x = 0;     ///< 2 s_x, the intrinsic spin of x
  Overlap projectile;  ///< x in a (core s)
  Overlap residual;    ///< x in F (core A)
  std::vector<InteractionTerm> interaction = PostForm();
};

/// The numerical settings of the finite-range amplitude.
struct TransferGrid {
  double step_fm = 0.1;           ///< step of the r_aA, r_sF grid: a multiple of both waves' steps
  double r_max_fm = 20.0;         ///< largest r_aA, r_sF (<= the matching radius of both waves)
  double vertex_range_fm = 10.0;  ///< the projectile overlap (V_sx u_a, or u_a in a prior form or with a remnant) taken as 0 beyond this r_sx
  int angle_nodes = 48;           ///< Gauss-Legendre nodes in cos(r_aA, r_sF)
};

/// What a transfer calculation gives at each angle.
struct TransferObservables {
  std::vector<double> cross_section;         ///< dsigma/dOmega, in the unit of the residual overlap's kind
  std::vector<double> analysing_power_iT11;  ///< iT11 of the projectile a (Madison frame, Eq. (19.8)); 0 for J_a = 0
};

/// The finite-range dsigma/dOmega at the c.m. angles theta (deg) of s relative to a, in the unit of the residual
/// overlap's kind (mb/sr for a bound state or a bin, mb/(sr MeV) for a scattering state: OverlapKinds()).
std::vector<double> FiniteRangeCrossSection(const TransferReaction &reaction, const DistortedWaves &entrance,
                                            const DistortedWaves &exit, const TransferGrid &grid,
                                            const std::vector<double> &theta_deg);

/// The cross section and the projectile's vector analysing power, finite range.
TransferObservables FiniteRangeObservables(const TransferReaction &reaction, const DistortedWaves &entrance,
                                           const DistortedWaves &exit, const TransferGrid &grid,
                                           const std::vector<double> &theta_deg);

/// The zero-range dsigma/dOmega (unit as above) with the constant D0 (MeV fm^3/2); needs l_a = 0 and the post form
/// without remnant (the interaction is not read).  The radial integral runs on the residual overlap's grid.
std::vector<double> ZeroRangeCrossSection(const TransferReaction &reaction, double d0_MeV_fm32,
                                          const DistortedWaves &entrance, const DistortedWaves &exit,
                                          const std::vector<double> &theta_deg);
/// The same with the analysing power.
TransferObservables ZeroRangeObservables(const TransferReaction &reaction, double d0_MeV_fm32,
                                         const DistortedWaves &entrance, const DistortedWaves &exit,
                                         const std::vector<double> &theta_deg);

/*!
 * The zero-range orbital integral of the note's Eq. (6.7),
 *
 *   O_m(theta) = int d^3r chi^(-)*_f(c_f r) (u(r)/r) Y*_lm(r^) chi^(+)_i(c_i r)   (fm^3/2 for u in fm^-1/2),
 *
 * for m = -l ... l (index m + l), k_i along z and k_f at theta (deg) in the xz plane; the post-form zero range is
 * c_i = 1, c_f = alpha; the engine's prior-form THM weight M(E) (src/ThmDistortion.h) is c_i = beta, c_f = 1 with
 * l = 0 and M = sqrt(4 pi) O_0 (its phi has no Y_00).  Spinless: waves without a spin-orbit term.
 */
std::vector<std::complex<double>> ZeroRangeOrbitalIntegral(const RadialOverlap &state, const DistortedWaves &entrance,
                                                           const DistortedWaves &exit, double c_i, double c_f,
                                                           double theta_deg);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
