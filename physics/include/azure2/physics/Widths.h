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
 *                                              Rev. Mod. Phys. 30 (1958) 257, section XIII);
 *   Gamma_W(XL)                                the Weisskopf single-particle estimate of a gamma width (eV), below.
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

/// The Weisskopf single-particle estimate of the radiative width of an electric (magnetic = false) or magnetic
/// 2^L-pole transition of energy E_gamma (MeV) in a nucleus of mass number A, in **eV** (the total width hbar * lambda
/// of the transition, not a reduced strength; no frame: E_gamma is the transition energy, recoil neglected), for the
/// nuclear radius R = r0 A^(1/3), r0 = 1.2 fm:
///
///   Gamma_W(EL) = 2 (L+1) / (L [(2L+1)!!]^2) (3/(L+3))^2 alpha hbar c (E_gamma/hbar c)^(2L+1) R^(2L),
///   Gamma_W(ML) = 80 (L+1) / (L [(2L+1)!!]^2) (3/(L+3))^2 alpha hbar c (hbar c / (2 m_p c^2))^2
///                 (E_gamma/hbar c)^(2L+1) R^(2L-2)                 (Bohr and Mottelson vol. 1, eqs. 3C-38, 3C-39),
///
/// i.e. hbar lambda(XL) = 8 pi (L+1) / (L [(2L+1)!!]^2) (E/hbar c)^(2L+1) B_W(XL) with the Weisskopf reduced strengths
/// B_W(EL) = (1/4 pi) (3/(L+3))^2 R^(2L) e^2 and B_W(ML) = (10/pi) (3/(L+3))^2 R^(2L-2) mu_N^2; alpha, hbar c and m_p
/// are Constants.h's.  Until 2026-11-18 (the owner's answer) the function used the two-digit coefficients of the
/// textbook tables (Blatt and Weisskopf 1952; Krane eqs. 10.13, 10.15): c_E1..E4 = 6.8e-2, 4.9e-8, 2.3e-14, 6.8e-21,
/// c_M1..M4 = 2.1e-2, 1.5e-8, 6.8e-15, 2.0e-21 eV MeV^-(2L+1), which are above the closed form by +0.79 % (E1),
/// +2.28 % (E2), +3.05 % (E3), -3.11 % (E4), +1.34 % (M1), +1.94 % (M2), -0.80 % (M3), -7.23 % (M4) (table / closed
/// form - 1, independent of E_gamma and A); a strength in W.u. moved by that factor (new / old = table / closed form).
/// NaN for L outside 1..4, E_gamma <= 0 or A <= 0.
double weisskopf_width_eV(bool magnetic, int L, double E_gamma_MeV, double mass_number);

}  // namespace physics
}  // namespace azure2

#endif
