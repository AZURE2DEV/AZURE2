#ifndef AZURE2_PHYSICS_CHANNELFUNCTIONS_H
#define AZURE2_PHYSICS_CHANNELFUNCTIONS_H

#include <complex>
#include <functional>

#include "azure2/physics/coulomb/Coulomb.h"

namespace azure2 {
namespace physics {

/*!
 * Channel functions of an open particle channel (docs/upgrade/coulomb.md section 2): from the Coulomb functions
 * F_l, G_l at rho = k a_c (Coulomb.h),
 *
 *   P_c = rho / (F^2 + G^2),   S_c = rho (F F' + G G') / (F^2 + G^2),   exp(-i phi_c) = (G - i F) / (F^2 + G^2)^1/2,
 *
 * all dimensionless, and S'_c = dS_c/dE_c (MeV^-1) by a Richardson-extrapolated central difference.  The order of
 * operations is the engine's (CoulFunc.cpp before this file existed), so the values are the same bits, including
 * the two guards just above a charged threshold (eta >~ 110), where F^2 + G^2 overflows and the formula is used with
 * F, G scaled by hypot(F, G) -- only then.  Kinematics: eta and rho at a channel energy E_c (MeV, c.m., above the
 * threshold of the pair in its states), non-relativistic or relativistic (Haeck, Hale, Paris and Thompson,
 * LA-UR-19-31084 (2019)), as the engine's Kinematics.h, which calls these.
 *
 * Valid for E_c > 0.  Below threshold (Whittaker functions), at threshold and for eta > 100 (P = 0 in double
 * precision) the engine's ChannelFunc applies its rules; those move here with the Whittaker functions (coulomb.md
 * section 1.4).
 */

/// The kinematic constants of a partition: charge numbers, reduced mass (u), masses (MeV/c^2), relativistic or not.
struct ChannelKinematics {
  int z1 = 0, z2 = 0;
  double mu_u = 0.0;               ///< reduced mass (u)
  double m1_MeV = 0.0, m2_MeV = 0.0;  ///< masses (MeV/c^2), used by the relativistic forms only
  bool relativistic = false;
};

/// Relativistic two-body kinematics at c.m. kinetic energy T (MeV; negative below threshold), masses in MeV/c^2:
/// k^2 (fm^-2), |k| (fm^-1), and eta (or eta_b = Z1 Z2 alpha mu_eff / kappa below threshold).
double relativistic_k2(double m1_MeV, double m2_MeV, double T_MeV);
double relativistic_k(double m1_MeV, double m2_MeV, double T_MeV);
double relativistic_eta(double z1z2, double m1_MeV, double m2_MeV, double T_MeV);

/// Sommerfeld parameter eta at E_c > 0 (MeV, c.m.).
double eta(const ChannelKinematics &kin, double E_MeV);
/// rho = k a at E_c > 0, radius a (fm).
double rho(const ChannelKinematics &kin, double a_fm, double E_MeV);
/// eta_b at binding energy B > 0 (MeV; E_c = -B).
double eta_bound(const ChannelKinematics &kin, double B_MeV);
/// kappa a at binding energy B > 0.
double rho_bound(const ChannelKinematics &kin, double a_fm, double B_MeV);

/// P = rho / (F^2 + G^2) from the Coulomb functions at rho.
double penetrability(const coulomb::CoulombWaves &w, double rho);
/// S = rho (F F' + G G') / (F^2 + G^2), with the overflow guard (rho G'/G when F^2 and F F' are both 0).
double shift(const coulomb::CoulombWaves &w, double rho);
/// exp(-i phi) = (G - i F) / (F^2 + G^2)^1/2, with the overflow guard (hypot).
std::complex<double> exp_minus_i_phi(const coulomb::CoulombWaves &w);
/// dS/dE at E > 0 (MeV^-1) from S(E) (any route: the Coulomb functions, or a table): Richardson's extrapolation
/// (4 D(h/2) - D(h)) / 3 of the central difference D(h) = [S(E + h) - S(E - h)] / 2h, h = 1e-3 E.  Error about 1e-8
/// relative with S from coulsm (S is good to 1e-10).
double shift_derivative(const std::function<double(double)> &S_of_E, double E_MeV);

/// Everything at one channel energy.
struct OpenChannelValues {
  double P = 0.0, S = 0.0, dS_dE = 0.0;
  std::complex<double> exp_minus_i_phi{1.0, 0.0};
  double sigma = 0.0;  ///< Coulomb phase sigma_l (rad)
  coulomb::Status status = coulomb::Status::failed;
};
/// P, S, S', exp(-i phi), sigma_l of channel (kin, l, a) at E_c > 0 (MeV, c.m.).  On a failure of the Coulomb
/// functions (at E_c or a point of the S' stencil), status says so and the numbers are NaN.
OpenChannelValues channel_functions(const ChannelKinematics &kin, int l, double a_fm, double E_MeV,
                                    coulomb::Method method = coulomb::Method::coulsm);

/// e^x K_n(x), integer n >= 0, x > 0: the trapezoid rule on int_0^inf e^(-x (cosh t - 1)) cosh(n t) dt in long
/// double (docs/upgrade/coulomb.md section 1.4, "Bound channel functions"); NaN for invalid arguments.
double bessel_k_scaled(int n, double x);

/*!
 * S_l(0), the limit of the shift function at threshold from either side: with c = eta k = Z1 Z2 alpha mu c^2 / hbar c
 * (fm^-1) and x = (8 c a)^(1/2), S_l(0) = -l - (x/2) K_2l(x) / K_2l+1(x) (the decaying solution sqrt(r) K_2l+1((8 c
 * r)^(1/2)) of the Coulomb equation at E = 0); -l without Coulomb (c = 0) or if K_2l+1 underflows.
 */
double zero_energy_shift(double c_fm, double a_fm, int l);

/// The bound-channel shift function of a neutral channel, S = rho k_l'(rho)/k_l(rho) + 1 by its recurrence
/// (S_0 = -rho, S_l = -rho^2 / (l - S_l-1) - l).
double neutral_bound_shift(int l, double rho);

/// What coulomb_bound_shift needs of a charged bound channel: its eta_b and rho = kappa a at a binding B (MeV), as
/// the caller's kinematics computes them, and c = eta k (fm^-1) and mu c^2 (MeV) for the threshold band.
struct BoundShiftArgs {
  std::function<void(double B_MeV, double &eta_b, double &rho)> eta_rho;
  double c_fm = 0.0;
  double mu_MeV = 0.0;
  double a_fm = 0.0;
};

/*!
 * S of a charged channel at binding B > 0: rho W'/W from whittaker_scaled where eta_b <= 1000 (B >= B_min, the binding
 * at eta_b = 1000), and closer to threshold the linear interpolation S(0) + (S(B_min) - S(0)) B / B_min, so S is
 * smooth through threshold.  S(0) where whittaker_scaled fails.
 */
double coulomb_bound_shift(int l, double B_MeV, const BoundShiftArgs &args);

}  // namespace physics
}  // namespace azure2

#endif
