#ifndef AZURE2_PHYSICS_OBSERVABLES_ADJOINT_H
#define AZURE2_PHYSICS_OBSERVABLES_ADJOINT_H

// The reverse mode of Observables.h (docs/upgrade/physics_functions.md, stage 6): from the cotangent `bar` of an
// observable of a reaction at one energy (dL/d(observable), L a real scalar) to the cotangents of the T-matrix
// elements of its pathways.  Convention: Tbar = dL/dRe T + i dL/dIm T, so that dL = Re(conj(Tbar) dT); an
// observable Re(c conj(T)) therefore gives Tbar = bar c, and w |T|^2 gives Tbar = bar 2 w T.
//
// As in Observables.h a caller walks its own pathways and adds what these return into the cotangent of each;
// `geom` is pi/k^2 in fm^2 and `i1i2` the spin factor of the entrance pair.  Each expression is the one of the
// engine's AMatrixFunc::PointTBar it was moved from, with its order of operations: the gradients of every project
// depend on it to the last digit.

#include "azure2/physics/Observables.h"

#include <cmath>
#include <complex>
#include <cstddef>

namespace azure2 {
namespace physics {

/// Angle-integrated cross section (add_partial_wave_cross_section, cross_section_barns, times `symmetry`, the
/// (1 + delta_12) of an identical entrance pair where the value carries it): the cotangent of every pathway of
/// one (J, l, l'), whose coherent sum is T_sum.
inline std::complex<double> partial_wave_cotangent(double bar, double geom, double J, double i1i2, double symmetry,
                                                   std::complex<double> T_sum) {
  const double w = geom * (2.0 * J + 1.0) * i1i2 * symmetry;
  return bar * (2.0 / 100.0) * w * T_sum;
}

/// Total cross section by the optical theorem (add_optical_theorem_term, total_cross_section): the cotangent of a
/// diagonal pathway of total angular momentum J.  Real: the value depends on Re T only.
inline double optical_theorem_cotangent(double bar, double geom, double i1i2, double J) {
  return bar * geom * i1i2 / 100. * 2. * (2. * J + 1.);
}

/// One pathway of a beta-delayed entrance (add_beta_delayed_term, cross_section_barns).
inline std::complex<double> beta_delayed_cotangent(double bar, std::complex<double> T) {
  return bar * (2.0 * 25.0 / 100.0) * T;
}

/// Phase shift (phase_shift_degrees of U): Ubar = bar (90 / pi) (-Im U + i Re U) / |U|^2.  U must not vanish.
inline std::complex<double> phase_shift_cotangent(double bar, std::complex<double> U) {
  const double absU2 = std::norm(U);
  return bar * (90.0 / kObservablePi) * std::complex<double>(-std::imag(U), std::real(U)) / absU2;
}

/// U = (e^{2 i omega} - T) / e^{2 i omega} (collision_from_t_matrix): dU/dT = -1 / e^{2 i omega}, so
/// Tbar = conj(-1 / e^{2 i omega}) Ubar.  `exp_coulomb_phase` = e^{i omega}.
inline std::complex<double> collision_cotangent_to_t_matrix(std::complex<double> exp_coulomb_phase,
                                                            std::complex<double> Ubar) {
  const std::complex<double> expCP2 = exp_coulomb_phase * exp_coulomb_phase;
  const std::complex<double> invExpCP2 = 1.0 / expCP2;
  return std::conj(-invExpCP2) * Ubar;
}

/// The factor of the differential sum: dsigma/dOmega = c_R Re sum z T1 conj(T2) P_L in b/sr
/// (scale_legendre_coefficients), c_R = geom i1i2 symmetry / (100 pi).
inline double legendre_sum_factor(double geom, double i1i2, double symmetry) {
  return geom * i1i2 * symmetry / (kObservablePi * 100.0);
}

/// What stands for P_L(cos theta) when the observable is the Legendre coefficient of order L itself: the
/// coefficient is the differential one times 1 / symmetry (absolute, A_L) or 4 pi / (symmetry sigma_int)
/// (normalised, a_L = A_L 4 pi / sigma_int; sigma_int in b).
inline double legendre_coefficient_factor(bool absolute, double symmetry, double sigma_int) {
  return absolute ? 1.0 / symmetry : 4.0 * kObservablePi / (symmetry * sigma_int);
}

/// The weight W of one interference term z T1 conj(T2) P_L of the differential sum: the value gains
/// W Re(T1 conj(T2)).
inline double interference_weight(double sum_factor, double z, double legendre) { return sum_factor * z * legendre; }

/// The same for an unobserved primary with an observed secondary (upos_term), of order `l` with its R_l.
inline double upos_interference_weight(double sum_factor, double z_upos, int l, double R_l, double legendre) {
  return sum_factor * z_upos * std::pow(2.0 * l + 1.0, 0.5) / 4.0 * R_l * legendre;
}

/// value += W Re(T1 conj(T2)) gives T1bar += bar W T2 and T2bar += bar W T1: one of the two, for the partner
/// `T_other`.
inline std::complex<double> interference_pair_cotangent(double bar, double W, std::complex<double> T_other) {
  return bar * W * T_other;
}

/// The interference term z T1 conj(T2) of order L into the angle-independent coefficients of its two pathways:
/// q1 (pathway 1, order L) += z T2, q2 (pathway 2, order L) += z T1.  The cotangent of a pathway is then
/// bar c_R sum_L q_L P_L (legendre_contraction).  q1 and q2 may be the same cell.
inline void add_interference_cotangents(double z, std::complex<double> T1, std::complex<double> T2,
                                        std::complex<double> &q1, std::complex<double> &q2) {
  q1 += z * T2;
  q2 += z * T1;
}

/// sum_L q_L P_L over n orders.
inline std::complex<double> legendre_contraction(const std::complex<double> *q, const double *legendre, int n) {
  std::complex<double> s(0.0, 0.0);
  for (int l = 0; l < n; l++) s += q[l] * legendre[l];
  return s;
}

/// The factor of the Coulomb-nuclear interference of elastic scattering (elastic_interference): i geom f_it /
/// (100 sqrt(pi)), f_it = 2 for identical particles.
inline std::complex<double> coulomb_interference_factor(double geom, double it_factor) {
  return std::complex<double>(0.0, 1.0) * geom * it_factor / (std::sqrt(kObservablePi) * 100.0);
}

/// The Coulomb amplitude of channel spin s: f_d + (-1)^s f_x (f_x = 0 but for identical particles with spin).
inline std::complex<double> channel_spin_coulomb_amplitude(double s, std::complex<double> f_direct,
                                                           std::complex<double> f_exchange) {
  return (((std::lround(2.0 * s) / 2) % 2 == 0) ? 1.0 : -1.0) * f_exchange + f_direct;
}

/// value += Re(B conj(T)), B = factor stat_spin f_C(s) P_l, for a diagonal pathway of orbital momentum l:
/// Tbar += bar B.
inline std::complex<double> coulomb_interference_cotangent(double bar, std::complex<double> factor, double stat_spin,
                                                           std::complex<double> f_coulomb, double legendre) {
  const std::complex<double> B = factor * stat_spin * f_coulomb * legendre;
  return bar * B;
}

/// One (t, t') pair of a capture analysing power into the sums at one angle: den += a_k Re(T1 conj T2) P_k,
/// num += b_k Re(i T1 conj T2) P_k^1 (add_capture_ay_term with the angular functions); A_y = num / den.
inline void add_capture_ay_sums(double ak, double bk, double legendre, double assoc_legendre, std::complex<double> T1,
                                std::complex<double> T2, double &den, double &num) {
  const std::complex<double> prod = T1 * std::conj(T2);
  if (ak != 0.0) den += ak * std::real(prod) * legendre;
  if (bk != 0.0) num += bk * (-std::imag(prod)) * assoc_legendre;
}

/// The cotangents of the two elements of one pair for A_y = num / den (bar = 1):
///   d1 = (-i b_k P_k^1 T2 - A_y a_k P_k T2) / den,   d2 = (+i b_k P_k^1 T1 - A_y a_k P_k T1) / den.
inline void capture_ay_cotangents(double ak, double bk, double legendre, double assoc_legendre, double ay, double den,
                                  std::complex<double> T1, std::complex<double> T2, std::complex<double> &d1,
                                  std::complex<double> &d2) {
  const std::complex<double> I2(0.0, 1.0);
  d1 = (-I2 * bk * assoc_legendre * T2 - ay * ak * legendre * T2) / den;
  d2 = (I2 * bk * assoc_legendre * T1 - ay * ak * legendre * T1) / den;
}

}  // namespace physics
}  // namespace azure2

#endif
