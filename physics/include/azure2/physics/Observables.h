#ifndef AZURE2_PHYSICS_OBSERVABLES_H
#define AZURE2_PHYSICS_OBSERVABLES_H

// The observables of a reaction from the T-matrix elements of its pathways at one energy
// (docs/source/developer/level_matrix.rst, "Cross Section Calculation"; docs/source/theory/conventions.rst):
//
//   angle-integrated   sigma = (pi/k^2) g sum_{J,l,l'} (2J+1) |sum_pathways T|^2          (coherent by J, l, l')
//   total              sigma = (pi/k^2) g sum_J (2J+1) sum_c 2 Re T_cc                    (optical theorem)
//   phase shift        delta = arg(U)/2,  U = (e^{2i omega} - T)/e^{2i omega}
//   Legendre           B_L = (pi/k^2) g f_sym / pi  Re sum Z1Z2 T1 conj(T2)               (per order L)
//   elastic            the Coulomb-nuclear interference per (channel spin, l), its direct and exchange coefficients
//   A_y of capture     sum_k b_k P_k^1 / sum_k a_k P_k, a_k with Re T T'*, b_k with Re (i T T'*)
//
// A caller walks its own list of pathways and hands each element here; `geom` is the caller's geometrical factor
// pi/k^2 in fm^2 (a cross section in barns is the sum / 100), `i1i2` the spin factor 1/((2 I_1 + 1)(2 I_2 + 1)) of
// the entrance pair.  Every function is a few operations per pathway, so all are inline; each adds into the
// caller's sum, which keeps the order of the additions the caller's.

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace azure2 {
namespace physics {

/// pi, the digits of the engine's constant.
constexpr double kObservablePi = 3.14159265358979323846;

/// The coherent sum of the T-matrix elements of the pathways of one (J, l, l').
struct PartialWaveSum {
  double J = 0.0;
  int l = 0;   ///< entrance orbital angular momentum
  int lp = 0;  ///< exit orbital angular momentum, or the multipolarity of an external-capture pathway
  std::complex<double> T;
};

/// The coherent sums of one channel-spin group, in order of first appearance.  Internal and external pathways of
/// one (J, l, l') add before the modulus is taken: that is their interference in the angle-integrated cross section.
struct PartialWaveSums {
  std::vector<PartialWaveSum> rows;
  void clear() { rows.clear(); }
  /// Adds T to the sum of (J, l, lp), starting it when it is the first.
  void add(double J, int l, int lp, std::complex<double> T) {
    for (PartialWaveSum &r : rows) {
      if (J == r.J && l == r.l && lp == r.lp) {
        r.T += T;
        return;
      }
    }
    rows.push_back(PartialWaveSum{J, l, lp, T});
  }
};

/// sum += geom (2J+1) i1i2 |T|^2 over the coherent sums: the angle-integrated cross section in fm^2.
inline void add_partial_wave_cross_section(const PartialWaveSums &sums, double geom, double i1i2,
                                           std::complex<double> &sum) {
  for (const PartialWaveSum &r : sums.rows) sum += geom * (2. * r.J + 1.) * i1i2 * (r.T) * std::conj(r.T);
}

/// sum -= geom (2J+1) i1i2 |T|^2: the particle exits taken off the R-matrix-capture cross section.
inline void subtract_partial_wave_cross_section(const PartialWaveSums &sums, double geom, double i1i2,
                                                std::complex<double> &sum) {
  for (const PartialWaveSum &r : sums.rows) sum -= geom * (2. * r.J + 1.) * i1i2 * (r.T) * std::conj(r.T);
}

/// The term of one elastic pathway (entrance = exit channel) of the R-matrix-capture cross section,
/// 2 geom (2J+1) i1i2 Im T, in fm^2.
inline double rmc_elastic_term(double geom, double J, double i1i2, std::complex<double> T) {
  return 2. * geom * (2. * J + 1.) * i1i2 * std::imag(T);
}

/// sum += 25 |T|^2: one pathway of a beta-delayed entrance (no geometrical factor).
inline void add_beta_delayed_term(std::complex<double> T, std::complex<double> &sum) {
  sum += 25. * T * std::conj(T);
}

/// A sum in fm^2 as a cross section in barns.
inline double cross_section_barns(std::complex<double> sum) { return std::real(sum) / 100.; }

/// total += (2J+1) 2 Re T_cc: one diagonal pathway of the optical theorem.
inline void add_optical_theorem_term(double J, std::complex<double> T, double &total) {
  total += (2. * J + 1.) * 2. * std::real(T);
}

/// The total cross section in barns from the sum of add_optical_theorem_term.
inline double total_cross_section(double geom, double i1i2, double total) { return geom * i1i2 * total / 100.; }

/// U of a diagonal pathway from its T: U = (e^{2i omega} - T) / e^{2i omega}, `exp_coulomb_phase` = e^{i omega}.
inline std::complex<double> collision_from_t_matrix(std::complex<double> exp_coulomb_phase, std::complex<double> T) {
  std::complex<double> expCoulPhaseSquared = exp_coulomb_phase * exp_coulomb_phase;
  return (expCoulPhaseSquared - T) / expCoulPhaseSquared;
}

/// The nuclear phase shift in degrees, half the argument of U = exp(2 i delta): in (-90, 90], and in [0, 180) for a
/// pair of identical particles (the convention of partial-wave analyses).
inline double phase_shift_degrees(std::complex<double> U, bool identical) {
  double phase = 180.0 / kObservablePi / 2.0 * std::atan2(std::imag(U), std::real(U));
  if (identical && phase < 0.0) phase += 180.0;
  return phase;
}

/// b += z T1 conj(T2): one interference term of the differential sum into its Legendre order.
inline void add_interference_product(double z, std::complex<double> T1, std::complex<double> T2,
                                     std::complex<double> &b) {
  b += z * T1 * std::conj(T2);
}

/// The sums of an order (the real parts of add_interference_product, or of add_upos_term) as the Legendre
/// coefficients B_L in b/sr: b / pi * geom * i1i2 * symmetry / 100.
inline void scale_legendre_coefficients(std::vector<double> &c, double geom, double i1i2, double symmetry) {
  for (double &b : c) b = b / kObservablePi * geom * i1i2 * symmetry / 100.;
}

/// R_k of a secondary decay of mixed multipolarity L, L + 1 with mixing ratio delta:
/// (R_LL + 2 delta R_LL' + delta^2 R_L'L') / (1 + delta^2).
inline double upos_mixed_rk(double R_L, double R_LLp, double R_LpLp, double delta) {
  return (R_L + 2. * delta * R_LLp + delta * delta * R_LpLp) / (1. + delta * delta);
}

/// One interference term of an unobserved primary with an observed secondary, of order l, with its R_l.
inline double upos_term(double z_upos, std::complex<double> T1, std::complex<double> T2, int l, double R_l) {
  return std::real(z_upos * T1 * std::conj(T2)) * std::pow(2. * l + 1., 0.5) / (4.) * R_l;
}

/// The Coulomb-nuclear interference of one (channel spin, l): c = sum statSpin conj(T) over the diagonal pathways.
struct CoulombInterferenceTerm {
  double s = 0.0;
  int l = 0;
  std::complex<double> c;
};

/// The terms of an elastic exit, in order of first appearance, without the angle-dependent f_C(s) P_l.
struct CoulombInterference {
  std::vector<CoulombInterferenceTerm> terms;
  void clear() { terms.clear(); }
  /// Adds the diagonal pathway of channel spin s and orbital momentum l.
  void add(double s, int l, double stat_spin, std::complex<double> T) {
    std::size_t i = 0;
    while (i < terms.size() && !(terms[i].s == s && terms[i].l == l)) i++;
    if (i == terms.size()) terms.push_back(CoulombInterferenceTerm{s, l, std::complex<double>(0.0, 0.0)});
    terms[i].c += stat_spin * std::conj(T);
  }
};

/*!
 * The coefficients of the Coulomb part of elastic scattering,
 *   dsigma/dOmega = |F_C|^2 + |F_N|^2 + 2 Re[F_C* F_N],
 *   2 Re[F_C* F_N] / 100 = Re sum_l (f_d D_l + f_x X_l) P_l,
 *   D_l = (i geom f_it / (100 sqrt(pi))) sum_s c_{s,l},  X_l = the same with (-1)^s,
 * f_it = 2 for identical particles.  `coulomb` multiplies |f_C|^2 and `cross` the exchange cross term of
 * identical particles with spin (of spin `spin`), set where `spin_dependent`.
 */
struct ElasticInterference {
  double coulomb = 0.0;
  double cross = 0.0;
  std::vector<std::complex<double>> direct;    ///< D_l
  std::vector<std::complex<double>> exchange;  ///< X_l (spin-dependent exchange only)
};

/// Fills `e` from the terms.  Returns -1, or the first l above `max_order` (the fill stops there).
inline int elastic_interference(const CoulombInterference &terms, double geom, bool identical, bool spin_dependent,
                                double spin, int max_order, ElasticInterference &e) {
  e.coulomb = geom / 100.;
  if (spin_dependent) {
    const double j = spin;
    e.cross = 2.0 * e.coulomb * ((std::lround(2.0 * j) % 2 == 0) ? 1.0 : -1.0) / (2.0 * j + 1.0);
  }
  const double itFactor = identical ? 2.0 : 1.0;
  const std::complex<double> cI = std::complex<double>(0., 1.) / std::sqrt(kObservablePi) * geom * itFactor / 100.;
  for (const CoulombInterferenceTerm &t : terms.terms) {
    if (t.l > max_order) return t.l;
    if ((int)e.direct.size() <= t.l) e.direct.resize(t.l + 1, std::complex<double>(0.0, 0.0));
    e.direct[t.l] += cI * t.c;
    if (!spin_dependent) continue;
    if ((int)e.exchange.size() <= t.l) e.exchange.resize(t.l + 1, std::complex<double>(0.0, 0.0));
    // channel spins of two identical particles are integers
    e.exchange[t.l] += (((std::lround(2.0 * t.s) / 2) % 2 == 0) ? 1.0 : -1.0) * cI * t.c;
  }
  return -1;
}

/// One (t, t') pair of a capture analysing power (Seyler and Weller 1979) into the sums of its order: the
/// denominator takes a_k Re T T'*, the numerator b_k Re (i T T'*).
inline void add_capture_ay_term(double ak, double bk, std::complex<double> T1, std::complex<double> T2,
                                double &denominator, double &numerator) {
  const std::complex<double> prod = T1 * std::conj(T2);
  if (ak != 0.0) denominator += ak * std::real(prod);
  if (bk != 0.0) numerator += bk * (-std::imag(prod));
}

}  // namespace physics
}  // namespace azure2

#endif
