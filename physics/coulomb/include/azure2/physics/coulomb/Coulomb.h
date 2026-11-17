#ifndef AZURE2_PHYSICS_COULOMB_COULOMB_H
#define AZURE2_PHYSICS_COULOMB_COULOMB_H

#include <string>

namespace azure2 {
namespace physics {
namespace coulomb {

/*!
 * The Coulomb functions of the radial equation u'' + [1 - 2 eta/rho - l(l+1)/rho^2] u = 0 for real Sommerfeld
 * parameter eta (> 0 repulsive, 0 neutral), real rho = k r > 0 and integer l >= 0 (docs/upgrade/coulomb.md
 * section 1).  F_l is regular at the origin, G_l irregular; F ~ sin(theta), G ~ cos(theta) for large rho with
 * theta = rho - eta ln(2 rho) - l pi/2 + sigma_l (Abramowitz and Stegun ch. 14, DLMF 33.2); the Wronskian is
 * F' G - F G' = 1, primes d/drho.  sigma_l(eta) = arg Gamma(l + 1 + i eta), the continuous branch.
 *
 * Every function is pure: no object, no static or global mutable state, no cache, no output.  The same arguments
 * give the same bits in any thread and in any order of calls.  Arguments are dimensionless; eta and rho from the
 * channel energy, radius and masses are the caller's (ChannelFunctions.h).
 *
 * The sources are those of pycoulomb (github.com/skowrons94/pycoulomb; physics/coulomb/README.md records the
 * commits), one adaptor per algorithm behind one signature.
 */

/// The algorithm for F, F', G, G'.  They compute the same functions and differ in accuracy and speed (README.md).
enum class Method {
  coulsm,  ///< N. Michel's cwfcomp (CPC 176 (2007) 232) with the SAMMY-side patch: AZURE2's routine since 2025-09.
  coul,    ///< cwfcomp as in AZURE2 before 2025-09.
  coulfg,  ///< A. R. Barnett's COULFG (CPC 27 (1982) 147), C translation from LLNL.
  coulcc   ///< COULCC (Thompson and Barnett, CPC 36 (1985) 363), C++ translation bit-identical to the Fortran.
};

/// Outcome of one evaluation.  On anything but ok every number of the result is NaN.
enum class Status {
  ok,
  invalid_argument,  ///< l < 0 or l > 10000, eta or rho not finite, rho <= 0
  failed,            ///< the algorithm threw, reported failure or returned a non-finite value
  out_of_range       ///< the value exists but is outside the range of a double (Whittaker.h: whittaker)
};

/// Coulomb functions at one (l, eta, rho).
struct CoulombWaves {
  double F, Fp, G, Gp;  ///< F_l, dF_l/drho, G_l, dG_l/drho
  double sigma;         ///< Coulomb phase sigma_l(eta) (rad)
  Status status;
};

/// F_l, F'_l, G_l, G'_l and sigma_l at (l, eta, rho) by the given algorithm (coulsm by default).
CoulombWaves coulomb(int l, double eta, double rho, Method method = Method::coulsm) noexcept;

/// sigma_l(eta) = Im ln Gamma(l + 1 + i eta), not reduced modulo 2 pi.  NaN for an invalid l or eta.
double coulomb_phase(int l, double eta) noexcept;

/*!
 * The regular Riccati-Bessel functions u_l(x) = x j_l(x) = F_l(0, x) for l = 0 ... l_max at x into u[0 ... l_max]
 * (0 for x <= 0).  Each by a recurrence in its stable direction: upward from j_0, j_1 where l_max < x, else Miller's
 * downward recurrence from far above both, normalised with sum (2l+1) j_l^2 = 1 and signed by j_0 or j_1.  Within
 * 3e-15 of mpmath for l <= 8, x = 1e-6 ... 5e4 (tests/coulomb_neutral_test.cpp through coulsm).  The derivative is
 * u_l' = u_{l-1} - (l/x) u_l.  Pure, as the functions above.
 */
void riccati_bessel(int l_max, double x, double *u);

/// "coulsm", "coul", "coulfg", "coulcc".
const char *method_name(Method method) noexcept;
/// Inverse of method_name (case-sensitive); false, method unchanged, for an unknown name.
bool method_from_name(const std::string &name, Method &method) noexcept;

}  // namespace coulomb
}  // namespace physics
}  // namespace azure2

#endif
