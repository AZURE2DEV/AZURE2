#ifndef AZURE2_PHYSICS_TRANSFER_RADIALNUMERICS_H
#define AZURE2_PHYSICS_TRANSFER_RADIALNUMERICS_H

#include <complex>
#include <vector>

#include "azure2/physics/coulomb/Coulomb.h"

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The numerical tools shared by the bound states, the distorted waves and the transfer amplitude
 * (docs/upgrade/transfer.md section 7): a radial grid r_i = i h, Numerov's method outward from the origin, Lagrange
 * interpolation on the grid, Simpson's rule, Gauss-Legendre nodes and the spherical harmonics Y_lm(theta, 0).
 * Every function is pure; the caller owns every vector.
 */

/*!
 * Numerov's method for u'' + q(r) u = 0 on r_i = i h (fm), i = 0 ... q.size() - 1, outward from the regular
 * solution u = r^{l+1} (1 + a_1 r + a_2 r^2 + a_3 r^3) at r = h, 2h (divided by h^{l+1}, so that it does not
 * underflow at large l), u_0 = 0, where near the origin q(r) = -l(l+1)/r^2 + q_m1/r + q_0 + O(r):
 * a_1 = -q_m1/(2l + 2), a_2 = -(q_m1 a_1 + q_0)/(2 (2l + 3)), a_3 = -(q_m1 a_2 + q_0 a_1)/(3 (2l + 4)) (the series
 * of the regular solution, n (n + 2l + 1) a_n = -q_m1 a_{n-1} - q_0 a_{n-2}; q_m1 (fm^-1) is -2 eta k for a point
 * charge, 0 otherwise; q_0 in fm^-2).  q (fm^-2) is k^2 - (2 mu/hbar^2) U(r) - l(l+1)/r^2 on the
 * grid (q[0] is not read).  Local error O(h^6), global O(h^4).  A solution that grows past 1e200 is rescaled (every
 * value so far divided by 1e200), so only its shape is meaningful; the caller normalises.  Returns u (q.size()
 * values).
 */
std::vector<std::complex<double>> NumerovOutward(const std::vector<std::complex<double>> &q, double h_fm, int l,
                                                 std::complex<double> q0_fm2, double q_m1_fm = 0.0);

/*!
 * Numerov's method inward: u given at the last two nodes (u_last at r_{n-1}, u_before at r_{n-2}), continued down to
 * node `lowest`.  Same conventions; a solution growing past 1e200 inward is rescaled.  Returns all n values (those
 * below `lowest` are 0).
 */
std::vector<std::complex<double>> NumerovInward(const std::vector<std::complex<double>> &q, double h_fm,
                                                std::complex<double> u_before, std::complex<double> u_last,
                                                int lowest);

/// The six-point Lagrange stencil at r on r_i = i h, i = 0 ... n - 1 (n >= 6): the first node (floor(r/h) - 2, kept
/// inside the grid) and the weights of the value and of d/dr (fm^-1), f(r) = sum_i value[i] f_{first + i}.
struct SixPointStencil {
  int first = 0;
  double value[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double derivative[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
};
SixPointStencil SixPointLagrange(double r_fm, double h_fm, int n);

/// f(r) from its values f_i on r_i = i h by six-point Lagrange interpolation (SixPointLagrange); 0 beyond the last
/// node.  Error O(h^6 f^(6)).
double InterpolateOnGrid(const std::vector<double> &f, double h_fm, double r_fm);
std::complex<double> InterpolateOnGrid(const std::vector<std::complex<double>> &f, double h_fm, double r_fm);

/// int_0^{r_{n-1}} f dr on r_i = i h by Simpson's rule (and the 3/8 rule on the last three intervals when their
/// number is odd).  Needs at least 4 nodes.
double Simpson(const std::vector<double> &f, double h_fm);
std::complex<double> Simpson(const std::vector<std::complex<double>> &f, double h_fm);

/*!
 * F_l, F_l', G_l, G_l' and sigma_l at (l, eta, rho) for the matching of waves: COULFG at and beyond the turning point
 * rho_tp = eta + sqrt(eta^2 + l(l+1)) (6e-14 against mpmath there, physics/coulomb/README.md), COULCC inside it
 * (1.5e-12).  COULCC is not used beyond the turning point: there it returns F and G with the wrong sign at isolated
 * points (l = 1, eta = 0.493007384468157, rho = 8.548514656481954: one point in 2.6e5 of a scan; reported in
 * docs/upgrade/transfer.md section 10).  Throws std::runtime_error when the functions fail.
 */
coulomb::CoulombWaves CoulombForMatching(int l, double eta, double rho);

/// The n-point Gauss-Legendre nodes x and weights w on [-1, 1] (Newton on P_n, to 1e-15).
void GaussLegendre(int n, std::vector<double> &x, std::vector<double> &w);

/*!
 * Y_lm(theta, phi = 0) for 0 <= l <= l_max, 0 <= m <= min(l, m_max) at x = cos(theta), Condon-Shortley phase
 * (Y_l,-m = (-1)^m Y_lm at phi = 0), by the stable recurrences of the normalised associated Legendre functions;
 * stored at [l (m_max + 1) + m] (0 where m > l).  Real at phi = 0.
 */
std::vector<double> SphericalHarmonicsAtPhiZero(int l_max, int m_max, double x);

/*!
 * The same at the polar angle theta (rad, any sign: a negative theta is the direction at phi = pi) with
 * dY_lm/dtheta = [sqrt((l-m)(l+m+1)) Y_l,m+1 - sqrt((l+m)(l-m+1)) Y_l,m-1]/2 (Y_l,-1 = -Y_l1), both stored at
 * [l (m_max + 1) + m] in the caller's vectors (resized; the same recurrence as SphericalHarmonicsAtPhiZero).
 */
void SphericalHarmonicsWithDerivative(int l_max, int m_max, double theta, std::vector<double> &y,
                                      std::vector<double> &dy);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
