#ifndef AZURE2_PHYSICS_TRANSFER_THMVERTEX_H
#define AZURE2_PHYSICS_TRANSFER_THMVERTEX_H

#include <complex>
#include <vector>

namespace azure2 {
namespace physics {
namespace transfer {

/*!
 * The entrance vertex of a Trojan-horse (THM) level, its first part (docs/upgrade/transfer.md section 8;
 * docs/upgrade/THM.md section 1).  For a level lambda of the x + A system with level function u(r) inside the
 * channel radius a (u in fm^-1/2, (T + V) u = E_lambda u for r < a, B = a u'(a)/u(a) its boundary value) and the
 * wave w(r) of the transferred x relative to A (w = r j_l(pr) in the plane-wave limit, p in fm^-1, the half-off-shell
 * energy eps = hbar^2 p^2 / 2 mu), the vertex is the volume integral
 *
 *   W = int_0^a w(r) V(r) u(r) dr                                           (MeV fm^3/2 for w in fm),      (8.1)
 *
 * and, when w solves the free radial equation at eps, Green's theorem splits it exactly into a surface and an
 * interior term,
 *
 *   W = (hbar^2/2 mu) u(a) [B w(a)/a - w'(a)] + (E_lambda - eps) int_0^a w u dr,                            (8.2)
 *
 * whose surface bracket is, for w = r j_l(pr), the plane-wave vertex of the engine (THM.md section 1)
 *
 *   M_l(p) = (B - 1) j_l(pa) - pa j_l'(pa).                                                                  (8.3)
 *
 * With distorted waves (the next part) w is the projection of the DWBA source on the channel and Eq. (8.2) gains
 * the term of the distorting potential; this part implements Eqs. (8.1)-(8.3) for any w the caller gives.
 */

/// M_l(p) of Eq. (8.3) (dimensionless) for p > 0 (fm^-1), a > 0 (fm) and the boundary value B; j_l from the
/// Coulomb functions at eta = 0 (F_l(0, x) = x j_l(x); RadialNumerics.h's CoulombForMatching).
double PlaneWaveSurfaceVertex(int l, double p_fm, double a_fm, double B);

/// The level function inside the channel radius: u (fm^-1/2) and V (MeV) on r_i = i h, i = 0 ... n - 1, a = r_{n-1}.
struct LevelFunction {
  double step_fm = 0.0;
  std::vector<double> u;
  std::vector<double> potential_MeV;
  double energy_MeV = 0.0;  ///< E_lambda (c.m. of x + A)
  double mu_MeV = 0.0;      ///< reduced mass of x + A (MeV)
};

/// B = a u'(a)/u(a) of a level function (u'(a) by the five-point one-sided difference, O(h^4)).
double BoundaryValue(const LevelFunction &level);

/// The two forms of the vertex: Eq. (8.1), and the surface and interior terms of Eq. (8.2) (all MeV fm^3/2).
struct EntranceVertexValue {
  std::complex<double> volume;
  std::complex<double> surface;
  std::complex<double> interior;
};

/*!
 * Eqs. (8.1)-(8.2) for a wave w on the level function's grid (w in fm), its derivative w'(a) (dimensionless), the
 * half-off-shell energy eps (MeV) and the boundary value B (the level's own: BoundaryValue).  Simpson's rule.
 */
EntranceVertexValue EntranceVertex(const LevelFunction &level, const std::vector<std::complex<double>> &w,
                                   std::complex<double> dw_da, double epsilon_MeV, double B);

/// The plane-wave case: w = r j_l(pr) on the level's grid, eps = hbar^2 p^2 / 2 mu, B = BoundaryValue(level).
EntranceVertexValue PlaneWaveEntranceVertex(const LevelFunction &level, int l, double p_fm);

}  // namespace transfer
}  // namespace physics
}  // namespace azure2

#endif
