#ifndef AZURE2_PHYSICS_POTENTIALTUNING_H
#define AZURE2_PHYSICS_POTENTIALTUNING_H

#include <string>

#include "azure2/physics/HybridPotential.h"

namespace azure2 {
namespace physics {

/*!
 * What a pair's potential gives by itself at low energy, for tuning its depth (docs/upgrade/potential_background.md
 * section 13): its own phase shift, its scattering length and effective range from the effective-range expansion of
 * those phase shifts, the number of its bound states, and the depth that gives a chosen scattering length.
 *
 * Every function takes the potential of one partial wave as ForPartialWave makes it (the wave's depth and Gaussian
 * width in the central term), the charge numbers z1, z2, the reduced mass mu (MeV), l and <2 l.i>, as IntegrateOutward.
 * The model is the engine's: the potential acts from 0 to the matching radius r_m (p.matching_radius_fm), the point
 * Coulomb (or free) waves beyond; the radial equation by IntegrateOutward's Numerov grid (O(h^4)).  Energies E are
 * c.m. (MeV), k = sqrt(2 mu E)/(hbar c) (fm^-1), eta = Z1 Z2 alpha mu c^2/(hbar c k), non-relativistic.
 */

/// h(eta) = Re psi(i eta) - ln eta (dimensionless) for eta > 0, psi the digamma function: by psi(z + N) = psi(z) +
/// sum_{n<N} 1/(z + n) with N = 12 and the asymptotic series of psi at z + N to the term in z^-14 (relative error below
/// 1e-15 for eta > 0; checked against mpmath's digamma).
double CoulombH(double eta);

/// C_0^2(eta) = 2 pi eta / (exp(2 pi eta) - 1) (dimensionless; 1 at eta = 0).
double CoulombC0Squared(double eta);

/*!
 * The nuclear phase shift delta_l (rad, in (-pi/2, pi/2]) of the potential alone at E > 0, relative to the point-Coulomb
 * waves (to the free waves for z1 z2 = 0): the regular solution u integrated out to r_m (IntegrateOutward with a = r_m)
 * and matched there to u ~ cos(delta) F_l + sin(delta) G_l, tan delta = (u' F - k u F')/(k u G' - u' G) at rho = k r_m
 * (F, G of coulomb::coulomb).  NaN where the Coulomb functions fail (eta > ~100).
 */
double PotentialPhaseShift(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                           double E_MeV);

/*!
 * The number of bound states of the potential in the partial wave (E < 0, the potential from 0 to infinity with its own
 * exterior): by Sturm's oscillation theorem, the nodes of the zero-energy regular solution u_0 on (0, r_m) plus one when
 * its log-derivative r_m u_0'/u_0 lies below that of the zero-energy decaying solution at r_m, S_l(0) (zero_energy_shift;
 * -l without Coulomb) -- the eigenvalues E_n < 0 are the zeros of r u'/u - S(E) at r_m, and r u'/u decreases and S(E)
 * increases with E between the poles of r u'/u.  Deep states (Pauli-forbidden ones of a cluster potential) are counted.
 */
struct BoundStateCount {
  int bound = 0;            ///< bound states with E < 0
  int nodes = 0;            ///< sign changes of u_0 on (0, r_m)
  double logDerivative = 0; ///< r_m u_0'(r_m)/u_0(r_m)
  double exterior = 0;      ///< S_l(0) at r_m
};
BoundStateCount PotentialBoundStates(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i);

/*!
 * The effective-range expansion of the potential's own phase shifts (potential_background.md section 13.1):
 *
 *   charged, l = 0:  K(E) = C_0^2(eta) k cot delta_0 + 2 k eta h(eta) = -1/a_0 + r_0 k^2/2 + P k^4 + ...
 *   neutral, any l:  K(E) = k^(2l+1) cot delta_l               = -1/a_l + r_l k^2/2 + ...
 *
 * (Jackson and Blatt's Coulomb-modified function; a_0, r_0 in fm, a_l in fm^(2l+1), r_l in fm^(1-2l); a_0 > 0 for
 * 3He + 4He).  K is evaluated at 24 energies equally spaced in k^2 over [k_lo^2, k_hi^2] with k_hi = 0.5/R_N (R_N the
 * nuclear range: the Gaussian width r0, or R + 2a of a Woods-Saxon term) and k_lo = k_hi/4 (raised to where eta = 30 for
 * a charged pair), and a cubic in k^2 is fitted by least squares; a_0 = -1/K(0), r_0 = 2 dK/dk^2(0).  `aSpread`,
 * `rSpread`: the largest difference to two other fits (a quadratic on the same points, a cubic on the lower two thirds
 * of the window), the stated accuracy of the extrapolation (about 1e-5 relative for a_0 of 3He + 4He, section 13.1).
 * Refused (ok = false, `why`) for a charged pair with l > 0 (the Coulomb-modified function of l > 0 is not built), a
 * window that eta = 30 leaves empty, a phase shift that fails.
 */
struct EffectiveRange {
  bool ok = false;
  std::string why;
  double a = 0.0, r = 0.0;              ///< scattering length, effective range
  double aSpread = 0.0, rSpread = 0.0;  ///< the differences to the two other fits
  double eLow = 0.0, eHigh = 0.0;       ///< the window (MeV, c.m.)
  double kLow = 0.0, kHigh = 0.0;       ///< the window (fm^-1)
  int points = 0;
  int boundStates = 0;                  ///< PotentialBoundStates' count
};
EffectiveRange PotentialEffectiveRange(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i);

/*!
 * The depth V0 (MeV, of the central term of p, i.e. of the partial wave) for which PotentialEffectiveRange gives the
 * scattering length `a_target` with the same number of bound states as p has (the branch of p): within a branch a_l
 * runs monotonically from +inf (just after a state is bound at threshold) through 0 to -inf (just before the next one
 * is), so the branch holds exactly one such depth.  Scanned from p's depth in steps of 0.5 % towards it (deeper when
 * a_l(p) > a_target) up to 4 V0 or down to 0, then bisected to 1e-10 relative.  found = false with `why` when the scan
 * leaves the branch or the range first.
 */
struct DepthForLength {
  bool found = false;
  std::string why;
  double depth_MeV = 0.0;
  double a = 0.0, r = 0.0;   ///< a_l and r_l at that depth
  int boundStates = 0;       ///< the branch: the number of bound states, as p's
  int nodes = 0;             ///< nodes of the zero-energy solution on (0, r_m) at that depth
};
DepthForLength DepthForScatteringLength(const HybridPotential &p, int z1, int z2, double mu_MeV, int l,
                                        double two_l_dot_i, double a_target);

}  // namespace physics
}  // namespace azure2

#endif
