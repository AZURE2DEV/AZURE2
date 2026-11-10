#ifndef AZURE2_PHYSICS_HYBRIDPOTENTIAL_H
#define AZURE2_PHYSICS_HYBRIDPOTENTIAL_H

#include <complex>
#include <functional>
#include <optional>
#include <vector>

namespace azure2 {
namespace physics {

/*!
 * The nuclear potential of a hybrid particle pair and the radial equation it enters (docs/upgrade/
 * spin_orbit_potential.md section 1.1).  For a pair with reduced mass mu (MeV) and charge numbers Z1, Z2, at the
 * channel energy E (MeV, c.m.) and orbital angular momentum l,
 *
 *   u''(r) = (2 mu / hbar^2) [V(r) - E] u(r),   V(r) = V_C(r) + V_N(r) + hbar^2 l(l+1) / (2 mu r^2),
 *
 * with V_C the potential of a uniformly charged sphere of radius R_C inside R_C and Z1 Z2 e^2 / r outside, and V_N a
 * Woods-Saxon well -V0 / (1 + exp((r - R)/a)) or a Gaussian -V0 exp(-(r/r0)^2), V0 > 0 attractive.  The equation is
 * integrated inward by Numerov's method from the matching radius r_m, where the two outermost nodes carry the
 * solution of the point-Coulomb problem the caller names; the potential therefore acts between the channel radius
 * a_c and r_m.  A value: compared by value, never changed after construction, held by the pair that uses it.
 */

/// The radial shape of a central term.
enum class FormShape {
  WoodsSaxon,  ///< 1 / (1 + exp((r - R)/a)): radius_fm = R, diffuseness_fm = a
  Gaussian     ///< exp(-(r/r0)^2): width_fm = r0
};

/// A form factor f(r) (dimensionless, f -> 1 at the centre); the fields its shape does not use are 0.
struct FormFactor {
  FormShape shape = FormShape::WoodsSaxon;
  double radius_fm = 0.0;       ///< Woods-Saxon R (fm)
  double diffuseness_fm = 0.0;  ///< Woods-Saxon a (fm)
  double width_fm = 0.0;        ///< Gaussian r0 (fm)
};
bool operator==(const FormFactor &a, const FormFactor &b);
inline bool operator!=(const FormFactor &a, const FormFactor &b) { return !(a == b); }

/// The potential of one hybrid pair, with the numerical settings of its radial integration.
struct HybridPotential {
  double coulomb_radius_fm = 3.5;     ///< R_C of the uniformly charged sphere (fm)
  double central_depth_MeV = 0.0;     ///< V0 > 0 (MeV), attractive
  FormFactor central;                 ///< shape of the central term
  double matching_radius_fm = 15.0;   ///< r_m (fm): beyond it the solution is the point-Coulomb one
  double step_fm = 0.0275;            ///< Numerov step (fm), adjusted so that a_c is a node
};
bool operator==(const HybridPotential &a, const HybridPotential &b);
inline bool operator!=(const HybridPotential &a, const HybridPotential &b) { return !(a == b); }

/// A Woods-Saxon potential of depth V0 (MeV), radius R and diffuseness a (fm), with the default R_C, r_m and step.
/// Throws std::invalid_argument unless V0, R, a > 0.
HybridPotential WoodsSaxonPotential(double depth_MeV, double radius_fm, double diffuseness_fm);
/// A Gaussian potential of depth V0 (MeV) and width r0 (fm), with the default R_C, r_m and step.
/// Throws std::invalid_argument unless V0, r0 > 0.
HybridPotential GaussianPotential(double depth_MeV, double width_fm);

/// The nuclear term V_N(r) (MeV) at r (fm).
double NuclearTerm(const HybridPotential &p, double r_fm);
/// V(r) = V_C + V_N + centrifugal (MeV) at r (fm) for charges z1, z2, reduced mass mu (MeV) and l.
double RadialPotential(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double r_fm);

/*!
 * Every value the radial solution depends on, as numbers, for keys of tables and memos that must not mix
 * potentials: {0} without a potential, else {1, shape, V0, R, a, r0, R_C, r_m, step}.
 */
std::vector<double> PotentialKey(const std::optional<HybridPotential> &p);

/// The solution of the radial equation and its r-derivative (fm^-1) at the channel radius.
struct RadialValue {
  std::complex<double> u;
  std::complex<double> du_dr;
};

/*!
 * Integrates the radial equation inward from r_m to the channel radius a (fm, 0 < a <= r_m) at energy E (MeV, c.m.)
 * for charges z1, z2, reduced mass mu (MeV) and l.  The step is the one closest to p.step_fm that puts a on a node;
 * outer(r) gives the solution beyond r_m (H^+_l = G + iF for E > 0) at the two nodes past r_m.  Returns u(a) (the
 * node value) and du/dr(a) by the five-point central difference, both O(step^4).  The loop stops at the last node the
 * difference reads.
 */
RadialValue IntegrateInward(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double a_fm, double E_MeV,
                            const std::function<std::complex<double>(double r_fm)> &outer);

}  // namespace physics
}  // namespace azure2

#endif
