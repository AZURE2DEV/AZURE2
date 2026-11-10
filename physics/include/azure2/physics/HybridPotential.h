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

/// The solution of the radial equation and its r-derivative (fm^-1) at the channel radius: u exp(log_scale).
struct RadialValue {
  std::complex<double> u;
  std::complex<double> du_dr;
  double log_scale = 0.0;  ///< ln of the factor u and du_dr were divided by (0 unless |u| passed 1e250 inward)
};

/*!
 * Integrates the radial equation inward from r_m to the channel radius a (fm, 0 < a <= r_m) at energy E (MeV, c.m.)
 * for charges z1, z2, reduced mass mu (MeV) and l.  The step is the one closest to p.step_fm that puts a on a node;
 * outer(r) gives the solution beyond r_m (H^+_l = G + iF for E > 0) at the two nodes past r_m.  Returns u(a) (the
 * node value) and du/dr(a) by the five-point central difference, both O(step^4).  The loop stops at the last node the
 * difference reads.  A
 * solution that grows past 1e250 inward (a decaying one far below the barrier) is rescaled, the factor kept in
 * log_scale; the outgoing wave of an open channel (eta <= 100: |G| of order e^(pi eta) <~ 1e137) is not.
 */
RadialValue IntegrateInward(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double a_fm, double E_MeV,
                            const std::function<std::complex<double>(double r_fm)> &outer);

/// The decaying solution W~_l of the radial equation at E = -B <= 0 at radius r, and its logarithmic derivative.
struct DecayingValue {
  double log_W = 0.0;     ///< ln |W~(r)|; W~ = W_{-eta_b,l+1/2}(2 kappa r) beyond r_m (unit ANC) for B > 0
  double sign = 1.0;      ///< sign of W~(r)
  double rho_dlog = 0.0;  ///< r W~'(r)/W~(r) = rho W~'/W~: the shift function S at E = -B when r = a
};

/*!
 * W~_l(r): the solution of the radial equation with the potential at E = -B (B >= 0, MeV) that is the point-Coulomb
 * decaying solution beyond r_m -- the Whittaker function W_{-eta_b,l+1/2}(2 kappa r) for B > 0 (unit ANC, so the ANC
 * of a bound state keeps its meaning), sqrt(r) K_{2l+1}(sqrt(8 c r)), c = Z1 Z2 alpha mu c^2/(hbar c), at B = 0
 * (r^-l without Coulomb; only rho_dlog has a meaning there) -- by IntegrateInward (stable: it grows inward).
 * Beyond r_m the point-Coulomb values.  Non-relativistic kinematics; mu in MeV.
 */
DecayingValue DecayingSolution(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double r_fm, double B_MeV);

/*!
 * The shift function S_l(-B) = a W~'(a)/W~(a) of a hybrid pair at E = -B <= 0, continuous down to and at threshold:
 * the exact ratio for eta_b <= 1000, and closer to threshold the linear interpolation between S(-B_min), eta_b(B_min)
 * = 1000, and S(0) of the zero-energy solution, as physics/ coulomb_bound_shift does for the point-Coulomb functions
 * (S - S(0) is linear in B there to about 1e-11 of S).
 */
double HybridBoundShift(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double a_fm, double B_MeV);

}  // namespace physics
}  // namespace azure2

#endif
