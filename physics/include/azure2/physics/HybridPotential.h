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

/// The particle of the pair whose spin i the spin-orbit term couples to l.
enum class SpinOrbitParticle {
  Light,  ///< particle 1 of the pair
  Heavy   ///< particle 2 of the pair
};

/*!
 * The spin-orbit term of a pair with one spinning particle (spin_orbit_potential.md section 2.1, Eq. (3); FRESCO's
 * convention): v_so(r) <2 l.i> added to V(r), with
 *
 *   v_so(r) = V_so (hbar / m_pi c)^2 (1/r) df/dr,   (hbar / m_pi c)^2 = 2.000 fm^2 (kSpinOrbitLengthSquared_fm2),
 *
 * f the Woods-Saxon form factor 1 / (1 + exp((r - R_so)/a_so)) and <2 l.i> = j(j+1) - l(l+1) - i(i+1) (TwoLDotI),
 * i the spin of the particle the term acts on; with one spinning particle j = J.  df/dr < 0, so V_so > 0 lowers
 * j = l + i.  For i = 1/2, 2 l.i = l.sigma: a strength quoted for l.s without the 2 is twice this V_so.
 */
struct SpinOrbitTerm {
  double depth_MeV = 0.0;  ///< V_so (MeV)
  FormFactor form;         ///< the Woods-Saxon form factor (R_so, a_so)
  SpinOrbitParticle acts_on = SpinOrbitParticle::Light;
};
bool operator==(const SpinOrbitTerm &a, const SpinOrbitTerm &b);
inline bool operator!=(const SpinOrbitTerm &a, const SpinOrbitTerm &b) { return !(a == b); }

/// (hbar / m_pi c)^2 of the spin-orbit term, 2.000 fm^2 exactly (FRESCO's CONLS).
constexpr double kSpinOrbitLengthSquared_fm2 = 2.0;

/// The 2J of the keys of channel functions for a channel without a spin-orbit term (every channel of a pair whose
/// potential has none): its functions do not depend on J.
constexpr int kNoSpinOrbit = -1;

/// <2 l.i> = j(j+1) - l(l+1) - i(i+1) for 2j = two_j, 2i = two_i (Eq. (4) of the note); l for j = l + 1/2 and
/// -(l + 1) for j = l - 1/2 when i = 1/2.  Throws std::invalid_argument unless |l - i| <= j <= l + i.
double TwoLDotI(int l, int two_j, int two_i);

/// The potential of one hybrid pair, with the numerical settings of its radial integration.
struct HybridPotential {
  double coulomb_radius_fm = 3.5;     ///< R_C of the uniformly charged sphere (fm)
  double central_depth_MeV = 0.0;     ///< V0 (MeV): > 0 attractive, < 0 repulsive
  FormFactor central;                 ///< shape of the central term
  double matching_radius_fm = 15.0;   ///< r_m (fm): beyond it the solution is the point-Coulomb one
  double step_fm = 0.0275;            ///< Numerov step (fm), adjusted so that a_c is a node
  std::optional<SpinOrbitTerm> spin_orbit;  ///< none: a central potential
  /// The boundary constant b of the potential's interior R-function R_pot = 1/(f - b) (interior_potential.md); none:
  /// the potential acts outside the channel radius only.
  std::optional<double> interior_boundary;
};
bool operator==(const HybridPotential &a, const HybridPotential &b);
inline bool operator!=(const HybridPotential &a, const HybridPotential &b) { return !(a == b); }

/// A Woods-Saxon potential of depth V0 (MeV), radius R and diffuseness a (fm), with the default R_C, r_m and step.
/// Throws std::invalid_argument unless V0, R, a > 0.
HybridPotential WoodsSaxonPotential(double depth_MeV, double radius_fm, double diffuseness_fm);
/// A Gaussian potential of depth V0 (MeV) and width r0 (fm), with the default R_C, r_m and step.
/// Throws std::invalid_argument unless V0, r0 > 0.
HybridPotential GaussianPotential(double depth_MeV, double width_fm);
/// A Woods-Saxon spin-orbit term of depth V_so (MeV), radius R_so and diffuseness a_so (fm) acting on `on`.
/// Throws std::invalid_argument unless V_so, R_so, a_so > 0.
SpinOrbitTerm WoodsSaxonSpinOrbit(double depth_MeV, double radius_fm, double diffuseness_fm, SpinOrbitParticle on);

/// V_C(r) (MeV) of a uniformly charged sphere of radius R_C (fm) for the charges z1, z2: z1 z2 e^2 (3 - r^2/R_C^2)/(2 R_C)
/// inside R_C (3 z1 z2 e^2/(2 R_C) at r = 0), z1 z2 e^2/r outside; R_C = 0 is a point charge.
double CoulombSphereTerm(int z1, int z2, double rc_fm, double r_fm);
/// The nuclear term V_N(r) (MeV) at r (fm).
double NuclearTerm(const HybridPotential &p, double r_fm);
/// v_so(r) = V_so (hbar/m_pi c)^2 (1/r) df/dr (MeV) of the term at r > 0 (fm); the factor <2 l.i> is the caller's.
double SpinOrbitRadialTerm(const SpinOrbitTerm &t, double r_fm);
/*!
 * V(r) = V_C + V_N + centrifugal + v_so(r) <2 l.i> (MeV) at r (fm) for charges z1, z2, reduced mass mu (MeV), l and
 * two_l_dot_i = <2 l.i> (TwoLDotI of the channel).  The spin-orbit term is added last, and not at all when the
 * potential has none or <2 l.i> = 0, so that a central potential is evaluated with the same operations as without it.
 */
double RadialPotential(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double r_fm);

/*!
 * Every value the radial solution depends on, as numbers, for keys of tables and memos that must not mix
 * potentials: {0} without a potential, else {1, shape, V0, R, a, r0, R_C, r_m, step}, followed by {V_so, R_so, a_so,
 * acts_on (0 light, 1 heavy)} when the potential has a spin-orbit term, and by {b} when it has an interior boundary
 * constant (so the length tells the four cases apart).
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
 * for charges z1, z2, reduced mass mu (MeV), l and <2 l.i> = two_l_dot_i (0 for a channel without spin-orbit term).  The step is the one closest to p.step_fm that puts a on a node;
 * outer(r) gives the solution beyond r_m (H^+_l = G + iF for E > 0) at the two nodes past r_m.  Returns u(a) (the
 * node value) and du/dr(a) by the five-point central difference, both O(step^4).  The loop stops at the last node the
 * difference reads.  A
 * solution that grows past 1e250 inward (a decaying one far below the barrier) is rescaled, the factor kept in
 * log_scale; the outgoing wave of an open channel (eta <= 100: |G| of order e^(pi eta) <~ 1e137) is not.
 */
RadialValue IntegrateInward(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm, double E_MeV,
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
DecayingValue DecayingSolution(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double r_fm,
                               double B_MeV);

/*!
 * The shift function S_l(-B) = a W~'(a)/W~(a) of a hybrid pair at E = -B <= 0, continuous down to and at threshold:
 * the exact ratio for eta_b <= 1000, and closer to threshold the linear interpolation between S(-B_min), eta_b(B_min)
 * = 1000, and S(0) of the zero-energy solution, as physics/ coulomb_bound_shift does for the point-Coulomb functions
 * (S - S(0) is linear in B there to about 1e-11 of S).
 */
double HybridBoundShift(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                        double B_MeV);

/// The regular solution u of the radial equation (u(0) = 0) at the channel radius, as the interior R-function of the
/// potential needs it (docs/upgrade/interior_potential.md).
struct RegularValue {
  double log_derivative = 0.0;  ///< f = a u'(a)/u(a), real at any E (infinite where u(a) = 0)
  double norm_fm = 0.0;         ///< int_0^a u^2 dr / u(a)^2 (fm): df/dE = -(2 mu a/hbar^2) norm_fm (note Eq. (6))
};

/*!
 * Integrates the radial equation outward (the stable direction of the regular solution) from the origin to the channel
 * radius a (fm) at energy E (MeV, c.m., any sign) for charges z1, z2, reduced mass mu (MeV), l and <2 l.i> =
 * two_l_dot_i.  Numerov on the grid of IntegrateInward (a on a node), started at the first two nodes from the series
 * u = r^(l+1) [1 + x/(2(2l+3)) + x^2/(8(2l+3)(2l+5)) + ...], x = (2 mu/hbar^2)(V - E) r^2 with V the potential without
 * its centrifugal term at the first node; u'(a) by the five-point central difference, the integral by Simpson's rule
 * on the nodes; all O(step^4).
 */
RegularValue RegularSolution(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                             double E_MeV);

/// S, P and exp(-i phi) of a channel, as the level and collision matrices read them.
struct ChannelTriple {
  double S = 0.0;
  double P = 0.0;
  std::complex<double> exp_minus_i_phi = 1.0;
};

/*!
 * The channel functions of a channel whose potential also acts inside the channel radius: with R = R_pot + the poles,
 * R_pot = 1/(f - b), the collision matrix keeps its form with (note Eqs. (4), (5))
 *
 *   S' + i P' = b + (L - b)(f - b)/(f - L),   L = S + i P,
 *   exp(-i phi') = exp(-i phi) conj(f - L)/|f - L| sgn(f - b)        (sgn(0) = +1),
 *
 * f the logarithmic derivative of the regular solution (RegularSolution), b the boundary constant of R_pot.  Without
 * levels U = exp(-2 i phi') is the scattering of the potential for any a and b.  f = +-infinity gives the unprimed
 * functions.  Below threshold (P = 0) S' is real and has a pole where f = S, a bound state of the potential.
 */
ChannelTriple InteriorChannelFunctions(const ChannelTriple &outer, double f, double b);

}  // namespace physics
}  // namespace azure2

#endif
