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

/*!
 * The interior background of a partial wave (docs/upgrade/potential_background.md sections 1-3): the R-function of
 * the potential inside the channel radius, R0(E) = 1/(l(E) - B_bg) - sum_p gamma_p^2/(E_p - E), added to the
 * diagonal R-matrix element of the channel, with the box states p named by their node numbers subtracted.
 */
struct InteriorBackground {
  double boundary = 0.0;               ///< B_bg (dimensionless): the boundary condition the background is defined with
  std::vector<int> subtract_nodes;     ///< node numbers of the box states whose poles are subtracted (rule 1, section 3)
  /// The interior direct capture of the background (potential_background.md section 5): the potential's interior wave,
  /// driven by the surface term of the channel's full wave at a, captured into the potential's bound state inside a
  /// (InteriorCaptureIntegral).  false: the background acts on the channel functions only.
  bool capture = false;
  /// The interior charge of the interior capture (potential_background.md section 14): a charge number q (units of e)
  /// added to the pair's own in the E1 (E2) operator inside the channel radius, q r (q r^2).  The operator of the
  /// relative motion carries Z1 (m2/M)^L + Z2 (-m1/M)^L, which nearly vanishes for E1 between nuclei of equal Z/A;
  /// the strength inside the radius is then the interior's (isospin mixing, the cluster's structure), a number the
  /// masses do not give.  0: the pair's operator alone.  Only with `capture`.
  double charge_e1 = 0.0;
  double charge_e2 = 0.0;
};
bool operator==(const InteriorBackground &a, const InteriorBackground &b);

/// What one partial wave of a pair has of its own: a depth in place of the central one and an interior background.
struct PartialWave {
  int l = 0;
  int two_j = -1;                          ///< 2J, or -1 (kNoSpinOrbit): every J of this l
  std::optional<double> depth_MeV;         ///< V0 of this partial wave (MeV, > 0 attractive); none: the central depth
  std::optional<double> width_fm;          ///< the Gaussian r0 of this partial wave (fm); none: the central one
  std::optional<InteriorBackground> background;  ///< none: no background (the channel functions of the exterior alone)
  /// `only bound_state` (potential_background.md section 5.9): the depth and width of this entry give the bound-state
  /// shape of the interior capture only; the channel functions of the wave are the point-Coulomb ones, as if the pair
  /// had no potential in it (BoundStateOnly).  false: the potential acts on the wave's channel functions.
  bool bound_state_only = false;
};
bool operator==(const PartialWave &a, const PartialWave &b);

/// The potential of one hybrid pair, with the numerical settings of its radial integration.
struct HybridPotential {
  double coulomb_radius_fm = 3.5;     ///< R_C of the uniformly charged sphere (fm)
  double central_depth_MeV = 0.0;     ///< V0 > 0 (MeV), attractive
  FormFactor central;                 ///< shape of the central term
  double matching_radius_fm = 15.0;   ///< r_m (fm): beyond it the solution is the point-Coulomb one
  double step_fm = 0.0275;            ///< Numerov step (fm), adjusted so that a_c is a node
  std::optional<SpinOrbitTerm> spin_orbit;  ///< none: a central potential
  std::vector<PartialWave> waves;     ///< per-partial-wave depths and backgrounds (empty: none; potential_background.md)
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
 * acts_on (0 light, 1 heavy)} when the potential has a spin-orbit term.
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

// ---- Partial waves and the interior background (docs/upgrade/potential_background.md) ----

/// The entry of p.waves for the channel (l, 2J): the one with this l and two_j, else the one with this l and every J
/// (two_j = -1); nullptr when there is none.
const PartialWave *FindWave(const HybridPotential &p, int l, int two_j);
/// Does a channel function of the pair depend on J (a spin-orbit term, or an entry of waves with a J)?
bool DependsOnJ(const HybridPotential &p);
/// The potential the radial equation of the channel (l, 2J) uses: p with the entry's depth and Gaussian width in place
/// of the central ones (p itself, field for field, when the entry has neither).
HybridPotential ForPartialWave(const HybridPotential &p, int l, int two_j);
/// ForPartialWave for the decaying solution W~ at or below threshold (the ANC <-> gamma factor, the bound final state of
/// external capture): throws std::runtime_error when the partial wave has an interior background, whose W~ is not
/// built (potential_background.md section 6; S(E <= 0) is, ChannelBoundShift).
HybridPotential BelowThresholdPotential(const HybridPotential &p, int l, int two_j);
/// The background of the channel (l, 2J), or nullptr.
const InteriorBackground *BackgroundOf(const HybridPotential &p, int l, int two_j);
/// Is the potential of the channel (l, 2J) for the bound-state shape only (PartialWave::bound_state_only)?  Then its
/// channel functions (waves, shift function, decaying solution) are the point-Coulomb ones.
bool BoundStateOnly(const HybridPotential &p, int l, int two_j);

/// The regular solution of the radial equation at the channel radius: u(a) and du/dr(a) (fm^-1), times exp(-log_scale).
struct InteriorValue {
  double u = 0.0;
  double du_dr = 0.0;
  double log_scale = 0.0;  ///< ln of the factor u and du_dr were divided by
};

/*!
 * Integrates the radial equation of p (the depth p.central_depth_MeV as given: call ForPartialWave first) outward from
 * the origin to a (fm) at energy E (MeV, c.m., any sign) for charges z1, z2, reduced mass mu (MeV), l and <2 l.i>, by
 * Numerov's method on the grid of IntegrateInward (the step closest to p.step_fm with a on a node), with u(0) = 0 and
 * u(h) = h^(l+1) (l = 1: the limit (g u)(0) = -2 u(h)/h^2 of the singular centrifugal term).  The derivative at a by
 * the same five-point difference (the loop runs two nodes past a).  `nodes`, if given, receives the number of sign
 * changes of u on the open interval (0, a).
 */
InteriorValue IntegrateOutward(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                               double a_fm, double E_MeV, int *nodes = nullptr);

/// l(E) = a u'(a)/u(a) of the regular solution (dimensionless), Eq. (2) of potential_background.md.
double InteriorLogDerivative(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                             double a_fm, double E_MeV);

/// R_pot(E; a, B) = u(a) / (a u'(a) - B u(a)) = 1/(l(E) - B) (dimensionless), Eq. (3).
double InteriorRFunction(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                         double E_MeV, double boundary);

/// A pole of R_pot: an eigenvalue of the potential in the box [0, a] with a u'(a) = B u(a), Eq. (4).
struct BoxState {
  int nodes = 0;            ///< sign changes of u_p on (0, a)
  double E_MeV = 0.0;       ///< E_p (MeV, c.m.)
  double gamma2_MeV = 0.0;  ///< gamma_p^2 = -1/(dl/dE) at E_p (MeV): the residue, R_pot ~ gamma_p^2/(E_p - E)
};

/*!
 * The box states of the channel with E_p in [E_lo, E_hi] (MeV), by a scan of f(E) = a u'(a) - B u(a) with the step
 * dE (MeV) and bisection to 1e-11 MeV; E_lo below the minimum of V on [0, a] is raised to it (no state lies below).
 * gamma_p^2 from a Richardson central difference of l(E) at E_p.  Ascending in E, so nodes = 0, 1, ... for the
 * states the window holds.
 */
std::vector<BoxState> InteriorBoxStates(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                                        double a_fm, double boundary, double E_lo_MeV, double E_hi_MeV, double dE_MeV = 0.25);

/*!
 * The background R0(E) of the channel (dimensionless): R_pot(E; a, B_bg) minus the poles gamma_p^2/(E_p - E) of the box
 * states of bg.subtract_nodes (Eqs. (3), (4)); p as for IntegrateOutward.  Throws std::runtime_error when a named
 * state is not found below 200 MeV.
 */
double BackgroundR(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                   double E_MeV, const InteriorBackground &bg);

/// The shift function with the background at E <= 0 (the same map as DressWaves on the decaying solution W~, real):
/// S~ = B + (S - B) / (1 - R0 (S - B)), S = rho W~'/W~.
double DressShift(double S, double R0, double boundary);

/*!
 * S(E = -B_MeV <= 0) of the channel (l, 2J) of the pair: HybridBoundShift with the partial wave's own potential
 * (ForPartialWave), changed by the partial wave's interior background when it has one (DressShift with R0 at E =
 * -B_MeV).  Continuous through threshold with the open-channel S~ of the dressed waves (potential_background.md
 * section 6).
 */
double ChannelBoundShift(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, int two_j, double two_l_dot_i,
                         double a_fm, double B_MeV);

/// The factor that takes the decaying solution W~(a) of the partial wave (l, 2J) at E <= 0 to the dressed one:
/// W~_d / W~ = 1 - R0 (S - B), S = a W~'/W~ undressed (HybridBoundShift), Eq. (9) on a real wave; 1 for a wave without
/// background.  The ANC <-> gamma factor of a level below threshold and the closed-channel capture read it.
double BackgroundWFactor(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, int two_j, double two_l_dot_i,
                         double a_fm, double E_MeV);

/// Half-width (MeV) of the window around a subtracted box pole inside which R0 and the capture kernel are the line
/// between their values at the window's edges (potential_background.md section 5.10).
constexpr double kSubtractedPoleWindow_MeV = 1e-3;

/// The box states of bg.subtract_nodes (InteriorBoxStates, searched upward in windows of 25 MeV up to 200 MeV), in the
/// order named; p as for IntegrateOutward.  Throws std::runtime_error when a named state is not found below 200 MeV.
std::vector<BoxState> SubtractedBoxStates(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                                          double a_fm, const InteriorBackground &bg);

/// Does some partial wave of p have a background with interior capture (InteriorBackground::capture)?
bool HasInteriorCapture(const HybridPotential &p);
/// The background of the channel (l, 2J) when it asks for interior capture, else nullptr.
const InteriorBackground *InteriorCaptureOf(const HybridPotential &p, int l, int two_j);

/*!
 * The interior direct-capture kernel of a background partial wave (potential_background.md section 5.2, Eq. (C5)):
 *
 *   K(E) = int_0^a phi_f(r) w(r) g_E(r) dr   (fm, times the unit of w),
 *   g_E(r) = R_pot(E) phi_E(r) - sum_{p in subtract} gamma_p^2 phi_p(r) / (E_p - E),
 *
 * with phi_E = u_E/u_E(a) the regular solution of the entrance partial wave's potential `pi` (ForPartialWave) at the
 * channel energy E (MeV, c.m., any sign), R_pot = 1/(l(E) - B_bg) (Eq. (3)), the subtracted box states p of bg with
 * phi_p their regular solutions at E_p normalised at a, and phi_f = u_f/u_f(a) the regular solution of the final
 * partial wave's potential `pf` at E_f < 0.  The interior wave of the entrance channel carried by the background is
 * g_E(r) D, D = (rho u' - B_bg u)(a) of the channel's full exterior wave (g_E(a) = R0, the background); the bound
 * state inside a is phi_f(r) times its exterior value at a.  w(r) is the multipole's radial weight (r^L, times the
 * finite-wavelength charge function).  Both regular solutions on the grid of IntegrateOutward (a on a node); the
 * integral by Simpson's rule (the 3/8 rule on the last three intervals when their number is odd): O(h^4).
 */
double InteriorCaptureIntegral(const HybridPotential &pi, int li, double two_l_dot_i_i, double E_MeV,
                               const InteriorBackground &bg, const HybridPotential &pf, int lf, double two_l_dot_i_f,
                               double Ef_MeV, int z1, int z2, double mu_MeV, double a_fm,
                               const std::function<double(double r_fm)> &weight);

/// InteriorCaptureIntegral at E without the window around a subtracted pole (its value at the window's edges).
double InteriorCaptureKernelAt(const HybridPotential &pi, int li, double two_l_dot_i_i, double E_MeV,
                               const InteriorBackground &bg, const HybridPotential &pf, int lf, double two_l_dot_i_f,
                               double Ef_MeV, int z1, int z2, double mu_MeV, double a_fm,
                               const std::function<double(double r_fm)> &weight);

/*!
 * The matching function of a bound state at the channel radius, continuous in E and in the depth: f(E) = sign W~(a)
 * (S u(a) - a u'(a)), u the regular solution (IntegrateOutward, scaled), S = rho W~'/W~ of the decaying solution of the
 * same potential (DecayingSolution).  f = 0 at an eigenvalue of the potential (Wronskian of u and W~ at a).
 */
double BoundMatching(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                     double E_MeV);

/// An eigenvalue of the potential below threshold (BoundMatching = 0), searched near E (MeV < 0).
struct PotentialEigenvalue {
  bool found = false;
  double E_MeV = 0.0;  ///< the eigenvalue nearest to the requested energy within the window (MeV, c.m.)
};
/// The eigenvalue of p nearest to E_MeV < 0 in [E_MeV - window, min(E_MeV + window, -1e-6)] (scan in steps of
/// window/100, bisection to 1e-10 MeV).
PotentialEigenvalue NearestEigenvalue(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i,
                                      double a_fm, double E_MeV, double window_MeV);
/// The central depth V0 (MeV) nearest to p.central_depth_MeV within +-30 % for which p binds at E_MeV < 0
/// (BoundMatching(E) = 0; scan in steps of 0.25 %, bisection to 1e-9 relative); found = false when none.
struct DepthAtEnergy {
  bool found = false;
  double depth_MeV = 0.0;
};
DepthAtEnergy DepthBindingAt(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                             double E_MeV);
/*!
 * The ANC (fm^-1/2) of the potential's own bound state at E < 0 (an eigenvalue of p, NearestEigenvalue), normalised
 * on [0, inf): 1 / (|W~(a)| (int_0^a (u/u(a))^2 dr + int_a^inf (W~/W~(a))^2 dr)^1/2), u the regular solution and W~ the
 * decaying solution of p (unit ANC; DecayingSolution).  The exterior integral by Simpson's rule in steps of 0.02 fm to
 * a + 25/kappa (W~^2 below e^-50 of its value at a).  Independent of the channel radius a when E is an eigenvalue.
 */
double PotentialBoundStateANC(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                              double E_MeV);
/// int_0^a (u(r)/u(a))^2 dr (fm) of the regular solution at E (Simpson on the grid of IntegrateOutward).
double InteriorNormIntegral(const HybridPotential &p, int z1, int z2, double mu_MeV, int l, double two_l_dot_i, double a_fm,
                            double E_MeV);

/// Waves at the channel radius: F, G and their rho-derivatives.
struct WavesAtRadius {
  double F = 0.0, Fp = 0.0, G = 0.0, Gp = 0.0;
};
/*!
 * The background as a change of the exterior waves at a (potential_background.md Eq. (9)): for u = F and u = G,
 *   u~ = u - R0 (rho u' - B u),   rho u~' = rho u' - B R0 (rho u' - B u),
 * a real map of determinant 1 (the Wronskian is kept), after which P~ = rho/|H~|^2, phi~ = arg H~ and
 * S~ = Re rho H~'/H~ are the channel functions of the levels on top of the background (Eq. (8)).
 */
WavesAtRadius DressWaves(const WavesAtRadius &w, double rho, double R0, double boundary);

}  // namespace physics
}  // namespace azure2

#endif
