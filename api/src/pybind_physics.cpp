// The table of physics functions of the _azure2 module (docs/upgrade/pyazr3.md section 2.7): one row per function,
// with the Python name, the arguments (unit, domain), the returned fields (unit), the docstring and an adapter that
// calls one library function of physics/ or experiment/ on one element.  One template (Evaluate) applies a row
// element-wise to equal-length arrays and returns the values and a status per element; pyazr/physics.py builds one
// Python function per row from physics_functions() (signature, docstring, broadcasting, DomainError).
//
// The row types (arguments, returned fields, the domain check, their description) are pybind_rows.h, shared with
// the stopping table of pybind_stopping.cpp.
//
// A new function is one row here and a test naming it (tests/consistency check 41).  No physics in this file: an
// adapter converts arguments (doubles to integers, masses to a ChannelKinematics) and calls the library.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/Frames.h"
#include "azure2/physics/AstrophysicalFactor.h"
#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/CoulombScattering.h"
#include "azure2/physics/Elements.h"
#include "azure2/physics/ResonanceQuantities.h"
#include "azure2/physics/Widths.h"
#include "azure2/physics/coulomb/Coulomb.h"
#include "azure2/physics/coulomb/Whittaker.h"
#include "azure2/physics/transfer/ThmVertex.h"
#include "pybind_rows.h"

namespace py = pybind11;
namespace phys = azure2::physics;
namespace coul = azure2::physics::coulomb;
namespace frames = azure2::experiment;
using namespace azure2_rows;

namespace {

/// The adapter of a row: reads the row's arguments a[0 .. n_args), writes its fields o[0 .. n_outs), returns a status.
using Adapter = int (*)(const double *a, double *o);

struct PhysicsRow {
  std::string name;
  std::string summary;   ///< one line: what it computes, with units
  std::string doc;       ///< the formula, conventions, validity
  std::string library;   ///< the library function called (header: function)
  std::string axis;      ///< the argument a plot puts on its horizontal axis
  std::vector<ArgSpec> args;
  std::vector<OutSpec> outs;
  Adapter call;
  std::string group = "";  ///< the row's section of the table (set by BuildTable from its section lines)
};

// ---------------------------------------------------------------- argument conversions (no physics)

/// The kinematics of a pair from its charges and atomic masses (u), as the engine builds it (Kinematics::Channel).
phys::ChannelKinematics Pair(double z1, double z2, double m1_u, double m2_u, double relativistic) {
  phys::ChannelKinematics k;
  k.z1 = static_cast<int>(z1);
  k.z2 = static_cast<int>(z2);
  k.mu_u = m1_u * m2_u / (m1_u + m2_u);
  k.m1_MeV = m1_u * phys::kAtomicMassUnit_MeV;
  k.m2_MeV = m2_u * phys::kAtomicMassUnit_MeV;
  k.relativistic = relativistic != 0.0;
  return k;
}
int Int(double x) { return static_cast<int>(x); }
int FromCoulomb(coul::Status s) {
  if (s == coul::Status::ok) return kOk;
  return s == coul::Status::invalid_argument ? kJointly : kFailed;
}
/// The Coulomb waves of channel (l, pair, a) at E_c: a[0] = l, a[1..4] = z1, z2, m1, m2, a[5] = a, a[6] = E, a[7] = rel.
coul::CoulombWaves ChannelWaves(const double *a, double &rho) {
  const phys::ChannelKinematics k = Pair(a[1], a[2], a[3], a[4], a[7]);
  rho = phys::rho(k, a[5], a[6]);
  return coul::coulomb(Int(a[0]), phys::eta(k, a[6]), rho, coul::Method::coulsm);
}

// ---------------------------------------------------------------- the rows

const std::vector<std::string> kMethods = {"coulsm", "coul", "coulfg", "coulcc"};  // coul::Method in order

// Examples: p + 12C (1H and 12C atomic masses, AME2020) unless a row says otherwise.
constexpr double kMassH1 = 1.00782503, kMassC12 = 12.0, kMassHe4 = 4.00260325, kMassN15 = 15.0001089;

ArgSpec L() {
  return Integer("l", "orbital angular momentum").ge(0).le(10000).example(0).sweep(0, 10, Spacing::linear);
}
ArgSpec Z(const std::string &n, const std::string &which, double example) {
  return Integer(n, "charge number of " + which).ge(0).example(example);
}
ArgSpec M(const std::string &n, const std::string &which, double example) {
  return Real(n, "u", "atomic mass of " + which).gt(0).example(example);
}
ArgSpec Ecm() {
  return Real("e_cm", "MeV", "c.m. kinetic energy of the pair").gt(0).example(1).sweep(0.01, 10, Spacing::log);
}
ArgSpec Radius() { return Real("radius", "fm", "channel radius a").gt(0).example(5).sweep(1, 10, Spacing::linear); }
ArgSpec Rel() {
  return Flag("relativistic", "relativistic kinematics (the engine's option relativistic; Haeck et al., LA-UR-19-31084)",
              false);
}
std::vector<ArgSpec> Channel() {
  return {L(), Z("z1", "particle 1", 1), Z("z2", "particle 2", 6), M("m1_u", "particle 1", kMassH1),
          M("m2_u", "particle 2", kMassC12), Radius(), Ecm(), Rel()};
}
const char *const kChannelDoc =
    "The channel is the pair (z1, m1_u) + (z2, m2_u) in orbital l at radius a and c.m. energy E > 0 (an open "
    "channel); eta and rho = k a are the engine's (non-relativistic unless relativistic=1), the Coulomb functions "
    "F_l, G_l at rho those of coulomb_waves with method coulsm, the engine's routine (a session with a Coulomb table "
    "interpolates them to about 1e-7).";

std::vector<PhysicsRow> BuildTable() {
  std::vector<PhysicsRow> t;
  std::vector<std::pair<std::size_t, std::string>> sections;  // (first row, group) of each section of the table
  const auto section = [&](const char *group) { sections.emplace_back(t.size(), group); };
  section("kinematics");
  t.push_back({"sommerfeld_eta", "Sommerfeld parameter eta = Z1 Z2 e^2 / (hbar v) (dimensionless).",
               "eta = Z1 Z2 alpha (mu c^2 / 2E)^(1/2) with mu = m1 m2 / (m1 + m2), or the relativistic form.",
               "physics/ChannelFunctions.h: eta", "e_cm",
               {Z("z1", "particle 1", 1), Z("z2", "particle 2", 6), M("m1_u", "particle 1", kMassH1),
                M("m2_u", "particle 2", kMassC12), Ecm(), Rel()},
               {{"eta", "", "Sommerfeld parameter"}},
               [](const double *a, double *o) {
                 o[0] = phys::eta(Pair(a[0], a[1], a[2], a[3], a[5]), a[4]);
                 return kOk;
               }});
  t.push_back({"wave_number", "Wave number k of the relative motion (fm^-1).",
               "k = (2 mu c^2 E)^(1/2) / hbar c, or the relativistic form; rho = k a.",
               "physics/ChannelFunctions.h: rho (at a = 1 fm)", "e_cm",
               {M("m1_u", "particle 1", kMassH1), M("m2_u", "particle 2", kMassC12), Ecm(), Rel()},
               {{"k", "fm^-1", "wave number"}},
               [](const double *a, double *o) {
                 o[0] = phys::rho(Pair(0, 0, a[0], a[1], a[3]), 1.0, a[2]);
                 return kOk;
               }});
  t.push_back({"lab_to_cm_energy", "c.m. kinetic energy of a beam on a target at rest (MeV).",
               "E_cm = E_lab m2 / (m1 + m2), or sqrt(s) - m1 - m2 with s = (m1 + m2)^2 + 2 m2 E_lab (relativistic): the "
               "engine's conversion of a lab data point.",
               "experiment/Frames.h: BeamToCM", "e_lab",
               {M("m1_u", "the beam particle", kMassH1), M("m2_u", "the target particle", kMassC12),
                Real("e_lab", "MeV", "lab kinetic energy of the beam").ge(0).example(1).sweep(0, 10, Spacing::linear),
                Rel()},
               {{"e_cm", "MeV", "c.m. kinetic energy"}},
               [](const double *a, double *o) {
                 frames::Partition p;
                 p.m1 = a[0], p.m2 = a[1];
                 o[0] = frames::BeamToCM(p, a[2], a[2], a[2], a[3] != 0.0).eCM;
                 return kOk;
               }});
  t.push_back({"cm_to_lab_energy", "Lab kinetic energy of the beam for a c.m. energy (MeV).",
               "The inverse of lab_to_cm_energy: E_lab = E_cm (m1 + m2) / m2, or the relativistic inverse.",
               "experiment/Frames.h: BeamLabFromCM", "e_cm",
               {M("m1_u", "the beam particle", kMassH1), M("m2_u", "the target particle", kMassC12),
                Real("e_cm", "MeV", "c.m. kinetic energy of the pair").ge(0).example(1).sweep(0, 10, Spacing::linear),
                Rel()},
               {{"e_lab", "MeV", "lab kinetic energy of the beam"}},
               [](const double *a, double *o) {
                 frames::Partition p;
                 p.m1 = a[0], p.m2 = a[1];
                 o[0] = frames::BeamLabFromCM(p, a[2], a[3] != 0.0);
                 return kOk;
               }});
  t.push_back({"elastic_lab_to_cm_angle",
               "c.m. angle and Jacobian dOmega_lab/dOmega_cm of elastic scattering at a lab angle (deg).",
               "theta_cm = theta_lab + asin((m1/m2) sin theta_lab), non-relativistic, forward branch; J = "
               "dOmega_lab / dOmega_cm takes a lab differential cross section to the c.m. one, dsigma/dOmega_cm = J "
               "dsigma/dOmega_lab.  Status -2 beyond the maximum lab angle asin(m2/m1) of a beam heavier than the "
               "target.",
               "experiment/Frames.h: ElasticLabToCM, ParticleJacobian", "theta_lab",
               {M("m1_u", "the beam (detected) particle", kMassH1), M("m2_u", "the target particle", kMassC12),
                Real("theta_lab", "deg", "lab angle of the detected beam particle").ge(0).le(180).example(90)
                    .sweep(0, 180, Spacing::linear)},
               {{"theta_cm", "deg", "c.m. angle"}, {"jacobian", "", "dOmega_lab / dOmega_cm"}},
               [](const double *a, double *o) {
                 frames::Partition p;
                 p.m1 = a[0], p.m2 = a[1];
                 frames::AngleLimit lim;
                 if (!frames::ElasticLabToCM(p, a[2], o[0], lim)) return kJointly;
                 // the Jacobian of an elastic map does not depend on the beam energy (Q = 0): 1 MeV
                 o[1] = frames::ParticleJacobian(p, p, 1.0, a[2], o[0]);
                 return kOk;
               }});
  t.push_back({"reaction_lab_to_cm_angle",
               "c.m. angle and Jacobian dOmega_lab/dOmega_cm of the detected particle of a reaction (deg).",
               "1 + 2 -> 3 + 4 with 3 detected at theta_lab, beam energy E_lab, Q value Q (Marion and Young, Nuclear "
               "Reaction Analysis, 1968), non-relativistic, forward branch.  Status -2 below the exit threshold and "
               "beyond the maximum lab angle.",
               "experiment/Frames.h: ReactionLabToCM, ParticleJacobian", "theta_lab",
               {M("m1_u", "the beam particle", kMassH1), M("m2_u", "the target particle", kMassN15),
                M("m3_u", "the detected particle", kMassHe4), M("m4_u", "the residual particle", kMassC12),
                Real("q_value", "MeV", "Q value of the reaction").example(4.965),
                Real("e_lab", "MeV", "lab kinetic energy of the beam").gt(0).example(1).sweep(0.1, 10, Spacing::log),
                Real("theta_lab", "deg", "lab angle of the detected particle").ge(0).le(180).example(90)
                    .sweep(0, 180, Spacing::linear)},
               {{"theta_cm", "deg", "c.m. angle"}, {"jacobian", "", "dOmega_lab / dOmega_cm"}},
               [](const double *a, double *o) {
                 frames::Partition in, out;
                 in.m1 = a[0], in.m2 = a[1], in.sepE = a[4];
                 out.m1 = a[2], out.m2 = a[3];
                 frames::AngleLimit lim;
                 if (!frames::ReactionLabToCM(in, out, a[5], a[6], o[0], lim)) return kJointly;
                 o[1] = frames::ParticleJacobian(in, out, a[5], a[6], o[0]);
                 return kOk;
               }});
  section("Coulomb functions");
  t.push_back({"coulomb_waves", "Coulomb functions F_l, G_l, their rho-derivatives and sigma_l at (l, eta, rho).",
               "Solutions of u'' + [1 - 2 eta/rho - l(l+1)/rho^2] u = 0, F regular, G irregular, F ~ sin(theta), G ~ "
               "cos(theta), theta = rho - eta ln(2 rho) - l pi/2 + sigma_l (DLMF 33.2); Wronskian F' G - F G' = 1.  "
               "method names the algorithm of pycoulomb (coulsm: the engine's).",
               "physics/coulomb/Coulomb.h: coulomb", "rho",
               {L(), Real("eta", "", "Sommerfeld parameter (0 neutral, > 0 repulsive)").ge(0).example(0.5)
                         .sweep(0, 10, Spacing::linear),
                Real("rho", "", "k r").gt(0).example(5).sweep(0.01, 20, Spacing::log),
                Choice("method", "algorithm", kMethods, 0)},
               {{"F", "", "regular function F_l"}, {"Fp", "", "dF_l/drho"}, {"G", "", "irregular function G_l"},
                {"Gp", "", "dG_l/drho"}, {"sigma", "rad", "Coulomb phase sigma_l"}},
               [](const double *a, double *o) {
                 const coul::CoulombWaves w = coul::coulomb(Int(a[0]), a[1], a[2], static_cast<coul::Method>(Int(a[3])));
                 o[0] = w.F, o[1] = w.Fp, o[2] = w.G, o[3] = w.Gp, o[4] = w.sigma;
                 return FromCoulomb(w.status);
               }});
  t.push_back({"coulomb_phase", "Coulomb phase shift sigma_l(eta) = arg Gamma(l + 1 + i eta) (rad).",
               "The continuous branch, not reduced modulo 2 pi.", "physics/coulomb/Coulomb.h: coulomb_phase", "eta",
               {L(), Real("eta", "", "Sommerfeld parameter").example(1).sweep(-10, 10, Spacing::linear)},
               {{"sigma", "rad", "Coulomb phase sigma_l"}},
               [](const double *a, double *o) {
                 o[0] = coul::coulomb_phase(Int(a[0]), a[1]);
                 return std::isfinite(o[0]) ? kOk : kFailed;
               }});
  t.push_back({"whittaker_wave", "Whittaker function W_l(rho) = W_{-eta_b, l+1/2}(2 rho) of a bound channel.",
               "The decaying solution of u'' - [1 + 2 eta_b/rho + l(l+1)/rho^2] u = 0, W ~ exp(-rho) (2 rho)^(-eta_b) "
               "for large rho (unit ANC); no GSL (relative error below 1e-12 against mpmath on l = 0-10, eta_b = "
               "0-1000, rho = 1e-4-200).  W = mantissa 10^log10_scale; the field W is 0 or inf where W leaves the "
               "double range (mantissa and log10_scale keep the value).",
               "physics/coulomb/Whittaker.h: whittaker_scaled", "rho",
               {L(), Real("eta_b", "", "bound-state Sommerfeld parameter Z1 Z2 e^2 mu / (hbar^2 kappa)").ge(0)
                         .example(1).sweep(0, 10, Spacing::linear),
                Real("rho", "", "kappa r").gt(0).example(5).sweep(0.01, 20, Spacing::log)},
               {{"mantissa", "", "W / 10^log10_scale, in [1, 10)"}, {"log10_scale", "", "decimal exponent of W"},
                {"log10_W", "", "log10 W"}, {"W", "", "W_l(rho)"}, {"rho_dlog", "", "rho W'/W, the shift function"}},
               [](const double *a, double *o) {
                 const coul::WhittakerScaled w = coul::whittaker_scaled(Int(a[0]), a[1], a[2]);
                 o[0] = w.mantissa, o[1] = static_cast<double>(w.log10_scale), o[4] = w.rho_dlog;
                 o[2] = std::log10(w.mantissa) + o[1];
                 o[3] = std::pow(10.0, o[2]);
                 return FromCoulomb(w.status);
               }});
  section("channel functions");
  t.push_back({"penetrability", "Penetrability P_l = rho / (F_l^2 + G_l^2) of an open channel (dimensionless).",
               kChannelDoc, "physics/ChannelFunctions.h: penetrability", "e_cm", Channel(),
               {{"P", "", "penetrability"}},
               [](const double *a, double *o) {
                 double r = 0.0;
                 const coul::CoulombWaves w = ChannelWaves(a, r);
                 o[0] = phys::penetrability(w, r);
                 return FromCoulomb(w.status);
               }});
  t.push_back({"shift", "Shift function S_l = rho (F F' + G G') / (F^2 + G^2) of an open channel (dimensionless).",
               kChannelDoc, "physics/ChannelFunctions.h: shift", "e_cm", Channel(),
               {{"S", "", "shift function"}},
               [](const double *a, double *o) {
                 double r = 0.0;
                 const coul::CoulombWaves w = ChannelWaves(a, r);
                 o[0] = phys::shift(w, r);
                 return FromCoulomb(w.status);
               }});
  t.push_back({"hard_sphere_phase", "Hard-sphere phase phi_l = arg(G_l + i F_l) of an open channel (rad).",
               std::string(kChannelDoc) + "  exp(-i phi) = (G - i F) / (F^2 + G^2)^(1/2) (Lane and Thomas); the "
               "hard-sphere phase shift is -phi (the session's delta_hs).",
               "physics/ChannelFunctions.h: exp_minus_i_phi", "e_cm", Channel(),
               {{"phi", "rad", "hard-sphere phase, in (-pi, pi]"}},
               [](const double *a, double *o) {
                 double r = 0.0;
                 const coul::CoulombWaves w = ChannelWaves(a, r);
                 o[0] = -std::arg(phys::exp_minus_i_phi(w));
                 return FromCoulomb(w.status);
               }});
  t.push_back({"shift_derivative", "Energy derivative dS_l/dE of the shift function of an open channel (MeV^-1).",
               std::string(kChannelDoc) + "  Richardson's extrapolation of central differences of S with h = 1e-3 E "
               "(about 1e-8 relative), the engine's.",
               "physics/ChannelFunctions.h: channel_functions (dS_dE)", "e_cm", Channel(),
               {{"dS_dE", "MeV^-1", "dS/dE"}},
               [](const double *a, double *o) {
                 const phys::OpenChannelValues v =
                     phys::channel_functions(Pair(a[1], a[2], a[3], a[4], a[7]), Int(a[0]), a[5], a[6]);
                 o[0] = v.dS_dE;
                 return FromCoulomb(v.status);
               }});
  section("THM");
  t.push_back({"thm_vertex_function", "Plane-wave surface vertex M_l(p) of a Trojan-horse level (dimensionless).",
               "M_l(p) = (B - 1) j_l(pa) - pa j_l'(pa) (docs/upgrade/transfer.md eq. 8.3, THM.md section 1): the "
               "surface term of the vertex int_0^a w V u dr for the plane wave w = r j_l(pr) of the transferred "
               "particle, p its momentum, a the channel radius, B the level's boundary value.",
               "physics/transfer/ThmVertex.h: PlaneWaveSurfaceVertex", "p",
               {L(), Real("p", "fm^-1", "momentum of the transferred particle").gt(0).example(0.5)
                         .sweep(0.01, 3, Spacing::linear),
                Radius(),
                Real("boundary", "", "boundary value B = a u'(a)/u(a) of the level").example(0)
                    .sweep(-5, 5, Spacing::linear)},
               {{"M", "", "plane-wave surface vertex M_l(p)"}},
               [](const double *a, double *o) {
                 o[0] = azure2::physics::transfer::PlaneWaveSurfaceVertex(Int(a[0]), a[1], a[2], a[3]);
                 return std::isfinite(o[0]) ? kOk : kFailed;
               }});
  section("widths");
  t.push_back({"wigner_limit", "Wigner limit gamma_W^2 = 3 hbar^2 / (2 mu a^2) of a reduced width (MeV).",
               "The engine's convention (Parameter.wigner_limit; the Teichmann-Wigner sum-rule unit); hbar^2 / (mu a^2) "
               "is 2/3 of it.  theta^2 = gamma^2 / gamma_W^2.",
               "physics/Widths.h: wigner_limit", "radius",
               {M("m1_u", "particle 1", kMassH1), M("m2_u", "particle 2", kMassC12), Radius()},
               {{"gamma2_W", "MeV", "Wigner limit of gamma^2"}},
               [](const double *a, double *o) {
                 o[0] = phys::wigner_limit(a[0] * a[1] / (a[0] + a[1]), a[2]);
                 return kOk;
               }});
  t.push_back({"single_particle_width", "Single-particle width Gamma_sp = 2 P_l(E) gamma_W^2 of an open channel (MeV).",
               std::string(kChannelDoc) + "  gamma_W^2 is wigner_limit; Gamma / Gamma_sp is the dimensionless reduced "
               "width theta^2 of a narrow level.",
               "physics/Widths.h: single_particle_width", "e_cm", Channel(),
               {{"gamma_sp", "MeV", "single-particle width"}},
               [](const double *a, double *o) {
                 o[0] = phys::single_particle_width(Pair(a[1], a[2], a[3], a[4], a[7]), Int(a[0]), a[5], a[6]);
                 return std::isfinite(o[0]) ? kOk : kFailed;
               }});
  t.push_back({"photon_penetrability", "Photon penetrability (E_gamma / hbar c)^(2L+1) of a multipole L (fm^-(2L+1)).",
               "The factor of a photon channel's width Gamma = 2 gamma^2 P (Lane and Thomas section XIII), the engine's.",
               "physics/Widths.h: photon_penetrability", "e_gamma",
               {Real("e_gamma", "MeV", "photon energy").gt(0).example(1).sweep(0.1, 20, Spacing::log),
                Integer("multipolarity", "multipolarity L").ge(1).le(10).example(1).sweep(1, 5, Spacing::linear)},
               {{"P", "fm^-(2L+1)", "photon penetrability"}},
               [](const double *a, double *o) {
                 o[0] = phys::photon_penetrability(a[0], Int(a[1]));
                 return kOk;
               }});
  section("resonance quantities");
  t.push_back({"statistical_factor", "Statistical factor omega = (2J+1)(1+delta_aA)/((2j_a+1)(2j_A+1)).",
               "Of a resonance strength (data_kinds.md eq. 8); spins doubled.", "physics/ResonanceQuantities.h: "
               "statistical_factor", "two_J",
               {Integer("two_J", "2J of the level").ge(0).example(1).sweep(0, 10, Spacing::linear),
                Integer("two_ja", "2j of entrance particle a").ge(0).example(1),
                Integer("two_jA", "2j of entrance particle A").ge(0).example(0),
                Flag("identical", "the entrance particles are identical", false)},
               {{"omega", "", "statistical factor"}},
               [](const double *a, double *o) {
                 o[0] = phys::statistical_factor(Int(a[0]), Int(a[1]), Int(a[2]), a[3] != 0.0);
                 return kOk;
               }});
  t.push_back({"resonance_strength", "Resonance strength omega gamma = omega Gamma_a Gamma_b / Gamma, with derivatives.",
               "From observed widths in one energy unit (the result carries it), the three widths independent "
               "(data_kinds.md eq. 8).  Status 4 for Gamma = 0.", "physics/ResonanceQuantities.h: resonance_strength",
               "gamma_a",
               {Real("omega", "", "statistical factor").example(1),
                Real("gamma_a", "width unit", "entrance width Gamma_a").example(0.1).sweep(0.001, 1, Spacing::log),
                Real("gamma_b", "width unit", "exit width Gamma_b").example(0.01).sweep(0.001, 1, Spacing::log),
                Real("gamma", "width unit", "total width Gamma").example(1).sweep(0.1, 10, Spacing::log)},
               {{"value", "width unit", "omega gamma"}, {"d_gamma_a", "", "d/dGamma_a"}, {"d_gamma_b", "", "d/dGamma_b"},
                {"d_gamma", "", "d/dGamma"}},
               [](const double *a, double *o) {
                 if (a[3] == 0.0) return 4;
                 const phys::ResonanceStrength r = phys::resonance_strength(a[0], a[1], a[2], a[3]);
                 o[0] = r.value, o[1] = r.d_entrance, o[2] = r.d_exit, o[3] = r.d_total;
                 return kOk;
               }});
  t.push_back({"branching_ratio", "Branching ratio Gamma_b / Gamma, with derivatives.", "Status 2 for Gamma = 0.",
               "physics/ResonanceQuantities.h: branching_ratio", "gamma_b",
               {Real("gamma_b", "width unit", "exit width Gamma_b").example(0.5).sweep(0, 1, Spacing::linear),
                Real("gamma", "width unit", "total width Gamma").example(1).sweep(0.1, 10, Spacing::log)},
               {{"value", "", "Gamma_b / Gamma"}, {"d_gamma_b", "width unit^-1", "d/dGamma_b"},
                {"d_gamma", "width unit^-1", "d/dGamma"}},
               [](const double *a, double *o) {
                 if (a[1] == 0.0) return 2;
                 const phys::BranchingRatio r = phys::branching_ratio(a[0], a[1]);
                 o[0] = r.value, o[1] = r.d_exit, o[2] = r.d_total;
                 return kOk;
               }});
  section("Coulomb scattering and S factor");
  t.push_back({"rutherford_cross_section", "Rutherford cross section (b/sr, c.m.) and its energy derivative.",
               "With the engine's eta and k (the Coulomb part of its elastic cross section).",
               "physics/CoulombScattering.h: rutherford_cross_section", "theta_cm",
               {Z("z1", "particle 1", 2), Z("z2", "particle 2", 6), M("m1_u", "particle 1", kMassHe4),
                M("m2_u", "particle 2", kMassC12), Ecm().example(2),
                Real("theta_cm", "deg", "c.m. angle").gt(0).lt(180).example(90).sweep(5, 175, Spacing::linear), Rel()},
               {{"value", "b/sr", "dsigma/dOmega_cm"}, {"d_dE", "b/sr/MeV", "d/dE at fixed angle"}},
               [](const double *a, double *o) {
                 const phys::CoulombCrossSection c =
                     phys::rutherford_cross_section(Pair(a[0], a[1], a[2], a[3], a[6]), a[4], a[5]);
                 o[0] = c.value, o[1] = c.d_dE;
                 return kOk;
               }});
  t.push_back({"mott_cross_section", "Mott cross section of two identical particles (b/sr, c.m.), with d/dE.",
               "The symmetrised Coulomb amplitude with the interference weight (-1)^(2s)/(2s+1).",
               "physics/CoulombScattering.h: mott_cross_section", "theta_cm",
               {Z("z", "each particle", 6), M("m_u", "each particle", kMassC12),
                Integer("two_s", "2s, s the spin of each particle").ge(0).example(0), Ecm().example(2),
                Real("theta_cm", "deg", "c.m. angle").gt(0).lt(180).example(90).sweep(5, 175, Spacing::linear), Rel()},
               {{"value", "b/sr", "dsigma/dOmega_cm"}, {"d_dE", "b/sr/MeV", "d/dE at fixed angle"}},
               [](const double *a, double *o) {
                 const phys::CoulombCrossSection c =
                     phys::mott_cross_section(Pair(a[0], a[0], a[1], a[1], a[5]), Int(a[2]), a[3], a[4]);
                 o[0] = c.value, o[1] = c.d_dE;
                 return kOk;
               }});
  t.push_back({"sfactor_conversion", "E exp(2 pi eta(E)) (MeV): S = sigma x this (S in MeV b, sigma in b).",
               "The engine's quantity sfactor; eta non-relativistic.  Status -2 for an uncharged pair.",
               "physics/AstrophysicalFactor.h: sfactor_conversion", "e_cm",
               {Z("z1", "particle 1", 1), Z("z2", "particle 2", 6), M("m1_u", "particle 1", kMassH1),
                M("m2_u", "particle 2", kMassC12), Ecm()},
               {{"factor", "MeV", "E exp(2 pi eta)"}},
               [](const double *a, double *o) {
                 o[0] = phys::sfactor_conversion(Pair(a[0], a[1], a[2], a[3], 0.0), a[4]);
                 return std::isfinite(o[0]) ? kOk : kJointly;
               }});
  t.push_back({"modified_sfactor_conversion", "E exp(2 pi eta(E) + g E) (MeV): S* = sigma x this.",
               "The modified S factor of heavy-ion fusion (Patterson, Winkler and Zaidins, Astrophys. J. 157 (1969) "
               "367), the engine's quantity modified_sfactor; g = 0.46 MeV^-1 for 12C + 12C.  Status -2 for an "
               "uncharged pair.",
               "physics/AstrophysicalFactor.h: modified_sfactor_conversion", "e_cm",
               {Z("z1", "particle 1", 6), Z("z2", "particle 2", 6), M("m1_u", "particle 1", kMassC12),
                M("m2_u", "particle 2", kMassC12), Real("g_per_MeV", "MeV^-1", "the constant g").example(0.46),
                Ecm().example(3).sweep(1, 8, Spacing::linear)},
               {{"factor", "MeV", "E exp(2 pi eta + g E)"}},
               [](const double *a, double *o) {
                 o[0] = phys::modified_sfactor_conversion(Pair(a[0], a[1], a[2], a[3], 0.0), a[5], a[4]);
                 return std::isfinite(o[0]) ? kOk : kJointly;
               }});
  for (std::size_t k = 0; k < sections.size(); k++) {
    const std::size_t end = k + 1 < sections.size() ? sections[k + 1].first : t.size();
    for (std::size_t i = sections[k].first; i < end; i++) t[i].group = sections[k].second;
  }
  return t;
}

const std::vector<PhysicsRow> &Table() {
  static const std::vector<PhysicsRow> table = BuildTable();  // immutable after the first call
  return table;
}

const PhysicsRow &Row(const std::string &name) {
  for (const PhysicsRow &r : Table())
    if (r.name == name) return r;
  throw std::invalid_argument("no physics function named '" + name + "' (pyazr.kinds()['physics'] lists them)");
}

// ---------------------------------------------------------------- the template

using Column = py::array_t<double, py::array::c_style | py::array::forcecast>;

/// The row `name` on equal-length 1-d arrays, one per argument in the row's order: (values of shape (n_outs, n),
/// status of shape (n,)), NaN where the status is not 0.
py::tuple Evaluate(const std::string &name, const std::vector<Column> &columns) {
  const PhysicsRow &row = Row(name);
  const std::size_t na = row.args.size(), no = row.outs.size();
  if (columns.size() != na)
    throw std::invalid_argument(name + ": " + std::to_string(na) + " argument arrays expected, " +
                                std::to_string(columns.size()) + " given");
  const py::ssize_t n = na ? columns[0].size() : 0;
  std::vector<const double *> in(na);
  for (std::size_t k = 0; k < na; k++) {
    if (columns[k].ndim() != 1 || columns[k].size() != n)
      throw std::invalid_argument(name + ": the argument arrays must be 1-d of one length (broadcast them first)");
    in[k] = columns[k].data();
  }
  py::array_t<double> values({static_cast<py::ssize_t>(no), n});
  py::array_t<int> status(n);
  double *v = values.mutable_data();
  int *st = status.mutable_data();
  {
    py::gil_scoped_release release;
    std::vector<double> a(na), o(no);
    for (py::ssize_t i = 0; i < n; i++) {
      for (std::size_t k = 0; k < na; k++) a[k] = in[k][i];
      int s = CheckDomain(row.args, a.data());
      if (s == kOk) s = row.call(a.data(), o.data());
      for (std::size_t j = 0; j < no; j++) v[j * n + i] = s == kOk ? o[j] : kNaN;
      st[i] = s;
    }
  }
  return py::make_tuple(values, status);
}

py::list Describe() {
  py::list rows;
  for (const PhysicsRow &r : Table()) {
    py::list args, outs;
    for (const ArgSpec &s : r.args) args.append(DescribeArgument(s));
    for (const OutSpec &s : r.outs) outs.append(DescribeReturn(s));
    py::dict row;
    row["name"] = r.name;
    row["summary"] = r.summary;
    row["doc"] = r.doc;
    row["library"] = r.library;
    row["axis"] = r.axis;
    row["group"] = r.group;
    row["arguments"] = args;
    row["returns"] = outs;
    rows.append(row);
  }
  return rows;
}

}  // namespace

void BindStopping(py::module_ &parent);  // pybind_stopping.cpp: stopping/, the submodule _azure2.stopping
void BindLevelDoors(py::module_ &m);     // pybind_level_doors.cpp: one block's formal <-> observed levels (pyazr.endf)
void BindCollisionMatrix(py::module_ &m);  // pybind_collision_matrix.cpp: one block's A, U and T at one energy

void BindPhysics(py::module_ &m) {
  BindStopping(m);
  BindLevelDoors(m);
  BindCollisionMatrix(m);
  m.def("physics_functions", &Describe,
        "The table of physics functions, one dict per row: name, summary, doc, library (the function called), axis "
        "(the argument a plot puts on its horizontal axis), group (its section of the table), arguments (name, unit, "
        "doc, kind real / integer / flag / choice, domain in words, lower and upper (None: unbounded), lower_open, "
        "upper_open, choices, default or None, example (a typical value) or None, sweep {lower, upper, spacing "
        "linear / log} (the range of a first plot over it) or None) and returns (name, unit, doc); plain JSON data.  pyazr.physics builds one Python function per row from it.");
  m.def("physics_evaluate", &Evaluate, py::arg("name"), py::arg("arguments"),
        "Row `name` element-wise on equal-length 1-d float arrays, one per argument in the row's order (a choice as "
        "its index, a flag as 0 or 1): (values, status), values of shape (n_returns, n), NaN where status != 0; "
        "status 0 computed, k > 0 argument k (1-based) outside its domain, -2 the arguments together outside the "
        "domain, -1 not converged or failed.");
  m.def(
      "element_symbols",
      []() {
        std::vector<std::string> out;
        for (int z = 0; z <= phys::kMaxElementZ; z++) out.push_back(phys::ElementSymbol(z));
        return out;
      },
      "Element symbols by charge number, index 0 (the neutron, \"n\") to 118 (physics/Elements.h, the one table).");
}
