#include "azure2/physics/transfer/GlobalOptical.h"

#include <algorithm>
#include <cmath>

/*
 * The formulas, transcribed from the papers as compiled in RIPL-3 (om-parameter-u.dat entries 6200, 2405/5405,
 * 7100/8100, 9100, 9600; Capote et al., Nucl. Data Sheets 110 (2009) 3107) and the FRONT front end of TWOFNR
 * (J.A. Tostevin, Surrey: Daehnick, Koning-Delaroche, Liang).  The independent check of every row is
 * physics/tests/transfer_global_optical_test.cpp, whose numbers come from tests/reference/thm_optical_reference.py
 * (it reads the RIPL-3 file itself, and FRONT21 for Daehnick and Liang).  E is the projectile lab energy (MeV),
 * A, Z, N those of the target, eta = (N - Z)/A.
 */

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

/// The potential of reduced radii r0 (fm) under R = r0 A_t^{1/3}: depths V, W, W_D (MeV), a W or W_D below 0 set
/// to 0 and a term of depth 0 without geometry.
OpticalPotential Central(int a_t, double v, double rv, double av, double w, double rw, double aw, double wd,
                         double rd, double ad, double rc) {
  OpticalPotential p;
  p.convention = RadiusConvention::TargetMass;
  p.target_mass_number = a_t;
  p.coulomb_r0_fm = rc;
  p.real_volume = {v, rv, av};
  w = std::max(0.0, w);
  wd = std::max(0.0, wd);
  p.imaginary_volume = w == 0.0 ? WoodsSaxonTerm() : WoodsSaxonTerm{w, rw, aw};
  p.imaginary_surface = wd == 0.0 ? WoodsSaxonTerm() : WoodsSaxonTerm{wd, rd, ad};
  return p;
}

// An & Cai 2006 (spin-orbit 3.557, 0.972, 1.011 not in the row).
OpticalPotential AnCai06(int, int, int z_t, int a_t, double E) {
  const double A = a_t, Z = z_t, a13 = std::cbrt(A);
  return Central(a_t, 91.85 - 0.249 * E + 0.000116 * E * E + 0.642 * Z / a13, 1.152 - 0.00776 / a13,
                 0.719 + 0.0126 * a13, 1.104 + 0.0622 * E, 1.305 + 0.0997 / a13, 0.855 - 0.100 * a13,
                 10.83 - 0.0306 * E, 1.334 + 0.152 / a13, 0.531 + 0.062 * a13, 1.303);
}

// Daehnick et al. 1980, the global set as FRONT21 codes it (spin-orbit not in the row).
OpticalPotential Daehnick80(int, int, int z_t, int a_t, double E) {
  const double A = a_t, Z = z_t, N = a_t - z_t, a13 = std::cbrt(A);
  const double g = std::exp(-(E / 100.0) * (E / 100.0));
  double ai = 0.53 + 0.07 * a13;
  for (double magic : {8.0, 20.0, 28.0, 50.0, 82.0, 126.0}) ai -= 0.04 * std::exp(-std::pow((magic - N) / 2.0, 2));
  return Central(a_t, 88.5 - 0.26 * E + 0.88 * Z / a13, 1.17, 0.709 + 0.0017 * E, (12.2 + 0.026 * E) * (1.0 - g),
                 1.325, ai, (12.2 + 0.026 * E) * g, 1.325, ai, 1.30);
}

// Koning & Delaroche 2003, global (spin-orbit not in the row); the projectile enters by its charge.
OpticalPotential KoningDelaroche03(int z_p, int, int z_t, int a_t, double E) {
  const double A = a_t, Z = z_t, N = a_t - z_t, a13 = std::cbrt(A), eta = (N - Z) / A;
  const bool proton = z_p == 1;
  const double ef = proton ? -8.4075 + 0.01378 * A : -11.2814 + 0.02646 * A;
  const double v1 = proton ? 59.30 + 21.0 * eta - 0.024 * A : 59.30 - 21.0 * eta - 0.024 * A;
  const double v2 = proton ? 0.007067 + 4.23e-6 * A : 0.007228 - 1.48e-6 * A;
  const double v3 = proton ? 1.729e-5 + 1.136e-8 * A : 1.994e-5 - 2.0e-8 * A;
  const double v4 = 7.0e-9;
  const double w1 = proton ? 14.667 + 0.009629 * A : 12.195 + 0.0167 * A;
  const double w2 = 73.55 + 0.0795 * A;
  const double d1 = proton ? 16.0 + 16.0 * eta : 16.0 - 16.0 * eta;
  const double d2 = 0.0180 + 0.003802 / (1.0 + std::exp((A - 156.0) / 8.0));
  const double d3 = 11.5;
  const double f = E - ef;
  double v = v1 * (1.0 - v2 * f + v3 * f * f - v4 * f * f * f), rc = 0.0;
  if (proton) {
    rc = 1.198 + 0.697 / (a13 * a13) + 12.994 / std::pow(a13, 5);
    const double vc = 1.73 / rc * Z / a13;
    v += vc * v1 * (v2 - 2.0 * v3 * f + 3.0 * v4 * f * f);
  }
  const double r = 1.3039 - 0.4054 / a13, a = 0.6778 - 1.487e-4 * A;
  return Central(a_t, v, r, a, w1 * f * f / (f * f + w2 * w2), r, a, d1 * f * f * std::exp(-d2 * f) / (f * f + d3 * d3),
                 1.3424 - 0.01585 * a13, proton ? 0.5187 + 5.205e-4 * A : 0.5446 - 1.656e-4 * A, rc);
}

// Becchetti & Greenlees 1971 (spin-orbit 2.5 not in the row).
OpticalPotential BecchettiGreenlees71(int z_p, int, int z_t, int a_t, double E) {
  const double A = a_t, Z = z_t, N = a_t - z_t, eta = (N - Z) / A;
  const bool helion = z_p == 2;
  return Central(a_t, helion ? 151.9 - 0.17 * E + 50.0 * eta : 165.0 - 0.17 * E - 6.4 * eta, 1.20, 0.72,
                 helion ? 41.7 - 0.33 * E + 44.0 * eta : 46.0 - 0.33 * E - 110.0 * eta, 1.40, helion ? 0.88 : 0.84,
                 0.0, 0.0, 0.0, 1.30);
}

// Liang, Li & Cai 2009 (spin-orbit not in the row).
OpticalPotential Liang09(int, int, int z_t, int a_t, double E) {
  const double A = a_t, Z = z_t, N = a_t - z_t, a13 = std::cbrt(A), eta = (N - Z) / A;
  return Central(a_t, 118.36 - 0.2071 * E + 6.3961e-5 * E * E + 26.001 * eta + 0.5668 * Z / a13,
                 1.1657 + 0.0401 / a13, 0.6641 + 0.0305 * a13, -6.8871 + 0.3115 * E - 6.8096e-4 * E * E,
                 1.4022 + 0.0418 / a13, 0.7732 + 0.0219 * a13, 20.119 - 0.1626 * E - 5.4067 * eta + 1.2087 * a13,
                 1.1802 + 0.0587 / a13, 0.6292 + 0.0657 * a13, 1.289);
}

// McFadden & Satchler 1966, energy independent.
OpticalPotential McFadden66(int, int, int, int a_t, double) {
  return Central(a_t, 185.0, 1.4, 0.52, 25.0, 1.4, 0.52, 0.0, 0.0, 0.0, 1.4);
}

// Avrigeanu, Hodgson & Avrigeanu 1994.
OpticalPotential Avrigeanu94(int, int, int z_t, int a_t, double E) {
  const double A = a_t, Z = z_t, a13 = std::cbrt(A);
  return Central(a_t, 101.1 - 0.248 * E + 6.052 * Z / a13, 1.245, 0.817 - 0.0085 * a13,
                 12.64 + 0.2 * E - 1.706 * a13, 1.57, 0.692 - 0.02 * a13, 0.0, 0.0, 0.0, 1.245);
}

const LightIon kNeutron{0, 1}, kProton{1, 1}, kDeuteron{1, 2}, kTriton{1, 3}, kHelion{2, 3}, kAlpha{2, 4};

const std::vector<GlobalOpticalSet> kSets = {
    {"ancai06", "d", {kDeuteron}, "An & Cai, PRC 73 (2006) 054605", 12.0, 238.0, 0.0, 183.0, AnCai06},
    {"daehnick80", "d", {kDeuteron}, "Daehnick, Childs & Vrcelj, PRC 21 (1980) 2253 (global set)", 27.0, 238.0, 11.8,
     90.0, Daehnick80},
    {"kd03", "n p", {kNeutron, kProton}, "Koning & Delaroche, NPA 713 (2003) 231 (global)", 24.0, 209.0, 0.001, 200.0,
     KoningDelaroche03},
    {"bg71", "t 3He", {kTriton, kHelion}, "Becchetti & Greenlees (1971)", 40.0, 208.0, 1.0, 40.0,
     BecchettiGreenlees71},
    {"liang09", "3He", {kHelion}, "Liang, Li & Cai, J. Phys. G 36 (2009) 085104", 9.0, 208.0, 0.0, 270.0, Liang09},
    {"mcfadden66", "4He", {kAlpha}, "McFadden & Satchler, NPA 84 (1966) 177", 16.0, 208.0, 1.0, 25.0, McFadden66},
    {"avrigeanu94", "4He", {kAlpha}, "Avrigeanu, Hodgson & Avrigeanu, PRC 49 (1994) 2136", 16.0, 250.0, 1.0, 73.0,
     Avrigeanu94},
};

}  // namespace

const std::vector<GlobalOpticalSet> &GlobalOpticalSets() { return kSets; }

int GlobalOpticalIndex(const std::string &name) {
  for (std::size_t i = 0; i < kSets.size(); ++i)
    if (name == kSets[i].name) return static_cast<int>(i);
  return -1;
}

std::string GlobalOpticalNames() {
  std::string out;
  for (const GlobalOpticalSet &s : kSets) out += (out.empty() ? "" : ", ") + std::string(s.name);
  return out;
}

bool Describes(const GlobalOpticalSet &set, int z_p, int a_p) {
  for (const LightIon &ion : set.describes)
    if (ion.z == z_p && ion.a == a_p) return true;
  return false;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
