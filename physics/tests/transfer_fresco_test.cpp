// transfer_fresco (docs/upgrade/transfer.md section 9): physics/transfer against FRESCO 3.5, from the reference files
// of tests/data/fresco/<case>/ (make_fresco_benchmarks.py; FRESCO is not run here).  The inputs below are those of
// each case's job.nin; the inputs are compared before the results, so that a disagreement is never one of inputs:
//
//  1. Constants and kinematics: FRESCO's 2u/hbar^2 and e^2 (set through its unitmass and finec) and, per partition,
//     E_cm, k and eta against ours, to the digits FRESCO prints.
//  2. Potentials: FRESCO's listing (fort.34) of every component at r = 1 ... 10 fm against TermsAt (the Coulomb term
//     per unit Z1 Z2, the spin-orbit term before 2 l.s), to 1e-6 relative.
//  3. Bound states: n-p (deuteron), p-15N (p1/2, spin-orbit, Coulomb), alpha-d (6Li, one node), alpha-12C (two nodes):
//     the depth (FRESCO's scale factor times the input depth, 4 decimals), rms radius, D0 and the ANC (FRESCO's
//     plateau of u/W at 12 fm, fort.46), to the digits printed.
//  4. Elastic scattering (p + 15N 5 MeV with spin-orbit, d + 15N 10 MeV spin 1 with spin-orbit, 6Li + 12C 20 MeV):
//     every S_lj to 2e-5 (l <= 39 for spin 1: the coupling coefficients stop at j = 40), sigma/sigma_R at every angle
//     to 1e-4 relative, the reaction cross section to 1e-5; and
//     FRESCO's own reference against its refined run to the same tolerances (its convergence).
//  5. Transfer: 15N(d,n)16O(g.s.) and 12C(6Li,d)16O(g.s.), zero range (D0 given, recoil: FRESCO's INH = 2) and
//     finite range (post, no remnant; for 6Li the nuclear part of V_alpha-d), dsigma/dOmega from 0 to 90 deg where it
//     is above 1 % of its maximum, to 1 %; FRESCO's reference against its refined run to 0.2 % (its convergence).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/transfer/Transfer.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;

namespace {

int failures = 0;
std::string data_dir;

void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what.c_str());
    ++failures;
  }
}

// ---------------------------------------------------------------- the reference files
using Line = std::vector<std::string>;
/// The lines of each run ("reference", "refined", ...) of a case.
std::map<std::string, std::vector<Line>> ReadCase(const std::string &name) {
  std::ifstream in(data_dir + "/" + name + "/reference.txt");
  if (!in) {
    Check(false, "cannot read " + data_dir + "/" + name + "/reference.txt");
    return {};
  }
  std::map<std::string, std::vector<Line>> runs;
  std::string run, text;
  while (std::getline(in, text)) {
    if (text.empty() || text[0] == '#') continue;
    std::istringstream s(text);
    Line l;
    for (std::string t; s >> t;) l.push_back(t);
    if (l[0] == "run") {
      run = l[1];
      continue;
    }
    runs[run].push_back(l);
  }
  return runs;
}

std::vector<Line> Select(const std::vector<Line> &lines, const std::string &key) {
  std::vector<Line> out;
  for (const Line &l : lines)
    if (l[0] == key) out.push_back(l);
  return out;
}

double D(const std::string &s) { return std::stod(s); }

// ---------------------------------------------------------------- the inputs of the jobs
OpticalPotential Potential(RadiusConvention c, double At, double Ap, double rc) {
  OpticalPotential p;
  p.convention = c;
  p.target_mass_number = At;
  p.projectile_mass_number = Ap;
  p.coulomb_r0_fm = rc;
  return p;
}

OpticalPotential P15N() {  // p + 15N (p15N_elastic kp 1)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 15, 0, 1.25);
  p.real_volume = {50.0, 1.2, 0.65};
  p.imaginary_volume = {1.0, 1.3, 0.5};
  p.imaginary_surface = {6.0, 1.3, 0.5};
  p.real_spin_orbit = {6.0, 1.1, 0.6};
  return p;
}
OpticalPotential D15N(bool spin_orbit) {  // d + 15N (d15N_elastic kp 1; dn_* kp 1 without spin-orbit)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 15, 0, 1.3);
  p.real_volume = {90.0, 1.15, 0.75};
  p.imaginary_volume = {2.0, 1.35, 0.6};
  p.imaginary_surface = {10.0, 1.35, 0.6};
  if (spin_orbit) p.real_spin_orbit = {6.0, 1.07, 0.66};
  return p;
}
OpticalPotential N16O() {  // n + 16O (dn_* kp 2)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 16, 0, 1.3);
  p.real_volume = {48.0, 1.2, 0.65};
  p.imaginary_volume = {1.0, 1.3, 0.5};
  p.imaginary_surface = {6.0, 1.3, 0.5};
  return p;
}
OpticalPotential NP() {  // n-p binding (dn_* kp 3)
  OpticalPotential p = Potential(RadiusConvention::Absolute, 1, 0, 1.0);
  p.real_volume = {60.0, 1.5, 0.5};
  return p;
}
OpticalPotential P15NBound() {  // p-15N binding (dn_* kp 4)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 15, 0, 1.25);
  p.real_volume = {50.0, 1.25, 0.65};
  p.real_spin_orbit = {6.0, 1.25, 0.65};
  return p;
}
OpticalPotential Li6C12() {  // 6Li + 12C (li6* kp 1)
  OpticalPotential p = Potential(RadiusConvention::BothMasses, 12, 6, 0.73);
  p.real_volume = {120.0, 0.80, 0.85};
  p.imaginary_volume = {35.0, 0.95, 0.85};
  return p;
}
OpticalPotential D16O() {  // d + 16O (li6d_* kp 2)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 16, 0, 1.3);
  p.real_volume = {88.0, 1.17, 0.75};
  p.imaginary_volume = {1.0, 1.32, 0.65};
  p.imaginary_surface = {12.0, 1.32, 0.65};
  return p;
}
OpticalPotential AlphaD() {  // alpha-d binding (li6d_* kp 3)
  OpticalPotential p = Potential(RadiusConvention::Absolute, 1, 0, 1.9);
  p.real_volume = {70.0, 1.9, 0.65};
  return p;
}
OpticalPotential AlphaC12() {  // alpha-12C binding (li6d_* kp 4)
  OpticalPotential p = Potential(RadiusConvention::BothMasses, 12, 4, 0.9);
  p.real_volume = {100.0, 0.9, 0.65};
  return p;
}

ScatteringRequest Scattering(int z1, int z2, double m1, double m2, int two_s, double E_cm_in, double Q, int l_max) {
  ScatteringRequest q;
  q.pair.z1 = z1;
  q.pair.z2 = z2;
  q.pair.mu_u = m1 * m2 / (m1 + m2);
  q.two_s = two_s;
  q.energy_MeV = E_cm_in + Q;
  q.step_fm = 0.05;
  q.matching_radius_fm = 30.0;
  q.l_max = l_max;
  return q;
}

BoundStateRequest Bound(int z1, int z2, double m_core, double m_x, int l, int two_j, int two_s, double B, int nodes) {
  BoundStateRequest q;
  q.pair.z1 = z1;
  q.pair.z2 = z2;
  q.pair.mu_u = m_core * m_x / (m_core + m_x);
  q.l = l;
  q.two_j = two_j;
  q.two_s = two_s;
  q.binding_MeV = B;
  q.nodes = nodes;
  q.step_fm = 0.05;
  q.r_max_fm = 30.0;  // FRESCO's rmatch: the rms radius depends on it (the deuteron tail)
  return q;
}

// ---------------------------------------------------------------- 1, 2. constants, kinematics, potentials
void CheckConstants(const std::string &name, const std::vector<Line> &run) {
  for (const Line &l : Select(run, "constants")) {
    const double two_u = 2.0 * kAtomicMassUnit_MeV / (kHbarC_MeV_fm * kHbarC_MeV_fm), e2 = kFineStructure * kHbarC_MeV_fm;
    Check(std::fabs(D(l[1]) - two_u) < 6e-8 && std::fabs(D(l[2]) - e2) < 6e-6, name + ": FRESCO's 2u/hbar^2, e^2");
  }
}

void CheckKinematics(const std::string &name, const std::vector<Line> &run, const std::string &projectile,
                     const ScatteringRequest &q) {
  bool found = false;
  for (const Line &l : Select(run, "kin")) {
    if (l[1] != projectile) continue;
    found = true;
    const double k = rho(q.pair, 1.0, q.energy_MeV), et = eta(q.pair, q.energy_MeV);
    char buf[200];
    std::snprintf(buf, sizeof buf, "%s, %s: E %.5f k %.5f eta %.5f against FRESCO's %s %s %s", name.c_str(),
                  projectile.c_str(), q.energy_MeV, k, et, l[2].c_str(), l[3].c_str(), l[4].c_str());
    Check(std::fabs(D(l[2]) - q.energy_MeV) < 6e-6 && std::fabs(D(l[3]) - k) < 6e-6 && std::fabs(D(l[4]) - et) < 6e-6,
          buf);
  }
  Check(found, name + ": no kinematics line for " + projectile);
}

void CheckPotentials(const std::string &name, const std::vector<Line> &run, const std::map<int, OpticalPotential> &kp) {
  int compared = 0;
  double worst = 0.0;
  for (const Line &l : Select(run, "pot")) {
    const int k = std::stoi(l[1]), type = std::stoi(l[2]);
    const double r = D(l[3]);
    const std::complex<double> fresco(D(l[4]), D(l[5]));
    if (!kp.count(k) || type > 3) continue;
    const PotentialTerms t = TermsAt(kp.at(k), 1, r);
    const std::complex<double> ours = type == 0   ? std::complex<double>(t.coulomb_MeV)
                                      : type == 1 ? t.volume_MeV
                                      : type == 2 ? t.surface_MeV
                                                  : t.spin_orbit_form_MeV;
    worst = std::max(worst, std::abs(ours - fresco) / std::max(std::abs(fresco), 1e-6));
    ++compared;
  }
  char buf[160];
  std::snprintf(buf, sizeof buf, "%s: %d potential values against FRESCO's listing, worst %.2e", name.c_str(),
                compared, worst);
  Check(compared > 0 && worst < 1e-6, buf);
}

// ---------------------------------------------------------------- 3. bound states
void CheckBound(const std::string &name, const std::vector<Line> &run, int kn, const OpticalPotential &well,
                const BoundState &s, bool d0_listed) {
  for (const Line &l : Select(run, "bound")) {
    if (std::stoi(l[1]) != kn) continue;
    const double depth = D(l[2]) * well.real_volume.depth_MeV;
    double anc_plateau = 0.0;
    for (const Line &a : Select(run, "ancratio"))
      if (std::stoi(a[1]) == kn && D(a[2]) == 12.0) anc_plateau = D(a[3]);
    char buf[300];
    std::snprintf(buf, sizeof buf,
                  "%s, bound state %d: depth %.5f (FRESCO %.5f), rms %.4f (%s), D0 %.2f (%s), ANC %.6f (plateau %.6f)",
                  name.c_str(), kn, s.depth_MeV, depth, s.rms_radius_fm, l[3].c_str(),
                  s.l == 0 ? ZeroRangeConstant(s) : 0.0, l[4].c_str(), s.anc_fm_half, anc_plateau);
    bool ok = std::fabs(s.depth_MeV - depth) < 0.51e-4 * well.real_volume.depth_MeV &&
              std::fabs(s.rms_radius_fm - D(l[3])) < 0.51e-3 &&
              std::fabs(s.anc_fm_half - anc_plateau) < 2e-5 * std::fabs(anc_plateau);
    if (s.l == 0 && d0_listed) ok = ok && std::fabs(ZeroRangeConstant(s) - D(l[4])) < 0.051;
    Check(ok, buf);
    std::printf("%s\n", buf);
    return;
  }
  Check(false, name + ": no bound-state line " + std::to_string(kn));
}

// ---------------------------------------------------------------- 4. elastic
void CheckElastic(const std::string &name, const OpticalPotential &p, const ScatteringRequest &q,
                  const std::string &projectile) {
  const auto runs = ReadCase(name);
  if (runs.empty()) return;
  const std::vector<Line> &ref = runs.at("reference"), &fine = runs.at("refined");
  CheckConstants(name, ref);
  CheckKinematics(name, ref, projectile, q);
  CheckPotentials(name, ref, {{1, p}});
  const DistortedWaves w = SolveDistortedWaves(p, q);
  // S: ours against FRESCO's reference, and FRESCO's reference against its refined run
  double worst = 0.0, worst_fresco = 0.0;
  const std::vector<Line> s_ref = Select(ref, "S"), s_fine = Select(fine, "S");
  for (std::size_t i = 0; i < s_ref.size(); ++i) {
    const Line &l = s_ref[i];
    const std::complex<double> S(D(l[4]), D(l[5]));
    if (std::stoi(l[1]) > q.l_max) continue;  // FRESCO's S beyond our cut (|S - 1| is checked below)
    worst = std::max(worst, std::abs(WaveOf(w, std::stoi(l[1]), std::stoi(l[2])).S - S));
    if (i < s_fine.size())
      worst_fresco = std::max(worst_fresco, std::abs(std::complex<double>(D(s_fine[i][4]), D(s_fine[i][5])) - S));
  }
  char buf[200];
  std::snprintf(buf, sizeof buf, "%s: %zu S-matrix elements, worst |dS| %.2e (FRESCO reference - refined: %.2e)",
                name.c_str(), s_ref.size(), worst, worst_fresco);
  Check(!s_ref.empty() && worst < 2e-5 && worst_fresco < 2e-5, buf);
  std::printf("%s\n", buf);
  double beyond = 0.0;  // what our cut l_max leaves out
  for (const Line &l : s_ref)
    if (std::stoi(l[1]) > q.l_max) beyond = std::max(beyond, std::abs(std::complex<double>(D(l[4]), D(l[5])) - 1.0));
  std::snprintf(buf, sizeof buf, "%s: FRESCO's |S - 1| beyond l = %d is %.2e", name.c_str(), q.l_max, beyond);
  Check(beyond < 1e-6, buf);
  worst = worst_fresco = 0.0;
  double at = 0.0;
  const std::vector<Line> r_ref = Select(ref, "ratio"), r_fine = Select(fine, "ratio");
  for (std::size_t i = 0; i < r_ref.size(); ++i) {
    const double th = D(r_ref[i][1]), fresco = D(r_ref[i][2]);
    const double d = std::fabs(ElasticOverRutherford(w, std::min(th, 180.0)) / fresco - 1.0);
    if (d > worst) worst = d, at = th;
    if (i < r_fine.size()) worst_fresco = std::max(worst_fresco, std::fabs(D(r_fine[i][2]) / fresco - 1.0));
  }
  std::snprintf(buf, sizeof buf, "%s: sigma/sigma_R at %zu angles, worst %.2e at %.0f deg (FRESCO: %.2e)",
                name.c_str(), r_ref.size(), worst, at, worst_fresco);
  Check(!r_ref.empty() && worst < 1e-4 && worst_fresco < 1e-4, buf);
  std::printf("%s\n", buf);
  const double reaction = D(Select(ref, "reaction")[0][1]);
  std::snprintf(buf, sizeof buf, "%s: reaction cross section %.5f mb against %.5f", name.c_str(),
                ReactionCrossSection(w), reaction);
  Check(std::fabs(ReactionCrossSection(w) / reaction - 1.0) < 1e-5, buf);
}

// ---------------------------------------------------------------- 5. transfer
void CheckTransfer(const std::string &name, const std::vector<double> &ours, double tolerance) {
  const auto runs = ReadCase(name);
  if (runs.empty()) return;
  const std::vector<Line> ref = Select(runs.at("reference"), "xs"), fine = Select(runs.at("refined"), "xs");
  double peak = 0.0;
  for (const Line &l : ref)
    if (D(l[1]) <= 90.0) peak = std::max(peak, D(l[2]));
  double worst = 0.0, worst_fresco = 0.0, at = 0.0;
  int compared = 0;
  for (std::size_t i = 0; i < ref.size() && i < ours.size(); ++i) {
    const double th = D(ref[i][1]), fresco = D(ref[i][2]);
    if (th > 90.0 || fresco < 0.01 * peak) continue;
    ++compared;
    const double d = std::fabs(ours[i] / fresco - 1.0);
    if (d > worst) worst = d, at = th;
    worst_fresco = std::max(worst_fresco, std::fabs(D(fine[i][2]) / fresco - 1.0));
  }
  char buf[200];
  std::snprintf(buf, sizeof buf,
                "%s: %d angles in 0-90 deg above 1%% of the peak, worst %.2e at %.0f deg (FRESCO reference - "
                "refined: %.2e)",
                name.c_str(), compared, worst, at, worst_fresco);
  Check(compared > 5 && worst < tolerance && worst_fresco < 2e-3, buf);
  std::printf("%s\n", buf);
}

std::vector<double> Angles() {
  std::vector<double> th;
  for (int i = 0; i <= 36; ++i) th.push_back(i == 0 ? 1e-3 : 5.0 * i);
  return th;
}

void TestDN() {
  const auto runs = ReadCase("dn_finite_range");
  if (runs.empty()) return;
  const std::vector<Line> &ref = runs.at("reference");
  const ScatteringRequest in = Scattering(1, 7, 2.0, 15.0, 0, 10.0 * 15.0 / 17.0, 0.0, 30);
  const ScatteringRequest out = Scattering(0, 8, 1.0, 16.0, 0, 10.0 * 15.0 / 17.0, 9.9028, 30);
  CheckConstants("dn_finite_range", ref);
  CheckKinematics("dn_finite_range", ref, "d", in);
  CheckKinematics("dn_finite_range", ref, "n", out);
  CheckPotentials("dn_finite_range", ref, {{1, D15N(false)}, {2, N16O()}, {3, NP()}, {4, P15NBound()}});
  TransferReaction r;
  r.mass_a_u = 2.0, r.mass_A_u = 15.0, r.mass_s_u = 1.0, r.mass_F_u = 16.0, r.two_s_x = 1;
  r.projectile.state = SolveBoundState(NP(), Bound(0, 1, 1.0, 1.0, 0, 1, 1, 2.2246, 0));
  r.projectile.two_J_core = 1, r.projectile.two_J_composite = 2;
  r.residual.state = SolveBoundState(P15NBound(), Bound(1, 7, 15.0, 1.0, 1, 1, 1, 12.1274, 0));
  r.residual.two_J_core = 1, r.residual.two_J_composite = 0;
  CheckBound("dn_finite_range", ref, 1, NP(), r.projectile.state, true);
  CheckBound("dn_finite_range", ref, 2, P15NBound(), r.residual.state, false);
  const DistortedWaves wi = SolveDistortedWaves(D15N(false), in), wf = SolveDistortedWaves(N16O(), out);
  TransferGrid g;
  g.step_fm = 0.1;
  g.r_max_fm = 25.0;
  g.vertex_range_fm = 10.0;
  g.angle_nodes = 48;
  CheckTransfer("dn_finite_range", FiniteRangeCrossSection(r, wi, wf, g, Angles()), 0.01);
  CheckTransfer("dn_zero_range", ZeroRangeCrossSection(r, -130.31, wi, wf, Angles()), 0.01);
}

void TestLi6D() {
  const auto runs = ReadCase("li6d_finite_range");
  if (runs.empty()) return;
  const std::vector<Line> &ref = runs.at("reference");
  const ScatteringRequest in = Scattering(3, 6, 6.0, 12.0, 0, 20.0 * 12.0 / 18.0, 0.0, 38);
  const ScatteringRequest out = Scattering(1, 8, 2.0, 16.0, 0, 20.0 * 12.0 / 18.0, 5.6883, 38);
  CheckConstants("li6d_finite_range", ref);
  CheckKinematics("li6d_finite_range", ref, "6Li", in);
  CheckKinematics("li6d_finite_range", ref, "d", out);
  CheckPotentials("li6d_finite_range", ref, {{1, Li6C12()}, {2, D16O()}, {3, AlphaD()}, {4, AlphaC12()}});
  TransferReaction r;
  r.mass_a_u = 6.0, r.mass_A_u = 12.0, r.mass_s_u = 2.0, r.mass_F_u = 16.0, r.two_s_x = 0;
  r.projectile.state = SolveBoundState(AlphaD(), Bound(2, 1, 2.0, 4.0, 0, 0, 0, 1.4743, 1));
  r.projectile.two_J_core = 2, r.projectile.two_J_composite = 2;
  r.residual.state = SolveBoundState(AlphaC12(), Bound(2, 6, 12.0, 4.0, 0, 0, 0, 7.162, 2));
  r.residual.two_J_core = 0, r.residual.two_J_composite = 0;
  CheckBound("li6d_finite_range", ref, 1, AlphaD(), r.projectile.state, false);  // its D0 is that of KRPOT
  CheckBound("li6d_finite_range", ref, 2, AlphaC12(), r.residual.state, false);
  const auto zr = ReadCase("li6d_zero_range");
  if (!zr.empty()) CheckBound("li6d_zero_range", zr.at("reference"), 1, AlphaD(), r.projectile.state, true);
  const DistortedWaves wi = SolveDistortedWaves(Li6C12(), in), wf = SolveDistortedWaves(D16O(), out);
  TransferGrid g;
  g.step_fm = 0.1;
  g.r_max_fm = 25.0;
  g.vertex_range_fm = 12.0;
  g.angle_nodes = 48;
  r.interaction = PostInteraction::NuclearPart;
  CheckTransfer("li6d_finite_range", FiniteRangeCrossSection(r, wi, wf, g, Angles()), 0.01);
  CheckTransfer("li6d_zero_range", ZeroRangeCrossSection(r, 68.3056, wi, wf, Angles()), 0.01);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: %s <tests/data/fresco>\n", argv[0]);
    return 2;
  }
  data_dir = argv[1];
  CheckElastic("p15N_elastic", P15N(), Scattering(1, 7, 1.0, 15.0, 1, 5.0 * 15.0 / 16.0, 0.0, 22), "p");
  CheckElastic("d15N_elastic", D15N(true), Scattering(1, 7, 2.0, 15.0, 2, 10.0 * 15.0 / 17.0, 0.0, 32), "d");
  ScatteringRequest li6 = Scattering(3, 6, 6.0, 12.0, 2, 20.0 * 12.0 / 18.0, 0.0, 39);
  li6.step_fm = 0.025;  // k h = 0.04: sigma/sigma_R in the 165 deg minimum (1e-4 of Rutherford) to 1e-5
  CheckElastic("li6c12_elastic", Li6C12(), li6, "6Li");
  TestDN();
  TestLi6D();
  if (failures) {
    std::printf("transfer_fresco: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_fresco: all checks passed\n");
  return 0;
}
