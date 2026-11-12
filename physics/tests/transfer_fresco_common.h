// Shared by the FRESCO tests of physics/transfer (transfer_fresco_test, transfer_continuum_fresco_test): the reader
// of tests/data/fresco/<case>/reference.txt (make_fresco_benchmarks.py) and the inputs of the jobs there.
#ifndef AZURE2_PHYSICS_TESTS_TRANSFER_FRESCO_COMMON_H
#define AZURE2_PHYSICS_TESTS_TRANSFER_FRESCO_COMMON_H

#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/transfer/BoundState.h"
#include "azure2/physics/transfer/DistortedWaves.h"

namespace fresco_reference {

using namespace azure2::physics;
using namespace azure2::physics::transfer;

inline int failures = 0;
inline std::string data_dir;

inline void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what.c_str());
    ++failures;
  }
}

// ---------------------------------------------------------------- the reference files
using Line = std::vector<std::string>;
/// The lines of each run ("reference", "refined", ...) of a case.
inline std::map<std::string, std::vector<Line>> ReadCase(const std::string &name) {
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

inline std::vector<Line> Select(const std::vector<Line> &lines, const std::string &key) {
  std::vector<Line> out;
  for (const Line &l : lines)
    if (l[0] == key) out.push_back(l);
  return out;
}

inline double D(const std::string &s) { return std::stod(s); }

// ---------------------------------------------------------------- the inputs of the jobs
inline OpticalPotential Potential(RadiusConvention c, double At, double Ap, double rc) {
  OpticalPotential p;
  p.convention = c;
  p.target_mass_number = At;
  p.projectile_mass_number = Ap;
  p.coulomb_r0_fm = rc;
  return p;
}

inline OpticalPotential P15N() {  // p + 15N (p15N_elastic kp 1)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 15, 0, 1.25);
  p.real_volume = {50.0, 1.2, 0.65};
  p.imaginary_volume = {1.0, 1.3, 0.5};
  p.imaginary_surface = {6.0, 1.3, 0.5};
  p.real_spin_orbit = {6.0, 1.1, 0.6};
  return p;
}
inline OpticalPotential D15N(bool spin_orbit) {  // d + 15N (d15N_elastic kp 1; dn_* kp 1 without spin-orbit)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 15, 0, 1.3);
  p.real_volume = {90.0, 1.15, 0.75};
  p.imaginary_volume = {2.0, 1.35, 0.6};
  p.imaginary_surface = {10.0, 1.35, 0.6};
  if (spin_orbit) p.real_spin_orbit = {6.0, 1.07, 0.66};
  return p;
}
inline OpticalPotential N16O() {  // n + 16O (dn_* kp 2)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 16, 0, 1.3);
  p.real_volume = {48.0, 1.2, 0.65};
  p.imaginary_volume = {1.0, 1.3, 0.5};
  p.imaginary_surface = {6.0, 1.3, 0.5};
  return p;
}
inline OpticalPotential N15Core() {  // n + 15N, the core-core potential of the remnant (dn_*remnant kp 5)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 15, 0, 1.3);
  p.real_volume = {48.0, 1.2, 0.65};
  p.imaginary_volume = {1.0, 1.3, 0.5};
  p.imaginary_surface = {6.0, 1.3, 0.5};
  return p;
}
inline OpticalPotential NP() {  // n-p binding (dn_* kp 3)
  OpticalPotential p = Potential(RadiusConvention::Absolute, 1, 0, 1.0);
  p.real_volume = {60.0, 1.5, 0.5};
  return p;
}
inline OpticalPotential P15NBound() {  // p-15N binding (dn_* kp 4)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 15, 0, 1.25);
  p.real_volume = {50.0, 1.25, 0.65};
  p.real_spin_orbit = {6.0, 1.25, 0.65};
  return p;
}
inline OpticalPotential Li6C12() {  // 6Li + 12C (li6* kp 1)
  OpticalPotential p = Potential(RadiusConvention::BothMasses, 12, 6, 0.73);
  p.real_volume = {120.0, 0.80, 0.85};
  p.imaginary_volume = {35.0, 0.95, 0.85};
  return p;
}
inline OpticalPotential D16O() {  // d + 16O (li6d_* kp 2)
  OpticalPotential p = Potential(RadiusConvention::TargetMass, 16, 0, 1.3);
  p.real_volume = {88.0, 1.17, 0.75};
  p.imaginary_volume = {1.0, 1.32, 0.65};
  p.imaginary_surface = {12.0, 1.32, 0.65};
  return p;
}
inline OpticalPotential AlphaD() {  // alpha-d binding (li6d_* kp 3)
  OpticalPotential p = Potential(RadiusConvention::Absolute, 1, 0, 1.9);
  p.real_volume = {70.0, 1.9, 0.65};
  return p;
}
inline OpticalPotential AlphaC12() {  // alpha-12C binding (li6d_* kp 4)
  OpticalPotential p = Potential(RadiusConvention::BothMasses, 12, 4, 0.9);
  p.real_volume = {100.0, 0.9, 0.65};
  return p;
}

inline ScatteringRequest Scattering(int z1, int z2, double m1, double m2, int two_s, double E_cm_in, double Q, int l_max) {
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

inline BoundStateRequest Bound(int z1, int z2, double m_core, double m_x, int l, int two_j, int two_s, double B, int nodes) {
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

inline void CheckKinematics(const std::string &name, const std::vector<Line> &run, const std::string &projectile,
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

}  // namespace fresco_reference

#endif
