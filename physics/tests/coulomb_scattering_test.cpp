// CoulombScattering.h against mpmath (tests/data/coulomb_scattering.csv, make_coulomb_scattering.py: the cross
// sections from the Coulomb amplitude, Mott as the channel-spin average of the symmetrised amplitude, d/dE by mpmath's
// derivative), and against closed forms:
//  - non-relativistic Rutherford = (Z1 Z2 alpha hbar c / 4E)^2 / sin^4(theta/2), d/dE = -2 sigma_R / E;
//  - Mott at 90 deg: 4 Rutherford(90) for spin 0, Rutherford(90) for spin 1/2 (the arithmetic below);
//  - Mott symmetric about 90 deg; Mott -> Rutherford(theta) + Rutherford(pi - theta) without the interference term
//    (two_s large: w -> 0);
//  - invalid arguments give NaN.
//   physics_coulomb_scattering_test tests/data/coulomb_scattering.csv
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include "azure2/physics/Constants.h"
#include "azure2/physics/CoulombScattering.h"

using namespace azure2::physics;

namespace {
int failures = 0;
void check(const std::string &what, bool ok, double got = 0.0, double want = 0.0) {
  std::printf("  %s  %s", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) std::printf("  (got %.17g, want %.17g)", got, want);
  std::printf("\n");
  if (!ok) failures++;
}
const double kPi = 3.14159265358979323846;
bool close(double a, double b, double rel) { return std::fabs(a - b) <= rel * std::fabs(b); }

ChannelKinematics Pair(int z1, int z2, double m1u, double m2u, bool relativistic) {
  ChannelKinematics k;
  k.z1 = z1;
  k.z2 = z2;
  k.mu_u = m1u * m2u / (m1u + m2u);
  k.m1_MeV = m1u * kAtomicMassUnit_MeV;
  k.m2_MeV = m2u * kAtomicMassUnit_MeV;
  k.relativistic = relativistic;
  return k;
}
}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: physics_coulomb_scattering_test coulomb_scattering.csv\n");
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::fprintf(stderr, "cannot read %s\n", argv[1]);
    return 2;
  }
  // mpmath table: value to 1e-11, d/dE to 1e-9 (relative); the largest Mott phase eta ln tan^2(theta/2) of the table
  // is about 130 rad (12C + 12C at 0.3 MeV, 10 deg), whose cosine double precision gives to about 1e-13 absolute.
  std::string line;
  int rows = 0, badValue = 0, badDerivative = 0;
  double worstValue = 0.0, worstDerivative = 0.0;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream row(line);
    int z1, z2, twoS, rel;
    double m1, m2, e, theta, value, dE;
    if (!(row >> z1 >> z2 >> m1 >> m2 >> twoS >> rel >> e >> theta >> value >> dE)) {
      check("parse " + line, false);
      continue;
    }
    const ChannelKinematics kin = Pair(z1, z2, m1, m2, rel != 0);
    const CoulombCrossSection got = twoS < 0 ? rutherford_cross_section(kin, e, theta) : mott_cross_section(kin, twoS, e, theta);
    const double rv = std::fabs(got.value - value) / std::fabs(value);
    const double rd = std::fabs(got.d_dE - dE) / std::fabs(dE);
    worstValue = std::max(worstValue, rv);
    worstDerivative = std::max(worstDerivative, rd);
    if (!(rv <= 1e-11)) {
      badValue++;
      std::printf("    value  %s: got %.17g\n", line.c_str(), got.value);
    }
    if (!(rd <= 1e-9)) {
      badDerivative++;
      std::printf("    d_dE   %s: got %.17g\n", line.c_str(), got.d_dE);
    }
    rows++;
  }
  std::printf("  %d rows; largest relative difference: value %.2e, d/dE %.2e\n", rows, worstValue, worstDerivative);
  check("mpmath table read (120 rows)", rows == 120, rows, 120);
  check("value against mpmath to 1e-11", badValue == 0, badValue, 0);
  check("d/dE against mpmath to 1e-9", badDerivative == 0, badDerivative, 0);

  // Closed forms, non-relativistic.
  const double e2 = kFineStructure * kHbarC_MeV_fm;  // MeV fm
  const ChannelKinematics aC = Pair(2, 6, 4.00260325, 12.0, false);
  for (double theta : {15.0, 60.0, 120.0, 175.0}) {
    const double E = 3.2;
    const double s = std::sin(theta * kPi / 360.0);
    const double want = std::pow(12.0 * e2 / (4.0 * E), 2.0) / std::pow(s, 4.0) * 0.01;
    const CoulombCrossSection r = rutherford_cross_section(aC, E, theta);
    check("Rutherford alpha + 12C at " + std::to_string(theta) + " deg = (Z1 Z2 e^2/4E)^2/sin^4", close(r.value, want, 1e-13),
          r.value, want);
    check("  d/dE = -2 sigma_R / E", close(r.d_dE, -2.0 * want / E, 1e-12), r.d_dE, -2.0 * want / E);
  }
  // Mott of spin-0 bosons at 90 deg: sin = cos = 1/sqrt 2, tan = 1, so the bracket is 4 + 4 + 2 x 1 x 1 x 4 = 16 and
  // sigma_M(90) = 16 eta^2/(4 k^2) = 4 sigma_R(90) (sigma_R(90) = eta^2/(4 k^2) x 4).
  const ChannelKinematics aa = Pair(2, 2, 4.00260325, 4.00260325, false);
  for (double E : {0.1, 1.0, 10.0}) {
    const double m = mott_cross_section(aa, 0, E, 90.0).value, r = rutherford_cross_section(aa, E, 90.0).value;
    check("Mott (spin 0) at 90 deg = 4 Rutherford(90) at E = " + std::to_string(E), close(m, 4.0 * r, 1e-14), m, 4.0 * r);
    // fermions of spin 1/2: w = -1/2, bracket 4 + 4 - 4 = 4: sigma_M(90) = Rutherford(90)
    const double mf = mott_cross_section(aa, 1, E, 90.0).value;
    check("Mott (spin 1/2) at 90 deg = Rutherford(90)", close(mf, r, 1e-14), mf, r);
  }
  // Symmetry about 90 deg and the incoherent limit.
  for (int twoS : {0, 1, 2, 3}) {
    const double a = mott_cross_section(aa, twoS, 2.0, 35.0).value, b = mott_cross_section(aa, twoS, 2.0, 145.0).value;
    check("Mott 2s = " + std::to_string(twoS) + " symmetric about 90 deg", close(a, b, 1e-12), a, b);
  }
  {
    const double big = mott_cross_section(aa, 2000001, 2.0, 35.0).value;
    const double sum = rutherford_cross_section(aa, 2.0, 35.0).value + rutherford_cross_section(aa, 2.0, 145.0).value;
    check("Mott -> Rutherford(theta) + Rutherford(pi - theta) as w -> 0", close(big, sum, 1e-6), big, sum);
  }
  // Invalid arguments.
  check("Rutherford at theta = 0: NaN", std::isnan(rutherford_cross_section(aC, 1.0, 0.0).value));
  check("Rutherford at E = 0: NaN", std::isnan(rutherford_cross_section(aC, 0.0, 30.0).value));
  check("Mott at theta = 180: NaN", std::isnan(mott_cross_section(aa, 0, 1.0, 180.0).value));
  check("Mott with 2s < 0: NaN", std::isnan(mott_cross_section(aa, -1, 1.0, 60.0).value));

  std::printf("%s\n", failures ? "FAILED" : "PASSED");
  return failures ? 1 : 0;
}
