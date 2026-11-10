// The Rutherford and Mott cross sections (CoulombScattering.h; docs/upgrade/data_kinds.md section 4.3), with eta and k
// from ChannelFunctions.h.
#include "azure2/physics/CoulombScattering.h"

#include <cmath>
#include <limits>

namespace azure2 {
namespace physics {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kFm2ToBarn = 0.01;
const double kNaN = std::numeric_limits<double>::quiet_NaN();

/// eta and k (fm^-1) at E, and their logarithmic derivatives d ln eta / dE, d ln k / dE (MeV^-1).
struct EtaK {
  double eta, k, dlnEta, dlnK;
};

EtaK EtaAndK(const ChannelKinematics &kin, double E) {
  EtaK r;
  r.eta = eta(kin, E);
  r.k = rho(kin, 1.0, E);  // rho = k a at a = 1 fm
  if (!kin.relativistic) {
    // eta ~ E^-1/2, k ~ E^1/2
    r.dlnEta = -0.5 / E;
    r.dlnK = 0.5 / E;
    return r;
  }
  // k^2 = T (T + 2m)(T + 2 m1)(T + 2 m2) / (4 (m + T)^2) / (hbar c)^2,
  // eta = Z1 Z2 alpha (2 m1 m2 + 2 m T + T^2) / (2 hbar c k (m + T))   (relativistic_k2, relativistic_eta)
  const double m1 = kin.m1_MeV, m2 = kin.m2_MeV, m = m1 + m2;
  r.dlnK = 0.5 * (1.0 / E + 1.0 / (E + 2.0 * m) + 1.0 / (E + 2.0 * m1) + 1.0 / (E + 2.0 * m2) - 2.0 / (m + E));
  const double a = 2.0 * m1 * m2 + 2.0 * m * E + E * E;
  r.dlnEta = 2.0 * (m + E) / a - r.dlnK - 1.0 / (m + E);
  return r;
}

bool Valid(double E, double theta_deg) { return E > 0.0 && theta_deg > 0.0 && theta_deg < 180.0; }

}  // namespace

CoulombCrossSection rutherford_cross_section(const ChannelKinematics &kin, double E_MeV, double theta_deg) {
  if (!Valid(E_MeV, theta_deg)) return {kNaN, kNaN};
  const EtaK ek = EtaAndK(kin, E_MeV);
  const double s = std::sin(theta_deg * kPi / 360.0);
  const double s2 = s * s;
  CoulombCrossSection out;
  out.value = ek.eta * ek.eta / (4.0 * ek.k * ek.k * s2 * s2) * kFm2ToBarn;
  out.d_dE = out.value * 2.0 * (ek.dlnEta - ek.dlnK);
  return out;
}

CoulombCrossSection mott_cross_section(const ChannelKinematics &kin, int two_s, double E_MeV, double theta_deg) {
  if (!Valid(E_MeV, theta_deg) || two_s < 0) return {kNaN, kNaN};
  const EtaK ek = EtaAndK(kin, E_MeV);
  const double half = theta_deg * kPi / 360.0;
  const double s = std::sin(half), c = std::cos(half);
  const double s2 = s * s, c2 = c * c;
  const double w = ((two_s % 2 == 0) ? 1.0 : -1.0) / (two_s + 1.0);  // (-1)^(2s) / (2s + 1)
  const double logTan2 = std::log(s2 / c2);
  const double phase = ek.eta * logTan2;
  const double bracket = 1.0 / (s2 * s2) + 1.0 / (c2 * c2) + 2.0 * w * std::cos(phase) / (s2 * c2);
  const double scale = ek.eta * ek.eta / (4.0 * ek.k * ek.k) * kFm2ToBarn;
  CoulombCrossSection out;
  out.value = scale * bracket;
  // d/dE: the prefactor eta^2/k^2 and the phase eta ln tan^2(theta/2) of the interference term
  const double dPhase = ek.eta * ek.dlnEta * logTan2;
  out.d_dE = out.value * 2.0 * (ek.dlnEta - ek.dlnK) - scale * 2.0 * w * std::sin(phase) * dPhase / (s2 * c2);
  return out;
}

}  // namespace physics
}  // namespace azure2
