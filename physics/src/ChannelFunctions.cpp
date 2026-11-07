// Channel functions of an open particle channel and the channel kinematics (ChannelFunctions.h,
// docs/upgrade/coulomb.md section 2).  Every expression keeps the engine's order of operations (Kinematics.h,
// CoulFunc.cpp, CoulFunc.h as they were before they called this file), so that the engine's numbers do not move.
#include "azure2/physics/ChannelFunctions.h"

#include <cmath>
#include <limits>

#include "azure2/physics/Constants.h"

namespace azure2 {
namespace physics {

namespace {
constexpr double uconv = kAtomicMassUnit_MeV;
constexpr double hbarc = kHbarC_MeV_fm;
constexpr double fstruc = kFineStructure;
}  // namespace

double relativistic_k2(double m1, double m2, double T) {
  const double m = m1 + m2;
  return T * (T + 2.0 * m) * (T + 2.0 * m1) * (T + 2.0 * m2) / (4.0 * (m + T) * (m + T)) / (hbarc * hbarc);
}

double relativistic_k(double m1, double m2, double T) { return std::sqrt(std::fabs(relativistic_k2(m1, m2, T))); }

double relativistic_eta(double z1z2, double m1, double m2, double T) {
  const double m = m1 + m2;
  const double p = relativistic_k(m1, m2, T) * hbarc;  // MeV
  return z1z2 * fstruc * (2.0 * m1 * m2 + 2.0 * m * T + T * T) / (2.0 * p * (m + T));
}

double eta(const ChannelKinematics &kin, double energy) {
  if (!kin.relativistic) return std::sqrt(uconv / 2.) * fstruc * kin.z1 * kin.z2 * std::sqrt(kin.mu_u / energy);
  return relativistic_eta((double)kin.z1 * kin.z2, kin.m1_MeV, kin.m2_MeV, energy);
}

double rho(const ChannelKinematics &kin, double radius, double energy) {
  if (!kin.relativistic) return std::sqrt(2. * uconv) / hbarc * radius * std::sqrt(kin.mu_u * energy);
  return relativistic_k(kin.m1_MeV, kin.m2_MeV, energy) * radius;
}

double eta_bound(const ChannelKinematics &kin, double binding) {
  if (!kin.relativistic) return std::sqrt(uconv / 2.) * fstruc * kin.z1 * kin.z2 * std::sqrt(kin.mu_u / binding);
  return relativistic_eta((double)kin.z1 * kin.z2, kin.m1_MeV, kin.m2_MeV, -binding);
}

double rho_bound(const ChannelKinematics &kin, double radius, double binding) {
  if (!kin.relativistic) return std::sqrt(2. * uconv) / hbarc * radius * std::sqrt(kin.mu_u * binding);
  return relativistic_k(kin.m1_MeV, kin.m2_MeV, -binding) * radius;
}

double penetrability(const coulomb::CoulombWaves &c, double rho) {
  return rho / (std::pow(c.F, 2.0) + std::pow(c.G, 2.0));
}

double shift(const coulomb::CoulombWaves &c, double rho) {
  double s;
  if (std::pow(c.F, 2.0) == 0. && c.F * c.Fp == 0.)
    s = rho * (c.Gp / c.G);
  else
    s = rho / (std::pow(c.F, 2.0) + std::pow(c.G, 2.0)) * (c.F * c.Fp + c.G * c.Gp);
  // Just above a charged-particle threshold (eta >~ 110) G^2 or G G' overflows and, with F F' a denormal rather
  // than 0, the formula gives NaN or -inf.  Scaling by hypot(F, G) avoids the overflow; it is used only when the
  // plain formula failed, so every value that formula gets right is unchanged.
  if (!std::isfinite(s)) {
    const double h = std::hypot(c.F, c.G);
    s = rho * ((c.F / h) * (c.Fp / h) + (c.G / h) * (c.Gp / h));
  }
  return s;
}

std::complex<double> exp_minus_i_phi(const coulomb::CoulombWaves &c) {
  // The plain formula returns 0 instead of a unit phase where G^2 overflows (eta >~ 110); the scaled form only then.
  const double s = std::sqrt(std::pow(c.F, 2.0) + std::pow(c.G, 2.0));
  const double re = c.G / s, im = -c.F / s;
  if (std::isfinite(re) && std::isfinite(im) && std::isfinite(s) && s > 0.0) return std::complex<double>(re, im);
  const double h = std::hypot(c.F, c.G);
  return std::complex<double>(c.G / h, -c.F / h);
}

double shift_derivative(const std::function<double(double)> &S_of_E, double energy) {
  const double h = 1e-3 * energy;
  auto cd = [&](double hh) { return (S_of_E(energy + hh) - S_of_E(energy - hh)) / (2.0 * hh); };
  return (4.0 * cd(0.5 * h) - cd(h)) / 3.0;
}

OpenChannelValues channel_functions(const ChannelKinematics &kin, int l, double a_fm, double E_MeV,
                                    coulomb::Method method) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  OpenChannelValues v;
  v.P = v.S = v.dS_dE = v.sigma = nan;
  v.exp_minus_i_phi = std::complex<double>(nan, nan);
  if (!(E_MeV > 0.0) || !(a_fm > 0.0)) {
    v.status = coulomb::Status::invalid_argument;
    return v;
  }
  bool ok = true;
  auto waves = [&](double e) {
    const coulomb::CoulombWaves w = coulomb::coulomb(l, eta(kin, e), rho(kin, a_fm, e), method);
    if (w.status != coulomb::Status::ok) ok = false;
    return w;
  };
  const coulomb::CoulombWaves w = waves(E_MeV);
  if (!ok) {
    v.status = w.status;
    return v;
  }
  const double r = rho(kin, a_fm, E_MeV);
  const double dS = shift_derivative([&](double e) { return shift(waves(e), rho(kin, a_fm, e)); }, E_MeV);
  if (!ok) return v;
  v.P = penetrability(w, r);
  v.S = shift(w, r);
  v.dS_dE = dS;
  v.exp_minus_i_phi = exp_minus_i_phi(w);
  v.sigma = w.sigma;
  v.status = coulomb::Status::ok;
  return v;
}

}  // namespace physics
}  // namespace azure2
