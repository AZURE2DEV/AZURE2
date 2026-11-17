// Argument checks, dispatch to one adaptor per algorithm (methods.h), the Coulomb phase, method names
// (docs/upgrade/coulomb.md section 1.2).  From pycoulomb's libcoulomb/src/coulomb.cpp (README.md).
#include <cmath>
#include <complex>
#include <algorithm>
#include <limits>
#include <vector>

#include "azure2/physics/coulomb/Coulomb.h"
#include "methods.h"

namespace azure2 {
namespace physics {
namespace coulomb {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const int kLMax = 10000;

bool ValidL(int l) { return l >= 0 && l <= kLMax; }

CoulombWaves Failed(Status s) { return CoulombWaves{kNaN, kNaN, kNaN, kNaN, kNaN, s}; }

using Adaptor = Status (*)(int, double, double, CoulombWaves &);

Adaptor AdaptorOf(Method m) {
  switch (m) {
    case Method::coulsm: return detail::coulsm;
    case Method::coul: return detail::coul;
    case Method::coulfg: return detail::coulfg;
    case Method::coulcc: return detail::coulcc;
  }
  return nullptr;
}

}  // namespace

CoulombWaves coulomb(int l, double eta, double rho, Method method) noexcept {
  if (!ValidL(l) || !std::isfinite(eta) || !std::isfinite(rho) || !(rho > 0.0)) return Failed(Status::invalid_argument);
  const Adaptor f = AdaptorOf(method);
  if (f == nullptr) return Failed(Status::invalid_argument);
  CoulombWaves w = Failed(Status::failed);
  Status s = Status::failed;
  try {
    s = f(l, eta, rho, w);
  } catch (...) {
    s = Status::failed;
  }
  if (s != Status::ok || !std::isfinite(w.F) || !std::isfinite(w.Fp) || !std::isfinite(w.G) || !std::isfinite(w.Gp))
    return Failed(Status::failed);
  w.sigma = coulomb_phase(l, eta);
  w.status = Status::ok;
  return w;
}

// sigma_l = Im ln Gamma(l + 1 + i eta).  The argument is shifted to x = l + 1 + n >= 12 by
// ln Gamma(w) = ln Gamma(w + n) - sum_k ln(w + k), whose imaginary parts are atan2(eta, l + 1 + k) (continuous,
// Re > 0), and Im ln Gamma(x + i eta) is Stirling's series with 8 terms (truncation < 1e-17 for |w| >= 12).
double coulomb_phase(int l, double eta) noexcept {
  if (!ValidL(l) || !std::isfinite(eta)) return kNaN;
  if (eta == 0.0) return 0.0;
  const int n = (l + 1 < 12) ? 12 - (l + 1) : 0;
  double shift = 0.0;
  for (int k = 0; k < n; ++k) shift += std::atan2(eta, l + 1.0 + k);
  const double x = l + 1.0 + n;
  const std::complex<double> w(x, eta);
  // (w - 1/2) ln w - w + ln(2 pi)/2: imaginary part
  double im = (x - 0.5) * std::atan2(eta, x) + eta * std::log(std::abs(w)) - eta;
  constexpr double c[8] = {1.0 / 12.0, -1.0 / 360.0, 1.0 / 1260.0, -1.0 / 1680.0,
                           1.0 / 1188.0, -691.0 / 360360.0, 1.0 / 156.0, -3617.0 / 122400.0};
  const std::complex<double> inv = 1.0 / w, inv2 = inv * inv;
  std::complex<double> p = inv, series = 0.0;
  for (int i = 0; i < 8; ++i) {
    series += c[i] * p;
    p *= inv2;
  }
  im += series.imag();
  return im - shift;
}

// Moved unchanged from physics/src/transfer/DistortedWaves.cpp (its RiccatiBessel, the transfer code's plane waves) on
// 2026-11-16, so that coulsm's F at eta = 0 can use it (docs/upgrade/coulomb.md section 10).
void riccati_bessel(int l_max, double x, double *u) {
  if (!(x > 0.0)) {
    for (int l = 0; l <= l_max; ++l) u[l] = 0.0;
    return;
  }
  if (l_max < x) {
    double jm = std::sin(x) / x, j = std::sin(x) / (x * x) - std::cos(x) / x;
    u[0] = x * jm;
    for (int l = 1; l <= l_max; ++l) {
      u[l] = x * j;
      const double next = (2 * l + 1) / x * j - jm;
      jm = j;
      j = next;
    }
    return;
  }
  const double top = std::max<double>(l_max, x);
  const int start = static_cast<int>(top + 30.0 + std::sqrt(40.0 * top));
  std::vector<double> j(start + 2, 0.0);
  j[start] = 1.0;
  for (int l = start; l >= 1; --l) {
    j[l - 1] = (2 * l + 1) / x * j[l] - j[l + 1];
    if (std::fabs(j[l - 1]) > 1e100)  // kept below 1e100, so that the sum of squares cannot overflow
      for (int m = l - 1; m <= start; ++m) j[m] *= 1e-100;
  }
  double sum = 0.0;
  for (int l = start; l >= 0; --l) sum += (2 * l + 1) * j[l] * j[l];
  // The sign from the closed form of j_0 or, where it is small, j_1.
  const double j0 = std::sin(x) / x, j1 = std::sin(x) / (x * x) - std::cos(x) / x;
  const double sign = std::fabs(j0) > std::fabs(j1) ? (j0 * j[0] < 0.0 ? -1.0 : 1.0) : (j1 * j[1] < 0.0 ? -1.0 : 1.0);
  const double scale = sign / std::sqrt(sum);
  for (int l = 0; l <= l_max; ++l) u[l] = x * j[l] * scale;
}

const char *method_name(Method method) noexcept {
  switch (method) {
    case Method::coulsm: return "coulsm";
    case Method::coul: return "coul";
    case Method::coulfg: return "coulfg";
    case Method::coulcc: return "coulcc";
  }
  return "unknown";
}

bool method_from_name(const std::string &name, Method &method) noexcept {
  for (Method m : {Method::coulsm, Method::coul, Method::coulfg, Method::coulcc}) {
    if (name == method_name(m)) {
      method = m;
      return true;
    }
  }
  return false;
}

}  // namespace coulomb
}  // namespace physics
}  // namespace azure2
