#include "azure2/physics/transfer/RadialNumerics.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

constexpr double kRescale = 1e200;

template <class T>
T Lagrange6(const std::vector<T> &f, double h, double r) {
  const int n = static_cast<int>(f.size());
  if (n == 0 || r > (n - 1) * h) return T(0.0);
  if (n < 6) throw std::invalid_argument("InterpolateOnGrid: fewer than 6 nodes");
  const double s = r / h;
  int first = static_cast<int>(std::floor(s)) - 2;  // nodes first ... first + 5 around r
  first = std::max(0, std::min(first, n - 6));
  T value(0.0);
  for (int i = 0; i < 6; ++i) {
    double w = 1.0;
    for (int j = 0; j < 6; ++j)
      if (j != i) w *= (s - (first + j)) / static_cast<double>(i - j);
    value += w * f[first + i];
  }
  return value;
}

template <class T>
T SimpsonRule(const std::vector<T> &f, double h) {
  const int n = static_cast<int>(f.size());
  if (n < 4) throw std::invalid_argument("Simpson: fewer than 4 nodes");
  const int intervals = n - 1;
  const int simpson_end = (intervals % 2 == 0) ? n - 1 : n - 4;  // last node of the Simpson part
  T sum = f[0] + f[simpson_end];
  for (int i = 1; i < simpson_end; ++i) sum += (i % 2 == 1 ? 4.0 : 2.0) * f[i];
  sum *= h / 3.0;
  if (simpson_end != n - 1)  // Simpson's 3/8 rule on the last three intervals
    sum += 3.0 * h / 8.0 * (f[n - 4] + 3.0 * f[n - 3] + 3.0 * f[n - 2] + f[n - 1]);
  return sum;
}

}  // namespace

std::vector<std::complex<double>> NumerovOutward(const std::vector<std::complex<double>> &q, double h, int l,
                                                 std::complex<double> q0, double q_m1) {
  const int n = static_cast<int>(q.size());
  if (n < 3 || !(h > 0.0) || l < 0) throw std::invalid_argument("NumerovOutward: needs 3 nodes, h > 0, l >= 0");
  std::vector<std::complex<double>> u(n, 0.0);
  const double a1 = -q_m1 / (2.0 * l + 2.0);
  const std::complex<double> a2 = -(q_m1 * a1 + q0) / (2.0 * (2 * l + 3));
  const std::complex<double> a3 = -(q_m1 * a2 + q0 * a1) / (3.0 * (2 * l + 4));
  // u / h^{l+1}: (r/h)^{l+1} times the series, which does not underflow for large l.
  for (int i = 1; i <= 2; ++i) {
    const double r = i * h;
    u[i] = std::pow(static_cast<double>(i), l + 1) * (1.0 + r * (a1 + r * (a2 + r * a3)));
  }
  const double h2 = h * h / 12.0;
  for (int i = 2; i < n - 1; ++i) {
    u[i + 1] = (2.0 * (1.0 - 5.0 * h2 * q[i]) * u[i] - (1.0 + h2 * q[i - 1]) * u[i - 1]) / (1.0 + h2 * q[i + 1]);
    if (std::abs(u[i + 1]) > kRescale)
      for (int j = 0; j <= i + 1; ++j) u[j] /= kRescale;
  }
  return u;
}

std::vector<std::complex<double>> NumerovInward(const std::vector<std::complex<double>> &q, double h,
                                                std::complex<double> u_before, std::complex<double> u_last,
                                                int lowest) {
  const int n = static_cast<int>(q.size());
  if (n < 3 || lowest < 0 || lowest > n - 2) throw std::invalid_argument("NumerovInward: bad node range");
  std::vector<std::complex<double>> u(n, 0.0);
  u[n - 1] = u_last;
  u[n - 2] = u_before;
  const double h2 = h * h / 12.0;
  for (int i = n - 2; i > lowest; --i) {
    u[i - 1] = (2.0 * (1.0 - 5.0 * h2 * q[i]) * u[i] - (1.0 + h2 * q[i + 1]) * u[i + 1]) / (1.0 + h2 * q[i - 1]);
    if (std::abs(u[i - 1]) > kRescale)
      for (int j = i - 1; j < n; ++j) u[j] /= kRescale;
  }
  return u;
}

coulomb::CoulombWaves CoulombForMatching(int l, double eta, double rho) {
  const double turning = eta + std::sqrt(eta * eta + l * (l + 1.0));
  const coulomb::CoulombWaves c =
      coulomb::coulomb(l, eta, rho, rho >= turning ? coulomb::Method::coulfg : coulomb::Method::coulcc);
  if (c.status != coulomb::Status::ok)
    throw std::runtime_error("Coulomb functions failed at l = " + std::to_string(l) + ", eta = " + std::to_string(eta) +
                             ", rho = " + std::to_string(rho));
  return c;
}

double InterpolateOnGrid(const std::vector<double> &f, double h, double r) { return Lagrange6(f, h, r); }
std::complex<double> InterpolateOnGrid(const std::vector<std::complex<double>> &f, double h, double r) {
  return Lagrange6(f, h, r);
}

double Simpson(const std::vector<double> &f, double h) { return SimpsonRule(f, h); }
std::complex<double> Simpson(const std::vector<std::complex<double>> &f, double h) { return SimpsonRule(f, h); }

void GaussLegendre(int n, std::vector<double> &x, std::vector<double> &w) {
  if (n < 1) throw std::invalid_argument("GaussLegendre: n >= 1");
  x.assign(n, 0.0);
  w.assign(n, 0.0);
  const double pi = std::acos(-1.0);
  for (int i = 0; i < (n + 1) / 2; ++i) {
    double z = std::cos(pi * (i + 0.75) / (n + 0.5));
    double dp = 0.0;
    for (int iter = 0; iter < 100; ++iter) {
      double p0 = 1.0, p1 = z;
      for (int k = 2; k <= n; ++k) {
        const double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (z * p1 - p0) / (z * z - 1.0);
      const double dz = p1 / dp;
      z -= dz;
      if (std::fabs(dz) < 1e-16) break;
    }
    x[i] = -z;
    x[n - 1 - i] = z;
    w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * dp * dp);
  }
}

std::vector<double> SphericalHarmonicsAtPhiZero(int l_max, int m_max, double x) {
  const int stride = m_max + 1;
  std::vector<double> y((l_max + 1) * stride, 0.0);
  const double pi = std::acos(-1.0);
  const double s = std::sqrt(std::max(0.0, 1.0 - x * x));
  double pmm = 1.0 / std::sqrt(4.0 * pi);  // normalised P_m^m, Condon-Shortley phase
  for (int m = 0; m <= std::min(m_max, l_max); ++m) {
    if (m > 0) pmm *= -std::sqrt((2.0 * m + 1.0) / (2.0 * m)) * s;
    y[m * stride + m] = pmm;
    if (m + 1 > l_max) continue;
    double p_lm2 = pmm, p_lm1 = x * std::sqrt(2.0 * m + 3.0) * pmm;
    y[(m + 1) * stride + m] = p_lm1;
    for (int l = m + 2; l <= l_max; ++l) {
      const double a = std::sqrt((4.0 * l * l - 1.0) / (static_cast<double>(l) * l - static_cast<double>(m) * m));
      const double b = std::sqrt((static_cast<double>(l - 1) * (l - 1) - static_cast<double>(m) * m) /
                                 (4.0 * (l - 1) * (l - 1) - 1.0));
      const double p = a * (x * p_lm1 - b * p_lm2);
      y[l * stride + m] = p;
      p_lm2 = p_lm1;
      p_lm1 = p;
    }
  }
  return y;
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
