#include "azure2/physics/ExternalCapture.h"

#include "azure2/physics/Constants.h"

#include <cmath>
#include <stdexcept>

namespace azure2 {
namespace physics {

double bound_state_coefficient(double radius_fm, double reduced_mass_u) {
  return 2.0 * radius_fm * reduced_mass_u * kAtomicMassUnit_MeV / std::pow(kHbarC_MeV_fm, 2.0);
}

double bound_state_normalization(double n0, const std::vector<BoundChannel> &channels) {
  // Eq. (1), in the engine's order of operations (CNuc.cpp, the six copies before step 4).
  double n = n0;
  for (const BoundChannel &c : channels)
    n += bound_state_coefficient(c.radius_fm, c.reduced_mass_u) * std::pow(c.gamma, 2.0) * c.nf_integral_fm;
  return n;
}

double inverse_sqrt_normalization(double n0, const std::vector<BoundChannel> &channels) {
  const double n = bound_state_normalization(n0, channels);
  if (!(n > 0.0)) throw std::domain_error("inverse_sqrt_normalization: N_f is not positive");
  return 1.0 / std::sqrt(n);
}

bool couples(Multipole X, int L, int l, double s, int l_f, double s_f, double J_f) {
  const bool orbital = std::abs(l - L) <= l_f && l_f <= l + L && std::fabs(s - l_f) <= J_f && J_f <= s + l_f &&
                       s == s_f;
  const bool spin = X == Multipole::M && std::fabs(s - L) <= s_f && s_f <= s + L && std::fabs(l - s_f) <= J_f &&
                    J_f <= l + s_f && l == l_f;
  return orbital || spin;
}

std::complex<double> external_reduced_width(Multipole X, int L, double J_f, const std::vector<CaptureChannel> &level,
                                            const std::vector<FinalChannel> &final_state,
                                            const RadialIntegral &integral) {
  std::complex<double> g(0.0, 0.0);
  if (!((X == Multipole::E && (L == 1 || L == 2)) || (X == Multipole::M && L == 1))) return g;
  for (std::size_t c = 0; c < level.size(); c++)
    for (std::size_t cf = 0; cf < final_state.size(); cf++) {
      const CaptureChannel &i = level[c];
      const FinalChannel &f = final_state[cf];
      if (i.partition != f.partition || !couples(X, L, i.l, i.s, f.l, f.s, J_f)) continue;
      g -= f.scale * i.gamma * integral(c, cf);
    }
  return g;
}

std::vector<double> capture_mask(const std::vector<double> &T, std::size_t n, const std::vector<unsigned char> &masked) {
  if (T.size() != n * n || masked.size() != n) throw std::invalid_argument("capture_mask: T is n x n, one flag per level");
  std::vector<double> Pi(n * n, 0.0);
  for (std::size_t i = 0; i < n; i++) Pi[i * n + i] = 1.0;
  // LU factors of T with partial pivoting, then x_m = T^-1 e_m and Pi -= x_m T_m for each masked m.
  std::vector<double> LU = T;
  std::vector<std::size_t> piv(n);
  for (std::size_t i = 0; i < n; i++) piv[i] = i;
  for (std::size_t k = 0; k < n; k++) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; i++)
      if (std::fabs(LU[i * n + k]) > std::fabs(LU[p * n + k])) p = i;
    if (LU[p * n + k] == 0.0) throw std::domain_error("capture_mask: the transform is singular");
    if (p != k) {
      for (std::size_t j = 0; j < n; j++) std::swap(LU[p * n + j], LU[k * n + j]);
      std::swap(piv[p], piv[k]);
    }
    for (std::size_t i = k + 1; i < n; i++) {
      LU[i * n + k] /= LU[k * n + k];
      for (std::size_t j = k + 1; j < n; j++) LU[i * n + j] -= LU[i * n + k] * LU[k * n + j];
    }
  }
  for (std::size_t m = 0; m < n; m++) {
    if (!masked[m]) continue;
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; i++) x[i] = (piv[i] == m) ? 1.0 : 0.0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < i; j++) x[i] -= LU[i * n + j] * x[j];
    for (std::size_t i = n; i-- > 0;) {
      for (std::size_t j = i + 1; j < n; j++) x[i] -= LU[i * n + j] * x[j];
      x[i] /= LU[i * n + i];
    }
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) Pi[i * n + j] -= x[i] * T[m * n + j];
  }
  return Pi;
}

}  // namespace physics
}  // namespace azure2
