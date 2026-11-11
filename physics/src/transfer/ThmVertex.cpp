#include "azure2/physics/transfer/ThmVertex.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "azure2/physics/Constants.h"
#include "azure2/physics/coulomb/Coulomb.h"
#include "azure2/physics/transfer/RadialNumerics.h"

namespace azure2 {
namespace physics {
namespace transfer {

namespace {

/// x j_l(x) and d/dx [x j_l(x)]: the Coulomb functions at eta = 0.
coulomb::CoulombWaves Riccati(int l, double x) { return CoulombForMatching(l, 0.0, x); }

void CheckLevel(const LevelFunction &v) {
  if (!(v.step_fm > 0.0) || v.u.size() < 6 || v.potential_MeV.size() != v.u.size() || !(v.mu_MeV > 0.0))
    throw std::invalid_argument("THM vertex: the level function needs h > 0, at least 6 nodes, V on the same grid, mu > 0");
}

}  // namespace

double PlaneWaveSurfaceVertex(int l, double p, double a, double B) {
  if (l < 0 || !(p > 0.0) || !(a > 0.0)) throw std::invalid_argument("M_l: l >= 0, p > 0, a > 0");
  const double x = p * a;
  const coulomb::CoulombWaves c = Riccati(l, x);
  const double jl = c.F / x;                // j_l(x)
  const double x_djl = c.Fp - jl;           // x j_l'(x) = d(x j_l)/dx - j_l
  return (B - 1.0) * jl - x_djl;
}

double BoundaryValue(const LevelFunction &v) {
  CheckLevel(v);
  const int n = static_cast<int>(v.u.size()) - 1;
  const double h = v.step_fm, a = n * h;
  const std::vector<double> &u = v.u;
  const double du = (25.0 * u[n] - 48.0 * u[n - 1] + 36.0 * u[n - 2] - 16.0 * u[n - 3] + 3.0 * u[n - 4]) / (12.0 * h);
  return a * du / u[n];
}

EntranceVertexValue EntranceVertex(const LevelFunction &v, const std::vector<std::complex<double>> &w,
                                   std::complex<double> dw_da, double eps, double B) {
  CheckLevel(v);
  if (w.size() != v.u.size()) throw std::invalid_argument("THM vertex: the wave must be on the level's grid");
  const std::size_t n = v.u.size();
  const double a = (n - 1) * v.step_fm;
  std::vector<std::complex<double>> wvu(n), wu(n);
  for (std::size_t i = 0; i < n; ++i) {
    wvu[i] = w[i] * v.potential_MeV[i] * v.u[i];
    wu[i] = w[i] * v.u[i];
  }
  const double hbar2_2mu = kHbarC_MeV_fm * kHbarC_MeV_fm / (2.0 * v.mu_MeV);
  EntranceVertexValue r;
  r.volume = Simpson(wvu, v.step_fm);
  r.surface = hbar2_2mu * v.u[n - 1] * (B * w[n - 1] / a - dw_da);
  r.interior = (v.energy_MeV - eps) * Simpson(wu, v.step_fm);
  return r;
}

EntranceVertexValue PlaneWaveEntranceVertex(const LevelFunction &v, int l, double p) {
  CheckLevel(v);
  if (l < 0 || !(p > 0.0)) throw std::invalid_argument("THM vertex: l >= 0, p > 0");
  const std::size_t n = v.u.size();
  std::vector<std::complex<double>> w(n, 0.0);
  for (std::size_t i = 1; i < n; ++i) w[i] = Riccati(l, p * i * v.step_fm).F / p;  // r j_l(pr)
  const std::complex<double> dw = Riccati(l, p * (n - 1) * v.step_fm).Fp;    // d(r j_l(pr))/dr
  const double eps = kHbarC_MeV_fm * kHbarC_MeV_fm * p * p / (2.0 * v.mu_MeV);
  return EntranceVertex(v, w, dw, eps, BoundaryValue(v));
}

}  // namespace transfer
}  // namespace physics
}  // namespace azure2
