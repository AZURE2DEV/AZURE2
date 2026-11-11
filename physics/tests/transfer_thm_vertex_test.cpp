// transfer_thm_vertex (docs/upgrade/transfer.md sections 8 and 9): the THM entrance vertex in the plane-wave limit
// against closed forms written out here (not the engine's M_l, not the library's Bessel functions).
//
//  1. M_l(p) = (B - 1) j_l(pa) - pa j_l'(pa) for l = 0, 1, 2, p = 0.1, 0.5, 1.3 fm^-1, a = 5 fm, B = -0.7 and 2.3:
//     PlaneWaveSurfaceVertex against j_0, j_1, j_2 and their derivatives in sines and cosines, 1e-6 absolute (|M| is
//     of order 1).
//  2. The full form factor of a square-well level (V = -40 MeV inside a = 5 fm, mu of p + 18O, E_lambda = 1.2 MeV,
//     u = K r j_l(Kr)), for the same l and p: the volume integral of Eq. (8.1) against the Lommel integral
//     -V0 K a^2 [K j_l(pa) j_l'(Ka) - p j_l'(pa) j_l(Ka)] / (p^2 - K^2), the sum surface + interior of Eq. (8.2)
//     against the same, the surface term against (hbar^2/2mu) u(a) M_l(p; B) with B = 1 + Ka j_l'(Ka)/j_l(Ka) and
//     BoundaryValue against that B; all 1e-6 relative.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/Constants.h"
#include "azure2/physics/transfer/ThmVertex.h"

using namespace azure2::physics;
using namespace azure2::physics::transfer;

namespace {

int failures = 0;

void Check(bool ok, const std::string &what) {
  if (!ok) {
    std::printf("FAIL: %s\n", what.c_str());
    ++failures;
  }
}

/// j_l(x) and j_l'(x) for l = 0, 1, 2 in closed form.
void Bessel(int l, double x, double &j, double &dj) {
  const double s = std::sin(x), c = std::cos(x);
  if (l == 0) {
    j = s / x;
    dj = c / x - s / (x * x);
  } else if (l == 1) {
    j = s / (x * x) - c / x;
    dj = 2.0 * c / (x * x) + s / x - 2.0 * s / (x * x * x);
  } else {
    j = (3.0 / (x * x * x) - 1.0 / x) * s - 3.0 * c / (x * x);
    dj = (-9.0 / (x * x * x * x) + 4.0 / (x * x)) * s + (9.0 / (x * x * x) - 1.0 / x) * c;
  }
}

double ClosedM(int l, double p, double a, double B) {
  double j, dj;
  Bessel(l, p * a, j, dj);
  return (B - 1.0) * j - p * a * dj;
}

void TestSurfaceVertex() {
  double worst = 0.0;
  for (int l = 0; l <= 2; ++l)
    for (double p : {0.1, 0.5, 1.3})
      for (double B : {-0.7, 2.3})
        worst = std::max(worst, std::fabs(PlaneWaveSurfaceVertex(l, p, 5.0, B) - ClosedM(l, p, 5.0, B)));
  char buf[120];
  std::snprintf(buf, sizeof buf, "M_l(p) against the closed forms: worst %.2e", worst);
  Check(worst < 1e-6, buf);
}

void TestFullFormFactor() {
  const double a = 5.0, V0 = 40.0, E_level = 1.2, h = 0.001;
  const double mu = 1.0 * 18.0 / 19.0 * kAtomicMassUnit_MeV;
  const double hb = kHbarC_MeV_fm * kHbarC_MeV_fm / (2.0 * mu);
  const double K = std::sqrt((E_level + V0) / hb);
  const int n = static_cast<int>(std::lround(a / h)) + 1;
  for (int l = 0; l <= 2; ++l) {
    LevelFunction level;
    level.step_fm = h;
    level.energy_MeV = E_level;
    level.mu_MeV = mu;
    level.u.assign(n, 0.0);
    level.potential_MeV.assign(n, -V0);
    for (int i = 1; i < n; ++i) {
      double j, dj;
      Bessel(l, K * i * h, j, dj);
      level.u[i] = K * i * h * j;
    }
    double jKa, djKa;
    Bessel(l, K * a, jKa, djKa);
    const double B = 1.0 + K * a * djKa / jKa;
    char buf[200];
    std::snprintf(buf, sizeof buf, "l = %d: BoundaryValue %.10f against %.10f", l, BoundaryValue(level), B);
    Check(std::fabs(BoundaryValue(level) - B) < 1e-6 * std::max(1.0, std::fabs(B)), buf);
    for (double p : {0.1, 0.5, 1.3}) {
      double jpa, djpa;
      Bessel(l, p * a, jpa, djpa);
      const double lommel = -V0 * K * a * a * (K * jpa * djKa - p * djpa * jKa) / (p * p - K * K);
      const EntranceVertexValue v = PlaneWaveEntranceVertex(level, l, p);
      const double surface = hb * level.u[n - 1] * ClosedM(l, p, a, B);
      const double scale = std::max({std::fabs(lommel), std::fabs(surface), 1e-3});
      std::snprintf(buf, sizeof buf, "l = %d, p = %.1f: volume %.10e, surface + interior %.10e, Lommel %.10e", l, p,
                    v.volume.real(), (v.surface + v.interior).real(), lommel);
      Check(std::fabs(v.volume.real() - lommel) < 1e-6 * scale &&
                std::fabs((v.surface + v.interior).real() - lommel) < 1e-6 * scale,
            buf);
      std::snprintf(buf, sizeof buf, "l = %d, p = %.1f: surface %.10e against (hbar^2/2mu) u(a) M_l = %.10e", l, p,
                    v.surface.real(), surface);
      Check(std::fabs(v.surface.real() - surface) < 1e-6 * scale, buf);
    }
  }
}

}  // namespace

int main() {
  TestSurfaceVertex();
  TestFullFormFactor();
  if (failures) {
    std::printf("transfer_thm_vertex: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("transfer_thm_vertex: all checks passed\n");
  return 0;
}
