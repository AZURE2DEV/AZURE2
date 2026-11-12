// AstrophysicalFactor.h: the modified S factor's conversion E exp(2 pi eta + g E) against mpmath (30 digits, the
// constants of Constants.h, 12C + 12C with masses 12 u), its ratio to the S factor's conversion (exp(g E)), g = 0
// bit for bit the S factor's, and NaN for an uncharged pair or E <= 0.
#include <cmath>
#include <cstdio>
#include <string>

#include "azure2/physics/AstrophysicalFactor.h"
#include "azure2/physics/Constants.h"

using namespace azure2::physics;

namespace {
int failures = 0;
void check(const std::string &what, bool ok, double got = 0.0, double want = 0.0) {
  std::printf("  %s  %s", ok ? "ok  " : "FAIL", what.c_str());
  if (!ok) std::printf("  (got %.17g, want %.17g)", got, want);
  std::printf("\n");
  if (!ok) failures++;
}
bool close(double a, double b, double rel) { return std::fabs(a - b) <= rel * std::fabs(b); }

ChannelKinematics Pair(int z1, int z2, double m1u, double m2u) {
  ChannelKinematics k;
  k.z1 = z1;
  k.z2 = z2;
  k.mu_u = m1u * m2u / (m1u + m2u);
  k.m1_MeV = m1u * kAtomicMassUnit_MeV;
  k.m2_MeV = m2u * kAtomicMassUnit_MeV;
  k.relativistic = false;
  return k;
}
}  // namespace

int main() {
  const ChannelKinematics cc = Pair(6, 6, 12.0, 12.0);
  // mpmath: E exp(2 pi eta + g E), eta = 36 alpha sqrt(mu u / (2 E)), mu = 6 u
  struct Row {
    double E, g, want;
  };
  const Row rows[] = {{2.0, 0.46, 3.1364658505264113648e+27},
                      {5.0, 0.46, 4416080488780248533.7},
                      {1.0, -0.3, 5.8180812229774902391e+37}};
  for (const Row &r : rows) {
    const double got = modified_sfactor_conversion(cc, r.E, r.g);
    check("S* conversion of 12C+12C at E = " + std::to_string(r.E) + " MeV, g = " + std::to_string(r.g) +
              " /MeV against mpmath to 1e-12",
          close(got, r.want, 1e-12), got, r.want);
    const double ratio = got / sfactor_conversion(cc, r.E);
    check("S* / S conversion = exp(g E)", close(ratio, std::exp(r.g * r.E), 1e-15), ratio, std::exp(r.g * r.E));
  }
  for (double E : {0.5, 2.0, 7.3}) {
    const double a = modified_sfactor_conversion(cc, E, 0.0), b = sfactor_conversion(cc, E);
    check("g = 0 gives the S factor's conversion bit for bit at E = " + std::to_string(E), a == b, a, b);
  }
  check("NaN for an uncharged pair", std::isnan(modified_sfactor_conversion(Pair(0, 6, 1.0087, 12.0), 1.0, 0.46)));
  check("NaN for E <= 0", std::isnan(modified_sfactor_conversion(cc, 0.0, 0.46)));
  std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "passed", failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
