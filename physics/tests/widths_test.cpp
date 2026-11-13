// Widths.h against closed forms computed here:
//  - gamma_W^2 = 3 hbar^2 / (2 mu a^2) for p + 12C at a = 5 fm and alpha + alpha at a = 4 fm, written out with the
//    constants of Constants.h; NaN for mu <= 0 or a <= 0;
//  - Gamma_sp = 2 P gamma_W^2 for a neutral s wave, where P_0 = rho = k a exactly (Lane and Thomas, Rev. Mod. Phys. 30
//    (1958) 257, table III), and for a charged channel against the penetrability of channel_functions;
//  - (|E| / hbar c)^(2L+1) for E1 and E2 photons, with the sign of E dropped.
//   physics_widths_test
#include <cmath>
#include <cstdio>
#include <string>

#include "azure2/physics/Constants.h"
#include "azure2/physics/Widths.h"

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
ChannelKinematics Pair(int z1, int z2, double m1_u, double m2_u) {
  ChannelKinematics k;
  k.z1 = z1;
  k.z2 = z2;
  k.mu_u = m1_u * m2_u / (m1_u + m2_u);
  k.m1_MeV = m1_u * kAtomicMassUnit_MeV;
  k.m2_MeV = m2_u * kAtomicMassUnit_MeV;
  return k;
}
}  // namespace

int main() {
  const double hc = 197.3269804, u = 931.49410242;
  // p + 12C, a = 5 fm: mu = 1.00727646688 x 12 / 13.00727646688 u.
  const double mu_p12C = 1.00727646688 * 12.0 / 13.00727646688;
  const double want_p = 3.0 * hc * hc / (2.0 * mu_p12C * u * 25.0);
  check("gamma_W^2 of p + 12C at 5 fm = 3 (hbar c)^2 / (2 mu c^2 a^2) = 2.6990 MeV",
        close(wigner_limit(mu_p12C, 5.0), want_p, 1e-15) && close(want_p, 2.69899, 1e-5), wigner_limit(mu_p12C, 5.0),
        want_p);
  const double mu_aa = 4.00260325413 / 2.0;
  const double want_a = 3.0 * hc * hc / (2.0 * mu_aa * u * 16.0);
  check("gamma_W^2 of alpha + alpha at 4 fm", close(wigner_limit(mu_aa, 4.0), want_a, 1e-15), wigner_limit(mu_aa, 4.0),
        want_a);
  check("gamma_W^2 is NaN for mu = 0, a = 0, a < 0",
        std::isnan(wigner_limit(0.0, 5.0)) && std::isnan(wigner_limit(1.0, 0.0)) && std::isnan(wigner_limit(1.0, -1.0)));

  // Neutral s wave: P_0 = rho, so Gamma_sp = 2 k a gamma_W^2.  n + 16O at 1 MeV, a = 4 fm.
  const ChannelKinematics n16O = Pair(0, 8, 1.00866491595, 15.99491461957);
  const double k = std::sqrt(2.0 * u * n16O.mu_u * 1.0) / hc;  // fm^-1
  const double want_n = 2.0 * k * 4.0 * wigner_limit(n16O.mu_u, 4.0);
  check("Gamma_sp of n + 16O, l = 0, 1 MeV, 4 fm = 2 k a gamma_W^2",
        close(single_particle_width(n16O, 0, 4.0, 1.0), want_n, 1e-13), single_particle_width(n16O, 0, 4.0, 1.0), want_n);
  // Charged: p + 12C, l = 1, 0.5 MeV, 5 fm, against channel_functions' P.
  const ChannelKinematics p12C = Pair(1, 6, 1.00727646688, 12.0);
  const double P = channel_functions(p12C, 1, 5.0, 0.5).P;
  const double want_c = 2.0 * P * wigner_limit(p12C.mu_u, 5.0);
  check("Gamma_sp of p + 12C, l = 1, 0.5 MeV, 5 fm = 2 P gamma_W^2 (P of channel_functions)",
        single_particle_width(p12C, 1, 5.0, 0.5) == want_c, single_particle_width(p12C, 1, 5.0, 0.5), want_c);
  check("Gamma_sp is NaN below threshold and for a <= 0",
        std::isnan(single_particle_width(p12C, 0, 5.0, -0.1)) && std::isnan(single_particle_width(p12C, 0, 0.0, 1.0)));

  // Photons.
  check("P_gamma of an E1 photon of 1 MeV = (1 / hbar c)^3",
        close(photon_penetrability(1.0, 1), std::pow(1.0 / hc, 3), 1e-15), photon_penetrability(1.0, 1),
        std::pow(1.0 / hc, 3));
  check("P_gamma of an E2 photon of 4.439 MeV = (4.439 / hbar c)^5",
        close(photon_penetrability(4.439, 2), std::pow(4.439 / hc, 5), 1e-15));
  check("P_gamma takes |E_gamma|", photon_penetrability(-2.0, 1) == photon_penetrability(2.0, 1));

  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "ok", failures);
  return failures ? 1 : 0;
}
