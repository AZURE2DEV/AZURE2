// Widths.h against closed forms computed here:
//  - gamma_W^2 = 3 hbar^2 / (2 mu a^2) for p + 12C at a = 5 fm and alpha + alpha at a = 4 fm, written out with the
//    constants of Constants.h; NaN for mu <= 0 or a <= 0;
//  - Gamma_sp = 2 P gamma_W^2 for a neutral s wave, where P_0 = rho = k a exactly (Lane and Thomas, Rev. Mod. Phys. 30
//    (1958) 257, table III), and for a charged channel against the penetrability of channel_functions;
//  - (|E| / hbar c)^(2L+1) for E1 and E2 photons, with the sign of E dropped;
//  - the Weisskopf estimate Gamma_W(XL) in eV against the closed forms of Blatt and Weisskopf (as Bohr and Mottelson
//    vol. 1 write them, eqs. 3C-38 and 3C-39) with R = 1.2 A^(1/3) fm, written out here from alpha, hbar c and the
//    proton mass: equal to 1e-13 for E1-E4 and M1-M4 (until 2026-11-18 the function held the two-digit textbook
//    coefficients, E1 +0.79 %, E2 +2.28 %, E3 +3.05 %, E4 -3.11 %, M1 +1.34 %, M2 +1.94 %, M3 -0.80 %, M4 -7.23 % from
//    the closed form; the table stays here, its offsets held to 0.01 %); two numbers by hand (an E1 of 1 MeV in A = 125,
//    an M1 of 2 MeV).  NaN outside E1-E4, M1-M4, for E <= 0 and A <= 0.
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

  // The Weisskopf estimate.  Closed form in eV for R = r0 A^(1/3):
  //   electric: 2 (L+1) / (L [(2L+1)!!]^2) (3/(L+3))^2 alpha hbar c (E/hbar c)^(2L+1) R^(2L)
  //   magnetic: the same with R^(2L-2) in place of R^(2L), times 40 (hbar c / (2 m_p c^2))^2
  //             (B_W(ML) = (10/pi) (3/(L+3))^2 R^(2L-2) mu_N^2 against B_W(EL) = (1/4pi) (3/(L+3))^2 R^(2L) e^2).
  const double alpha = 1.0 / 137.035999084, mp = 938.27208816, r0 = 1.2;
  const auto closed = [&](bool magnetic, int L, double E, double A) {
    double dfact = 1.0;
    for (int k = 2 * L + 1; k > 0; k -= 2) dfact *= k;
    const double R = r0 * std::cbrt(A);
    const double common = 2.0 * (L + 1) / (L * dfact * dfact) * std::pow(3.0 / (L + 3), 2) * alpha * hc *
                          std::pow(E / hc, 2 * L + 1);
    const double mev = magnetic ? common * 40.0 * std::pow(hc / (2.0 * mp), 2) * std::pow(R, 2 * L - 2)
                                : common * std::pow(R, 2 * L);
    return mev * 1e6;
  };
  // The function is the closed form: to 1e-13 over A = 4-208, E = 0.05-20 MeV, every multipole.  The two-digit table
  // it replaced on 2026-11-18 (the owner's answer) stands beside it, so that the move of every W.u. stays on record:
  // table / closed form - 1 per multipole, in percent, held to 0.01 %.
  const double tableE[4] = {6.8e-2, 4.9e-8, 2.3e-14, 6.8e-21}, tableM[4] = {2.1e-2, 1.5e-8, 6.8e-15, 2.0e-21};
  const double offE[4] = {0.789, 2.285, 3.054, -3.114}, offM[4] = {1.338, 1.942, -0.805, -7.225};
  for (int magnetic = 0; magnetic < 2; magnetic++) {
    for (int L = 1; L <= 4; L++) {
      const double want = (magnetic ? offM : offE)[L - 1];
      const double coefficient = (magnetic ? tableM : tableE)[L - 1];
      const int thirds = magnetic ? 2 * L - 2 : 2 * L;
      double worst = 0.0, tableWorst = 0.0;
      for (double A : {4.0, 13.0, 27.0, 100.0, 208.0})
        for (double E : {0.05, 0.5, 2.0, 8.0, 20.0}) {
          const double c = closed(magnetic, L, E, A);
          worst = std::fmax(worst, std::fabs(weisskopf_width_eV(magnetic, L, E, A) / c - 1.0));
          const double table = coefficient * std::pow(A, thirds / 3.0) * std::pow(E, 2 * L + 1);
          tableWorst = std::fmax(tableWorst, std::fabs((table / c - 1.0) * 100.0 - want));
        }
      const std::string name = std::string(magnetic ? "M" : "E") + std::to_string(L);
      check("Weisskopf " + name + ": the closed form (r0 = 1.2 fm) to 1e-13, A = 4-208, E = 0.05-20 MeV", worst <= 1e-13,
            worst, 1e-13);
      check("Weisskopf " + name + ": the former table / closed form - 1 = " + std::to_string(want) + " % to 0.01 %",
            tableWorst <= 0.01, tableWorst, 0.01);
    }
  }
  check("Weisskopf E1, 1 MeV, A = 125: 1.68669 eV (the former table: 6.8e-2 x 125^(2/3) = 1.7 eV)",
        close(weisskopf_width_eV(false, 1, 1.0, 125.0), 1.6866885576579, 1e-12), weisskopf_width_eV(false, 1, 1.0, 125.0),
        1.6866885576579);
  check("Weisskopf M1 does not depend on A: 0.165782 eV at 2 MeV (the former table: 2.1e-2 E^3 = 0.168 eV)",
        weisskopf_width_eV(true, 1, 2.0, 13.0) == weisskopf_width_eV(true, 1, 2.0, 208.0) &&
            close(weisskopf_width_eV(true, 1, 2.0, 13.0), 0.16578213605135, 1e-12),
        weisskopf_width_eV(true, 1, 2.0, 13.0), 0.16578213605135);
  check("Weisskopf is NaN for L = 0, L = 5, E = 0, E < 0, A = 0",
        std::isnan(weisskopf_width_eV(false, 0, 1.0, 13.0)) && std::isnan(weisskopf_width_eV(true, 5, 1.0, 13.0)) &&
            std::isnan(weisskopf_width_eV(false, 1, 0.0, 13.0)) && std::isnan(weisskopf_width_eV(false, 1, -1.0, 13.0)) &&
            std::isnan(weisskopf_width_eV(true, 2, 1.0, 0.0)));

  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "ok", failures);
  return failures ? 1 : 0;
}
