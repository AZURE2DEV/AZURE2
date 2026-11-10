// ResonanceQuantities.h against numbers computed by hand (written here, with their arithmetic) and the derivatives
// against central differences:
//  - T2's level 5 (tests/upgrade/T2/README.md section 1): J = 3-, p + 15N (spins 1/2, 1/2), omega = 7 / (2 x 2) = 7/4;
//    Gamma_p = 1.702 keV, Gamma_a0 = 12.92 keV, Gamma_a1 = 0, so Gamma = 14.622 keV and
//    omega gamma = 1.75 x 1.702 x 12.92 / 14.622 = 2.6318028 keV (docs/upgrade/project_format3.md section 9.3);
//  - identical bosons 12C + 12C (spins 0), J = 2: omega = 5 x 2 / (1 x 1) = 10; alpha + alpha, J = 0: omega = 2;
//  - branching ratio 12.92 / 14.622 = 0.88360006...; total width of signed widths (-689.6, 36.05, 20) = 745.65.
//   physics_resonance_quantities_test
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/physics/ResonanceQuantities.h"

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
}  // namespace

int main() {
  // Statistical factors.
  check("omega of T2 level 5 (3-, p + 15N) = 7/4", statistical_factor(6, 1, 1, false) == 1.75, statistical_factor(6, 1, 1, false), 1.75);
  check("omega of 12C + 12C, J = 2 (identical bosons) = 10", statistical_factor(4, 0, 0, true) == 10.0);
  check("omega of a + a, J = 0 = 2", statistical_factor(0, 0, 0, true) == 2.0);
  check("omega of a + 12C, J = 1 = 3", statistical_factor(2, 0, 0, false) == 3.0);

  // T2 level 5 by hand, keV.
  const double Gp = 1.702, Ga0 = 12.92, Ga1 = 0.0;
  const double G = total_width({Gp, Ga0, Ga1});
  check("Gamma of T2 level 5 = 14.622 keV", close(G, 14.622, 1e-15), G, 14.622);
  const ResonanceStrength wg = resonance_strength(1.75, Gp, Ga0, G);
  const double hand = 1.75 * 1.702 * 12.92 / 14.622;
  check("omega gamma of T2 level 5 = 2.6318028 keV", close(wg.value, 2.6318028, 5e-8), wg.value, 2.6318028);
  check("omega gamma = the expression written out", wg.value == hand, wg.value, hand);
  const BranchingRatio br = branching_ratio(Ga0, G);
  check("branching ratio Gamma_a0 / Gamma = 0.88360006", close(br.value, 0.88360006, 1e-8), br.value, 0.88360006);
  std::vector<double> dG;
  const double Gs = total_width({-689.6, 36.05, 20.0}, &dG);
  check("total width of signed widths = 745.65", close(Gs, 745.65, 1e-15), Gs, 745.65);
  check("d Gamma / d Gamma_c = sign", dG.size() == 3 && dG[0] == -1.0 && dG[1] == 1.0 && dG[2] == 1.0);
  check("no width, no strength", resonance_strength(1.75, 0.0, 0.0, 0.0).value == 0.0 && branching_ratio(0.0, 0.0).value == 0.0);

  // Derivatives against central differences (each argument alone, the others held).
  const double h = 1e-6;
  auto wgOf = [](double a, double b, double g) { return resonance_strength(1.75, a, b, g).value; };
  const double dA = (wgOf(Gp + h, Ga0, G) - wgOf(Gp - h, Ga0, G)) / (2 * h);
  const double dB = (wgOf(Gp, Ga0 + h, G) - wgOf(Gp, Ga0 - h, G)) / (2 * h);
  const double dT = (wgOf(Gp, Ga0, G + h) - wgOf(Gp, Ga0, G - h)) / (2 * h);
  check("d omega gamma / d Gamma_a = central difference", close(wg.d_entrance, dA, 1e-8), wg.d_entrance, dA);
  check("d omega gamma / d Gamma_b = central difference", close(wg.d_exit, dB, 1e-8), wg.d_exit, dB);
  check("d omega gamma / d Gamma = central difference", close(wg.d_total, dT, 1e-8), wg.d_total, dT);
  // The chain rule through Gamma = Gamma_a + Gamma_b + Gamma_other: d/dGamma_a = d_entrance + d_total.
  auto wgSum = [](double a, double b) { return resonance_strength(1.75, a, b, total_width({a, b, 0.0})).value; };
  const double dChain = (wgSum(Gp + h, Ga0) - wgSum(Gp - h, Ga0)) / (2 * h);
  check("d omega gamma / d Gamma_a with Gamma the sum = d_entrance + d_total", close(wg.d_entrance + wg.d_total, dChain, 1e-8),
        wg.d_entrance + wg.d_total, dChain);
  auto brOf = [](double b, double g) { return branching_ratio(b, g).value; };
  const double dBb = (brOf(Ga0 + h, G) - brOf(Ga0 - h, G)) / (2 * h);
  const double dBt = (brOf(Ga0, G + h) - brOf(Ga0, G - h)) / (2 * h);
  check("d BR / d Gamma_b = central difference", close(br.d_exit, dBb, 1e-8), br.d_exit, dBb);
  check("d BR / d Gamma = central difference", close(br.d_total, dBt, 1e-8), br.d_total, dBt);

  if (failures) std::printf("%d FAILED\n", failures);
  else std::printf("all passed\n");
  return failures ? 1 : 0;
}
