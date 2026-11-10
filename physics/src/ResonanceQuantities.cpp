// The quantities of an isolated resonance from observed widths (ResonanceQuantities.h; data_kinds.md section 4.7).
#include "azure2/physics/ResonanceQuantities.h"

#include <cmath>

namespace azure2 {
namespace physics {

double statistical_factor(int two_J, int two_ja, int two_jA, bool identical) {
  return (two_J + 1.0) * (identical ? 2.0 : 1.0) / ((two_ja + 1.0) * (two_jA + 1.0));
}

ResonanceStrength resonance_strength(double omega, double Gamma_a, double Gamma_b, double Gamma) {
  ResonanceStrength out;
  if (Gamma == 0.0) return out;
  out.value = omega * Gamma_a * Gamma_b / Gamma;
  out.d_entrance = omega * Gamma_b / Gamma;
  out.d_exit = omega * Gamma_a / Gamma;
  out.d_total = -out.value / Gamma;
  return out;
}

BranchingRatio branching_ratio(double Gamma_b, double Gamma) {
  BranchingRatio out;
  if (Gamma == 0.0) return out;
  out.value = Gamma_b / Gamma;
  out.d_exit = 1.0 / Gamma;
  out.d_total = -Gamma_b / (Gamma * Gamma);
  return out;
}

double total_width(const std::vector<double> &widths, std::vector<double> *d_widths) {
  double sum = 0.0;
  if (d_widths) d_widths->assign(widths.size(), 0.0);
  for (std::size_t c = 0; c < widths.size(); c++) {
    sum += std::fabs(widths[c]);
    if (d_widths) (*d_widths)[c] = widths[c] > 0.0 ? 1.0 : (widths[c] < 0.0 ? -1.0 : 0.0);
  }
  return sum;
}

}  // namespace physics
}  // namespace azure2
