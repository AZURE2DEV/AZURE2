// inference_diagonal_priors: the objective's diagonal priors as one list (azure2/inference/DiagonalPriors.h).  Values
// by hand; the running sum and the subtotals bit for bit the loops of GaussianPriorChi2 they replace (a sum formed in
// another order, or one entry left out of a subtotal, fails); the gradient and the curvature against central
// differences of the value at p, not at the entry's value (objective_arrays.md R3); an empty list adds nothing.
#include "azure2/inference/DiagonalPriors.h"

#include <cmath>

#include "check.h"

using namespace azure2::inference;

int main() {
  // Five priors of three sources: two shifts of segment 2 (MeV), a norm of segment 2, a screening U_e (eV), and an
  // apparatus parameter of row 1.
  DiagonalPriors list;
  list.push_back({3, 1.06, 1.00, 0.05, PriorSource::kNorm, 2});           // pull 1.2, chi2 1.44
  list.push_back({4, 0.003, 0.0, 0.002, PriorSource::kShift, 2});         // pull 1.5, chi2 2.25
  list.push_back({5, 0.1, 0.0, 0.05, PriorSource::kShiftSqrt, 2});        // pull 2, chi2 4
  list.push_back({7, 300.0, 340.0, 20.0, PriorSource::kScreening, 0});    // pull -2, chi2 4
  list.push_back({9, 1.25, 1.0, 0.5, PriorSource::kApparatus, 1});        // pull 0.5, chi2 0.25
  const double byHand = 1.44 + 2.25 + 4.0 + 4.0 + 0.25;

  check("sum by hand: 1.44 + 2.25 + 4 + 4 + 0.25 = 11.94", std::fabs(PriorsChi2(list) - byHand) < 1e-12);
  // The running sum of a term (SegmentTerms adds each prior to its chi^2): the loop it replaces, bit for bit.
  double running = 123.456, loop = 123.456;
  AddPriorsChi2(list, running);
  for (const DiagonalPrior &e : list) loop += GaussianPriorChi2(e.value, e.centre, e.sigma);
  check("AddPriorsChi2 = the running loop of GaussianPriorChi2, bit for bit", running == loop);
  check("AddPriorsChi2 adds into the sum, not from 0", running != PriorsChi2(list));
  // A subtotal from 0 (WorkingCopy's sums, the writer's totals): the loop from 0 over the source, bit for bit.
  double shifts = 0.0;
  shifts += GaussianPriorChi2(0.003, 0.0, 0.002);
  check("PriorsChi2(source) = the loop from 0, bit for bit", PriorsChi2(list, PriorSource::kShift) == shifts);
  check("PriorsChi2(source, group): the norm of segment 2 = 1.44",
        std::fabs(PriorsChi2(list, PriorSource::kNorm, 2) - 1.44) < 1e-13);
  check("PriorsChi2(source, another group) = 0", PriorsChi2(list, PriorSource::kNorm, 1) == 0.0);
  check("PriorsChi2(apparatus, row 1) = 0.25, row 0 = 0",
        std::fabs(PriorsChi2(list, PriorSource::kApparatus, 1) - 0.25) < 1e-15 &&
            PriorsChi2(list, PriorSource::kApparatus, 0) == 0.0);
  check("PriorsChi2 of a source the list lacks = 0", PriorsChi2(list, PriorSource::kDensity) == 0.0);
  const double sum = PriorsChi2(list, PriorSource::kNorm) + PriorsChi2(list, PriorSource::kShift) +
                     PriorsChi2(list, PriorSource::kShiftSqrt) + PriorsChi2(list, PriorSource::kScreening) +
                     PriorsChi2(list, PriorSource::kApparatus);
  check("the subtotals add up to the whole (1e-13)", std::fabs(sum - PriorsChi2(list)) < 1e-13);

  // The gradient and the curvature at a vector p whose components differ from the entries' values (the value is
  // the clamped one; the gradient is on p).  Index 0 carries no prior; index 8 is outside every entry.
  std::vector<double> p(10, 0.0);
  p[3] = 1.02, p[4] = -0.001, p[5] = 0.07, p[7] = 355.0, p[9] = 0.9;
  std::vector<double> grad(10, 0.5), g2(10, 0.25);
  AddPriorsGradient(list, p, grad);
  AddPriorsCurvature(list, g2);
  auto valueAt = [&](const std::vector<double> &q) {
    DiagonalPriors at = list;
    for (DiagonalPrior &e : at) e.value = q[e.index];
    return PriorsChi2(at);
  };
  bool gradient = true, curvature = true;
  for (const DiagonalPrior &e : list) {
    const int i = e.index;
    const double h = 1e-4 * e.sigma;
    std::vector<double> up = p, down = p;
    up[i] += h, down[i] -= h;
    const double fd = (valueAt(up) - valueAt(down)) / (2.0 * h);
    const double fd2 = (valueAt(up) - 2.0 * valueAt(p) + valueAt(down)) / (h * h);
    gradient = gradient && std::fabs((grad[i] - 0.5) - fd) <= 1e-7 * std::fabs(fd);
    curvature = curvature && std::fabs((g2[i] - 0.25) - fd2) <= 1e-4 * std::fabs(fd2);
  }
  check("gradient at p = central difference of the value (1e-7 relative), added to what was there", gradient);
  check("curvature = second difference of the value (1e-4 relative), added to what was there", curvature);
  check("gradient by hand at p[3]: 0.5 + 2 x 0.02 / 0.05^2 = 16.5", std::fabs(grad[3] - 16.5) < 1e-11);
  check("curvature by hand at index 7: 0.25 + 2/20^2 = 0.255", std::fabs(g2[7] - 0.255) < 1e-15);
  check("indices without a prior are untouched", grad[0] == 0.5 && grad[8] == 0.5 && g2[0] == 0.25 && g2[8] == 0.25);

  // An entry without an index (a gather without a map), or with one beyond the vectors, adds no gradient.
  DiagonalPriors noIndex = {{-1, 2.0, 1.0, 1.0, PriorSource::kDoppler, 0}, {12, 2.0, 1.0, 1.0, PriorSource::kDoppler, 0}};
  std::vector<double> gradNone(10, 0.0), g2None(10, 0.0);
  AddPriorsGradient(noIndex, p, gradNone);
  AddPriorsCurvature(noIndex, g2None);
  bool untouched = true;
  for (int i = 0; i < 10; i++) untouched = untouched && gradNone[i] == 0.0 && g2None[i] == 0.0;
  check("index -1 or beyond the vector: no gradient, no curvature", untouched);
  check("... but its value counts: ((2 - 1)/1)^2 twice = 2", PriorsChi2(noIndex) == 2.0);

  // An empty list adds nothing.
  const DiagonalPriors empty;
  double unchanged = 7.25;
  AddPriorsChi2(empty, unchanged);
  std::vector<double> g(3, 1.0);
  AddPriorsGradient(empty, g, g);
  check("an empty list: chi2 unchanged, total 0, no gradient",
        unchanged == 7.25 && PriorsChi2(empty) == 0.0 && g[0] == 1.0 && g[2] == 1.0);
  return result();
}
