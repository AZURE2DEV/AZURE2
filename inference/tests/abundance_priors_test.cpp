// inference_abundance_priors: the abundance term (T7) on arrays (azure2/inference/AbundancePriors.h).  Two samples of
// two nuclei by hand, one prior replaced (inactive), one sample without a sum term; the value bit for bit the
// written-out loop of the objective (priors then sum, sample by sample); the sampler's part = the value minus the
// priors; the gradient against central differences; the rows' squares = the sums' part, their derivative against
// central differences of the residual.
#include "azure2/inference/AbundancePriors.h"

#include <cmath>

#include "check.h"

using namespace azure2::inference;

int main() {
  AbundanceBlock b;
  b.numSamples = 2;
  b.numNuclei = 2;
  // Sample 1: w = (0.70, 0.32), the second prior replaced; the sum 1.02 against 1 with s = 0.01.
  // Sample 2: w = (0.90, 0.20), both priors active, no sum term (s = 0).
  b.w = {0.70, 0.32, 0.90, 0.20};
  b.centre = {0.75, 0.25, 0.92, 0.10};
  b.sigma = {0.05, 0.02, 0.01, 0.05};
  b.active = {1, 0, 1, 1};
  b.sumCentre = {1.0, 1.0};
  b.sumSigma = {0.01, 0.0};
  // By hand: priors (-0.05/0.05)^2 = 1, (-0.02/0.01)^2 = 4, (0.1/0.05)^2 = 4; sum of sample 1 (0.02/0.01)^2 = 4.
  const double value = AbundanceChi2(b, true), sums = AbundanceChi2(b, false);
  check("value by hand: 1 + 4 + 4 + 4 = 13 (1e-12)", std::fabs(value - 13.0) < 1e-12);
  check("sampler's part by hand: the sum of sample 1, 4 (1e-12)", std::fabs(sums - 4.0) < 1e-12);
  check("AbundanceSum: 1.02 and 1.10 (1e-15)",
        std::fabs(AbundanceSum(b, 1) - 1.02) < 1e-15 && std::fabs(AbundanceSum(b, 2) - 1.10) < 1e-15);

  // The objective's loop, written out (AbundanceTerm before the arrays): one sum from 0, sample by sample, the
  // active priors in nucleus order, then the sum term.  Bit for bit.
  double loop = 0.0;
  loop += GaussianPriorChi2(0.70, 0.75, 0.05);
  loop += GaussianPriorChi2(0.70 + 0.32, 1.0, 0.01);
  loop += GaussianPriorChi2(0.90, 0.92, 0.01);
  loop += GaussianPriorChi2(0.20, 0.10, 0.05);
  check("value = the objective's loop, bit for bit", value == loop);
  double priors = 0.0;
  for (int k = 0; k < 4; k++)
    if (b.active[k]) priors += GaussianPriorChi2(b.w[k], b.centre[k], b.sigma[k]);
  check("sampler's part = value - priors (1e-12)", std::fabs(sums - (value - priors)) < 1e-12);

  // The gradient: by hand 2*0.02/1e-4 - 2*0.05/0.0025 = 360, 400 (replaced prior), -400, 80; and against central
  // differences of the value.
  std::vector<double> dw;
  AbundanceGradient(b, dw);
  check("gradient by hand: 360, 400, -400, 80 (1e-9)",
        dw.size() == 4 && std::fabs(dw[0] - 360.0) < 1e-9 && std::fabs(dw[1] - 400.0) < 1e-9 &&
            std::fabs(dw[2] + 400.0) < 1e-9 && std::fabs(dw[3] - 80.0) < 1e-9);
  bool central = true;
  for (int k = 0; k < 4; k++) {
    const double h = 1e-6;
    AbundanceBlock up = b, down = b;
    up.w[k] += h;
    down.w[k] -= h;
    const double fd = (AbundanceChi2(up, true) - AbundanceChi2(down, true)) / (2.0 * h);
    if (!(std::fabs(fd - dw[k]) < 1e-5 * (1.0 + std::fabs(dw[k])))) central = false;
  }
  check("gradient = central differences of the value (1e-5 relative)", central);

  // The rows: one, for sample 1, r = 2 and dr/dw = 1/s = 100 on both abundances of the sample.
  const std::vector<AbundanceSumRow> rows = AbundanceSumRows(b);
  check("one row, of sample 1, r = 2, dr/dw = 100 (1e-12)",
        rows.size() == 1 && rows[0].sample == 1 && std::fabs(rows[0].residual - 2.0) < 1e-12 &&
            std::fabs(rows[0].derivative - 100.0) < 1e-12);
  double squares = 0.0;
  for (const AbundanceSumRow &r : rows) squares += r.residual * r.residual;
  check("sum of the rows' squares = the sums' part of the value (1e-12)", std::fabs(squares - sums) < 1e-12);
  bool rowCentral = true;
  for (int k = 0; k < 2; k++) {
    const double h = 1e-6;
    AbundanceBlock up = b, down = b;
    up.w[k] += h;
    down.w[k] -= h;
    const double fd = (AbundanceSumRows(up)[0].residual - AbundanceSumRows(down)[0].residual) / (2.0 * h);
    if (!(std::fabs(fd - rows[0].derivative) < 1e-6 * rows[0].derivative)) rowCentral = false;
  }
  check("dr/dw = central differences of the residual (1e-6 relative)", rowCentral);
  // Gauss-Newton: 2 r dr/dw of the rows is the sums' part of the gradient (sample 1: 2 * 2 * 100 = 400).
  check("2 r dr/dw = the sum term's gradient (1e-9)", std::fabs(2.0 * rows[0].residual * rows[0].derivative - 400.0) < 1e-9);

  // Nothing without a block.
  AbundanceBlock empty;
  std::vector<double> none;
  AbundanceGradient(empty, none);
  check("an empty block: value 0, no gradient, no rows",
        AbundanceChi2(empty, true) == 0.0 && none.empty() && AbundanceSumRows(empty).empty());
  return result();
}
