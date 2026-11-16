// inference_gaussian_prior: one Gaussian prior ((x - x0)/sigma)^2 (azure2/inference/GaussianPrior.h): values by hand,
// the gradient and the curvature against central differences of the value.  A factor 2 lost in the gradient, or a
// sigma in place of sigma^2, fails here.
#include "azure2/inference/GaussianPrior.h"

#include <cmath>

#include "check.h"

using namespace azure2::inference;

int main() {
  // A normalization 1.06 with a prior 1.00 +- 0.05: pull 1.2, chi2 1.44.
  check("chi2 by hand: ((1.06 - 1.00)/0.05)^2 = 1.44", std::fabs(GaussianPriorChi2(1.06, 1.00, 0.05) - 1.44) < 1e-13);
  check("zero at the centre", GaussianPriorChi2(2.5, 2.5, 0.1) == 0.0 && GaussianPriorGradient(2.5, 2.5, 0.1) == 0.0);
  check("symmetric about the centre", GaussianPriorChi2(0.75, 1.0, 0.5) == GaussianPriorChi2(1.25, 1.0, 0.5));
  check("one unit at one sigma", std::fabs(GaussianPriorChi2(300.0 + 20.0, 300.0, 20.0) - 1.0) < 1e-15);
  // The value is exactly the square of the pull, as the engine always wrote it.
  const double pull = (7.3 - 6.9) / 0.37;
  check("the square of the pull, bit for bit", GaussianPriorChi2(7.3, 6.9, 0.37) == pull * pull);

  bool gradient = true, curvature = true;
  const double cases[][3] = {{1.06, 1.0, 0.05}, {-3.0e-3, 0.0, 2.0e-3}, {2.1e15, 2.0e15, 3.0e14}, {250.0, 300.0, 20.0}};
  for (const auto &k : cases) {
    const double x = k[0], c = k[1], s = k[2], h = 1e-4 * s;
    const double fd = (GaussianPriorChi2(x + h, c, s) - GaussianPriorChi2(x - h, c, s)) / (2.0 * h);
    const double fd2 = (GaussianPriorChi2(x + h, c, s) - 2.0 * GaussianPriorChi2(x, c, s) + GaussianPriorChi2(x - h, c, s)) / (h * h);
    gradient = gradient && std::fabs(GaussianPriorGradient(x, c, s) - fd) <= 1e-8 * std::fabs(fd);
    curvature = curvature && std::fabs(GaussianPriorCurvature(s) - fd2) <= 1e-5 * std::fabs(fd2);
  }
  check("gradient = central difference of the value (1e-8 relative)", gradient);
  check("curvature = second difference of the value (1e-5 relative)", curvature);
  check("gradient by hand: 2 x 0.06 / 0.05^2 = 48", std::fabs(GaussianPriorGradient(1.06, 1.0, 0.05) - 48.0) < 1e-11);
  check("curvature by hand: 2 / 0.05^2 = 800", std::fabs(GaussianPriorCurvature(0.05) - 800.0) < 1e-10);
  return result();
}
