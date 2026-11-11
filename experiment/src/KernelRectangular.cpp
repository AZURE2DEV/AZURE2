// The rectangular energy resolution (Kernels.h): the average of the node values over a window of full width w about
// the point's energy, w a parameter of the detector's energy_resolution (format 3: energy_resolution rectangular width
// <w> <unit>).
#include <algorithm>

#include "azure2/experiment/KernelQuadrature.h"
#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LinearForm.h"

namespace azure2 {
namespace experiment {
namespace {

/// The c.m. half width of the window: w (lab, as the record gives it) times the point's c.m./lab factor, halved.
double HalfWidth(const KernelParameters &p, double cmFactor) { return 0.5 * cmFactor * p.values[0]; }

/*!
 *   Y(E0) = (1/w) int_{E0 - w/2}^{E0 + w/2} sigma(E) dE,
 * the node values linear between nodes, 2-point Gauss-Legendre on every interval (exact for that piecewise-linear
 * sigma), normalised by the integral of 1 over the same window, so a constant survives exactly.
 */
template <class V>
V FoldRectangular(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes, const V *values) {
  if (nodes.n < 2) return V(0.0);
  const double half = HalfWidth(p, point.cmFactor);
  auto one = [](double) { return 1.0; };
  V num = 0.0;
  double den = 0.0;
  GridQuadrature<V>(nodes, values, point.energy - half, point.energy + half, false, one, one, num, den);
  return (den > 0.0) ? num / den : V(0.0);
}

/// dY/dw = [sigma(E0 + w/2) + sigma(E0 - w/2)] / (2 w) - Y / w (lab w; the c.m. edges move by f/2 per unit of w).
double RectangularDerivative(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes,
                             const double *values, int parameter) {
  if (parameter != 0 || nodes.n < 2) return 0.0;
  const double w = p.values[0];
  const double half = HalfWidth(p, point.cmFactor);
  const double y = FoldRectangular<double>(p, point, nodes, values);
  const double edges = InterpolateNodes<double>(nodes, values, point.energy + half) +
                       InterpolateNodes<double>(nodes, values, point.energy - half);
  return edges / (2.0 * w) - y / w;
}

/// The window itself.
void RectangularSupport(const KernelParameters &p, const SupportPoint &s, double &high, double &low) {
  const double half = HalfWidth(p, s.cmFactor);
  high = s.energy + half;
  low = std::max(s.energy - half, s.minimumEnergy);
}

}  // namespace

Kernel RectangularKernel() {
  Kernel k;
  k.name = "rectangular";
  k.keyword = "detector.energy_resolution rectangular";
  k.tableForm = true;
  k.parameters = {{"width", "energy", "MeV", true, "full width w of the window, lab"}};
  k.fold = FoldRectangular<double>;
  k.weights = FoldRectangular<LinearForm>;
  k.derivative = RectangularDerivative;
  k.support = RectangularSupport;
  k.doc = "a rectangular energy resolution: the average of the cross section over [E0 - w/2, E0 + w/2]";
  return k;
}

}  // namespace experiment
}  // namespace azure2
