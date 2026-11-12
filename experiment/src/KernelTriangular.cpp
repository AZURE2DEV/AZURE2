// The triangular energy resolution (Kernels.h): the average of the node values weighted by a triangle of base width w
// (full width at half maximum w/2) about the point's energy, w a parameter of the detector's energy_resolution (format
// 3: energy_resolution triangular width <w> <unit> fixed|fit ...).
#include "azure2/experiment/KernelQuadrature.h"
#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LinearForm.h"

namespace azure2 {
namespace experiment {
namespace {

/// The c.m. half base h of the triangle: w (lab, as the record gives it) times the point's c.m./lab factor, halved.
double HalfBase(const KernelParameters &p, double cmFactor) { return 0.5 * cmFactor * p.values[0]; }

/*!
 *   Y(E0) = int t(E - E0) sigma(E) dE / int t,  t(x) = 1 - |x|/h on [-h, h],
 * the node values linear between nodes, 2-point Gauss-Legendre on every interval of each half of the triangle (the
 * weight is linear there: exact for that piecewise-linear sigma), normalised by the same quadrature of t.
 */
template <class V>
V FoldTriangular(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes, const V *values) {
  if (nodes.n < 2) return V(0.0);
  const double h = HalfBase(p, point.cmFactor), e0 = point.energy;
  auto t = [&](double e) { return 1.0 - (e > e0 ? e - e0 : e0 - e) / h; };
  V num = 0.0;
  double den = 0.0;
  GridQuadrature<V>(nodes, values, e0 - h, e0, false, t, t, num, den);
  GridQuadrature<V>(nodes, values, e0, e0 + h, false, t, t, num, den);
  return (den > 0.0) ? num / den : V(0.0);
}

/// dY/dw = 2 (R - Y)/w, R the rectangular average over [E0 - h, E0 + h] (t vanishes at its ends: no edge terms).
double TriangularDerivative(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes,
                            const double *values, int parameter) {
  if (parameter != 0 || nodes.n < 2) return 0.0;
  const double w = p.values[0], h = HalfBase(p, point.cmFactor);
  auto one = [](double) { return 1.0; };
  double num = 0.0, den = 0.0;
  GridQuadrature<double>(nodes, values, point.energy - h, point.energy + h, false, one, one, num, den);
  const double r = (den > 0.0) ? num / den : 0.0;
  return 2.0 * (r - FoldTriangular<double>(p, point, nodes, values)) / w;
}

/// The base of the triangle.
void TriangularSupport(const KernelParameters &p, const SupportPoint &s, double &high, double &low) {
  const double h = HalfBase(p, s.cmFactor);
  high = s.energy + h;
  low = (s.energy - h > s.minimumEnergy) ? s.energy - h : s.minimumEnergy;
}

}  // namespace

Kernel TriangularKernel() {
  Kernel k;
  k.name = "triangular";
  k.keyword = "detector.energy_resolution triangular";
  k.tableForm = true;
  k.parameters = {{"width", "energy", "MeV", true, "base width w of the triangle (FWHM w/2), lab"}};
  k.fold = FoldTriangular<double>;
  k.weights = FoldTriangular<LinearForm>;
  k.derivative = TriangularDerivative;
  k.thinLimit = AverageThinLimit;
  k.support = TriangularSupport;
  k.doc = "a triangular energy resolution: the cross section weighted by 1 - |E - E0|/(w/2) over [E0 - w/2, E0 + w/2]";
  return k;
}

}  // namespace experiment
}  // namespace azure2
