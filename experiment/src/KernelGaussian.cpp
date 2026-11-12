// The Gaussian beam-spread kernel (Kernels.h): the normalised Gaussian average of the node values about the point's
// energy, sigma from the beam record (a constant or a formula, KernelParameters::sigma).
#include "azure2/experiment/KernelQuadrature.h"
#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LinearForm.h"

namespace azure2 {
namespace experiment {
namespace {

/*!
 *   Y(E0) = int K(E'; E0) sigma(E') dE' / int K(E'; E0) dE',  K a Gaussian of width sigma_b(E0) about E0,
 * normalised by its own integral over the nodes' window, so a constant cross section is reproduced exactly however
 * the window is truncated.  Without a spread (sigma <= 0) the value at E0, interpolated on the nodes.
 */
template <class V>
V FoldGaussian(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes, const V *values) {
  const int numPoints = nodes.n;
  if (numPoints < 2) return V(0.0);
  const double eHigh = nodes.High();
  const double eLow = nodes.Low();
  const double centroid = point.energy;
  const double s = p.sigma(centroid, point.cmFactor);
  if (!(s > 0.)) return InterpolateNodes<V>(nodes, values, centroid);
  auto k = [&](double e) { return GaussianDensity(e, centroid, s); };
  V num = 0.0;
  double den = 0.0;
  GridQuadrature<V>(nodes, values, eLow, eHigh, false, k, k, num, den);
  return (den > 0.0) ? num / den : V(0.0);
}

/// convolutionRange sigma on each side (a THM (HOES) point: thmConvolutionRange sigma).
void GaussianSupport(const KernelParameters &p, const SupportPoint &s, double &high, double &low) {
  const double reach = (s.thm ? kThmConvolutionRange : kConvolutionRange) * p.sigma(s.energy, s.cmFactor);
  high = s.energy + reach;
  low = s.energy - reach;
  if (low < s.minimumEnergy) low = s.minimumEnergy;
}

}  // namespace

Kernel GaussianKernel() {
  Kernel k;
  k.name = "gaussian";
  k.keyword = "beam.spread";
  k.fold = FoldGaussian<double>;
  k.weights = FoldGaussian<LinearForm>;
  k.support = GaussianSupport;
  k.doc = "a Gaussian beam-energy spread, sigma constant or a formula of the energy, normalised over its window";
  return k;
}

}  // namespace experiment
}  // namespace azure2
