// The one quadrature of the kernels (Kernels.h): 2-point Gauss-Legendre on every interval of a point's nodes.
#ifndef AZURE2_EXPERIMENT_KERNELQUADRATURE_H
#define AZURE2_EXPERIMENT_KERNELQUADRATURE_H

#include <algorithm>
#include <cmath>

#include "azure2/experiment/Kernels.h"

namespace azure2 {
namespace experiment {

/*!
 * 2-point Gauss-Legendre over the nodes (energies descending from node 1), restricted to the c.m. window [lo, hi]:
 * the node values (values[j - 1]) and, for a yield, the stopping are interpolated linearly between the two nodes that
 * bracket each Gauss point.  Accumulates num += integral of wNum(E) value(E) [/ (stopping(E) 1e24)] and den +=
 * integral of wDen(E), so a kernel average is num / den and a yield is num alone.  Exact for cubics and indifferent
 * to the spacing of an adaptive grid.
 */
template <class V, class WNum, class WDen>
void GridQuadrature(const KernelNodes &nodes, const V *values, double lo, double hi, bool divideByStopping, WNum wNum,
                    WDen wDen, V &num, double &den) {
  const double x1 = -1.0 / sqrt(3.0), x2 = 1.0 / sqrt(3.0);
  const int n = nodes.n;
  for (int j = 1; j < n; j++) {
    const double Ea = nodes.Energy(j);
    const double Eb = nodes.Energy(j + 1);
    if (Eb >= hi || Ea <= lo) continue;
    const double a = std::min(Ea, hi), b = std::max(Eb, lo);
    const double mid = 0.5 * (a + b), half = 0.5 * (a - b);
    if (half < 1.0e-12) continue;
    for (int gp = 0; gp < 2; gp++) {
      const double e = mid + half * (gp == 0 ? x1 : x2);
      const double frac = (Ea - e) / (Ea - Eb);  // 0 at node j, 1 at node j + 1
      V sigma = (1.0 - frac) * values[j - 1] + frac * values[j];
      double f = 1.0;
      if (divideByStopping) {
        double stopping = (1.0 - frac) * nodes.Stopping(j) + frac * nodes.Stopping(j + 1);
        f = 1.0 / stopping / 1.0e24;
      }
      num += (half * wNum(e) * f) * sigma;
      den += half * wDen(e);
    }
  }
}

/// The node values interpolated linearly at energy e inside the nodes' span; 0 outside it.
template <class V>
V InterpolateNodes(const KernelNodes &nodes, const V *values, double e) {
  for (int j = 1; j < nodes.n; j++) {
    double Ea = nodes.Energy(j), Eb = nodes.Energy(j + 1);
    if (e <= Ea && e >= Eb) {
      double frac = (Ea - e) / (Ea - Eb);
      return (1.0 - frac) * values[j - 1] + frac * values[j];
    }
  }
  return V(0.0);
}

}  // namespace experiment
}  // namespace azure2

#endif
