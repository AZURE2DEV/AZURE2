// The one quadrature of the kernels (Kernels.h): 2-point Gauss-Legendre on every interval of a point's nodes.
#ifndef AZURE2_EXPERIMENT_KERNELQUADRATURE_H
#define AZURE2_EXPERIMENT_KERNELQUADRATURE_H

#include <algorithm>
#include <cmath>

#include "azure2/experiment/Kernels.h"

namespace azure2 {
namespace experiment {

/// The first interval j in [1, n) (nodes j, j + 1) that is not wholly above `hi`, i.e. the first with
/// Energy(j + 1) < hi; n if there is none.  The energies descend, so the intervals before it all lie above `hi`.
/// Written on the negated test so that a NaN `hi` selects every interval, as a walk over all of them would.
inline int FirstIntervalBelow(const KernelNodes &nodes, double hi) {
  int a = 1, b = nodes.n;
  while (a < b) {
    const int m = a + (b - a) / 2;
    if (nodes.Energy(m + 1) >= hi)
      a = m + 1;
    else
      b = m;
  }
  return a;
}

/// One past the last interval j in [1, n) that is not wholly below `lo`, i.e. the first j with Energy(j) <= lo; n
/// if there is none (a NaN `lo` selects every interval).
inline int EndIntervalAbove(const KernelNodes &nodes, double lo) {
  int a = 1, b = nodes.n;
  while (a < b) {
    const int m = a + (b - a) / 2;
    if (nodes.Energy(m) <= lo)
      b = m;
    else
      a = m + 1;
  }
  return a;
}

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
  // The intervals that meet the window are consecutive (descending energies): found by bisection, so that a narrow
  // window on many nodes (the Gaussian inside the target's depth integral) costs its own intervals, not all of them.
  // The test below stays: it decides, the bisection only skips intervals it would skip.
  const int first = FirstIntervalBelow(nodes, hi), end = EndIntervalAbove(nodes, lo);
  for (int j = first; j < end; j++) {
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
