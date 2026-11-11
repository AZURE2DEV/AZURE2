// The back of a target layer on a point's sub-point grid: the energy E_out at which the beam leaves a layer of areal
// density N, found on the same nodes and with the same quadrature as the yield it bounds (the target kernel,
// KernelTarget.cpp).
#ifndef AZURE2_EXPERIMENT_LAYEREXIT_H
#define AZURE2_EXPERIMENT_LAYEREXIT_H

#include <algorithm>
#include <cmath>

namespace azure2 {
namespace experiment {

/*!
 * E_out of a layer (docs/upgrade/energy_grid.md section 5.2): the energy, in [E_n, top], with
 *
 *   sum over the sub-point intervals [E_{j+1}, E_j] of the two-point Gauss-Legendre rule of 1 / (eps_lin(E) 1e24),
 *   each interval clipped to [E_out, top]  =  N 1e-24,
 *
 * eps_lin the stopping interpolated linearly between the interval's two nodes.  This is, term by term and in the same
 * order of operations, the weight GridQuadrature (KernelQuadrature.h) puts on the cross section of a target
 * integration, so a constant cross section sigma_0 folded over [E_out, top] is 1e-24 sigma_0 N to rounding: the
 * layer's thickness and its yield come from one epsilon on one set of nodes (an E_out from another integration of
 * epsilon, e.g. TargetEffect::TargetThickness, differs by that integration's quadrature error).
 *
 * Units: energies and stopping in the frame the sub-points carry them (c.m. MeV and c.m. MeV cm^2 per active nucleus,
 * the lab stopping times the point's c.m./lab factor; dE_cm / eps_cm = dE_lab / eps_lab, so the sum is the lab areal
 * density), \p density in atoms/cm^2 (active nuclei).  \p energy(j) and \p stopping(j), j = 1 .. n, are the nodes in
 * descending energy.  The stopping-scale nuisance of the stopping module will multiply eps here and in GridQuadrature
 * together (energy_grid.md section 6.2), never one of them alone.
 *
 * Returns false when the grid ends before the layer does (the areal density between E_n and top is less than N), or
 * when a node's stopping is not positive and finite; \p eOut is then E_n.  Inside the interval that holds it, E_out is
 * found by bisection to the last bit of the energy (the clipped rule is monotone in its lower end for a positive
 * epsilon).
 */
template <class Energy, class Stopping>
bool LayerExitEnergy(int n, Energy energy, Stopping stopping, double top, double density, double &eOut) {
  const double x1 = -1.0 / std::sqrt(3.0), x2 = 1.0 / std::sqrt(3.0);
  const double target = density * 1.0e-24;  // atoms per barn: the unit of the weights below
  eOut = n > 0 ? energy(n) : top;
  if (!(target > 0.0)) return false;
  // The clipped rule on interval j (nodes Ea > Eb) over [b, a], as GridQuadrature evaluates it.
  auto clipped = [&](int j, double Ea, double Eb, double a, double b) {
    const double mid = 0.5 * (a + b), half = 0.5 * (a - b);
    if (half < 1.0e-12) return 0.0;
    const double sa = stopping(j), sb = stopping(j + 1);
    double sum = 0.0;
    for (int gp = 0; gp < 2; gp++) {
      const double e = mid + half * (gp == 0 ? x1 : x2);
      const double frac = (Ea - e) / (Ea - Eb);
      const double s = (1.0 - frac) * sa + frac * sb;
      sum += half * (1.0 / s / 1.0e24);
    }
    return sum;
  };
  double accumulated = 0.0;
  for (int j = 1; j < n; j++) {
    const double Ea = energy(j), Eb = energy(j + 1);
    if (!(stopping(j) > 0.0) || !(stopping(j + 1) > 0.0) || !std::isfinite(stopping(j)) || !std::isfinite(stopping(j + 1))) return false;
    if (Eb >= top) continue;
    const double a = std::min(Ea, top);
    const double whole = clipped(j, Ea, Eb, a, Eb);
    if (accumulated + whole < target) {
      accumulated += whole;
      continue;
    }
    // E_out lies in [Eb, a]: the lower end b with accumulated + clipped(b) = target
    const double rest = target - accumulated;
    double lo = Eb, hi = a;  // clipped(lo) >= rest, clipped(hi) = 0 < rest
    for (int it = 0; it < 200; it++) {
      const double m = 0.5 * (lo + hi);
      if (!(m > lo && m < hi)) break;
      (clipped(j, Ea, Eb, a, m) >= rest ? lo : hi) = m;
    }
    eOut = lo;
    return true;
  }
  return false;
}

}  // namespace experiment
}  // namespace azure2

#endif
