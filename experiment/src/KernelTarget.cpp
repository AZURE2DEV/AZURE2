// The target kernel (Kernels.h): the integration of a point's node values over the depth of its target (thin, thick,
// layer yield), with the beam spread and the straggling as a Gaussian average at each depth.
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "azure2/experiment/KernelQuadrature.h"
#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LayerExit.h"
#include "azure2/experiment/LinearForm.h"
#include "azure2/physics/Constants.h"

namespace azure2 {
namespace experiment {
namespace {
constexpr double kPi = 3.14159265358979323846;
}

double GaussianDensity(double e, double centroid, double sigma) {
  const double pi = kPi;
  return pow(2. * pi, -0.5) / sigma * exp(-pow(e - centroid, 2.0) / 2.0 / pow(sigma, 2.0));
}

double StragglingSigma(const KernelParameters &p, double deltaE) {
  if (!p.straggling || deltaE <= 0.0) return 0.0;
  return p.stragglingCoefficient * std::sqrt(deltaE * 1000.0) / 1000.0;
}

/*!
 * Bohr's energy straggling after traversing an areal density N x (atoms/cm^2) of the target:
 * Omega_lab^2 = 4 pi Z1^2 e^4 Z2 N x (N. Bohr, Phil. Mag. 30 (1915) 581), a lab width converted to the c.m. with the
 * entrance pair's mass ratio.  e^2 = alpha hbar c in MeV fm, 1 fm = 1e-13 cm.  Z2 is the charge of the target nucleus:
 * for a compound target whose stopping cross section is per molecule or per active atom the electrons of the other
 * constituents are not counted, and an explicit coefficient should be entered instead.
 */
double BohrStragglingSigma(const KernelParameters &p, double arealDensity, double cmFactor, int z1, int z2) {
  if (!p.BohrStraggling() || !(arealDensity > 0.0) || z1 == 0 || z2 <= 0) return 0.0;
  const double pi = kPi;
  const double e2 = physics::kFineStructure * physics::kHbarC_MeV_fm * 1.0e-13;  // MeV cm
  const double omega2 = 4.0 * pi * z1 * z1 * z2 * e2 * e2 * arealDensity;  // lab MeV^2
  return cmFactor * std::sqrt(omega2);
}

namespace {

/*!
 * Target integration:  Y(E0) = (1/N) int_{E0-T}^{E0} <sigma>(E_d) / epsilon(E_d) dE_d  over the depth, where
 * <sigma>(E_d) is sigma itself, or its Gaussian average about E_d of width sqrt(sigma_b(E_d)^2 +
 * sigma_straggling(E0 - E_d)^2) when there is a beam spread and/or straggling, normalised by its own integral over the
 * window (the quadrature averages sigma/epsilon over that Gaussian, which is <sigma>/epsilon(E_d) to the curvature of
 * 1/epsilon across it).  Energies and epsilon are c.m.: dE_cm/eps_cm = dE_lab/eps_lab.  An infinitely thick target
 * (density 0) and a yield on a layer are not divided by N:
 *   Y(E0) = 1e-24 int_{E_out}^{E0} <sigma>(E_d) / epsilon(E_d) dE_d   reactions per incident particle (sigma in b,
 * epsilon in MeV cm^2 per active nucleus), with E_out = 1 keV for the thick target and, for the layer, the energy at
 * which the integral of 1/epsilon on the same nodes reaches the working copy's N (LayerExit.h).
 * A thin target inside one interval of the nodes is no special case: the quadrature integrates the linear interpolant
 * of sigma (and of 1/epsilon) between the two enclosing nodes over [E0 - T, E0], i.e. it gives the interpolated sigma at
 * the depth-averaged energy int E dE/epsilon / int dE/epsilon (exactly so for 1/epsilon constant across the target).
 */
template <class V>
V FoldTarget(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes, const V *values) {
  const int numPoints = nodes.n;
  if (numPoints < 2) return V(0.0);
  const double eHigh = nodes.High();
  const double eLow = nodes.Low();
  const double centroid = point.energy;
  const double range = kConvolutionRange;
  auto one = [](double) { return 1.0; };
  // sigma(E) / (epsilon(E) 1e24) interpolated at one energy, for a depth node with no spread to average over
  auto yieldAt = [&](double e) {
    if (const int j = BracketingInterval(nodes, e)) {
      double Ea = nodes.Energy(j), Eb = nodes.Energy(j + 1);
      double frac = (Ea - e) / (Ea - Eb);
      V sigma = (1.0 - frac) * values[j - 1] + frac * values[j];
      double stopping = (1.0 - frac) * nodes.Stopping(j) + frac * nodes.Stopping(j + 1);
      return (1.0 / stopping / 1.0e24) * sigma;
    }
    return V(0.0);
  };
  // Areal density (atoms/cm^2) of target traversed from the grid's top node down to c.m. energy e, int dE / epsilon
  // over the nodes (trapezoid in 1/epsilon).  For Bohr straggling.  A node with a non-positive stopping power
  // contributes nothing.
  std::vector<double> arealCum;
  auto invEps = [](double eps) { return (eps > 0.0) ? 1.0 / eps : 0.0; };
  if (p.BohrStraggling()) {
    arealCum.assign(numPoints + 1, 0.0);
    for (int j = 1; j < numPoints; j++) {
      double Ea = nodes.Energy(j), Eb = nodes.Energy(j + 1);
      double ia = invEps(nodes.Stopping(j)), ib = invEps(nodes.Stopping(j + 1));
      arealCum[j + 1] = arealCum[j] + 0.5 * (ia + ib) * (Ea - Eb);
    }
  }
  auto arealTo = [&](double e) {
    if (arealCum.empty() || e >= nodes.Energy(1)) return 0.0;
    if (const int j = BracketingInterval(nodes, e)) {
      double Ea = nodes.Energy(j), Eb = nodes.Energy(j + 1);
      double frac = (Ea - e) / (Ea - Eb);
      double ia = invEps(nodes.Stopping(j)), ib = invEps(nodes.Stopping(j + 1));
      double ie = (1.0 - frac) * ia + frac * ib;
      return arealCum[j] + 0.5 * (ia + ie) * (Ea - e);
    }
    return arealCum[numPoints];
  };

  const double surface = centroid;
  double back = centroid - point.targetThickness;
  V outer = 0.0;
  const bool noSpread = !p.gaussian && !p.straggling;
  // A yield on a layer (docs/upgrade/energy_grid.md section 5.2): the depth integral runs from the top of the
  // integration down to E_out, found on these nodes with this quadrature (LayerExit.h), so that
  // int_{E_out}^{top} dE/epsilon is the areal density the working copy holds.
  if (p.layerYield) {
    const double top = noSpread ? eHigh : surface;
    auto energy = [&](int j) { return nodes.Energy(j); };
    auto stopping = [&](int j) { return nodes.Stopping(j); };
    if (!LayerExitEnergy(numPoints, energy, stopping, top, point.density, back))
      return V(std::numeric_limits<double>::quiet_NaN());  // refused at set-up
  }
  if (noSpread) {
    // The grid was built to span the target (the support below), so it is integrated whole: clipping it to
    // [E0 - T, E0] here would drop a point whose segment energy shift moved E0 away from the grid's own centre.
    // A fitted density lays the grid out for the largest density the fit may reach, so the depth integral is cut at
    // the back of the target, E0 - T(N): smooth in N, grid unchanged.
    double unused = 0.0;
    const double lo = (p.densityFit || p.layerYield) ? std::max(eLow, back) : eLow;
    GridQuadrature<V>(nodes, values, lo, eHigh, true, one, one, outer, unused);
  } else {
    const double x1 = -1.0 / sqrt(3.0), x2 = 1.0 / sqrt(3.0);
    // the Gaussian's width at a depth: the beam spread and the straggling down to it
    auto sigmaAt = [&](double eDepth) {
      double beam = p.gaussian ? p.sigma(eDepth, point.cmFactor) : 0.0;
      if (!(beam > 0.)) beam = 0.0;
      double straggle = p.BohrStraggling()
                            ? BohrStragglingSigma(p, arealTo(eDepth) - arealTo(surface), point.cmFactor, point.z1, point.z2)
                            : StragglingSigma(p, surface - eDepth);
      return std::sqrt(beam * beam + straggle * straggle);
    };
    // one panel [b, a] of the depth integral: 2-point Gauss-Legendre, the Gaussian average over every node of its
    // window at each of the two depths
    auto panel = [&](double a, double b) {
      double mid = 0.5 * (a + b), half = 0.5 * (a - b);
      if (half < 1.0e-10) return;
      for (int gp = 0; gp < 2; gp++) {
        double eDepth = mid + half * (gp == 0 ? x1 : x2);
        double sigmaEff = sigmaAt(eDepth);
        V inner = 0.0;
        if (sigmaEff > 1.0e-10) {
          double den = 0.0;
          // the quadrature asks for the numerator's and the denominator's weight at each Gauss point: one evaluation
          double atEnergy = std::numeric_limits<double>::quiet_NaN(), density = 0.0;
          auto k = [&](double e) {
            if (e != atEnergy) {
              atEnergy = e;
              density = GaussianDensity(e, eDepth, sigmaEff);
            }
            return density;
          };
          GridQuadrature<V>(nodes, values, std::max(eLow, eDepth - range * sigmaEff), std::min(eHigh, eDepth + range * sigmaEff), true, k, k, inner, den);
          inner = (den > 0.0) ? inner / den : V(0.0);
        } else {
          inner = yieldAt(eDepth);
        }
        outer += half * inner;
      }
    };
    const int depthSteps = p.depthStepsPerSigma;
    if (depthSteps <= 0) {
      // the depth integral on the node intervals inside the target (--fold-depth-steps=0)
      for (int i = 1; i < numPoints; i++) {
        double Ea = nodes.Energy(i), Eb = nodes.Energy(i + 1);
        if (Eb >= surface || Ea <= back) continue;
        panel(std::min(Ea, surface), std::max(Eb, back));
      }
    } else {
      // The depth integral on its own panels (docs/upgrade/energy_grid.md section 3.4).  The Gaussian average is
      // smooth in the depth on the scale of its width, whatever the nodes resolve, so a panel is sigma / depthSteps
      // wide (sigma at the panel's top, the smallest in the panel when the straggling grows with the depth), and it
      // never ends above the next node: where sigma / depthSteps is shorter than the way to the next node (the first
      // depths under pure straggling, sigma -> 0 at the surface) the panels are the node intervals, as above.  Every
      // panel passes at least one node, so there are never more panels than node intervals.  A NaN sigma gives the
      // node interval (and the depth's value, below, the interpolated one, as on the node intervals).
      const double top = std::min(surface, nodes.Energy(1)), bottom = std::max(back, nodes.Energy(numPoints));
      double a = top;
      int i = FirstIntervalBelow(nodes, a);  // Energy(i) >= a > Energy(i + 1)
      while (a > bottom && i < numPoints) {
        const double node = nodes.Energy(i + 1), h = sigmaAt(a) / depthSteps;
        double b = a - h;
        if (!(h > 0.0) || b > node) b = node;
        // No sliver above the back: a remainder the panel would skip (shorter than 2e-10 MeV, below) is joined to
        // this panel, so that the panels tile [back, top] whole (a constant folds to itself).
        if (b < bottom + 2.0e-10) b = bottom;
        panel(a, b);
        a = b;
        while (i < numPoints && nodes.Energy(i + 1) >= a) i++;
      }
    }
  }
  // a thick target (density 0) and a layer yield are the yield per incident particle, 1e-24 int sigma/epsilon dE;
  // a thin target's quantity is the target-averaged cross section, that divided by 1e-24 N
  return (p.thickTarget || p.layerYield) ? outer : outer / (point.density * 1.0e-24);
}

/*!
 * The nodes reach forward above and backward below the point: the target thickness below, padded with a beam spread
 * by a further factor 5 of its convolutionRange sigma on each side (the nested integral), and the straggling's own
 * convolutionRange sigma at the back of the target.
 */
void TargetSupport(const KernelParameters &p, const SupportPoint &s, double &high, double &low) {
  const double beam = p.gaussian ? kConvolutionRange * p.sigma(s.energy, s.cmFactor) : 0.0;
  double backward = s.targetThickness + (p.gaussian ? 5.0 * beam : 0.0);
  const double forward = p.gaussian ? 5.0 * beam : 0.0;
  if (p.BohrStraggling()) {
    // the whole target: its areal density (the largest a fitted density is laid out for), or for a thick target the
    // thickness over the surface stopping power
    double areal = p.layoutDensity;
    if (p.thickTarget) {
      const double eps = s.cmFactor * p.stopping(s.energy / s.cmFactor, s.cmFactor);
      areal = (eps > 0.0) ? s.targetThickness / eps : 0.0;
    }
    backward += kConvolutionRange * BohrStragglingSigma(p, areal, s.cmFactor, s.z1, s.z2);
  } else {
    backward += kConvolutionRange * StragglingSigma(p, s.targetThickness);
  }
  high = s.energy + forward;
  low = s.energy - backward;
  if (low < s.minimumEnergy) low = s.minimumEnergy;
}

/*!
 * The thin limit (Kernel::thinLimit) of a target that loses less than MeshZeroWidth() (1 eV): sigma constant across
 * it.  A thin target's quantity is the target-averaged cross section, weight 1; a layer yield 1e-24 N (reactions per
 * incident particle, sigma in b, N atoms/cm^2); a thick target 1e-24 T / epsilon, T its c.m. thickness down to 1 keV
 * and epsilon the node's c.m. stopping (dE_cm / eps_cm = dE_lab / eps_lab).
 */
double TargetThinLimit(const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes) {
  if (p.layerYield) return 1.0e-24 * point.density;
  if (p.thickTarget) {
    const double eps = nodes.Stopping(1);
    return (eps > 0.0) ? 1.0e-24 * point.targetThickness / eps : std::numeric_limits<double>::quiet_NaN();
  }
  return 1.0;
}

}  // namespace

Kernel TargetKernel() {
  Kernel k;
  k.name = "target";
  k.keyword = "target (with dataset.target)";
  k.usesTarget = true;
  k.fold = FoldTarget<double>;
  k.weights = FoldTarget<LinearForm>;
  k.support = TargetSupport;
  k.thinLimit = TargetThinLimit;
  k.doc = "the depth integral over a thin target (the target-averaged cross section), a thick target or a layer "
          "(the yield per incident particle), with the beam spread and the straggling as a Gaussian at each depth";
  return k;
}

}  // namespace experiment
}  // namespace azure2
