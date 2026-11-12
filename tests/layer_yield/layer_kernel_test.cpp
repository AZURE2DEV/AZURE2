/*!
 * The layer kernel alone (TargetEffect::IsLayerYield, azure2/experiment/LayerExit.h, experiment/src/KernelTarget.cpp), on sub-point values
 * chosen by hand -- no R-matrix, no project (the pattern of tests/reference/beam_profile_reference_test.cpp).  The
 * project of this directory (check.sh) tests the same kernel through the engine against closed forms to the output's
 * 11 digits; here the identity of energy_grid.md section 5.3 is tested to rounding:
 *
 *   1. a constant cross section sigma_0 (b) folds to 1e-24 sigma_0 N (reactions per incident particle) to 1e-14, for
 *      epsilon = a0 x^a1 (lab MeV cm^2, x = E/MeV) on a uniform and on an irregular grid, and for a c.m./lab factor of
 *      0.9375 (sub-point energies and stopping both c.m.); the same point as a thin target is sigma_0 to 1e-5 only
 *      (its E_out is TargetThickness's, another integration of epsilon), which shows the identity is not trivial;
 *   2. the weights of the fold (FoldWeights, the analytic gradient's) add up to 1e-24 N to 1e-14;
 *   3. E_out of a constant epsilon is E_0 - N epsilon to 1e-14 (the rule is exact for it);
 *   4. a cross section linear in E folds to 1e-24 N sigma(E_bar_Q) with E_bar_Q the weighted mean energy of the nodes,
 *      checked against a constant epsilon, where E_bar = E_0 - N epsilon / 2 exactly: 1e-13;
 *   5. a grid that ends before the layer does (E_out below the lowest node) gives NaN, never a truncated yield;
 *   6. a single node (a support narrower than MeshZeroWidth(), a layer losing less than 1 eV) is the thin limit, to
 *      rounding: 1e-24 sigma N for the layer (weight 1e-24 N), sigma for the thin target (weight 1), and
 *      1e-24 sigma T / eps for the infinitely thick target (density 0, T its thickness) -- each was 0 before.
 *
 * Run:  <build>/tests/layer_yield/layer_kernel_test          (ctest: layer_kernel)
 */
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Config.h"
#include "PointFold.h"
#include "azure2/experiment/LayerExit.h"
#include "TargetEffect.h"
#include "Unit.h"
#include "WorkingCopy.h"

// Each consumer of the engine defines this (src/AZURE2.cpp, api/src/pybind_api.cpp); CoulFunc reads it.
Config *g_config = nullptr;

namespace {

int failures = 0;

void check(const std::string &name, double got, double want, double rtol) {
  const double rel = std::fabs(got - want) / std::max(std::fabs(want), 1.0e-300);
  const bool ok = rel <= rtol;
  std::printf("  %s  %s: %.15e vs %.15e (relative %.2e, tolerance %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), got,
              want, rel, rtol);
  if (!ok) failures++;
}

TargetEffect Layer(const Config &configure, double density, const std::string &formula, std::vector<double> params,
                   bool layerYield) {
  azure2::experiment::TargetEffectSpec s;
  s.active = true;
  s.targetName = "L";
  s.numIntegrationPoints = 20;
  s.targetIntegration = true;
  s.density = density;
  s.layerYield = layerYield;
  s.stoppingPowerEquation = formula;
  s.stoppingPowerParameters = std::move(params);
  s.hasStragglingFlag = true;
  s.straggling = false;
  return TargetEffect(s, configure);
}

struct Fold {
  double value = 0.0;
  double weightSum = 0.0;
};

// A point at c.m. energy e0 (lab e0 / f) with sub-points at the given c.m. energies (descending), each with the
// c.m. stopping f eps_lab(E_lab) and the cross section sigma(E): the yield and the sum of its weights.
template <class Sigma>
Fold Run(const Config &configure, const TargetEffect &effect, const std::vector<double> &nodes, double f, Sigma sigma,
         double thickness) {
  WorkingCopy data;
  data.AddTargetEffect(effect);
  const size_t rec = azure2::Results::RecordSize(0, 0);
  azure2::Results results;
  results.slot.assign(rec * (nodes.size() + 1), 0.0);
  const double e0 = nodes.front();
  Unit point(0.0, e0, 1, 2, false, false, false, 0.0, 0, 0);
  point.SetResultRecord(0, 0, 0);
  point.SetTargetEffectNum(1);
  point.SetLabEnergy(e0 / f);
  point.SetCMEnergy(e0);
  point.SetTargetKinematics(f, 1, 7);
  point.SetTargetThickness(thickness);
  const int unit = data.AppendUnit(point);
  const int first = data.NumSlots();
  for (size_t i = 0; i < nodes.size(); i++) {
    Site sub = data.UnitAt(unit).OwnSite();  // the point's settings, at the node
    sub.SetCMEnergy(nodes[i]);
    sub.SetLabEnergy(nodes[i]);
    sub.SetResultRecord(rec * (i + 1), 0, 0);
    sub.SetFitCrossSection(results, sigma(nodes[i]));
    sub.SetStoppingPower(f * effect.StoppingAt(nodes[i] / f, configure, 1.0));
    data.UnitAt(unit).MakeSubPoint(sub);
    data.AppendSlot(sub);
  }
  data.SetSlotRange(unit, first, (int)nodes.size());
  FoldPoint(data.UnitAt(unit), data, configure, results);
  Fold out;
  out.value = data.UnitAt(unit).GetFitCrossSection(results);
  std::vector<double> w;
  if (FoldWeights(data.UnitAt(unit), data, configure, w))
    for (double x : w) out.weightSum += x;
  return out;
}

std::vector<double> Uniform(double hi, double lo, int n) {
  std::vector<double> e;
  for (int i = 0; i <= n; i++) e.push_back(hi - (hi - lo) * i / n);
  return e;
}

}  // namespace

int main() {
  std::ostringstream sink;
  Config configure(sink);
  g_config = &configure;
  std::cout << "layer kernel" << std::endl;

  const double N = 2.0e18, sigma0 = 0.37;
  const double a0 = 4.5e-20, a1 = -0.6;  // eps = a0 x^a1, lab MeV cm^2 per active nucleus
  const double eLab = 0.8;
  for (double f : {1.0, 0.9375}) {
    const TargetEffect layer = Layer(configure, N, "a0*x^a1", {a0, a1}, true);
    const double delta = f * layer.TargetThickness(eLab, configure);  // c.m.
    const double e0 = f * eLab;
    std::vector<double> nodes = Uniform(e0, e0 - 1.02 * delta, 12);
    const char *frame = f == 1.0 ? "f = 1" : "f = 0.9375";
    Fold y = Run(configure, layer, nodes, f, [&](double) { return sigma0; }, delta);
    check(std::string("1. constant sigma, uniform grid, ") + frame + ": Y = 1e-24 sigma_0 N", y.value, 1e-24 * sigma0 * N, 1e-14);
    check(std::string("2. the weights add up to 1e-24 N, ") + frame, y.weightSum, 1e-24 * N, 1e-14);
    // an irregular grid: nodes pulled towards the surface
    std::vector<double> odd;
    for (int i = 0; i <= 15; i++) odd.push_back(e0 - 1.03 * delta * std::pow(i / 15.0, 1.7));
    y = Run(configure, layer, odd, f, [&](double) { return sigma0; }, delta);
    check(std::string("1. constant sigma, irregular grid, ") + frame, y.value, 1e-24 * sigma0 * N, 1e-14);
    // the same point as a thin target: the target-averaged sigma, from TargetThickness's E_out
    const TargetEffect thin = Layer(configure, N, "a0*x^a1", {a0, a1}, false);
    y = Run(configure, thin, Uniform(e0, e0 - delta, 12), f, [&](double) { return sigma0; }, delta);
    const double rel = std::fabs(y.value / sigma0 - 1);
    std::printf("  ..    the thin target on its own grid gives sigma_0 to %.2e (TargetThickness against the quadrature)\n", rel);
    if (!(rel > 1e-12 && rel < 1e-3)) {
      std::printf("  FAIL  the thin target should differ from sigma_0 by the quadrature (1e-12 .. 1e-3)\n");
      failures++;
    }
  }

  // 3, 4: a constant epsilon (the rule is exact for 1/eps constant): E_out = E0 - N eps, E_bar = E0 - N eps / 2
  {
    const double eps = 5.0e-21, e0 = 0.96;
    const TargetEffect layer = Layer(configure, N, "a0", {eps}, true);
    std::vector<double> nodes = Uniform(e0, e0 - 1.02 * N * eps, 9);
    double eOut = 0.0;
    const bool found = azure2::experiment::LayerExitEnergy(
        (int)nodes.size(), [&](int j) { return nodes[j - 1]; }, [&](int) { return eps; }, e0, N, eOut);
    check("3. constant epsilon: E_out = E0 - N eps", found ? eOut : 0.0, e0 - N * eps, 1e-14);
    const double s0 = 0.2, s1 = 3.0;  // sigma = s0 + s1 (E - E0), b
    Fold y = Run(configure, layer, nodes, 1.0, [&](double e) { return s0 + s1 * (e - e0); }, N * eps);
    check("4. linear sigma, constant epsilon: Y = 1e-24 N sigma(E0 - N eps/2)", y.value,
          1e-24 * N * (s0 - s1 * 0.5 * N * eps), 1e-13);
  }

  // 5. a grid that ends inside the layer
  {
    const TargetEffect layer = Layer(configure, N, "a0*x^a1", {a0, a1}, true);
    const double delta = layer.TargetThickness(eLab, configure);
    Fold y = Run(configure, layer, Uniform(eLab, eLab - 0.9 * delta, 12), 1.0, [&](double) { return sigma0; }, delta);
    const bool nan = std::isnan(y.value);
    std::printf("  %s  5. a grid that ends inside the layer gives NaN (got %g)\n", nan ? "ok  " : "FAIL", y.value);
    if (!nan) failures++;
  }

  // 6. a single node: the thin limit of each kind of target
  {
    const double e0 = 0.6, s1 = 0.41;
    const std::vector<double> one = {e0};
    const TargetEffect layer = Layer(configure, 3.0e14, "a0*x^a1", {a0, a1}, true);
    Fold y = Run(configure, layer, one, 1.0, [&](double) { return s1; }, 1.0e-7);
    check("6. one node, layer: Y = 1e-24 sigma N", y.value, 1e-24 * s1 * 3.0e14, 1e-15);
    check("6. one node, layer: the weight is 1e-24 N", y.weightSum, 1e-24 * 3.0e14, 1e-15);
    const TargetEffect thin = Layer(configure, 3.0e14, "a0*x^a1", {a0, a1}, false);
    y = Run(configure, thin, one, 1.0, [&](double) { return s1; }, 1.0e-7);
    check("6. one node, thin target: sigma", y.value, s1, 1e-15);
    check("6. one node, thin target: the weight is 1", y.weightSum, 1.0, 1e-15);
    const TargetEffect thick = Layer(configure, 0.0, "a0*x^a1", {a0, a1}, false);
    const double t = 4.0e-7, eps = thick.StoppingAt(e0, configure, 1.0);
    y = Run(configure, thick, one, 1.0, [&](double) { return s1; }, t);
    check("6. one node, thick target: 1e-24 sigma T / eps", y.value, 1e-24 * s1 * t / eps, 1e-15);
  }

  if (failures) {
    std::printf("FAIL layer_kernel (%d)\n", failures);
    return 1;
  }
  std::printf("PASS layer_kernel\n");
  return 0;
}
