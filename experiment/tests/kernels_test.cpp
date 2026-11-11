// The table of apparatus kernels (azure2/experiment/Kernels.h) on nodes laid out by hand: every row is complete, the
// weights the gradient reads (the kernel on LinearForm) reproduce the forward value, and each kernel against a closed
// form (a constant survives every normalised kernel; a thin target of constant stopping averages a linear cross
// section to its value at the mid-depth; a thick target integrates a constant to 1e-24 sigma dE/epsilon).
#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LinearForm.h"

using namespace azure2::experiment;

namespace {

int failures = 0;

void Check(const std::string &what, bool ok) {
  std::printf("%s %s\n", ok ? "  ok  " : "  FAIL", what.c_str());
  if (!ok) failures++;
}

void Close(const std::string &what, double got, double want, double tolerance) {
  const double scale = std::max(1.0, std::fabs(want));
  const bool ok = std::fabs(got - want) <= tolerance * scale;
  std::printf("%s %s: %.17g (want %.17g)\n", ok ? "  ok  " : "  FAIL", what.c_str(), got, want);
  if (!ok) failures++;
}

/// Nodes from hi down to lo, n intervals, with a stopping per node.
struct Grid {
  std::vector<double> energy, stopping;
  KernelNodes nodes;
  Grid(double hi, double lo, int n, double eps) {
    for (int i = 0; i <= n; i++) energy.push_back(hi - (hi - lo) * i / n);
    stopping.assign(energy.size(), eps);
    nodes.n = static_cast<int>(energy.size());
    nodes.energy = energy.data();
    nodes.stopping = stopping.data();
  }
};

double ConstantSigma(const void *owner, const void *, double, double cmFactor) {
  return cmFactor * *static_cast<const double *>(owner);
}

/// The fold of `values` and the dot product of its LinearForm weights with them.
void FoldBoth(const Kernel &k, const KernelParameters &p, const KernelPoint &point, const Grid &g,
              const std::vector<double> &values, double &value, double &fromWeights) {
  EffectKernels ek;
  ek.resolution = &k;
  value = Fold(ek, Broadening::kApply, p, point, g.nodes, values.data());
  std::vector<LinearForm> unit;
  for (int j = 1; j <= g.nodes.n; j++) unit.push_back(LinearForm::Unit(j));
  const LinearForm form = FoldForm(ek, Broadening::kApply, p, point, g.nodes, unit.data());
  fromWeights = 0.0;
  for (const auto &t : form.Terms()) fromWeights += t.second * values[t.first - 1];
}

}  // namespace

int main() {
  std::printf("the table\n");
  std::set<std::string> names;
  for (const Kernel &k : Kernels()) {
    Check(k.name + ": a unique name", names.insert(k.name).second);
    Check(k.name + ": a format-3 keyword", !k.keyword.empty());
    Check(k.name + ": a support", k.support != nullptr);
    if (k.stage == Kernel::Stage::kResolution)
      Check(k.name + ": fold and weights", k.fold && k.weights && !k.broaden);
    else
      Check(k.name + ": the broadening on both value types, d/dT and its weights",
            k.broaden && k.broadenForm && k.broadenDT && k.broadeningWeights && !k.fold);
    Check(k.name + ": KernelNamed finds it", KernelNamed(k.name) == &k);
    if (k.tableForm) {
      const std::string field = k.keyword.substr(0, k.keyword.find(' '));
      const std::string form = k.keyword.substr(k.keyword.find(' ') + 1);
      Check(k.name + ": KernelOfForm finds it", KernelOfForm(field.substr(0, field.find('.')),
                                                             field.substr(field.find('.') + 1), form) == &k);
      for (const KernelParameter &p : k.parameters)
        Check(k.name + "." + p.name + ": a dimension and a unit", !p.dimension.empty() && !p.unit.empty());
    }
  }
  Check("no kernel named 'none'", KernelNamed("none") == nullptr);

  const double sigmaLab = 0.004;  // MeV
  KernelParameters p;
  p.sigma.at = ConstantSigma;
  p.sigma.owner = &sigmaLab;
  KernelPoint point;
  point.energy = 1.0;
  point.labEnergy = 1.2;
  point.nominalLabEnergy = 1.2;
  point.cmFactor = 1.0 / 1.2;

  std::printf("Gaussian beam spread\n");
  {
    const Kernel &k = *KernelNamed("gaussian");
    p.gaussian = true;
    Grid g(1.0 + 5.0 * sigmaLab, 1.0 - 5.0 * sigmaLab, 40, 1.0);
    double v, w;
    FoldBoth(k, p, point, g, std::vector<double>(g.energy.size(), 3.0), v, w);
    Close("a constant survives", v, 3.0, 1e-14);
    Close("the weights give the same value", w, v, 1e-14);
    std::vector<double> linear;
    for (double e : g.energy) linear.push_back(2.0 + 5.0 * (e - 1.0));
    FoldBoth(k, p, point, g, linear, v, w);
    Close("a linear cross section about the centre is its value there", v, 2.0, 1e-13);
    Close("the weights give the same value", w, v, 1e-14);
  }

  std::printf("thin target, constant stopping\n");
  {
    const Kernel &k = *KernelNamed("target");
    p.gaussian = false;
    const double eps = 1e-20;      // c.m. MeV cm^2 (10 eV per 1e15 atoms/cm^2)
    const double density = 1e18;   // atoms/cm^2
    point.density = density;
    point.targetThickness = eps * density;  // 0.01 MeV
    Grid g(1.0, 1.0 - point.targetThickness, 20, eps);
    std::vector<double> linear;
    for (double e : g.energy) linear.push_back(2.0 + 5.0 * (e - 1.0));
    double v, w;
    FoldBoth(k, p, point, g, linear, v, w);
    Close("the depth average of a linear cross section is its value at mid-depth", v, 2.0 - 5.0 * 0.005, 1e-12);
    Close("the weights give the same value", w, v, 1e-14);
    p.thickTarget = true;
    FoldBoth(k, p, point, g, std::vector<double>(g.energy.size(), 3.0), v, w);
    Close("thick: 1e-24 sigma dE / epsilon", v, 1e-24 * 3.0 * point.targetThickness / eps, 1e-12);
    p.thickTarget = false;
  }

  std::printf("support\n");
  {
    SupportPoint s;
    s.energy = 1.0;
    s.cmFactor = 1.0;
    s.minimumEnergy = 0.001;
    double hi, lo;
    EffectKernels ek;
    ek.resolution = KernelNamed("gaussian");
    p.gaussian = true;
    EffectSupport(ek, p, s, hi, lo);
    Close("Gaussian: E0 + 5 sigma", hi, 1.0 + 5.0 * sigmaLab, 1e-15);
    Close("Gaussian: E0 - 5 sigma", lo, 1.0 - 5.0 * sigmaLab, 1e-15);
    ek.resolution = nullptr;
    EffectSupport(ek, p, s, hi, lo);
    Check("no kernel: the point itself", hi == 1.0 && lo == 1.0);
  }

  std::printf("%s\n", failures ? "FAILED" : "all passed");
  return failures ? 1 : 0;
}
