// The table of apparatus kernels (Kernels.h) and the fold of an effect: its broadening, then its resolution kernel.
// A new kernel is one function in a file of its own (KernelGaussian.cpp is the pattern) and one line below.
#include <vector>

#include "azure2/experiment/Kernels.h"
#include "azure2/experiment/LinearForm.h"

namespace azure2 {
namespace experiment {

// The rows, one per kernel file.
Kernel TargetKernel();
Kernel GaussianKernel();
Kernel BeamProfileKernel();
Kernel UdrKernelRow();
Kernel DopplerKernel();

const std::vector<Kernel> &Kernels() {
  static const std::vector<Kernel> table = {
      TargetKernel(),
      GaussianKernel(),
      BeamProfileKernel(),
      UdrKernelRow(),
      DopplerKernel(),
  };
  return table;
}

const Kernel *KernelNamed(const std::string &name) {
  for (const Kernel &k : Kernels())
    if (k.name == name) return &k;
  return nullptr;
}

const Kernel *KernelOfForm(const std::string &record, const std::string &field, const std::string &form) {
  const std::string keyword = record + "." + field + " " + form;
  for (const Kernel &k : Kernels())
    if (k.tableForm && k.keyword == keyword) return &k;
  return nullptr;
}

namespace {

template <class V>
V FoldOf(const Kernel &k, const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes, const V *values);
template <>
double FoldOf(const Kernel &k, const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes,
              const double *values) {
  return k.fold(p, point, nodes, values);
}
template <>
LinearForm FoldOf(const Kernel &k, const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes,
                  const LinearForm *values) {
  return k.weights(p, point, nodes, values);
}
template <class V>
void BroadenOf(const Kernel &k, const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes,
               const V *values, const double *e, int m, V *out);
template <>
void BroadenOf(const Kernel &k, const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes,
               const double *values, const double *e, int m, double *out) {
  k.broaden(p, point, nodes, values, e, m, out);
}
template <>
void BroadenOf(const Kernel &k, const KernelParameters &p, const KernelPoint &point, const KernelNodes &nodes,
               const LinearForm *values, const double *e, int m, LinearForm *out) {
  k.broadenForm(p, point, nodes, values, e, m, out);
}

/// The broadening (SAMMY's order: first) at every node, or at the point without a resolution kernel; then the
/// resolution kernel on the broadened values.  `broaden` is the broadening of V (its value or, d/dT, its derivative).
template <class V, class Broaden>
V Compose(const EffectKernels &k, Broadening broadening, const KernelParameters &p, const KernelPoint &point,
          const KernelNodes &nodes, const V *values, Broaden broaden) {
  if (!k.broadening || broadening == Broadening::kAlreadyApplied) {
    if (!k.resolution) return V(0.0);
    return FoldOf<V>(*k.resolution, p, point, nodes, values);
  }
  const int n = nodes.n;
  if (n < 2) return V(0.0);
  if (!k.resolution) {
    V at = 0.0;
    broaden(&point.energy, 1, &at);
    return at;
  }
  std::vector<V> b(n);
  broaden(nodes.energy, n, b.data());
  return FoldOf<V>(*k.resolution, p, point, nodes, b.data());
}

}  // namespace

double Fold(const EffectKernels &k, Broadening broadening, const KernelParameters &p, const KernelPoint &point,
            const KernelNodes &nodes, const double *values) {
  return Compose<double>(k, broadening, p, point, nodes, values, [&](const double *e, int m, double *out) {
    BroadenOf<double>(*k.broadening, p, point, nodes, values, e, m, out);
  });
}

LinearForm FoldForm(const EffectKernels &k, Broadening broadening, const KernelParameters &p, const KernelPoint &point,
                    const KernelNodes &nodes, const LinearForm *values) {
  return Compose<LinearForm>(k, broadening, p, point, nodes, values, [&](const double *e, int m, LinearForm *out) {
    BroadenOf<LinearForm>(*k.broadening, p, point, nodes, values, e, m, out);
  });
}

double FoldTemperatureDerivative(const EffectKernels &k, Broadening broadening, const KernelParameters &p,
                                 const KernelPoint &point, const KernelNodes &nodes, const double *values) {
  // The resolution kernel after the broadening is linear in its node values: apply it to d(broadened)/dT.
  if (!k.broadening || broadening == Broadening::kAlreadyApplied) return 0.0;
  return Compose<double>(k, broadening, p, point, nodes, values, [&](const double *e, int m, double *out) {
    k.broadening->broadenDT(p, point, nodes, values, e, m, out);
  });
}

void BroadeningWeights(const EffectKernels &k, const KernelParameters &p, const KernelPoint &point,
                       const KernelNodes &nodes, double e, std::vector<std::pair<int, double>> &w,
                       std::vector<std::pair<int, double>> *dwdT) {
  w.clear();
  if (dwdT) dwdT->clear();
  if (k.broadening) k.broadening->broadeningWeights(p, point, nodes, e, w, dwdT);
}

void EffectSupport(const EffectKernels &k, const KernelParameters &p, const SupportPoint &s, double &high, double &low) {
  if (k.resolution) {
    k.resolution->support(p, s, high, low);
  } else {
    high = s.energy;
    low = s.energy;
    if (low < s.minimumEnergy) low = s.minimumEnergy;
  }
  if (k.broadening) k.broadening->support(p, s, high, low);
}

}  // namespace experiment
}  // namespace azure2
