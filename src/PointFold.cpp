// The apparatus fold of a point (include/PointFold.h): what the kernels read is gathered from the point, its
// sub-points and the working copy, and the effect's kernels are called once.  The kernels themselves are in
// experiment/ (Kernels.h); nothing here depends on which kernel the effect has.
#include "PointFold.h"

#include <cstddef>

#include "Config.h"
#include "DataSet.h"
#include "Results.h"
#include "ShiftedFold.h"
#include "TargetEffect.h"
#include "Unit.h"
#include "WorkingCopy.h"

using azure2::experiment::Broadening;
using azure2::experiment::LinearForm;

namespace {

/// One point's view for its kernels: the effect's parameters, the point, and its sub-points as nodes: all of them,
/// or for a point with an energy shift on the shared mesh those of its slots [first, first + nodes.n) that enclose
/// its moved support (ShiftedFold.h), with the stopping translated with the shift.
struct PointKernelInput {
  const TargetEffect *effect;
  azure2::experiment::KernelParameters parameters;
  azure2::experiment::KernelPoint point;
  std::vector<double> energy, stopping;
  azure2::experiment::KernelNodes nodes;
  int first = 0;

  PointKernelInput(const Unit &unit, const WorkingCopy &data, const Config &configure)
      : effect(data.Data().GetTargetEffect(unit.GetTargetEffectNum())), parameters(effect->KernelParametersFor(configure)) {
    const TargetEffectValues &values = data.EffectValues(unit.GetTargetEffectNum());
    if (!values.kernelValues.empty()) parameters.values = values.kernelValues.data();  // a fitted kernel parameter's value
    const azure2::experiment::Kernel *resolution = effect->Kernels().resolution;
    const bool target = resolution && resolution->usesTarget;
    point.energy = unit.GetCMEnergy();
    point.labEnergy = unit.GetLabEnergy();
    point.nominalLabEnergy = unit.GetOriginalEnergy();
    point.cmFactor = unit.GetCMFactor();
    point.z1 = unit.GetBeamCharge();
    point.z2 = unit.GetTargetCharge();
    point.targetThickness = target ? unit.GetTargetThickness() : 0.0;
    point.density = values.density;
    point.binLow = unit.GetBinLowCM();
    point.binHigh = unit.GetBinHighCM();
    point.photoQ = unit.GetPhotoQ();
    point.photoMass = unit.GetPhotoCompoundMass();
    if (const DataSet::PointSetUp *setUp = data.Data().SetUpOf(unit)) {
      point.udr = setUp->udr.get();
      point.udrMass = setUp->udrMass;
    }
    point.temperature = values.dopplerTemperature;
    const Unit::SubGridSpec &grid = data.Data().GridSpecOf(unit);
    const bool shifted = azure2::IsShiftedOnMesh(unit, data);
    const azure2::ShiftedNodes moved = shifted ? azure2::ShiftedFoldNodes(unit, data) : azure2::ShiftedNodes();
    first = shifted ? moved.first : 0;
    const int n = shifted ? moved.count : unit.NumSubPoints();
    energy.resize(n > 0 ? n : 0);
    for (int j = 1; j <= n; j++) energy[j - 1] = unit.GetSubPoint(data, first + j)->GetCMEnergy();
    if (target && shifted) {
      azure2::TranslatedStopping(unit, data, configure, moved, stopping);
    } else if (target) {
      stopping.resize(energy.size());
      for (int j = 1; j <= n; j++) stopping[j - 1] = unit.GetSubPoint(data, j)->GetStoppingPower();
    }
    nodes.n = n > 0 ? n : 0;
    nodes.energy = energy.data();
    nodes.stopping = target ? stopping.data() : nullptr;
    // On the shared mesh the outer nodes lie beyond the support: the kernel integrates over the support alone, moved
    // with the point's shift.
    if (shifted) {
      nodes.supportLow = moved.low;
      nodes.supportHigh = moved.high;
    } else if (grid.mesh >= 0) {
      nodes.supportLow = grid.endEnergy;
      nodes.supportHigh = grid.startEnergy;
    }
  }
};

}  // namespace

double FoldValues(const Unit &unit, const WorkingCopy &data, const Config &configure, const std::vector<double> &values,
                  Broadening broadening) {
  const PointKernelInput in(unit, data, configure);
  return azure2::experiment::Fold(in.effect->Kernels(), broadening, in.parameters, in.point, in.nodes,
                                  values.data() + in.first);
}

void FoldPoint(Unit &unit, const WorkingCopy &data, const Config &configure, azure2::Results &results) {
  const int n = unit.NumSubPoints();
  std::vector<double> values(n > 0 ? n : 0);
  for (int j = 1; j <= n; j++) values[j - 1] = unit.SlotValue(*unit.GetSubPoint(data, j), j, results);
  unit.SetFitCrossSection(results, FoldValues(unit, data, configure, values, Broadening::kApply));
}

bool FoldWeights(const Unit &unit, const WorkingCopy &data, const Config &configure, std::vector<double> &c,
                 Broadening broadening) {
  const PointKernelInput in(unit, data, configure);
  const int n = in.nodes.n;
  c.assign((std::size_t)unit.NumSubPoints(), 0.0);
  std::vector<LinearForm> forms;
  forms.reserve(n);
  for (int j = 1; j <= n; j++) forms.push_back(LinearForm::Unit(j));
  const LinearForm y = azure2::experiment::FoldForm(in.effect->Kernels(), broadening, in.parameters, in.point, in.nodes,
                                                    forms.data());
  if (!y.Valid()) return false;
  for (const auto &t : y.Terms()) {
    if (t.first < 1 || t.first > n) return false;
    c[(std::size_t)(in.first + t.first - 1)] += t.second;
  }
  return true;
}

double FoldTemperatureDerivative(const Unit &unit, const WorkingCopy &data, const Config &configure,
                                 const std::vector<double> &values) {
  const PointKernelInput in(unit, data, configure);
  if (static_cast<int>(values.size()) < in.first + in.nodes.n) return 0.0;
  return azure2::experiment::FoldTemperatureDerivative(in.effect->Kernels(), Broadening::kApply, in.parameters,
                                                       in.point, in.nodes, values.data() + in.first);
}

void FoldBroadeningWeights(const Unit &unit, const WorkingCopy &data, const Config &configure, double eCM,
                           std::vector<std::pair<int, double>> &w, std::vector<std::pair<int, double>> *dwdT) {
  const PointKernelInput in(unit, data, configure);
  azure2::experiment::BroadeningWeights(in.effect->Kernels(), in.parameters, in.point, in.nodes, eCM, w, dwdT);
  for (auto &x : w) x.first += in.first;  // node numbers of the point's slots
  if (dwdT)
    for (auto &x : *dwdT) x.first += in.first;
}
