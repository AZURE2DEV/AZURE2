// The parameters of kind `apparatus` (ApparatusParameters.h, docs/upgrade/apparatus_parameters.md): the rows, the
// grouping of the effects' declarations into parameters, their slots in the vector, the hand-over, the priors, the
// labels and the API's listing.  The one engine file a fitted number of an apparatus record passes through.
#include "ApparatusParameters.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <ostream>

#include "AZUREGrad.h"
#include "AZUREParams.h"
#include "Config.h"
#include "WorkingCopy.h"

std::string ApparatusDeclaration::Name() const {
  std::string name;
  for (const auto &part : owner) name += (name.empty() ? "" : ".") + part.second;
  return name;
}

const std::vector<ApparatusRow> &ApparatusRows() {
  static const std::vector<ApparatusRow> rows = {StoppingScaleRow(), KernelParameterRow()};
  return rows;
}

double ApparatusParam::Penalty(double value) const {
  if (!declaration.vary || priorReplaced || !(declaration.sigma > 0.0)) return 0.0;
  return GaussianPriorChi2(value, declaration.centre, declaration.sigma);
}

// ---------------------------------------------------------------- the row of a table kernel's own parameters
namespace {

void DeclareKernelParameters(const TargetEffect &te, std::vector<ApparatusDeclaration> &out) {
  for (const azure2::experiment::ApparatusParameterSpec &a : te.GetKernelParameterFits()) {
    ApparatusDeclaration d;
    d.slot = a.slot;
    d.key = a.key;
    d.owner = a.owner;
    d.unit = a.unit;
    d.layoutBase = te.GetKernelLayoutBase(a.slot);
    d.centre = a.centre.Or(a.value);
    d.vary = a.vary;
    d.sigma = a.priorWidth;
    d.layout = te.GetKernelLayoutValue(a.slot);
    out.push_back(d);
  }
}

// The working copy's values of the kernel: a copy of the effect's, with the fitted slots written over.
bool WriteKernelParameter(const TargetEffect &te, int slot, double v, TargetEffectValues &values) {
  if (values.kernelValues.empty()) values.kernelValues = te.GetKernelValues();
  if (slot < 0 || slot >= (int)values.kernelValues.size() || values.kernelValues[slot] == v) return false;
  values.kernelValues[slot] = v;
  return true;
}

double ReadKernelParameter(const TargetEffect &te, int slot, const TargetEffectValues &values) {
  const std::vector<double> &v = values.kernelValues.empty() ? te.GetKernelValues() : values.kernelValues;
  return (slot >= 0 && slot < (int)v.size()) ? v[slot] : 0.0;
}

}  // namespace

ApparatusRow KernelParameterRow() {
  // The fold reads the working copy's values (PointFold.cpp), so nothing follows a change.
  return {"kernel", "kernel parameter", "Kernel-Parameter", "a parameter of a table kernel of the effect's resolution",
          DeclareKernelParameters, WriteKernelParameter, ReadKernelParameter, nullptr};
}

// ---------------------------------------------------------------- the groups
void WorkingCopy::SetUpApparatusParams(const Config &configure) {
  DataSet &d = Draft();
  d.apparatusParams_.clear();
  for (int row = 0; row < (int)ApparatusRows().size(); row++) {
    std::map<std::string, int> byKey;
    for (int n = 1; n <= d.NumTargetEffects(); n++) {
      std::vector<ApparatusDeclaration> declared;
      ApparatusRows()[row].declare(*d.GetTargetEffect(n), declared);
      for (const ApparatusDeclaration &decl : declared) {
        auto it = byKey.find(decl.key);
        if (it == byKey.end()) {
          ApparatusParam p;
          p.row = row;
          p.declaration = decl;
          p.lead = n;
          d.apparatusParams_.push_back(p);
          it = byKey.emplace(decl.key, (int)d.apparatusParams_.size()).first;
          if (!decl.note.empty()) configure.Note(decl.note);
        }
        d.apparatusParams_[it->second - 1].effects.push_back(n);
      }
    }
  }
}

void AddApparatusParams(WorkingCopy &data, AZUREParams &params) {
  ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  for (int k = 1; k <= data.Data().NumApparatusParams(); k++) {
    const ApparatusParam &a = data.Data().GetApparatusParam(k);
    const ApparatusDeclaration &d = a.declaration;
    const std::string name = d.Name();
    const double v = a.Row().read(*data.Data().GetTargetEffect(a.lead), d.slot, data.EffectValues(a.lead));
    params.Add(name, v, (d.sigma > 0.0) ? 0.1 * d.sigma : 0.01 * v,
               ParameterRecord::Apparatus(a.lead, a.Row().name, d.unit, d.owner));
    if (!d.vary) {
      p.Fix(name);
      continue;
    }
    p.SetLimits(name, 0.1 * d.layoutBase, d.layout);
  }
}

// ---------------------------------------------------------------- the hand-over
void WorkingCopy::FillApparatusFromParams(const vector_r &p, const Config &configure) {
  if (Data().NumApparatusParams() == 0 || set_->apparatusParamOffset_ < 0) return;
  const int numRows = (int)ApparatusRows().size();
  std::vector<std::vector<char>> changed(numRows, std::vector<char>(Data().NumTargetEffects(), 0));
  std::vector<char> anyOfRow(numRows, 0);
  for (int k = 1; k <= Data().NumApparatusParams(); k++) {
    const int idx = set_->apparatusParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the apparatus block
    const ApparatusParam &a = Data().GetApparatusParam(k);
    double v = p[idx];
    // The sub-point grids span the layout value (as the layout density): beyond it the fold would run off its nodes.
    const double vMax = a.declaration.layout, vMin = 0.1 * a.declaration.layoutBase;
    if (v > vMax || !(v >= vMin)) {
      if (configure.FirstTime("apparatus parameter range"))
        configure.outStream << "WARNING: " << a.declaration.Name() << " = " << v << " is outside [" << vMin << ", " << vMax
                            << "], the range its sub-point grids span; the calculation uses the nearest end."
                            << " (Printed once.)" << std::endl;
      v = (v > vMax) ? vMax : vMin;
    }
    for (int n : a.effects)
      if (a.Row().write(*Data().GetTargetEffect(n), a.declaration.slot, v, EffectValues(n)))
        changed[a.row][n - 1] = anyOfRow[a.row] = 1;
  }
  for (int row = 0; row < numRows; row++)
    if (anyOfRow[row] && ApparatusRows()[row].follow) ApparatusRows()[row].follow(*this, changed[row], configure);
}

double WorkingCopy::ApparatusPenalty(int row) const {
  double chi2 = 0.0;
  for (int k = 1; k <= Data().NumApparatusParams(); k++) {
    const ApparatusParam &a = Data().GetApparatusParam(k);
    if (row >= 0 && a.row != row) continue;
    chi2 += a.Penalty(a.Row().read(*Data().GetTargetEffect(a.lead), a.declaration.slot, EffectValues(a.lead)));
  }
  return chi2;
}

// ---------------------------------------------------------------- outputs, labels, the API, the index map
void WriteApparatusChiSquared(const WorkingCopy &data, std::ostream &out) {
  for (int row = 0; row < (int)ApparatusRows().size(); row++) {
    bool anyPrior = false;
    for (int k = 1; k <= data.Data().NumApparatusParams(); k++) {
      const ApparatusParam &a = data.Data().GetApparatusParam(k);
      anyPrior = anyPrior || (a.row == row && a.declaration.vary && a.declaration.sigma > 0.0);
    }
    if (anyPrior) out << " Total-" << ApparatusRows()[row].chiLabel << "-Chi-Squared: " << data.ApparatusPenalty(row);
  }
}

std::string ApparatusLabel(const DataSet &data, int k) {
  const ApparatusParam &a = data.GetApparatusParam(k);
  return std::string(a.Row().label) + " " + a.declaration.Name() + " (" + a.Row().doc + ")";
}

void AppendApparatusParams(ParamIndexMap &map, const WorkingCopy *data) {
  if (!data) return;
  for (int k = 1; k <= data->Data().NumApparatusParams(); k++) map.AppendApparatus(k);
}
