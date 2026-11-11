// The model pass and the residual (ModelPass.h; docs/upgrade/data_model.md steps S4 and S10).  A segment's model value
// is its expression tree (Predictions.h, experiment/Expression.h) evaluated on its terms' leaf values; every node
// computes in the order the hand-written combination before S10 did, so every model value keeps its bits.  It goes to
// Results::model, never into a leaf slot.
#include "ModelPass.h"

#include <algorithm>
#include <atomic>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "CNuc.h"
#include "Config.h"
#include "CoulombReference.h"
#include "EData.h"
#include "EPoint.h"
#include "ESegment.h"
#include "PPair.h"
#include "Predictions.h"
#include "TargetEffect.h"

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

/*!
 * Calculate a point; if the calculation throws (a singular matrix, a Coulomb function out of range, ...), mark it
 * failed (EPoint::MarkFailed: NaN, reported once) and return false, so that neither the point nor chi^2 can pass for
 * a result.  Before 2026-11-01 a failed base point gave 0, a failed or non-finite component was left out of a SUM and
 * a failed RATIO denominator gave 0, all without a word: the natural SUM segments of tests/natural_neon on the
 * R-matrix route came back 0 at 84 points.
 */
bool CalculatePoint(EData &data, EPoint *point, const char *role, CNuc *cnuc, const Config &configure,
                    azure2::Results &results) {
  try {
    point->Calculate(&data, cnuc, configure, results);
    return true;
  } catch (const std::exception &e) {
    point->MarkFailed(&data, results, role, configure, e.what());
  } catch (...) {
    point->MarkFailed(&data, results, role, configure, "unknown exception.");
  }
  return false;
}

// Whether the combination calculates `point` itself: the leaves are not in yet, and the point is its own host.
bool CalculatesItself(const EPoint *point, azure2::LeafValues leaves) {
  return leaves == azure2::LeafValues::Calculate && (!point->IsMapped() || ScopedPointMappingOff::Active());
}

// A sample segment's point whose combination fails: its per-nucleus terms are NaN too (written by this evaluation).
double FailedSample(const ESegment &segment, EPoint *basePoint, CNuc *cnuc, azure2::Results &results) {
  if (segment.IsNatural()) basePoint->SetNaturalTerms(results, std::vector<double>(cnuc->NumNuclei(), kNaN));
  return kNaN;
}


using azure2::experiment::FactorRef;
using azure2::experiment::LeafPart;
using azure2::experiment::ParamRef;
using azure2::experiment::Prediction;

/// Up to N values on the stack, more on the heap: the inputs of a tree at one point, so that an evaluation of a tree
/// of up to N leaves, factors and parameters allocates nothing.
template <class T, std::size_t N>
class LocalArray {
 public:
  explicit LocalArray(std::size_t n) : n_(n) {
    if (n > N) heap_.resize(n);
  }
  T *data() { return n_ > N ? heap_.data() : fixed_; }
  T &operator[](std::size_t k) { return data()[k]; }

 private:
  std::size_t n_;
  T fixed_[N];
  std::vector<T> heap_;
};

/// The point of term t (0: the segment's own, k: its k-th leaf's) at pointIndex; null if it has none.
EPoint *TermPoint(EData &data, const ESegment &segment, int term, int pointIndex) {
  if (term < 0 || term >= segment.NumTerms() || pointIndex >= segment.NumTermPoints(term)) return nullptr;
  return segment.GetTermPoint(data, term, pointIndex + 1);
}

double LeafValue(const EPoint *p, LeafPart part, const azure2::Results &results) {
  switch (part) {
    case LeafPart::E1: return p->GetFitE1CrossSection(results);
    case LeafPart::E2: return p->GetFitE2CrossSection(results);
    case LeafPart::Value: break;
  }
  return p->GetFitCrossSection(results);
}

/// What a term is called in a message when its calculation fails.
const char *TermRole(const Prediction &p, const ESegment &segment, int term) {
  if (term == 0) return "base term";
  if (p.kind == azure2::experiment::TreeKind::Ratio) return "RATIO denominator";
  return segment.IsYieldTotalLeaf(term) ? "yield total" : "SUM component";
}

/*!
 * The base point's kernel for a transmission or yield (experiment::SubPointFold): the layout of its sub-points, with
 * the free-gas Doppler broadening of the point's target effect as a linear map before the transform (SAMMY's order:
 * every term's sub-point values broadened on the common grid, the broadened values transformed, the kernel applied
 * without its own Doppler step), and EPoint::IntegrateSubPointValues on a local vector, never on the sub-points.
 * For the backward pass of a fitted temperature, the map's derivative d/dT goes beside it, its cotangent into
 * parameter slot `temperatureSlot`.  EData::CheckTransformedSegments makes sure every term shares the base point's
 * sub-point grid.
 */
class BaseKernel : public azure2::experiment::SubPointFold {
 public:
  BaseKernel(const EData &data, EPoint *base, const Config &configure, int temperatureSlot)
      : data_(data), base_(base), configure_(configure) {
    if (!ESegment::HasSubPointIntegration(data, base)) return;
    stage_.subPoints = true;
    stage_.nsp = static_cast<std::size_t>(base->NumSubPoints());
    const TargetEffect *effect = data.Data().GetTargetEffect(base->GetTargetEffectNum());
    if (!effect || !effect->IsDoppler()) return;
    const bool temperature = temperatureSlot >= 0 && data.Data().DopplerParamOfEffect(base->GetTargetEffectNum()) > 0;
    // one row at the point (no kernel after the Doppler step) or one per sub-point (the kernel's grid)
    std::vector<double> energies;
    stage_.fold = effect->HasResolutionStage();
    if (!stage_.fold)
      energies.push_back(base->GetCMEnergy());
    else
      for (int j = 1; j <= base->NumSubPoints(); j++) energies.push_back(base->GetSubPoint(data, j)->GetCMEnergy());
    std::vector<std::pair<int, double>> dw, dwdT;
    rowStart_.push_back(0);
    rowStartT_.push_back(0);
    for (double e : energies) {
      base->DopplerWeights(data, e, dw, temperature ? &dwdT : nullptr);
      for (const auto &t : dw) {
        col_.push_back(t.first - 1);
        w_.push_back(t.second);
      }
      rowStart_.push_back(col_.size());
      for (const auto &t : dwdT) {
        colT_.push_back(t.first - 1);
        wT_.push_back(t.second);
      }
      rowStartT_.push_back(colT_.size());
    }
    stage_.premap = {energies.size(), rowStart_.data(), col_.data(), w_.data()};
    if (temperature) {
      stage_.premapDT = {energies.size(), rowStartT_.data(), colT_.data(), wT_.data()};
      stage_.dopplerParam = temperatureSlot;
    }
  }
  azure2::experiment::TransformStage stage(std::size_t) const override { return stage_; }
  double fold(std::size_t, const double *f, std::size_t n) const override {
    EPoint::DopplerAppliedScope dopplerDone(base_);
    return base_->IntegrateSubPointValues(data_, configure_, std::vector<double>(f, f + n));
  }
  bool weights(std::size_t, double *c, std::size_t n) const override {
    std::vector<double> weights;
    EPoint::DopplerAppliedScope dopplerDone(base_);
    if (!base_->TargetEffectWeights(data_, configure_, weights) || weights.size() != n) return false;
    std::copy(weights.begin(), weights.end(), c);
    return true;
  }

 private:
  const EData &data_;
  EPoint *base_;
  const Config &configure_;
  azure2::experiment::TransformStage stage_;
  std::vector<std::size_t> rowStart_, rowStartT_;
  std::vector<int> col_, colT_;
  std::vector<double> w_, wT_;
};

/*!
 * The inputs of a segment's tree at one point, as the model pass and its adjoint read them: the leaf values from the
 * terms' points' records (the E1 / E2 part where the leaf says so), the frame factors, the abundances and thickness at
 * the point, and for a transmission or yield its terms' sub-point values and the base point's kernel.  On the stack
 * up to 16 leaves, factors and parameters (nothing allocated on the plain, SUM, RATIO and sample paths).  Every term
 * point must exist (ComposeModelValue checks).  One parameter slot more than the tree's is kept for the cotangent of a
 * fitted Doppler temperature.
 */
struct PointInputs {
  PointInputs(EData &data, const ESegment &segment, const Prediction &p, int pointIndex, CNuc *cnuc,
              const Config &configure, const azure2::Results &results)
      : nl(p.leaves.size()), nf(p.factors.size()), np(p.params.size()), point(nl), leafValue(nl), factorValue(nf),
        paramValue(np + 1), leafArray(nl), factorArray(nf) {
    EPoint *basePoint = segment.GetPoint(data, pointIndex + 1);
    for (std::size_t k = 0; k < nl; k++) {
      point[k] = TermPoint(data, segment, p.leaves[k].term, pointIndex);
      const LeafPart part = p.leaves[k].part;
      leafValue[k] = azure2::experiment::IsCoulombReference(part)
                         ? azure2::CoulombReferenceValue(data, point[k], part, cnuc, configure)
                         : LeafValue(point[k], part, results);
      leafArray[k] = {&leafValue[k], 1};
    }
    for (std::size_t f = 0; f < nf; f++) {
      const FactorRef &r = p.factors[f];
      EPoint *termPoint = TermPoint(data, segment, r.term, pointIndex);
      switch (r.kind) {
        case FactorRef::Kind::TermFrame:
          factorValue[f] = basePoint && termPoint ? ESegment::TermFrameFactor(data.State(segment), r.term, pointIndex) : 1.0;
          break;
        case FactorRef::Kind::RatioFrame: factorValue[f] = segment.RatioTermFrameFactor(r.term, termPoint, cnuc); break;
        case FactorRef::Kind::SFactor: factorValue[f] = azure2::SFactorConversionValue(termPoint, cnuc); break;
      }
      factorArray[f] = {&factorValue[f], 1};
    }
    for (std::size_t q = 0; q < np; q++) {
      const ParamRef &r = p.params[q];
      paramValue[q] = (r.kind == ParamRef::Kind::Abundance) ? cnuc->GetAbundance(segment.GetNaturalSample(), r.nucleus)
                                                            : segment.GetSampleThickness(data.State(segment));
    }
    paramValue[np] = 0.0;
    in.leaves = {leafArray.data(), nl};
    in.factors = {factorArray.data(), nf};
    in.params = paramValue.data();
    if (p.tree.Transforms().empty()) return;
    kernel = std::make_unique<BaseKernel>(data, basePoint, configure, static_cast<int>(np));
    folds[0] = kernel.get();
    in.folds = {folds, 1};
    if (!ESegment::HasSubPointIntegration(data, basePoint)) return;
    for (std::size_t k = 0; k < nl; k++) {
      subOffsets.push_back(subValues.size());
      for (int j = 1; j <= point[k]->NumSubPoints(); j++)
        subValues.push_back(LeafValue(point[k]->GetSubPoint(data, j), p.leaves[k].part, results));
      subOffsets.push_back(subValues.size());
    }
    for (std::size_t k = 0; k < nl; k++) subArray.push_back({subValues.data(), &subOffsets[2 * k]});
    in.subLeaves = subArray;
  }

  std::size_t nl, nf, np;
  LocalArray<EPoint *, 16> point;
  LocalArray<double, 16> leafValue, factorValue, paramValue;
  LocalArray<azure2::experiment::PointArray, 16> leafArray, factorArray;
  std::vector<double> subValues;
  std::vector<std::size_t> subOffsets;
  std::vector<azure2::experiment::SubPointArray> subArray;
  std::unique_ptr<BaseKernel> kernel;
  const azure2::experiment::SubPointFold *folds[1] = {nullptr};
  azure2::experiment::Inputs in;
};

const Prediction &PredictionOf(const ESegment &segment) {
  const Prediction *p = segment.GetPrediction();
  if (!p) throw std::logic_error("EvaluateModelPoint: segment " + std::to_string(segment.GetSegmentKey()) +
                                 " has no model expression (EData::BuildPredictions)");
  return *p;
}

/*!
 * The model value of a segment at a point: its expression tree (Predictions.h) evaluated on its terms' values
 * (PointInputs).  The per-nucleus terms of a sample segment go into the base point's leaf record (the decomposition
 * output).
 */
double ComposeModelValue(EData &data, const ESegment &segment, int pointIndex, CNuc *cnuc, const Config &configure,
                         azure2::Results &results, azure2::LeafValues leaves) {
  const Prediction &p = PredictionOf(segment);
  EPoint *basePoint = segment.GetPoint(data, pointIndex + 1);
  // The terms' points, calculated here when the leaves are not in yet (in term order, the base first).
  for (const azure2::experiment::LeafSelector &leaf : p.leaves) {
    EPoint *termPoint = TermPoint(data, segment, leaf.term, pointIndex);
    if (!termPoint) return FailedSample(segment, basePoint, cnuc, results);  // a leaf has a point at every row
    if (azure2::experiment::IsCoulombReference(leaf.part)) continue;  // analytic: nothing to calculate
    if (CalculatesItself(termPoint, leaves) &&
        !CalculatePoint(data, termPoint, TermRole(p, segment, leaf.term), cnuc, configure, results))
      return FailedSample(segment, basePoint, cnuc, results);
  }
  PointInputs inputs(data, segment, p, pointIndex, cnuc, configure, results);
  azure2::experiment::Workspace ws;
  double m;
  if (segment.IsNatural()) {
    std::vector<double> naturalTerms(cnuc->NumNuclei(), 0.0);
    m = p.tree.Value(inputs.in, 0, ws, naturalTerms.data());
    basePoint->SetNaturalTerms(results, naturalTerms);
  } else {
    m = p.tree.Value(inputs.in, 0, ws);
  }
  if (p.kind == azure2::experiment::TreeKind::Ratio && inputs.leafValue[1] == 0.0) {
    // No ratio exists (e.g. the denominator's exit channel is closed at this energy): NaN and a warning, once.
    // Before 2026-11-01 the point was 0, which a fit compared with the datum as if it were a calculated ratio.
    static std::atomic<bool> reported(false);
    if (!reported.exchange(true)) {
#pragma omp critical(azure_failed_point)
      configure.outStream << "WARNING: the RATIO denominator of segment " << basePoint->GetSegmentKey()
                          << " is 0 at E_lab " << basePoint->GetLabEnergy() << " MeV, angle "
                          << basePoint->GetLabAngle() << " deg: the ratio is NaN there (further zero"
                          << " denominators are not reported)." << std::endl;
    }
  }
  return m;
}

}  // namespace

namespace azure2 {

double EvaluateModelPoint(EData &data, const ESegment &segment, int pointIndex, CNuc *compound, const Config &config,
                          Results &results, LeafValues leaves) {
  EPoint *point = (pointIndex >= 0 && pointIndex < segment.NumPoints()) ? segment.GetPoint(data, pointIndex + 1) : nullptr;
  if (!point) throw std::logic_error("EvaluateModelPoint: segment " + std::to_string(segment.GetSegmentKey()) +
                                     " has no point " + std::to_string(pointIndex));
  const double m = ComposeModelValue(data, segment, pointIndex, compound, config, results, leaves);
  point->SetModelValue(results, m);
  point->CopyAngularDistsToModel(results);
  return m;
}

void EvaluateModelSegment(EData &data, const ESegment &segment, CNuc *compound, const Config &config, Results &results,
                          LeafValues leaves) {
  for (int pid = 0; pid < segment.NumPoints(); pid++)
    EvaluateModelPoint(data, segment, pid, compound, config, results, leaves);
}

void EvaluateModel(EData &data, CNuc *compound, const Config &config, Results &results, LeafValues leaves) {
  for (int i = 1; i <= data.Data().NumSegments(); i++)
    EvaluateModelSegment(data, *data.Data().GetSegment(i), compound, config, results, leaves);
}

bool ModelPointBackward(EData &data, const ESegment &segment, int pointIndex, CNuc *compound, const Config &config,
                        const Results &results, double g, ModelCotangents &out) {
  const Prediction &p = PredictionOf(segment);
  for (const experiment::LeafSelector &leaf : p.leaves)
    if (!TermPoint(data, segment, leaf.term, pointIndex)) return false;
  PointInputs inputs(data, segment, p, pointIndex, compound, config, results);
  const std::size_t nl = inputs.nl, np = inputs.np;
  out.leaves.assign(nl, ModelCotangents::Leaf());
  out.params.assign(np + 1, 0.0);
  out.subPoints = !inputs.subArray.empty();
  out.dopplerOrder = out.subPoints && inputs.kernel->stage(0).premap.rowStart != nullptr;
  out.refs = &p.params;
  // the sub-point cotangents laid out as the inputs' sub-point values (one run per leaf), then handed out per leaf
  std::vector<double> subCotangents(inputs.subValues.size(), 0.0);
  LocalArray<double *, 16> leafCot(nl), subCot(nl);
  for (std::size_t k = 0; k < nl; k++) {
    ModelCotangents::Leaf &l = out.leaves[k];
    l.point = inputs.point[k];
    l.xsComponent = static_cast<int>(p.leaves[k].part);
    leafCot[k] = &l.cotangent;
    subCot[k] = out.subPoints ? subCotangents.data() : nullptr;
  }
  experiment::Cotangents cot;
  cot.leaves = {leafCot.data(), nl};
  cot.subLeaves = {subCot.data(), nl};
  cot.params = out.params.data();
  experiment::Workspace ws;
  const bool ok = p.tree.Backward(inputs.in, 1, &g, cot, ws);
  if (out.subPoints)
    for (std::size_t k = 0; k < nl; k++)
      out.leaves[k].sub.assign(subCotangents.begin() + static_cast<std::ptrdiff_t>(inputs.subOffsets[2 * k]),
                               subCotangents.begin() + static_cast<std::ptrdiff_t>(inputs.subOffsets[2 * k + 1]));
  // A Coulomb leaf depends on no parameter the caller back-propagates (its energy shift is differentiated by
  // finite differences, as every leaf's): it leaves the list, whose entries go through the R-matrix adjoint.
  for (std::size_t k = nl; k-- > 0;)
    if (experiment::IsCoulombReference(p.leaves[k].part))
      out.leaves.erase(out.leaves.begin() + static_cast<std::ptrdiff_t>(k));
  return ok;
}

Residual PointResidual(const EData &data, const ESegment &segment, int i, const Results &results, double norm) {
  return DataResidual(segment.GetPoints(data)[i].GetModelValue(results), ESegment::CMValue(data.State(segment), i), ESegment::CMError(data.State(segment), i), norm);
}

}  // namespace azure2

void EPoint::NoModelRecord(const azure2::Results &r) const {
  throw std::logic_error(
      "EPoint: a model value was read or written on a point without a model record in the buffer it was given "
      "(segment " + std::to_string(segment_key_) + ", E_lab " + std::to_string(lab_energy_) + " MeV, offset " +
      (model_slot_ == kNoRecord ? std::string("none") : std::to_string(model_slot_)) + ", " +
      std::to_string(r.model.size()) + " model values): only the points of a data set's own segments have one "
      "(EData::LayoutResults), and EData::BeginEvaluation sizes the buffer");
}
