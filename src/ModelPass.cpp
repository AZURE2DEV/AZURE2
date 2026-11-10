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
bool CalculatePoint(EPoint *point, const char *role, CNuc *cnuc, const Config &configure, azure2::Results &results) {
  try {
    point->Calculate(cnuc, configure, results);
    return true;
  } catch (const std::exception &e) {
    point->MarkFailed(results, role, configure, e.what());
  } catch (...) {
    point->MarkFailed(results, role, configure, "unknown exception.");
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

/// The segment of term t: the segment itself (0) or its t-th component segment.
ESegment *TermSegment(ESegment &segment, int term) {
  if (term == 0) return &segment;
  const std::vector<ESegment *> &c = segment.GetComponentSegments();
  return (term - 1 < static_cast<int>(c.size())) ? c[term - 1] : nullptr;
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

/// The point of term t (0: the segment's own, k: its k-th component segment's) at pointIndex; null if it has none.
EPoint *TermPoint(ESegment &segment, int term, int pointIndex) {
  if (term == 0) return segment.GetPoint(pointIndex + 1);
  ESegment *c = TermSegment(segment, term);
  return (c && pointIndex < c->NumPoints()) ? c->GetPoint(pointIndex + 1) : nullptr;
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
  const ESegment *c = segment.GetComponentSegments()[term - 1];
  return c->IsYieldTotalTerm() ? "yield total" : "SUM component";
}

/*!
 * The base point's kernel for a transmission or yield (experiment::SubPointFold): the layout of its sub-points, with
 * the free-gas Doppler broadening of the point's target effect as a linear map before the transform (SAMMY's order:
 * every term's sub-point values broadened on the common grid, the broadened values transformed, the kernel applied
 * without its own Doppler step), and EPoint::IntegrateSubPointValues on a local vector, never on the sub-points.
 * EData::CheckTransformedSegments makes sure every term shares the base point's sub-point grid.
 */
class BaseKernel : public azure2::experiment::SubPointFold {
 public:
  BaseKernel(EPoint *base, const Config &configure) : base_(base), configure_(configure) {
    if (!ESegment::HasSubPointIntegration(base)) return;
    stage_.subPoints = true;
    stage_.nsp = static_cast<std::size_t>(base->NumSubPoints());
    const TargetEffect *effect = base->GetParentData()->GetTargetEffect(base->GetTargetEffectNum());
    if (!effect || !effect->IsDoppler()) return;
    // one row at the point (no kernel after the Doppler step) or one per sub-point (the kernel's grid)
    std::vector<double> energies;
    stage_.fold = effect->HasResolutionStage();
    if (!stage_.fold)
      energies.push_back(base->GetCMEnergy());
    else
      for (int j = 1; j <= base->NumSubPoints(); j++) energies.push_back(base->GetSubPoint(j)->GetCMEnergy());
    std::vector<std::pair<int, double>> dw;
    rowStart_.push_back(0);
    for (double e : energies) {
      base->DopplerWeights(e, dw);
      for (const auto &t : dw) {
        col_.push_back(t.first - 1);
        w_.push_back(t.second);
      }
      rowStart_.push_back(col_.size());
    }
    stage_.premap = {energies.size(), rowStart_.data(), col_.data(), w_.data()};
  }
  azure2::experiment::TransformStage stage(std::size_t) const override { return stage_; }
  double fold(std::size_t, const double *f, std::size_t n) const override {
    EPoint::DopplerAppliedScope dopplerDone(base_);
    return base_->IntegrateSubPointValues(configure_, std::vector<double>(f, f + n));
  }
  bool weights(std::size_t, double *c, std::size_t n) const override {
    std::vector<double> weights;
    EPoint::DopplerAppliedScope dopplerDone(base_);
    if (!base_->TargetEffectWeights(configure_, weights) || weights.size() != n) return false;
    std::copy(weights.begin(), weights.end(), c);
    return true;
  }

 private:
  EPoint *base_;
  const Config &configure_;
  azure2::experiment::TransformStage stage_;
  std::vector<std::size_t> rowStart_;
  std::vector<int> col_;
  std::vector<double> w_;
};

/*!
 * The model value of a segment at a point: its expression tree (Predictions.h) evaluated on the values of its terms'
 * points.  The leaves are read from the points' records (the E1 / E2 part where the leaf says so); the frame factors
 * and the abundances and thickness are taken at the point; a transmission or yield reads its terms' sub-point values
 * and the base point's kernel.  The per-nucleus terms of a sample segment go into the base point's leaf record (the
 * decomposition output).
 */
double ComposeModelValue(ESegment &segment, int pointIndex, CNuc *cnuc, const Config &configure,
                         azure2::Results &results, azure2::LeafValues leaves) {
  const Prediction *p = segment.GetPrediction();
  if (!p) throw std::logic_error("EvaluateModelPoint: segment " + std::to_string(segment.GetSegmentKey()) +
                                 " has no model expression (EData::BuildPredictions)");
  EPoint *basePoint = segment.GetPoint(pointIndex + 1);
  const std::size_t nl = p->leaves.size(), nf = p->factors.size(), np = p->params.size();

  // The terms' points, calculated here when the leaves are not in yet (in term order, the base first).
  LocalArray<EPoint *, 16> point(nl);
  for (std::size_t k = 0; k < nl; k++) {
    const int term = p->leaves[k].term;
    point[k] = TermPoint(segment, term, pointIndex);
    if (!point[k]) return FailedSample(segment, basePoint, cnuc, results);  // component segments copy the base's points
    if (CalculatesItself(point[k], leaves) &&
        !CalculatePoint(point[k], TermRole(*p, segment, term), cnuc, configure, results))
      return FailedSample(segment, basePoint, cnuc, results);
  }

  LocalArray<double, 16> leafValue(nl), factorValue(nf), paramValue(np);
  LocalArray<azure2::experiment::PointArray, 16> leafArray(nl), factorArray(nf);
  for (std::size_t k = 0; k < nl; k++) {
    leafValue[k] = LeafValue(point[k], p->leaves[k].part, results);
    leafArray[k] = {&leafValue[k], 1};
  }
  for (std::size_t f = 0; f < nf; f++) {
    const FactorRef &r = p->factors[f];
    EPoint *termPoint = TermPoint(segment, r.term, pointIndex);
    factorValue[f] = (r.kind == FactorRef::Kind::TermFrame)
                         ? ESegment::TermFrameFactor(basePoint, termPoint)
                         : ESegment::RatioTermFrameFactor(termPoint, TermSegment(segment, r.term), cnuc);
    factorArray[f] = {&factorValue[f], 1};
  }
  for (std::size_t q = 0; q < np; q++) {
    const ParamRef &r = p->params[q];
    paramValue[q] = (r.kind == ParamRef::Kind::Abundance) ? cnuc->GetAbundance(segment.GetNaturalSample(), r.nucleus)
                                                          : segment.GetSampleThickness();
  }
  azure2::experiment::Inputs in;
  in.leaves = {leafArray.data(), nl};
  in.factors = {factorArray.data(), nf};
  in.params = paramValue.data();

  // A transmission or yield under a sub-point effect: its terms' sub-point values and the base point's kernel.
  std::vector<double> subValues;
  std::vector<std::size_t> subOffsets;
  std::vector<azure2::experiment::SubPointArray> subArray;
  std::unique_ptr<BaseKernel> kernel;
  const azure2::experiment::SubPointFold *folds[1] = {nullptr};
  if (!p->tree.Transforms().empty()) {
    kernel = std::make_unique<BaseKernel>(basePoint, configure);
    folds[0] = kernel.get();
    in.folds = {folds, 1};
    if (ESegment::HasSubPointIntegration(basePoint)) {
      for (std::size_t k = 0; k < nl; k++) {
        subOffsets.push_back(subValues.size());
        for (int j = 1; j <= point[k]->NumSubPoints(); j++)
          subValues.push_back(LeafValue(point[k]->GetSubPoint(j), p->leaves[k].part, results));
        subOffsets.push_back(subValues.size());
      }
      for (std::size_t k = 0; k < nl; k++) subArray.push_back({subValues.data(), &subOffsets[2 * k]});
      in.subLeaves = subArray;
    }
  }

  azure2::experiment::Workspace ws;
  double m;
  if (segment.IsNatural()) {
    std::vector<double> naturalTerms(cnuc->NumNuclei(), 0.0);
    m = p->tree.Value(in, 0, ws, naturalTerms.data());
    basePoint->SetNaturalTerms(results, naturalTerms);
  } else {
    m = p->tree.Value(in, 0, ws);
  }
  if (p->kind == azure2::experiment::TreeKind::Ratio && leafValue[1] == 0.0) {
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

double EvaluateModelPoint(ESegment &segment, int pointIndex, CNuc *compound, const Config &config, Results &results,
                          LeafValues leaves) {
  EPoint *point = (pointIndex >= 0 && pointIndex < segment.NumPoints()) ? segment.GetPoint(pointIndex + 1) : nullptr;
  if (!point) throw std::logic_error("EvaluateModelPoint: segment " + std::to_string(segment.GetSegmentKey()) +
                                     " has no point " + std::to_string(pointIndex));
  const double m = ComposeModelValue(segment, pointIndex, compound, config, results, leaves);
  point->SetModelValue(results, m);
  point->CopyAngularDistsToModel(results);
  return m;
}

void EvaluateModelSegment(ESegment &segment, CNuc *compound, const Config &config, Results &results,
                          LeafValues leaves) {
  for (int pid = 0; pid < segment.NumPoints(); pid++)
    EvaluateModelPoint(segment, pid, compound, config, results, leaves);
}

void EvaluateModel(EData &data, CNuc *compound, const Config &config, Results &results, LeafValues leaves) {
  for (int i = 1; i <= data.NumSegments(); i++) EvaluateModelSegment(*data.GetSegment(i), compound, config, results, leaves);
}

Residual PointResidual(const EPoint &point, const Results &results, double norm) {
  return DataResidual(point.GetModelValue(results), point.GetCMCrossSection(), point.GetCMCrossSectionError(), norm);
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
