// The fit objective as one ordered list of terms (Objective.h; docs/upgrade/objective.md).  The expressions are those
// of AZURECalc::operator() before Phase 2, moved here unchanged, in its order.
#include "Objective.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>

#include "AChannel.h"
#include "ALevel.h"
#include "CompoundSystem.h"
#include "JGroup.h"
#include "PPair.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "Unit.h"
#include "ModelPass.h"
#include "InferenceTypes.h"
#include "ModelData.h"
#include "DiagonalPriors.h"
#include "ApparatusParameters.h"

namespace {

// T1-T4: per segment its data chi^2, its normalization prior and its energy-shift priors (the priors gathered by
// AppendSegmentPriors, DiagonalPriors.h).
class SegmentTerms : public ObjectiveTerm {
 public:
  std::string Name() const override { return "data, normalizations and energy shifts"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    WorkingCopy *ld = c.data;
    if (c.segmentChi2) c.segmentChi2->assign(ld->Data().NumSegments() + 1, 0.0);
    ObjectiveContext::SegmentPriors *out = c.segmentPriorsOut;
    if (out) {
      out->norm.assign(ld->Data().NumSegments() + 1, 0.0);
      out->shift.assign(ld->Data().NumSegments() + 1, 0.0);
      out->shiftSqrt.assign(ld->Data().NumSegments() + 1, 0.0);
      out->shiftFn.assign(ld->Data().NumSegments() + 1, 0.0);
    }
    for (int i = 1; i <= ld->Data().NumSegments(); i++) {
      const Segment *segment = ld->Data().GetSegment(i);
      if (!segment) continue;
      if (c.thmOnly && !segment->IsTHM()) continue;
      // A THM segment with a free norm has an arbitrary overall scale: during a fit it is profiled out
      // analytically (n* set to its optimum, no norm penalty) so it cannot trade against the direct-data
      // normalizations.  A THM experiment (<thm> experiment[...]) profiles one norm (and a background) over all of
      // its segments: once, at its last segment, when the models of all of them are in.  It sets their norms and
      // chi^2.
      double segmentChiSquared = 0.0;
      const int thmGroup = ld->Data().ThmGroupOf(i);
      if (thmGroup >= 0) {
        if (ld->Data().IsLastOfThmGroup(thmGroup, i)) {
          const double g = ld->ProfileThmGroup(thmGroup, c.results);
          chi2 += g;
          if (c.segmentChi2) (*c.segmentChi2)[i] = g;
        }
      } else if (segment->IsProfiledNorm()) {
        segmentChiSquared = segment->ProfileNormChiSquared(*c.data, c.results, c.results.objective.norm[i]);
      } else {
        for (int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++)
          segmentChiSquared += azure2::PointResidual(*c.data, *segment, pointIdx, c.results, c.data->State(*segment).norm).Chi2();
      }
      // The segment's priors (AppendSegmentPriors): its normalization prior (profiled norms have none; a fixed norm
      // pays none, objective.md D5), then its energy-shift priors, each added to the running sum.
      if (c.segmentPriors) {
        DiagonalPriors priors;
        AppendSegmentPriors(*c.data, i, nullptr, priors);
        AddPriorsChi2(priors, chi2);
        if (out) {
          out->norm[i] = PriorsChi2(priors, PriorSource::kNorm);
          out->shift[i] = PriorsChi2(priors, PriorSource::kShift);
          out->shiftSqrt[i] = PriorsChi2(priors, PriorSource::kShiftSqrt);
          out->shiftFn[i] = PriorsChi2(priors, PriorSource::kShiftFn);
        }
      }
      if (thmGroup < 0) {
        c.results.objective.segmentChi2[i] = segmentChiSquared;
        if (c.segmentChi2) (*c.segmentChi2)[i] = segmentChiSquared;
        chi2 += segmentChiSquared;
      }
    }
  }
  // The per-segment table, Total-Chi-Squared, Total-Norm-Chi-Squared and Total-N are the writer's own
  // (WorkingCopy::WriteChiSquaredFile); here the energy-shift priors the fit minimised (AppendShiftPriors), when any
  // shift is a fit parameter.  Without this the file reported less than the minimiser's own total by this amount.
  void WriteTotal(const WorkingCopy &data, const Config &, const azure2::Results &, std::ostream &out) const override {
    bool anyShiftVaried = false;
    DiagonalPriors shiftPriors;
    for (int i = 1; i <= data.Data().NumSegments(); i++) {
      const Segment *segment = data.Data().GetSegment(i);
      if (segment->IsVaryEnergyShift() || segment->IsVaryEnergyShiftSqrt()) anyShiftVaried = true;
      if (segment->OwnsShiftFn())
        for (int q = 0; q < segment->NumShiftFnParams(); q++)
          if (segment->IsVaryShiftFn(q)) anyShiftVaried = true;
      AppendShiftPriors(data, i, nullptr, shiftPriors);
    }
    if (anyShiftVaried) out << " Total-Shift-Chi-Squared: " << PriorsChi2(shiftPriors);
  }
};

// T2: r^T M^-1 r over each covariance group in place of its points' r_i^2.  (A THM segment carries no covariance
// block: WorkingCopy refuses one.)
class DataCovarianceTerm : public ObjectiveTerm {
 public:
  std::string Name() const override { return "data covariance"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (c.data->Data().HasCovariance()) chi2 += c.data->CovarianceChi2Correction(c.results, c.segmentChi2);
  }
};

class NuisanceTerm : public ObjectiveTerm {
 public:
  std::string Name() const override { return "nuisance priors"; }
  // The sampler: normalised Gaussian log-priors from the limits manager (AZURECalcMCMC::BuildAutoPriors).
  void AddForSampler(const ObjectiveContext &, double &) const override {}
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!c.nuisance) return;
    const double v = c.nuisance(c.p);
    chi2 += v;
    StoreTotal(c, v);
  }
  // Part of the objective in every path (objective.md D6); the total of the last evaluation with the limits manager.
  std::string TotalLabel() const override { return "Total-Nuisance-Chi-Squared"; }
  void WriteTotal(const WorkingCopy &, const Config &, const azure2::Results &results, std::ostream &out) const override {
    if (results.objective.Total(Slot()) != 0.0) WriteStoredTotal(results, out);
  }
};

// T6: Park formalism, every level's overlap J = 1 - sum_c gamma_c^2 dS_c/dE kept positive: the sum over the R-matrix
// levels of (J/kParkNormScale)^2 for J < 0 (zero for every valid parameter set).  A level with J <= 0 has no standard
// R-matrix counterpart (its Brune reduced widths would be imaginary), and nothing in the level matrix itself keeps a
// minimizer out of that region, so this is the wall.  The shift functions must be current (CompoundSystem::CalcShiftFunctions).
// The gather is ParkLevelsOf, the arithmetic inference/'s ParkWall.h, which the gradient (AddParkPenaltyGradient)
// and the rows (AppendModelRows) call.
class ParkTerm : public ObjectiveTerm {
 public:
  bool ModelOnly() const override { return true; }
  std::string Name() const override { return "Park J > 0"; }
  // The sampler: a support, J < 0 has no R-matrix meaning (decision D2); inside it the wall is 0.
  void AddForSampler(const ObjectiveContext &, double &) const override {}
  bool ExcludesForSampler(const ObjectiveContext &c) const override {
    return (c.config.paramMask & Config::USE_PARK_FORMALISM) && ParkWallChi2(ParkLevelsOf(c.compound, c.state, false), kParkNormScale) > 0.0;
  }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!(c.config.paramMask & Config::USE_PARK_FORMALISM)) return;
    const double parkPenalty = ParkWallChi2(ParkLevelsOf(c.compound, c.state, false), kParkNormScale);
    chi2 += parkPenalty;
    StoreTotal(c, parkPenalty);
  }
  // Park formalism: the J > 0 wall, when a level sits beyond it (0 otherwise, printed in the Park formalism).
  std::string TotalLabel() const override { return "Total-Park-Chi-Squared"; }
  void WriteTotal(const WorkingCopy &, const Config &config, const azure2::Results &results, std::ostream &out) const override {
    if (config.paramMask & Config::USE_PARK_FORMALISM) WriteStoredTotal(results, out);
  }
};

// T8: one term per kind of target-effect prior, in the engine's order: the subtotal from 0 of its priors
// (AppendTargetPriors, DiagonalPriors.h), added to the running sum.
class TargetPriorTerm : public ObjectiveTerm {
 public:
  enum Kind { Screening, Density, Doppler, Thickness, Apparatus };
  explicit TargetPriorTerm(Kind k) : kind_(k) {}
  // The sampler: normalised Gaussian log-priors from ForEachDiagonalPrior (DiagonalPriors.h).
  void AddForSampler(const ObjectiveContext &, double &) const override {}
  std::string Name() const override {
    static const char *names[] = {"screening priors", "density priors", "Doppler priors", "thickness priors",
                                  "apparatus-parameter priors"};
    return names[kind_];
  }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    DiagonalPriors priors;
    AppendTargetPriors(*c.data, Source(), nullptr, priors);
    chi2 += PriorsChi2(priors);
  }
  // The apparatus priors are printed one field per row (WriteApparatusChiSquared): no field holds the whole term.
  std::string TotalLabel() const override {
    switch (kind_) {
      case Screening: return "Total-Screening-Chi-Squared";
      case Density: return "Total-Density-Chi-Squared";
      case Doppler: return "Total-Doppler-Chi-Squared";
      case Thickness: return "Total-Thickness-Chi-Squared";
      case Apparatus: break;
    }
    return "";
  }
  // The subtotal of the kind at the working copy's present values, re-formed by the term's own gather (a stored total
  // could be stale after pyazr's write_outputs: objective_arrays.md section 4), when a parameter of the kind has a prior.
  void WriteTotal(const WorkingCopy &data, const Config &, const azure2::Results &, std::ostream &out) const override {
    const DataSet &set = data.Data();
    bool anyPrior = false;
    switch (kind_) {
      case Screening:
        for (int k = 1; k <= set.NumScreeningParams(); k++)
          if (set.GetScreeningTargetEffect(k)->GetScreeningSigma() > 0.0) anyPrior = true;
        break;
      case Density:
        for (int k = 1; k <= set.NumDensityParams(); k++)
          if (set.GetDensityTargetEffect(k)->GetDensitySigma() > 0.0) anyPrior = true;
        break;
      case Doppler:
        for (int k = 1; k <= set.NumDopplerParams(); k++)
          if (set.GetDopplerTargetEffect(k)->GetDopplerSigma() > 0.0) anyPrior = true;
        break;
      case Thickness:
        for (int k = 1; k <= set.NumThicknessParams(); k++) {
          const Segment *seg = set.GetSegment(set.GetThicknessSegment(k));
          if (seg->GetThicknessSigma() > 0.0 && !seg->IsThicknessGroupFollower()) anyPrior = true;
        }
        break;
      case Apparatus:  // one line per row with a fitted prior
        WriteApparatusChiSquared(data, out);
        return;
    }
    if (!anyPrior) return;
    DiagonalPriors priors;
    AppendTargetPriors(data, Source(), nullptr, priors);
    out << " " << TotalLabel() << ": " << PriorsChi2(priors);
  }

 private:
  PriorSource Source() const {
    switch (kind_) {
      case Screening: return PriorSource::kScreening;
      case Density: return PriorSource::kDensity;
      case Doppler: return PriorSource::kDoppler;
      case Thickness: return PriorSource::kThickness;
      case Apparatus: break;
    }
    return PriorSource::kApparatus;
  }
  Kind kind_;
};

// T7: isotopic abundances (the samples of DataSet::Abundances at the working copy's values), summed over samples:
// every entry with a prior adds ((w - w0)/sigma)^2, and each sample whose width s is > 0 adds ((sum w) - c)^2 / s^2,
// c the sum's centre (1 unless written).  The gather is AbundanceBlockOf (DiagonalPriors.h), the arithmetic
// inference/'s AbundancePriors.h, which the gradient (AddAbundancePenaltyGradient) and the rows (AppendModelRows) call.
class AbundanceTerm : public ObjectiveTerm {
 public:
  bool ModelOnly() const override { return true; }
  std::string Name() const override { return "abundances"; }
  // The sampler: the per-nucleus priors are normalised Gaussian log-priors (ForEachDiagonalPrior); the sum-to-one
  // term, which couples them, is here.
  // (Adding the full term there counted the per-nucleus priors twice: an abundance's posterior 1/sqrt(2) too narrow.)
  void AddForSampler(const ObjectiveContext &c, double &chi2) const override {
    if (c.data->Data().Abundances().Present()) chi2 += AbundanceChi2(AbundanceBlockOf(*c.data, c.config.abundanceSumSigma), false);
  }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!c.data->Data().Abundances().Present()) return;
    const double abundancePenalty = AbundanceChi2(AbundanceBlockOf(*c.data, c.config.abundanceSumSigma), true);
    chi2 += abundancePenalty;
    StoreTotal(c, abundancePenalty);
  }
  // Isotopic abundances: their priors and the sum-to-one term.
  std::string TotalLabel() const override { return "Total-Abundance-Chi-Squared"; }
  void WriteTotal(const WorkingCopy &, const Config &, const azure2::Results &results, std::ostream &out) const override {
    if (results.objective.Total(Slot()) != 0.0) WriteStoredTotal(results, out);
  }
};

// T9: the model data (ModelData.h): measured scattering lengths and the level
// quantities of <modelData> (resonance strengths, energies, widths, branching ratios, ANCs).
class ModelDataTerm : public ObjectiveTerm {
 public:
  bool ModelOnly() const override { return true; }
  std::string Name() const override { return "model data"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    const ModelDataSet &set = c.compound->GetModelData();
    if (set.NumResiduals() == 0) return;
    const double v = ModelData::Chi2(c.compound, c.state, c.inputs, c.config, set);
    chi2 += v;
    if (set.hasSection) StoreTotal(c, v);
  }
  // The <modelData> data (modelData.out), whenever the project has the section.
  std::string TotalLabel() const override { return "Total-Model-Data-Chi-Squared"; }
  void WriteTotal(const WorkingCopy &, const Config &, const azure2::Results &results, std::ostream &out) const override {
    if (results.objective.Stored(Slot())) WriteStoredTotal(results, out);
  }
};

// T11: the priors of the gamma_decay records and the wall on a negative remainder (GammaDecayMap.h,
// gamma_parametrisation.md section 2.6), from the records' values of the last fill.
class GammaDecayTerm : public ObjectiveTerm {
 public:
  bool ModelOnly() const override { return true; }
  std::string Name() const override { return "gamma decay priors"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (c.compound->GetGammaDecays().Empty()) return;
    const double v = c.compound->GetGammaDecays().Penalty(*c.compound, c.state);
    chi2 += v;
    StoreTotal(c, v);
  }
  // Whenever the project has a record.  Without this field the file's totals did not add up to the minimised
  // chi-squared for such a project.
  std::string TotalLabel() const override { return "Total-Gamma-Decay-Chi-Squared"; }
  void WriteTotal(const WorkingCopy &, const Config &, const azure2::Results &results, std::ostream &out) const override {
    if (results.objective.Stored(Slot())) WriteStoredTotal(results, out);
  }
};

// T10: the prior parameter covariance (<config> priorCovariance=), (p - mu)^T C^-1 (p - mu).
class PriorCovarianceTerm : public ObjectiveTerm {
 public:
  std::string Name() const override { return "prior covariance"; }
  // The sampler: a normalised multivariate Gaussian log-prior (AZURECalcMCMC::CalculateLogPrior).
  void AddForSampler(const ObjectiveContext &, double &) const override {}
  void Add(const ObjectiveContext &c, double &chi2) const override {
    const PriorCovariance *pc = c.data->Data().GetPriorCovariance();
    if (!pc) return;
    const double pcChi2 = pc->Chi2(c.p);
    chi2 += pcChi2;
    StoreTotal(c, pcChi2);
  }
  // Whenever the key is set; with 15 significant digits (a newer field: the older ones keep their 6).
  std::string TotalLabel() const override { return "Total-Prior-Covariance-Chi-Squared"; }
  void WriteTotal(const WorkingCopy &, const Config &config, const azure2::Results &results, std::ostream &out) const override {
    if (config.priorCovariance.empty()) return;
    const std::streamsize prec = out.precision(15);
    WriteStoredTotal(results, out);
    out.precision(prec);
  }
};

}  // namespace

void ObjectiveTerm::WriteStoredTotal(const azure2::Results &results, std::ostream &out) const {
  out << " " << TotalLabel() << ": " << results.objective.Total(Slot());
}

template <class T, class... A>
const ObjectiveTerm *Objective::Append(A &&...args) {
  std::unique_ptr<ObjectiveTerm> term = std::make_unique<T>(std::forward<A>(args)...);
  term->slot_ = (int)terms_.Terms().size();
  const ObjectiveTerm *t = term.get();
  terms_.Append(std::move(term));
  return t;
}

namespace {
// A walk's order: the terms named, in that order, then every other term of the list in the value's order.
std::vector<const ObjectiveTerm *> WalkOrder(const std::vector<std::unique_ptr<ObjectiveTerm>> &terms,
                                             std::vector<const ObjectiveTerm *> first) {
  for (const auto &t : terms)
    if (std::find(first.begin(), first.end(), t.get()) == first.end()) first.push_back(t.get());
  return first;
}
}  // namespace

Objective::Objective() {
  // The value: the order of the sum (part of its bits).  A new term is appended here.
  const ObjectiveTerm *segments = Append<SegmentTerms>();
  Append<DataCovarianceTerm>();
  const ObjectiveTerm *nuisance = Append<NuisanceTerm>();
  const ObjectiveTerm *park = Append<ParkTerm>();
  const ObjectiveTerm *screening = Append<TargetPriorTerm>(TargetPriorTerm::Screening);
  const ObjectiveTerm *density = Append<TargetPriorTerm>(TargetPriorTerm::Density);
  const ObjectiveTerm *doppler = Append<TargetPriorTerm>(TargetPriorTerm::Doppler);
  const ObjectiveTerm *thickness = Append<TargetPriorTerm>(TargetPriorTerm::Thickness);
  const ObjectiveTerm *apparatus = Append<TargetPriorTerm>(TargetPriorTerm::Apparatus);
  const ObjectiveTerm *abundance = Append<AbundanceTerm>();
  const ObjectiveTerm *modelData = Append<ModelDataTerm>();
  priorCovariance_ = Append<PriorCovarianceTerm>();
  const ObjectiveTerm *gammaDecay = Append<GammaDecayTerm>();
  // The other walks have orders of their own, older than the list (each is part of a result: the gradient's and the
  // rows' sums, the bytes of chiSquared.out).  Each names the terms that take part, in that order; every other term
  // follows in the value's order (its function does nothing), so a new term is visited last by every walk.
  gradientWalk_ = WalkOrder(Terms(), {gammaDecay, modelData, park, abundance});
  rowWalk_ = WalkOrder(Terms(), {park, gammaDecay, abundance});
  reportWalk_ = WalkOrder(Terms(), {park, segments, abundance, screening, density, apparatus, doppler, thickness,
                                    nuisance, modelData, gammaDecay, priorCovariance_});
}

void Objective::WriteTotals(const WorkingCopy &data, const Config &config, const azure2::Results &results,
                            std::ostream &out) const {
  for (const ObjectiveTerm *t : reportWalk_) t->WriteTotal(data, config, results, out);
}

void Objective::RecordPriorCovariance(const ObjectiveContext &c) const {
  double scratch = 0.0;
  priorCovariance_->Add(c, scratch);
}

const Objective &Objective::Engine() {
  static const Objective objective;
  return objective;
}

double Objective::DataChi2(ObjectiveContext c) const {
  c.segmentPriors = false;
  c.thmOnly = false;
  return terms_.SumOfFirst(c, 2);  // T1 (T3, T4 off) and T2
}

double Objective::SamplerChi2(ObjectiveContext c) const {
  c.segmentPriors = false;  // T3, T4: the sampler's normalised log-priors
  c.thmOnly = false;
  return terms_.SumForSampler(c);
}

bool Objective::ExcludesForSampler(const ObjectiveContext &c) const { return terms_.ExcludesForSampler(c); }

double Objective::ModelChi2(const ObjectiveContext &c) const { return terms_.SumOfModelOnly(c); }

double Objective::ModelDataChi2(const ObjectiveContext &c) const {
  return ModelData::Chi2(c.compound, c.state, c.inputs, c.config, c.compound->GetModelData());
}

void Objective::ModelDataResiduals(const ObjectiveContext &c, vector_r &residuals) const {
  ModelData::Residuals(c.compound, c.state, c.inputs, c.config, c.compound->GetModelData(), residuals);
}

bool Objective::ModelDataGradient(const ObjectiveContext &c, const ParamIndexMap &pmap,
                                  const vector_matrix_r *shiftDeriv, GradAccum &accum) const {
  // T11 (gamma_decay priors, the remainder's wall) depends on the parameters alone: its gradient goes straight to their
  // indices (GradAccum::direct).
  c.compound->GetGammaDecays().AddPenaltyGradient(*c.compound, c.state, accum.direct);
  return ModelData::Gradient(c.compound, c.state, c.inputs, c.config, c.compound->GetModelData(), pmap, shiftDeriv, accum);
}

bool Objective::AppendModelDataRows(const ObjectiveContext &c, const ParamIndexMap &pmap,
                                    const vector_matrix_r *shiftDeriv, int nCols, vector_r &residuals,
                                    vector_r &jac) const {
  return ModelData::AppendRows(c.compound, c.state, c.inputs, c.config, c.compound->GetModelData(), pmap, shiftDeriv, nCols, residuals,
                               jac);
}

double Objective::Chi2(const ObjectiveContext &c) const {
  if (c.thmOnly) return terms_.SumOfFirst(c, 1);
  if (!c.termsOut) return terms_.Sum(c);
  // With a breakdown: each term on its own (a small prior next to a large data chi^2 keeps its digits), summed in
  // order.  The total may then differ from the plain sum in its last bits.
  return terms_.SumWithParts(c, *c.termsOut);
}

std::vector<std::string> ReplaceDiagonalPriors(const CompoundSystem *compound, WorkingCopy *data, const PriorCovariance &pc) {
  std::vector<std::string> replaced;
  const ParamIndexMap pmap = BuildParamIndexMap(compound, data, std::vector<bool>());
  // The flags are settings of the data set: one revision of it (WorkingCopy::ReviseDataSet) for all of them.
  data->ReviseDataSet([&](DataSet &ds) {
    for (size_t a = 0; a < pc.FullIndex().size(); a++) {
      const int f = pc.FullIndex()[a];
      if (f < 0 || f >= pmap.NumFull()) continue;
      const ParamDesc &d = pmap.Desc(f);
      bool had = false;
      switch (d.kind) {
        case ParamKind::Norm: {
          Segment *seg = ds.GetSegment(d.segment);
          had = seg->GetNormError() != 0.0;
          seg->ReplaceNormPrior();
          break;
        }
        case ParamKind::EnergyShift: {
          Segment *seg = ds.GetSegment(d.segment);
          had = seg->GetEnergyShiftError() != 0.0;
          seg->ReplaceShiftPrior();
          break;
        }
        case ParamKind::EnergyShiftSqrt: {
          Segment *seg = ds.GetSegment(d.segment);
          had = seg->GetEnergyShiftSqrtError() != 0.0;
          seg->ReplaceShiftSqrtPrior();
          break;
        }
        case ParamKind::EnergyShiftFn: {
          Segment *seg = ds.GetSegment(d.segment);
          had = d.fnTerm >= 0 && d.fnTerm < (int)seg->GetShiftFnErrors().size() && seg->GetShiftFnErrors()[d.fnTerm] != 0.0;
          seg->ReplaceShiftFnPrior(d.fnTerm);
          break;
        }
        case ParamKind::Thickness: {
          Segment *seg = ds.GetSegment(ds.GetThicknessSegment(d.thickness));
          had = seg->GetThicknessSigma() > 0.0;
          seg->ReplaceThicknessPrior();
          break;
        }
        case ParamKind::Screening: {
          TargetEffect *te = ds.GetScreeningTargetEffect(d.screening);
          had = te->GetScreeningSigma() > 0.0;
          te->ReplaceScreeningPrior();
          break;
        }
        case ParamKind::Density: {
          TargetEffect *te = ds.GetDensityTargetEffect(d.density);
          had = te->GetDensitySigma() > 0.0;
          te->ReplaceDensityPrior();
          break;
        }
        case ParamKind::Doppler: {
          TargetEffect *te = ds.GetDopplerTargetEffect(d.doppler);
          had = te->GetDopplerSigma() > 0.0;
          te->ReplaceDopplerPrior();
          break;
        }
        case ParamKind::Apparatus: {
          ApparatusParam &a = ds.GetApparatusParam(d.apparatus);
          had = a.declaration.sigma > 0.0;
          a.priorReplaced = true;
          break;
        }
        case ParamKind::Abundance: {
          had = ds.Abundances().Sample(d.sample).entries[d.nucleus - 1].sigma > 0.0;
          ds.ReplaceAbundancePrior(d.sample, d.nucleus);
          break;
        }
        default:  // levels and THM coherent background: a nuisance prior, skipped by index
          break;
      }
      if (had) replaced.push_back(pc.Names()[a]);
    }
  });
  return replaced;
}

std::vector<ParkLevel> ParkLevelsOf(const CompoundSystem *compound, const LevelState &state, bool channels) {
  std::vector<ParkLevel> levels;
  for (int j = 1; j <= compound->NumJGroups(); j++) {
    const JGroup *jg = compound->GetJGroup(j);
    if (!jg->IsInRMatrix()) continue;
    for (int la = 1; la <= jg->NumLevels(); la++) {
      if (!jg->GetLevel(la)->IsInRMatrix()) continue;
      ParkLevel l;
      l.group = j;
      l.level = la;
      l.J = state.ParkNorm(j, la);
      for (int ch = 1; channels && azure2::inference::OnParkWall(l.J) && ch <= jg->NumChannels(); ch++) {
        const AChannel *channel = jg->GetChannel(ch);
        if (channel->GetRadType() != 'P' || compound->GetPair(channel->GetPairNum())->GetPType() != 0) continue;
        l.channels.push_back(ch);
        l.gamma.push_back(state.Gamma(j, la, ch));
        l.dS.push_back(state.ShiftDerivative(j, la, ch));
        l.d2S.push_back(state.ShiftSecondDerivative(j, la, ch));
      }
      levels.push_back(std::move(l));
    }
  }
  return levels;
}

void AppendModelRows(const CompoundSystem *compound, const LevelState &state, const WorkingCopy &data, const Config &config, const ParamIndexMap *pmap, int nCols,
                     vector_r &residuals, vector_r *jac) {
  auto appendRow = [&](double r, const GradAccum *acc, const std::vector<std::pair<int, double>> *direct) {
    residuals.push_back(r);
    if (!jac || !pmap) return;
    const size_t base = jac->size();
    jac->resize(base + (size_t)nCols, 0.0);
    vector_r rowFull(pmap->NumFull(), 0.0);
    if (acc) acc->Scatter(*pmap, rowFull);
    if (direct)
      for (const auto &d : *direct)
        if (d.first >= 0 && d.first < pmap->NumFull()) rowFull[d.first] += d.second;
    for (int f = 0; f < pmap->NumFull(); f++) {
      const int a = pmap->FullToPacked(f);
      if (a >= 0 && a < nCols) (*jac)[base + (size_t)a] = rowFull[f];
    }
  };
  // Park: one row per R-matrix level, r = J/s on the wall, else 0 with no Jacobian (ParkWallRow).
  if (config.paramMask & Config::USE_PARK_FORMALISM) {
    std::vector<double> dGamma;
    double r, dE;
    for (const ParkLevel &l : ParkLevelsOf(compound, state, true)) {
      if (!ParkWallRow(l, kParkNormScale, r, dGamma, dE)) {
        appendRow(0.0, nullptr, nullptr);
        continue;
      }
      GradAccum acc;
      acc.Init(compound, state);
      for (size_t c = 0; c < l.channels.size(); c++) acc.AddGamma(l.group, l.level, l.channels[c], dGamma[c]);
      acc.AddE(l.group, l.level, dE);
      appendRow(r, &acc, nullptr);
    }
  }
  // gamma_decay priors and the remainder's wall: r and dr/dq (GammaDecaySet::PenaltyRows).
  for (const GammaDecaySet::PenaltyRow &row : compound->GetGammaDecays().PenaltyRows(*compound, state))
    appendRow(row.r, nullptr, &row.d);
  // Abundances: ((sum w) - c)/s per sample (c the sum's centre, 1 unless written) with s > 0; d/dw_k = 1/s on every
  // abundance of the sample (AbundanceSumRows).
  const AbundanceSamples &samples = data.Data().Abundances();
  if (samples.Present()) {
    const AbundanceBlock block = AbundanceBlockOf(data, config.abundanceSumSigma);
    for (const AbundanceSumRow &row : AbundanceSumRows(block)) {
      std::vector<std::pair<int, double>> direct;
      if (pmap)
        for (int i = 1; i <= block.numNuclei; i++) direct.push_back({pmap->AbundanceIndex(samples.Ordinal(row.sample, i)), row.derivative});
      appendRow(row.residual, nullptr, &direct);
    }
  }
}
