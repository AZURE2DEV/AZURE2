// The fit objective as one ordered list of terms (Objective.h; docs/upgrade/objective.md).  The expressions are those
// of AZURECalc::operator() before Phase 2, moved here unchanged, in its order.
#include "Objective.h"

#include <cmath>

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
#include "PriorCovariance.h"
#include "ModelData.h"

namespace {

// T1-T4: per segment its data chi^2, its normalization prior and its energy-shift priors.
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
        // Normalization prior (profiled norms have none; a fixed norm pays none, objective.md D5).
        if (c.segmentPriors && segment->IsVaryNorm() && !segment->IsNormPriorReplaced()) {
          double dataNorm = c.data->State(*segment).norm;
          double dataNormNominal = segment->GetNominalNorm();
          double dataNormError = dataNormNominal / 100. * segment->GetNormError();
          if (dataNormError != 0.) {
            const double term = pow((dataNorm - dataNormNominal) / dataNormError, 2.0);
            chi2 += term;
            if (out) out->norm[i] = term;
          }
        }
      }
      // Energy-shift priors.
      if (c.segmentPriors) segment->AddShiftPenalty(c.data->State(*segment), chi2);
      if (out && c.segmentPriors) segment->ShiftPenaltyParts(c.data->State(*segment), out->shift[i], out->shiftSqrt[i], out->shiftFn[i]);
      if (thmGroup < 0) {
        c.results.objective.segmentChi2[i] = segmentChiSquared;
        if (c.segmentChi2) (*c.segmentChi2)[i] = segmentChiSquared;
        chi2 += segmentChiSquared;
      }
    }
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
    c.results.objective.nuisance = v;
  }
};

// T6: Park formalism, every level's overlap J = 1 - sum_c gamma_c^2 dS_c/dE kept positive: the sum over the R-matrix
// levels of (J/kParkNormScale)^2 for J < 0 (zero for every valid parameter set).  A level with J <= 0 has no standard
// R-matrix counterpart (its Brune reduced widths would be imaginary), and nothing in the level matrix itself keeps a
// minimizer out of that region, so this is the wall.  The shift functions must be current (CompoundSystem::CalcShiftFunctions).
class ParkTerm : public ObjectiveTerm {
 public:
  bool ModelOnly() const override { return true; }
  std::string Name() const override { return "Park J > 0"; }
  // The sampler: a support, J < 0 has no R-matrix meaning (decision D2); inside it the wall is 0.
  void AddForSampler(const ObjectiveContext &, double &) const override {}
  bool ExcludesForSampler(const ObjectiveContext &c) const override {
    return (c.config.paramMask & Config::USE_PARK_FORMALISM) && Wall(*c.compound, c.state) > 0.0;
  }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!(c.config.paramMask & Config::USE_PARK_FORMALISM)) return;
    const double parkPenalty = Wall(*c.compound, c.state);
    chi2 += parkPenalty;
    c.results.objective.park = parkPenalty;
  }

 private:
  static double Wall(const CompoundSystem &compound, const LevelState &state) {
    double penalty = 0.0;
    for (int j = 1; j <= compound.NumJGroups(); j++) {
      const JGroup *group = compound.GetJGroup(j);
      if (!group->IsInRMatrix()) continue;
      for (int la = 1; la <= group->NumLevels(); la++) {
        const ALevel *level = group->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        const double parkNorm = state.ParkNorm(j, la);
        if (parkNorm < 0.0) penalty += pow(parkNorm / kParkNormScale, 2.0);
      }
    }
    return penalty;
  }
};

// T8: one term per kind of target-effect prior, in the engine's order.
class TargetPriorTerm : public ObjectiveTerm {
 public:
  enum Kind { Screening, Density, Doppler, Thickness, Apparatus };
  explicit TargetPriorTerm(Kind k) : kind_(k) {}
  // The sampler: normalised Gaussian log-priors from ForEachDiagonalPrior.
  void AddForSampler(const ObjectiveContext &, double &) const override {}
  std::string Name() const override {
    static const char *names[] = {"screening priors", "density priors", "Doppler priors", "thickness priors",
                                  "apparatus-parameter priors"};
    return names[kind_];
  }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    switch (kind_) {
      case Screening: chi2 += c.data->ScreeningPenalty(); break;
      case Density: chi2 += c.data->DensityPenalty(); break;
      case Doppler: chi2 += c.data->DopplerPenalty(); break;
      case Thickness: chi2 += c.data->ThicknessPenalty(); break;
      case Apparatus: chi2 += c.data->ApparatusPenalty(); break;
    }
  }

 private:
  Kind kind_;
};

// T7: isotopic abundances (the samples of DataSet::Abundances at the working copy's values), summed over samples:
// every entry with a prior adds ((w - w0)/sigma)^2, and each sample whose width s is > 0 adds ((sum w) - c)^2 / s^2,
// c the sum's centre (1 unless written).
class AbundanceTerm : public ObjectiveTerm {
 public:
  bool ModelOnly() const override { return true; }
  std::string Name() const override { return "abundances"; }
  // The sampler: the per-nucleus priors are normalised Gaussian log-priors (ForEachDiagonalPrior); the sum-to-one
  // term, which couples them, is here.
  // (Adding the full term there counted the per-nucleus priors twice: an abundance's posterior 1/sqrt(2) too narrow.)
  void AddForSampler(const ObjectiveContext &c, double &chi2) const override {
    if (c.data->Data().Abundances().Present()) chi2 += Penalty(c, Parts::kSumsOnly);
  }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!c.data->Data().Abundances().Present()) return;
    const double abundancePenalty = Penalty(c, Parts::kPriorsAndSums);
    chi2 += abundancePenalty;
    c.results.objective.abundance = abundancePenalty;
  }

 private:
  enum class Parts { kSumsOnly, kPriorsAndSums };
  static double Penalty(const ObjectiveContext &c, Parts parts) {
    const AbundanceSamples &samples = c.data->Data().Abundances();
    const std::vector<double> &values = c.data->AbundanceValues();
    double penalty = 0.0;
    for (int s = 1; s <= samples.NumSamples(); s++) {
      const AbundanceSample &sample = samples.Sample(s);
      for (int i = 1; parts == Parts::kPriorsAndSums && i <= (int)sample.entries.size(); i++) {
        const SampleAbundance &e = sample.entries[i - 1];
        const double w = values[samples.Ordinal(s, i) - 1];
        if (e.sigma > 0.0 && !e.priorReplaced) penalty += pow((w - e.nominal) / e.sigma, 2.0);
      }
      const double sig = samples.SumSigma(s, c.config.abundanceSumSigma);
      if (sig > 0.0) penalty += pow((samples.Sum(values, s) - sample.sumCentre) / sig, 2.0);
    }
    return penalty;
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
    if (set.hasSection) {
      c.results.objective.modelData = v;
      c.results.objective.hasModelData = true;
    }
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
    chi2 += c.compound->GetGammaDecays().Penalty(*c.compound, c.state);
  }
};

// T12: the priors of the fitted isospin shifts and mixes and the soft ties of the isospin levels (IsospinTies.h), at
// the levels and isospin values of the last fill.
class IsospinTerm : public ObjectiveTerm {
 public:
  bool ModelOnly() const override { return true; }
  std::string Name() const override { return "isospin priors and soft ties"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    const IsospinSet &isospin = c.compound->GetIsospin();
    if (!isospin.HasPenalty()) return;
    const double v = isospin.Penalty(*c.compound, c.state);
    chi2 += v;
    c.results.objective.isospin = v;
    c.results.objective.hasIsospin = true;
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
    c.results.objective.priorCovariance = pcChi2;
  }
};

}  // namespace

Objective::Objective() {
  terms_.emplace_back(new SegmentTerms);
  terms_.emplace_back(new DataCovarianceTerm);
  terms_.emplace_back(new NuisanceTerm);
  terms_.emplace_back(new ParkTerm);
  terms_.emplace_back(new TargetPriorTerm(TargetPriorTerm::Screening));
  terms_.emplace_back(new TargetPriorTerm(TargetPriorTerm::Density));
  terms_.emplace_back(new TargetPriorTerm(TargetPriorTerm::Doppler));
  terms_.emplace_back(new TargetPriorTerm(TargetPriorTerm::Thickness));
  terms_.emplace_back(new TargetPriorTerm(TargetPriorTerm::Apparatus));
  terms_.emplace_back(new AbundanceTerm);
  terms_.emplace_back(new ModelDataTerm);
  terms_.emplace_back(new PriorCovarianceTerm);
  terms_.emplace_back(new GammaDecayTerm);
  terms_.emplace_back(new IsospinTerm);
}

const Objective &Objective::Engine() {
  static const Objective objective;
  return objective;
}

double Objective::DataChi2(ObjectiveContext c) const {
  c.segmentPriors = false;
  c.thmOnly = false;
  double chi2 = 0.0;
  terms_[0]->Add(c, chi2);  // T1 (T3, T4 off)
  terms_[1]->Add(c, chi2);  // T2
  return chi2;
}

double Objective::SamplerChi2(ObjectiveContext c) const {
  c.segmentPriors = false;  // T3, T4: the sampler's normalised log-priors
  c.thmOnly = false;
  double chi2 = 0.0;
  for (const auto &t : terms_) t->AddForSampler(c, chi2);
  return chi2;
}

bool Objective::ExcludesForSampler(const ObjectiveContext &c) const {
  for (const auto &t : terms_)
    if (t->ExcludesForSampler(c)) return true;
  return false;
}

double Objective::ModelChi2(const ObjectiveContext &c) const {
  double chi2 = 0.0;
  for (const auto &t : terms_)
    if (t->ModelOnly()) t->Add(c, chi2);
  return chi2;
}

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
  // T12 (isospin priors and soft ties) likewise.
  c.compound->GetIsospin().AddPenaltyGradient(*c.compound, c.state, accum.direct);
  return ModelData::Gradient(c.compound, c.state, c.inputs, c.config, c.compound->GetModelData(), pmap, shiftDeriv, accum);
}

bool Objective::AppendModelDataRows(const ObjectiveContext &c, const ParamIndexMap &pmap,
                                    const vector_matrix_r *shiftDeriv, int nCols, vector_r &residuals,
                                    vector_r &jac) const {
  return ModelData::AppendRows(c.compound, c.state, c.inputs, c.config, c.compound->GetModelData(), pmap, shiftDeriv, nCols, residuals,
                               jac);
}

double Objective::Chi2(const ObjectiveContext &c) const {
  double chi2 = 0.0;
  if (c.thmOnly) {
    terms_.front()->Add(c, chi2);
    return chi2;
  }
  if (!c.termsOut) {
    for (const auto &t : terms_) t->Add(c, chi2);
    return chi2;
  }
  // With a breakdown: each term on its own (a small prior next to a large data chi^2 keeps its digits), summed in
  // order.  The total may then differ from the plain sum in its last bits.
  c.termsOut->clear();
  for (const auto &t : terms_) {
    double v = 0.0;
    t->Add(c, v);
    c.termsOut->push_back(v);
    chi2 += v;
  }
  return chi2;
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
  // Park: r = J/s for J < 0; dJ/dgamma_c = -2 gamma_c S'_c, dJ/dE = -sum_c gamma_c^2 S''_c (AddParkPenaltyGradient).
  if (config.paramMask & Config::USE_PARK_FORMALISM) {
    const double s = kParkNormScale;
    for (int j = 1; j <= compound->NumJGroups(); j++) {
      const JGroup *jg = compound->GetJGroup(j);
      if (!jg->IsInRMatrix()) continue;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        const ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        const double J = state.ParkNorm(j, la);
        if (!(J < 0.0)) {
          appendRow(0.0, nullptr, nullptr);
          continue;
        }
        GradAccum acc;
        acc.Init(compound, state);
        double dE = 0.0;
        for (int ch = 1; ch <= jg->NumChannels(); ch++) {
          const AChannel *channel = jg->GetChannel(ch);
          if (channel->GetRadType() != 'P' || compound->GetPair(channel->GetPairNum())->GetPType() != 0) continue;
          const double g = state.Gamma(j, la, ch);
          acc.AddGamma(j, la, ch, -2.0 * g * state.ShiftDerivative(j, la, ch) / s);
          dE -= g * g * state.ShiftSecondDerivative(j, la, ch);
        }
        acc.AddE(j, la, dE / s);
        appendRow(J / s, &acc, nullptr);
      }
    }
  }
  // gamma_decay priors and the remainder's wall: r and dr/dq (GammaDecaySet::PenaltyRows).
  for (const GammaDecaySet::PenaltyRow &row : compound->GetGammaDecays().PenaltyRows(*compound, state))
    appendRow(row.r, nullptr, &row.d);
  // The isospin priors and soft ties: r and dr/dq (IsospinSet::PenaltyRows).
  for (const IsospinSet::PenaltyRow &row : compound->GetIsospin().PenaltyRows(*compound, state))
    appendRow(row.r, nullptr, &row.d);
  // Abundances: ((sum w) - c)/s per sample (c the sum's centre, 1 unless written) with s > 0; d/dw_k = 1/s on every abundance of the sample.
  const AbundanceSamples &samples = data.Data().Abundances();
  if (samples.Present()) {
    for (int smp = 1; smp <= samples.NumSamples(); smp++) {
      const double s = samples.SumSigma(smp, config.abundanceSumSigma);
      if (!(s > 0.0)) continue;
      std::vector<std::pair<int, double>> direct;
      if (pmap)
        for (int k = 1; k <= samples.NumParams(); k++) {
          int sample, nucleus;
          samples.Of(k, sample, nucleus);
          if (sample == smp) direct.push_back({pmap->AbundanceIndex(k), 1.0 / s});
        }
      appendRow((samples.Sum(data.AbundanceValues(), smp) - samples.Sample(smp).sumCentre) / s, nullptr, &direct);
    }
  }
}
