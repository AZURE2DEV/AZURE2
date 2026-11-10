// The fit objective as one ordered list of terms (Objective.h; docs/upgrade/objective.md).  The expressions are those
// of AZURECalc::operator() before Phase 2, moved here unchanged, in its order.
#include "Objective.h"

#include <cmath>

#include "AChannel.h"
#include "ALevel.h"
#include "CNuc.h"
#include "JGroup.h"
#include "PPair.h"
#include "Config.h"
#include "EData.h"
#include "ESegment.h"
#include "EPoint.h"
#include "ModelPass.h"
#include "PriorCovariance.h"
#include "ScatteringLength.h"

namespace {

// T1-T4: per segment its data chi^2, its normalization prior and its energy-shift priors.
class SegmentTerms : public ObjectiveTerm {
 public:
  std::string Name() const override { return "data, normalizations and energy shifts"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    EData *ld = c.data;
    if (c.segmentChi2) c.segmentChi2->assign(ld->NumSegments() + 1, 0.0);
    ObjectiveContext::SegmentPriors *out = c.segmentPriorsOut;
    if (out) {
      out->norm.assign(ld->NumSegments() + 1, 0.0);
      out->shift.assign(ld->NumSegments() + 1, 0.0);
      out->shiftSqrt.assign(ld->NumSegments() + 1, 0.0);
      out->shiftFn.assign(ld->NumSegments() + 1, 0.0);
    }
    for (int i = 1; i <= ld->NumSegments(); i++) {
      ESegment *segment = ld->GetSegment(i);
      if (!segment) continue;
      if (c.thmOnly && !segment->IsTHM()) continue;
      // A THM segment with a free norm has an arbitrary overall scale: during a fit it is profiled out
      // analytically (n* set to its optimum, no norm penalty) so it cannot trade against the direct-data
      // normalizations.  A THM experiment (<thm> experiment[...]) profiles one norm (and a background) over all of
      // its segments: once, at its last segment, when the models of all of them are in.  It sets their norms and
      // chi^2.
      double segmentChiSquared = 0.0;
      const int thmGroup = ld->ThmGroupOf(i);
      if (thmGroup >= 0) {
        if (ld->IsLastOfThmGroup(thmGroup, i)) {
          const double g = ld->ProfileThmGroup(thmGroup, c.results);
          chi2 += g;
          if (c.segmentChi2) (*c.segmentChi2)[i] = g;
        }
      } else if (segment->IsProfiledNorm()) {
        segmentChiSquared = segment->ProfileNormChiSquared(c.results, c.results.objective.norm[i]);
      } else {
        for (int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++)
          segmentChiSquared += azure2::PointResidual(*segment->GetPoint(pointIdx + 1), c.results, segment->GetNorm()).Chi2();
        // Normalization prior (profiled norms have none; a fixed norm pays none, objective.md D5).
        if (c.segmentPriors && segment->IsVaryNorm() && !segment->IsNormPriorReplaced()) {
          double dataNorm = segment->GetNorm();
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
      if (c.segmentPriors) segment->AddShiftPenalty(chi2);
      if (out && c.segmentPriors) segment->ShiftPenaltyParts(out->shift[i], out->shiftSqrt[i], out->shiftFn[i]);
      if (thmGroup < 0) {
        c.results.objective.segmentChi2[i] = segmentChiSquared;
        if (c.segmentChi2) (*c.segmentChi2)[i] = segmentChiSquared;
        chi2 += segmentChiSquared;
      }
    }
  }
};

// T2: r^T M^-1 r over each covariance group in place of its points' r_i^2.  (A THM segment carries no covariance
// block: EData refuses one.)
class DataCovarianceTerm : public ObjectiveTerm {
 public:
  std::string Name() const override { return "data covariance"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (c.data->HasCovariance()) chi2 += c.data->CovarianceChi2Correction(c.results, c.segmentChi2);
  }
};

class NuisanceTerm : public ObjectiveTerm {
 public:
  std::string Name() const override { return "nuisance priors"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!c.nuisance) return;
    const double v = c.nuisance(c.p);
    chi2 += v;
    c.results.objective.nuisance = v;
  }
};

// T6: Park formalism, every level's overlap J kept positive.
class ParkTerm : public ObjectiveTerm {
 public:
  bool ModelOnly() const override { return true; }
  std::string Name() const override { return "Park J > 0"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!(c.config.paramMask & Config::USE_PARK_FORMALISM)) return;
    const double parkPenalty = c.compound->ParkNormPenalty();
    chi2 += parkPenalty;
    c.results.objective.park = parkPenalty;
  }
};

// T8: one term per kind of target-effect prior, in the engine's order.
class TargetPriorTerm : public ObjectiveTerm {
 public:
  enum Kind { Screening, Density, Doppler, Thickness };
  explicit TargetPriorTerm(Kind k) : kind_(k) {}
  std::string Name() const override {
    static const char *names[] = {"screening priors", "density priors", "Doppler priors", "thickness priors"};
    return names[kind_];
  }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    switch (kind_) {
      case Screening: chi2 += c.data->ScreeningPenalty(); break;
      case Density: chi2 += c.data->DensityPenalty(); break;
      case Doppler: chi2 += c.data->DopplerPenalty(); break;
      case Thickness: chi2 += c.data->ThicknessPenalty(); break;
    }
  }

 private:
  Kind kind_;
};

// T7: isotopic abundances, their priors and the sum-to-one term.
class AbundanceTerm : public ObjectiveTerm {
 public:
  bool ModelOnly() const override { return true; }
  std::string Name() const override { return "abundances"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!c.compound->HasAbundanceBlock()) return;
    const double abundancePenalty = c.compound->AbundancePenalty(c.config.abundanceSumSigma);
    chi2 += abundancePenalty;
    c.results.objective.abundance = abundancePenalty;
  }
};

// T9: measured scattering lengths (<scatteringLengths>), evaluated at each pair's threshold.
class ScatteringLengthTerm : public ObjectiveTerm {
 public:
  bool ModelOnly() const override { return true; }
  std::string Name() const override { return "scattering lengths"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!c.compound->HasScatteringLengthConstraints()) return;
    const double slChi2 = ScatteringLength::Chi2(c.compound, c.config, c.compound->ScatteringLengthConstraints());
    chi2 += slChi2;
    c.results.objective.scatteringLength = slChi2;
  }
};

// T10: the prior parameter covariance (<config> priorCovariance=), (p - mu)^T C^-1 (p - mu).
class PriorCovarianceTerm : public ObjectiveTerm {
 public:
  std::string Name() const override { return "prior covariance"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    const PriorCovariance *pc = c.data->GetPriorCovariance();
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
  terms_.emplace_back(new AbundanceTerm);
  terms_.emplace_back(new ScatteringLengthTerm);
  terms_.emplace_back(new PriorCovarianceTerm);
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

double Objective::ModelChi2(const ObjectiveContext &c) const {
  double chi2 = 0.0;
  for (const auto &t : terms_)
    if (t->ModelOnly()) t->Add(c, chi2);
  return chi2;
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

std::vector<std::string> ReplaceDiagonalPriors(CNuc *compound, EData *data, const PriorCovariance &pc) {
  std::vector<std::string> replaced;
  const ParamIndexMap pmap = BuildParamIndexMap(compound, data, std::vector<bool>());
  for (size_t a = 0; a < pc.FullIndex().size(); a++) {
    const int f = pc.FullIndex()[a];
    if (f < 0 || f >= pmap.NumFull()) continue;
    const ParamDesc &d = pmap.Desc(f);
    bool had = false;
    switch (d.kind) {
      case ParamKind::Norm: {
        ESegment *seg = data->GetSegment(d.segment);
        had = seg->GetNormError() != 0.0;
        seg->ReplaceNormPrior();
        break;
      }
      case ParamKind::EnergyShift: {
        ESegment *seg = data->GetSegment(d.segment);
        had = seg->GetEnergyShiftError() != 0.0;
        seg->ReplaceShiftPrior();
        break;
      }
      case ParamKind::EnergyShiftSqrt: {
        ESegment *seg = data->GetSegment(d.segment);
        had = seg->GetEnergyShiftSqrtError() != 0.0;
        seg->ReplaceShiftSqrtPrior();
        break;
      }
      case ParamKind::EnergyShiftFn: {
        ESegment *seg = data->GetSegment(d.segment);
        had = d.fnTerm >= 0 && d.fnTerm < (int)seg->GetShiftFnErrors().size() && seg->GetShiftFnErrors()[d.fnTerm] != 0.0;
        seg->ReplaceShiftFnPrior(d.fnTerm);
        break;
      }
      case ParamKind::Thickness: {
        ESegment *seg = data->GetSegment(data->GetThicknessSegment(d.thickness));
        had = seg->GetThicknessSigma() > 0.0;
        seg->ReplaceThicknessPrior();
        break;
      }
      case ParamKind::Screening: {
        TargetEffect *te = data->GetScreeningTargetEffect(d.screening);
        had = te->GetScreeningSigma() > 0.0;
        te->ReplaceScreeningPrior();
        break;
      }
      case ParamKind::Density: {
        TargetEffect *te = data->GetDensityTargetEffect(d.density);
        had = te->GetDensitySigma() > 0.0;
        te->ReplaceDensityPrior();
        break;
      }
      case ParamKind::Doppler: {
        TargetEffect *te = data->GetDopplerTargetEffect(d.doppler);
        had = te->GetDopplerSigma() > 0.0;
        te->ReplaceDopplerPrior();
        break;
      }
      case ParamKind::Abundance: {
        SampleAbundance &e = compound->GetSample(d.sample)->entries[d.nucleus - 1];
        had = e.sigma > 0.0;
        e.priorReplaced = true;
        break;
      }
      default:  // levels and THM coherent background: a nuisance prior, skipped by index
        break;
    }
    if (had) replaced.push_back(pc.Names()[a]);
  }
  return replaced;
}

void AppendModelRows(CNuc *compound, const Config &config, const ParamIndexMap *pmap, int nCols, vector_r &residuals,
                     vector_r *jac) {
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
    const double s = CNuc::kParkNormScale;
    for (int j = 1; j <= compound->NumJGroups(); j++) {
      JGroup *jg = compound->GetJGroup(j);
      if (!jg->IsInRMatrix()) continue;
      for (int la = 1; la <= jg->NumLevels(); la++) {
        ALevel *level = jg->GetLevel(la);
        if (!level->IsInRMatrix()) continue;
        const double J = level->GetParkNorm();
        if (!(J < 0.0)) {
          appendRow(0.0, nullptr, nullptr);
          continue;
        }
        GradAccum acc;
        acc.Init(compound);
        double dE = 0.0;
        for (int ch = 1; ch <= jg->NumChannels(); ch++) {
          AChannel *channel = jg->GetChannel(ch);
          if (channel->GetRadType() != 'P' || compound->GetPair(channel->GetPairNum())->GetPType() != 0) continue;
          const double g = level->GetFitGamma(ch);
          acc.AddGamma(j, la, ch, -2.0 * g * level->GetShiftDerivative(ch) / s);
          dE -= g * g * level->GetShiftSecondDerivative(ch);
        }
        acc.AddE(j, la, dE / s);
        appendRow(J / s, &acc, nullptr);
      }
    }
  }
  // Abundances: ((sum w) - 1)/s per sample with s > 0; d/dw_k = 1/s on every abundance of the sample.
  if (compound->HasAbundanceBlock()) {
    for (int smp = 1; smp <= compound->NumSamples(); smp++) {
      const double s = compound->SampleSumSigma(smp, config.abundanceSumSigma);
      if (!(s > 0.0)) continue;
      std::vector<std::pair<int, double>> direct;
      if (pmap)
        for (int k = 1; k <= compound->NumAbundanceParams(); k++) {
          int sample, nucleus;
          compound->AbundanceOf(k, sample, nucleus);
          if (sample == smp) direct.push_back({pmap->AbundanceIndex(k), 1.0 / s});
        }
      appendRow((compound->AbundanceSum(smp) - 1.0) / s, nullptr, &direct);
    }
  }
}
