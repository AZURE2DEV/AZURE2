// The fit objective as one ordered list of terms (Objective.h; docs/upgrade/objective.md).  The expressions are those
// of AZURECalc::operator() before Phase 2, moved here unchanged, in its order.
#include "Objective.h"

#include <cmath>

#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "ESegment.h"
#include "EPoint.h"
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
          const double g = ld->ProfileThmGroup(thmGroup);
          chi2 += g;
          if (c.segmentChi2) (*c.segmentChi2)[i] = g;
        }
      } else if (segment->IsProfiledNorm()) {
        segmentChiSquared = segment->ProfileNormChiSquared();
      } else {
        for (int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
          EPoint *point = segment->GetPoint(pointIdx + 1);
          if (point) {
            double residual = point->GetFitCrossSection() - point->GetCMCrossSection() * segment->GetNorm();
            double error = point->GetCMCrossSectionError() * segment->GetNorm();
            if (error != 0.0) segmentChiSquared += (residual * residual) / (error * error);
          }
        }
        // Normalization prior (profiled norms have none).
        double dataNorm = segment->GetNorm();
        double dataNormNominal = segment->GetNominalNorm();
        double dataNormError = dataNormNominal / 100. * segment->GetNormError();
        if (dataNormError != 0.) chi2 += pow((dataNorm - dataNormNominal) / dataNormError, 2.0);
      }
      // Energy-shift priors.
      segment->AddShiftPenalty(chi2);
      if (thmGroup < 0) {
        if (c.storeTotals) segment->SetSegmentChiSquared(segmentChiSquared);
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
    if (c.data->HasCovariance()) chi2 += c.data->CovarianceChi2Correction(c.segmentChi2);
  }
};

class NuisanceTerm : public ObjectiveTerm {
 public:
  std::string Name() const override { return "nuisance priors"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (c.nuisance) chi2 += c.nuisance(c.p);
  }
};

// T6: Park formalism, every level's overlap J kept positive.
class ParkTerm : public ObjectiveTerm {
 public:
  std::string Name() const override { return "Park J > 0"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!(c.config.paramMask & Config::USE_PARK_FORMALISM)) return;
    const double parkPenalty = c.compound->ParkNormPenalty();
    chi2 += parkPenalty;
    if (c.storeTotals) c.data->SetParkPenalty(parkPenalty);
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
  std::string Name() const override { return "abundances"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!c.compound->HasAbundanceBlock()) return;
    const double abundancePenalty = c.compound->AbundancePenalty(c.config.abundanceSumSigma);
    chi2 += abundancePenalty;
    if (c.storeTotals) c.data->SetAbundancePenalty(abundancePenalty);
  }
};

// T9: measured scattering lengths (<scatteringLengths>), evaluated at each pair's threshold.
class ScatteringLengthTerm : public ObjectiveTerm {
 public:
  std::string Name() const override { return "scattering lengths"; }
  void Add(const ObjectiveContext &c, double &chi2) const override {
    if (!c.compound->HasScatteringLengthConstraints()) return;
    const double slChi2 = ScatteringLength::Chi2(c.compound, c.config, c.compound->ScatteringLengthConstraints());
    chi2 += slChi2;
    if (c.storeTotals) c.data->SetScatteringLengthChi2(slChi2);
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
    if (c.storeTotals) c.data->SetPriorCovarianceChi2(pcChi2);
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

double Objective::Chi2(const ObjectiveContext &c) const {
  double chi2 = 0.0;
  if (c.thmOnly) {
    terms_.front()->Add(c, chi2);
    return chi2;
  }
  for (const auto &t : terms_) t->Add(c, chi2);
  return chi2;
}
