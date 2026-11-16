// AZUREAPI, the objective's value: the session's objects filled at a full vector (FillFromFullRWA), the data
// chi-squared and its residuals on them (EvaluateFilledChi2), the terms beyond the data (FilledModelPenalties),
// the chi-squared in each basis (CalculateChi2RWA, CalculateChi2Physical, CalculateResidualsRWA), the diagonal
// priors of the limits manager (NuisancePriorOf, NuisanceChi2, AddNuisanceGradient), and the whole objective as
// the program's fit evaluates it (CalculateObjectiveRWA, ObjectiveTermNames).  One unit because the two bases and
// the two forms (a sum, a list of terms) must give one number; the derivatives are AZUREAPIGradient.cpp.
#include <algorithm>
#include "Objective.h"
#include "ParameterLimitsManager.h"
#include "AZUREAPI.h"
#include "CompoundSystemSetUp.h"
#include "AdaptiveExtrapolation.h"
#include "ChannelFunc.h"
#include "ModelPass.h"
#include "ScatteringLength.h"
#include "AZUREParams.h"

#include "GSLException.h"
#include <stdexcept>
#include "ECAmplitudeCache.h"

#include "Config.h"
#include "CompoundSystem.h"
#include "LevelReport.h"
#include "CaptureAyRequests.h"
#include "ProjectSetup.h"
#include "PPair.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "Unit.h"
#include "TargetEffect.h"
#include "AMatrixFunc.h"
#include "AZUREGrad.h"
#include "CovarianceBand.h"
#include "PriorCovarianceFile.h"
#include "AZURECalc.h"
#include "CoulFunc.h"
#include "ECIntegral.h"
#include "Profile.h"
#include "InputFiles.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <new>
#include <cmath>
#include "AZUREAPIInternal.h"

// Every segment of the data, as the chi-squared loops visit them.
static std::vector<const Segment *> AllSegments(WorkingCopy *ld) {
  std::vector<const Segment *> out;
  for (int i = 1; i <= ld->Data().NumSegments(); i++) out.push_back(ld->Data().GetSegment(i));
  return out;
}

void AZUREAPI::FillFromFullRWA(const vector_r &full) {
  const CompoundSystem *lc = compound();
  WorkingCopy *ld = data();
  lc->FillCompoundFromParams(full, state());
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, results(), lc, &state(), &inputs(), &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure(), state());
}


// The data term of AZURECalc::operator(), on the filled canonical objects.
double AZUREAPI::EvaluateFilledChi2(vector_r *res) {
  const CompoundSystem *lc = compound();
  WorkingCopy *ld = data();
  if (res) res->clear();
  double chiSquared = 0.0;
  // Rows of the segments of a THM experiment, filled once it is profiled.
  std::vector<std::vector<std::pair<int, size_t>>> groupRows(ld->Data().NumThmGroups());
  // Every point in one parallel pass first (dev 37cdd34), as the CLI's AZURECalc::PrecalculatePoints.
  AZUREProfile::ScopedReport profileReport(configure().profile, "api chi2");
  const bool pre = PrecalculateAll(*ld, AllSegments(ld), lc, state(), configure(), results());
  for (int i = 1; i <= ld->Data().NumSegments(); i++) {
    const Segment *segment = ld->Data().GetSegment(i);
    if (!segment) continue;
    azure2::EvaluateModelSegment(*ld, *segment, lc, state(), configure(), results(),
                                 pre ? azure2::LeafValues::Ready : azure2::LeafValues::Calculate);

    // A THM experiment: one profile (shared norm, background) over its
    // segments, at the last of them, as the CLI.
    const int g = ld->Data().ThmGroupOf(i);
    if (g >= 0) {
      if (res) {
        groupRows[g].push_back(std::make_pair(i, res->size()));
        for (int pid = 0; pid < segment->NumPoints(); pid++)
          if (segment->GetPoint(*ld, pid + 1)) res->push_back(0.0);
      }
      if (ld->Data().IsLastOfThmGroup(g, i)) {
        chiSquared += ld->ProfileThmGroup(g, results());
        if (res) {
          std::vector<double> r;
          for (const std::pair<int, size_t> &row : groupRows[g]) {
            ld->ThmGroupResiduals(g, row.first, results(), r);
            std::copy(r.begin(), r.end(), res->begin() + row.second);
          }
        }
      }
      continue;
    }

    // A free THM norm is the arbitrary HOES scale: profiled to its optimum
    // (n* stored in the buffer), exactly as the CLI does.
    double segChi = 0.0;
    const bool profiled = segment->IsProfiledNorm();
    if (profiled) segChi = segment->ProfileNormChiSquared(*ld, results(), results().objective.norm[i]);

    const double norm = ld->NormUsed(i, results());
    for (int pid = 0; pid < segment->NumPoints(); pid++) {
      const double r = azure2::PointResidual(*ld, *segment, pid, results(), norm).Pull();
      if (!profiled) segChi += r * r;
      if (res) res->push_back(r);
    }

    results().objective.segmentChi2[i] = segChi;
    chiSquared += segChi;
  }
  // Data covariance: r^T M^-1 r over each group in place of its points' r_i^2 (never a THM segment).
  if (ld->Data().HasCovariance()) chiSquared += ld->CovarianceChi2Correction(results());
  return chiSquared;
}

// The terms of the objective beyond the data (EvaluateFilledChi2), on the filled canonical compound: the same in
// both bases (CalculateChi2RWA and CalculateChi2Physical must give one objective).
double AZUREAPI::FilledModelPenalties() {
  // Park's J > 0 wall, the abundance priors and sum-to-one, measured scattering lengths (Objective T6, T7, T9).
  const vector_r none;
  return Objective::Engine().ModelChi2(ObjectiveContext{compound(), state(), inputs(), data(), results(), none, configure()});
}

double AZUREAPI::CalculateChi2RWA(const vector_r &rwaParams) {
  // One request at a time: operate on the canonical compound/data in place
  // (re-filled here from these parameters) instead of cloning.
  const vector_r full = MapPackedToFull(rwaParams, all_rwa_, fixed_);
  FillFromFullRWA(full);
  // Prior parameter covariance (dev b7e77292): not part of this chi2 (pyazr's objective() adds
  // penalties()["prior_covariance"]), recorded for chiSquared.out (WriteOutputFiles).
  if (const PriorCovariance *pc = data()->Data().GetPriorCovariance()) results().objective.priorCovariance = pc->Chi2(full);
  return EvaluateFilledChi2(nullptr) + FilledModelPenalties();
}

vector_r AZUREAPI::CalculateResidualsRWA(const vector_r &params) {
  FillFromFullRWA(MapPackedToFull(params, all_rwa_, fixed_));
  vector_r res;
  EvaluateFilledChi2(&res);
  // Data covariance: whitened data rows, w = W r per group, so that sum w^2 is the chi-squared.
  if (data()->Data().HasCovariance()) data()->WhitenResidualRows(res, nullptr, 0);
  return res;
}

int AZUREAPI::RefreezeCovariance() {
  WorkingCopy *ld = data();
  int n = 0;
  for (int g = 0; g < ld->Data().NumCovarianceGroups(); g++)
    if (ld->Data().GetCovarianceGroup(g).HasRelative()) n++;
  ld->RefreezeCovariance();
  return n;
}

double AZUREAPI::CalculateChi2Physical(const vector_r &physicalParams) {
  int k = 0;
  vector_r params_ = all_;
  for (int i = 0; i < all_.size(); ++i) {
    if (!fixed_[i]) {
      params_[i] = physicalParams[k];
      ++k;
    }
  }

  // One request at a time: operate on the canonical compound/data in place
  // (re-filled below from these parameters) instead of cloning.
  const CompoundSystem *localCompound = compound();
  WorkingCopy *localData = data();

  AZUREParams params;
  localCompound->FillCompoundFromParamsPhysical(params_, state(), inputs());
  localData->FillAbundancesFromVector(params_);
  bool isValid = localCompound->TransformIn(configure(), state(), inputs());

  // A physical vector the transform rejects (a width beyond what the level
  // matrix can hold) is not a model: it used to come back as chi2 = 0, the
  // global minimum of any minimiser driving this function.
  if (!isValid) return std::numeric_limits<double>::infinity();

  localCompound->FillMnParams(params, state(), inputs());
  localData->FillMnParams(params);
  localCompound->FillCompoundFromParams(params.GetMinuitParams().Params(), state());
  // The tail of the physical vector (normalizations, energy shifts, screening)
  // is the same in both bases (all_ is the transformed R-matrix block plus
  // all_rwa_'s tail); the norms were never copied out of it, so the segments
  // kept whatever norm the previous call had set (pyazr calculate_chi2 ignored
  // its norm entries, 2026-10-26; thm 783fef9, dev 7284d09).
  localData->FillNormsFromParams(params_);
  // After the compound: the sub-point grids are laid out from its fit values.
  localData->FillEnergyShiftsFromParams(params_, results(), localCompound, &state(), &inputs(), &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure(), state());

  if (const PriorCovariance *pc = localData->Data().GetPriorCovariance()) results().objective.priorCovariance = pc->Chi2(params_);
  // The same terms CalculateChi2RWA adds: the two bases must give one objective (dev 7284d09).
  return EvaluateFilledChi2(nullptr) + FilledModelPenalties();
}

bool AZUREAPI::NuisancePriorOf(int f, double &nominal, double &sigma) {
  ReadFitSettings();
  if (!nuisance_) return false;
  // A parameter of the prior covariance has no diagonal prior of its own (objective.md D3).
  if (const PriorCovariance *pc = data()->Data().GetPriorCovariance()) {
    const std::vector<int> &listed = pc->FullIndex();
    if (std::find(listed.begin(), listed.end(), f) != listed.end()) return false;
  }
  return nuisance_->NuisancePrior(f, nominal, sigma);
}

double AZUREAPI::NuisanceChi2(const vector_r &full) {
  double chi2 = 0.0;
  for (int f = 0; f < (int)full.size(); f++) {
    double nominal, sigma;
    if (!NuisancePriorOf(f, nominal, sigma)) continue;
    chi2 += GaussianPriorChi2(full[f], nominal, sigma);
  }
  return chi2;
}

void AZUREAPI::AddNuisanceGradient(const vector_r &full, vector_r &gradFull) {
  for (int f = 0; f < (int)full.size() && f < (int)gradFull.size(); f++) {
    double nominal, sigma;
    if (NuisancePriorOf(f, nominal, sigma)) gradFull[f] += GaussianPriorGradient(full[f], nominal, sigma);
  }
}

std::vector<std::string> AZUREAPI::ObjectiveTermNames() const {
  std::vector<std::string> names;
  for (const auto &t : Objective::Engine().Terms()) names.push_back(t->Name());
  return names;
}

double AZUREAPI::CalculateObjectiveRWA(const vector_r &params, std::vector<double> *terms,
                                       ObjectiveContext::SegmentPriors *segments) {
  const vector_r full = MapPackedToFull(params, all_rwa_, fixed_);
  FillFromFullRWA(full);
  EvaluateFilledChi2(nullptr);  // the model value of every point (the objective's data term recomputes the sum)
  ObjectiveContext context{compound(), state(), inputs(), data(), results(), full, configure()};
  context.termsOut = terms;
  context.segmentPriorsOut = segments;
  context.nuisance = [this](const vector_r &q) { return NuisanceChi2(q); };
  return Objective::Engine().Chi2(context);
}
