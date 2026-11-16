// AZURECalc, the priors of the limits manager: the diagonal Gaussian prior a `nuisance` flag of the parameter
// limits puts on a parameter (NuisancePriorOf, with the rule that a parameter of the prior covariance has none),
// its chi-squared (CalculateNuisanceChiSquared, which the objective calls through ObjectiveContext::nuisance), its
// gradient (AddNuisanceGradient), and whether any such prior exists (HasNuisanceParameters, which switches the
// analytic G2 off).  One unit because the four must index the parameters the same way.
#include "AZURECalc.h"
#include "azure2/support/Environment.h"
#include "Objective.h"
#include "InferenceTypes.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "Config.h"
#include "CompoundSystem.h"
#include "WorkingCopy.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "AZUREGrad.h"

namespace {

// A parameter whose prior is its own (normalization, energy shift, abundance, screening, density, Doppler,
// thickness): the limits manager's nuisance flag adds nothing to it (CalculateNuisanceChiSquared skips it).  By
// kind (ParamKindRow::ownPrior), not by name.
bool HasOwnPriorKind(const ParamIndexMap &pmap, int fullIndex) {
  if (fullIndex < 0 || fullIndex >= pmap.NumFull()) return false;
  return KindRow(pmap.Desc(fullIndex).kind).ownPrior;
}

}  // namespace

bool AZURECalc::HasNuisanceParameters() const {
  // The CLI always hands AZURECalc a ParameterLimitsManager, so "a manager
  // exists" is not the question; whether it declares any nuisance parameter is.
  if (!limitsManager_) return false;
  AZUREParams fp;
  compound()->FillMnParams(fp, state_, inputs_, &configure());
  data()->FillMnParams(fp);
  // Index the parameters exactly as CalculateNuisanceChiSquared and
  // AddNuisanceGradient do, and ignore the ones they skip: a nuisance flag on a
  // norm or shift adds nothing to chi^2 or the gradient (those penalties come
  // from the segment's own norm/shift error, which G2 already includes), so it
  // must not switch off the analytic G2 or the segment-local shift derivative.
  // The GUI writes nuisance = 1 on penalized norms, so without this filter a
  // model with norm uncertainties lost both (11B+a, 2026-10-10).
  const ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), std::vector<bool>());
  int nf = 0;
  for (int i = 0; i < (int)fp.GetMinuitParams().Params().size(); i++) {
    const ROOT::Minuit2::MinuitParameter &mp = fp.GetMinuitParams().Parameter(i);
    const std::string &name = mp.GetName();
    if (mp.IsFixed() && !fp.Record(i).Row().segmentField) continue;  // the limits manager's numbering
    if (!HasOwnPriorKind(pmap, i) && limitsManager_->IsNuisanceParameterByIndex(nf)) {
      if (azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kGradDebug))
        std::cerr << "[grad] analytic G2 off: nuisance prior on '" << name << "'" << std::endl;
      return true;
    }
    nf++;
  }
  return false;
}

bool AZURECalc::CoveredByPriorCovariance(int f) const {
  const PriorCovariance *pc = data()->Data().GetPriorCovariance();
  if (!pc) return false;
  const std::vector<int> &listed = pc->FullIndex();
  return std::find(listed.begin(), listed.end(), f) != listed.end();
}

bool AZURECalc::NuisancePriorOf(int f, double &nominal, double &sigma) const {
  if (!limitsManager_) return false;
  // A parameter of the prior covariance has no diagonal prior of its own: the multivariate Gaussian replaces it
  // (docs/upgrade/objective.md D3).
  if (CoveredByPriorCovariance(f)) return false;
  return limitsManager_->NuisancePrior(f, nominal, sigma);
}

void AZURECalc::AddNuisanceGradient(const vector_r &p, std::vector<double> &grad) const {
  if (!limitsManager_) return;
  // p and grad are the full parameter vector (fixed parameters included), and
  // the prior is looked up by that same full index.
  const int n = (int)std::min(p.size(), grad.size());
  for (int f = 0; f < n; f++) {
    double nominalValue, paramError;
    if (!NuisancePriorOf(f, nominalValue, paramError)) continue;
    // d/dp of ((p - nominal)/error)^2.
    grad[f] += GaussianPriorGradient(p[f], nominalValue, paramError);
  }
}

double AZURECalc::CalculateNuisanceChiSquared(const vector_r &p) const {
  if (!limitsManager_) return 0.0;
  // p is the full parameter vector (fixed parameters included), as
  // CompoundSystem::FillCompoundFromParams reads it; the prior is looked up by the same
  // full index.  (Indexing p by the limits manager's non-fixed index, as this
  // once did, evaluates a prior on whatever parameter sits at that position as
  // soon as a fixed parameter precedes the one carrying it.)
  double nuisanceChiSquared = 0.0;
  for (int f = 0; f < (int)p.size(); f++) {
    double nominalValue, paramError;
    if (!NuisancePriorOf(f, nominalValue, paramError)) continue;
    nuisanceChiSquared += GaussianPriorChi2(p[f], nominalValue, paramError);
  }
  return nuisanceChiSquared;
}
