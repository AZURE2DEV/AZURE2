// AZUREAPI: the objective as one evaluation with its breakdown, as least-squares rows, and as the ensemble sampler's
// log-probability (AZUREAPI.h; docs/upgrade/pyazr3.md section 7, docs/upgrade/objective.md section 6).
#include "AZUREAPI.h"

#include <cmath>
#include <limits>

#include "AZUREGrad.h"
#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "Objective.h"
#include "ParameterKind.h"
#include "PriorCovariance.h"
#include "WorkingCopy.h"

namespace {

// The full RWA vector: every parameter's value, the free ones from `free` in their order.
vector_r FullOf(const vector_r &free, const vector_r &all, const std::vector<bool> &fixed) {
  vector_r full = all;
  for (int i = 0, k = 0; i < (int)all.size(); ++i)
    if (!fixed[i] && k < (int)free.size()) full[i] = free[k++];
  return full;
}

}  // namespace

AZUREAPI::ObjectiveBreakdown AZUREAPI::EvaluateObjectiveRWA(const vector_r &params) {
  ObjectiveBreakdown out;
  const vector_r full = FullOf(params, all_rwa_, fixed_);
  FillFromFullRWA(full);
  EvaluateFilledChi2(nullptr);  // the model value of every point: the one model evaluation
  ObjectiveContext context{compound(), state(), inputs(), data(), results(), full, configure()};
  context.nuisance = [this](const vector_r &q) { return NuisanceChi2(q); };
  // The sum as CalculateObjectiveRWA forms it (the same steps, so the same bits), then the breakdown: each term
  // summed on its own, from the model values already in the buffer (no point is evaluated again).
  out.total = Objective::Engine().Chi2(context);
  context.termsOut = &out.terms;
  context.segmentPriorsOut = &out.segments;
  Objective::Engine().Chi2(context);
  out.excluded = OutsideSamplerSupport(params) || Objective::Engine().ExcludesForSampler(context);
  return out;
}

bool AZUREAPI::ObjectiveRowsRWA(const vector_r &params, vector_r &rows, vector_r *jacobianOut,
                                const std::vector<int> &skipColumns) {
  const vector_r full = FullOf(params, all_rwa_, fixed_);
  vector_r residuals, jacobian;
  int nCols = (int)params.size();
  if (jacobianOut) {
    // The rows and Jacobian of the residual path, as they are: data, model data, prior covariance.
    const vector_r packed = CalculateResidualJacobianRWA(params, skipColumns);
    if (packed.size() < 2) return false;
    const size_t nRes = (size_t)packed[0];
    nCols = (int)packed[1];
    residuals.assign(packed.begin() + 2, packed.begin() + 2 + nRes);
    jacobian.assign(packed.begin() + 2 + nRes, packed.end());
  } else {
    // The same rows from a forward pass (CalculateResidualsRWA's data rows, then the model data and the prior
    // covariance, in the order of CalculateResidualJacobianRWA).
    FillFromFullRWA(full);
    EvaluateFilledChi2(&residuals);
    if (data()->Data().HasCovariance()) data()->WhitenResidualRows(residuals, nullptr, 0);
    Objective::Engine().ModelDataResiduals(ObjectiveContext{compound(), state(), inputs(), data(), results(), full, configure()},
                                           residuals);
    if (const PriorCovariance *pc = data()->Data().GetPriorCovariance())
      pc->AppendRows(full, [](int) { return -1; }, 0, residuals, nullptr);
  }
  // The rows the program's least-squares drivers add (AZURECalc::ResidualJacobian, PrepareFreeParams): Park's wall,
  // the gamma_decay priors and wall, the abundance sums; then every diagonal prior and every nuisance prior.
  const ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), fixed_);
  AppendModelRows(compound(), state(), *data(), configure(), jacobianOut ? &pmap : nullptr, nCols, residuals,
                  jacobianOut ? &jacobian : nullptr);
  auto priorRow = [&](int f, double nominal, double sigma) {
    residuals.push_back((full[f] - nominal) / sigma);
    if (!jacobianOut) return;
    const size_t base = jacobian.size();
    jacobian.resize(base + (size_t)nCols, 0.0);
    const int a = pmap.FullToPacked(f);
    if (a >= 0 && a < nCols) jacobian[base + (size_t)a] = 1.0 / sigma;
  };
  ForEachDiagonalPrior(data(), pmap, [&](int f, double nominal, double sigma) {
    if (f >= 0 && f < (int)full.size()) priorRow(f, nominal, sigma);
  });
  for (int f = 0; f < (int)full.size(); f++) {
    double nominal, sigma;
    if (NuisancePriorOf(f, nominal, sigma)) priorRow(f, nominal, sigma);
  }
  rows = residuals;
  if (jacobianOut) *jacobianOut = jacobian;
  return true;
}

bool AZUREAPI::OutsideSamplerSupport(const vector_r &params) {
  ReadFitSettings();
  return !PosteriorSupport::OfParameters(*nuisanceParams_).Contains(FullOf(params, all_rwa_, fixed_));
}

double AZUREAPI::LogProbabilityRWA(const vector_r &params) {
  // The support first, as the sampler tests its prior before the likelihood: a point outside is not evaluated.
  if (OutsideSamplerSupport(params)) return -std::numeric_limits<double>::infinity();
  const double objective = CalculateObjectiveRWA(params);
  const vector_r full = FullOf(params, all_rwa_, fixed_);
  ObjectiveContext context{compound(), state(), inputs(), data(), results(), full, configure()};
  if (Objective::Engine().ExcludesForSampler(context)) return -std::numeric_limits<double>::infinity();
  return -0.5 * objective;
}
