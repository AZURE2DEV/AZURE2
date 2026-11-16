// The makers of the posterior's support from the fit's parameters (include/PosteriorSupport.h); the rule itself is
// azure2/inference/PosteriorSupport.h.
#include "PosteriorSupport.h"

#include <limits>

#include "AZUREParams.h"
#include "ParameterLimitsManager.h"

PosteriorSupport PosteriorSupportOfFit(const Config &configure, const CompoundSystem *compound, const LevelState &state,
                                       const LevelInputs &inputs, WorkingCopy *data, AZUREParams &params) {
  {
    ParameterLimitsManager settings(&configure, compound, state, inputs, data, &params);
    if (settings.ReadParameterSettings()) settings.ApplyAllParameterSettings(params);
  }
  return PosteriorSupportOfParameters(params);
}

PosteriorSupport PosteriorSupportOfParameters(const AZUREParams &params) {
  PosteriorSupport support;
  const ROOT::Minuit2::MnUserParameters &p = params.GetMinuitParams();
  for (int i = 0; i < (int)p.Params().size(); i++) {
    const ROOT::Minuit2::MinuitParameter &par = p.Parameter(i);
    if (par.IsFixed()) continue;
    support.limits.push_back({i, par.HasLowerLimit() ? par.LowerLimit() : -std::numeric_limits<double>::infinity(),
                              par.HasUpperLimit() ? par.UpperLimit() : std::numeric_limits<double>::infinity()});
  }
  return support;
}
