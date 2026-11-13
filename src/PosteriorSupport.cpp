// The support of the posterior: the fit's limits of every parameter kind (include/PosteriorSupport.h).
#include "PosteriorSupport.h"

#include <limits>

#include "AZUREParams.h"
#include "ParameterLimitsManager.h"

PosteriorSupport PosteriorSupport::OfFit(const Config &configure, CNuc *compound, WorkingCopy *data,
                                         AZUREParams &params) {
  {
    ParameterLimitsManager settings(&configure, compound, data, &params);
    if (settings.ReadParameterSettings()) settings.ApplyAllParameterSettings(params);
  }
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

bool PosteriorSupport::Contains(const vector_r &fitVector) const {
  for (const Limit &l : limits)
    if (l.index >= (int)fitVector.size() || !(fitVector[l.index] >= l.lower && fitVector[l.index] <= l.upper)) return false;
  return true;
}
