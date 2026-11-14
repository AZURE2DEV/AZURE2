// AZUREAPI: the program's fit on a session's objective (AZUREAPI.h, FitWithEngine; docs/upgrade/fit.md).  The same
// table of fitters and the same functions as the program's fit stage (src/ProgramFit.cpp), not a copy.
#include "AZUREAPI.h"

#include <sstream>
#include <stdexcept>

#include <Minuit2/FunctionMinimum.h>

#include "AZURECalc.h"
#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "FitErrors.h"
#include "FitEvaluations.h"
#include "Fitters.h"
#include "ParameterLimitsManager.h"
#include "WorkingCopy.h"

AZUREAPI::EngineFit AZUREAPI::FitWithEngine(const EngineFitRequest &request) {
  const FitterRow *row = FitterNamed(request.method);
  if (!row) {
    std::string names;
    for (const FitterRow &r : Fitters()) names += (names.empty() ? "" : ", ") + std::string(r.name);
    throw std::invalid_argument("unknown fit method '" + request.method + "'; the engine's fitters: " + names + ".");
  }
  if (request.minos && !row->minos)
    throw std::invalid_argument("MINOS follows a Minuit2 minimum only; '" + request.method + "' has none.");
  size_t nFree = 0;
  for (bool f : fixed_) nFree += f ? 0 : 1;
  if (!request.x0.empty() && request.x0.size() != nFree)
    throw std::invalid_argument("x0 has " + std::to_string(request.x0.size()) + " values; the session has " +
                                std::to_string(nFree) + " free parameters.");

  // The parameter vector as the program makes it (limits, fixed flags), at the session's values.
  AZUREParams params;
  compound()->FillMnParams(params, state(), inputs(), &configure());
  data()->FillMnParams(params);
  ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  for (size_t i = 0, k = 0; i < all_rwa_.size() && i < mp.Params().size(); i++) {
    double value = all_rwa_[i];
    if (!fixed_[i] && !request.x0.empty()) value = request.x0[k];
    if (!fixed_[i]) k++;
    mp.SetValue(i, value);
  }
  // <parameterSettings> as the program applies it at its start (limits; nuisance priors in the objective).
  ParameterLimitsManager limits(&configure(), compound(), state(), inputs(), data(), &params);
  limits.ReadParameterSettings();
  limits.ApplyAllParameterSettings(params);

  FitterStates states;  // the fit's pooled level states, copies of the session's
  AZURECalc objective(data(), results(), compound(), state(), inputs(), states, configure(), &limits);
  objective.SetErrorDef(1.0);
  std::vector<bool> fixedMask(mp.Params().size());
  for (size_t i = 0; i < fixedMask.size(); i++) fixedMask[i] = mp.Parameter(i).IsFixed();
  objective.SetFixedMask(fixedMask);

  std::ostringstream log;
  FitControl control;
  control.maxCalls = request.maxCalls;
  control.refreezePasses = configure().refreezeCovariance;
  control.progress = request.progress;
  control.stop = request.stop;
  control.log = &log;

  EngineFit fit;
  data()->SetFit(true);
  try {
    fit.result = RunFitter(*row, objective, params, control);
    if (row->minos && !fit.result.fellBackToMigrad && fit.result.status != "stopped") {
      // The error analysis reports its evaluations and obeys the stop request as the fit does.
      bool stopped = false;
      {
        FitEvaluations evaluations(objective, control, row->name, "hesse");
        if (HasFreeParameterForHesse(fit.result) && !RunHesse(objective, fit.result))
          log << "HESSE did not return a valid covariance; the covariance is MIGRAD's error matrix." << std::endl;
        stopped = evaluations.Stopped();
      }
      if (request.minos && !stopped) {
        FitEvaluations evaluations(objective, control, row->name, "minos");
        RunMinos(objective, params, configure().chiVariance, fit.result, log);
        objective.SetErrorDef(1.0);
        stopped = evaluations.Stopped();
      }
      if (stopped) {
        fit.result.status = "stopped";
        fit.result.converged = false;
        fit.result.covariance.clear();
        fit.result.covarianceFromHesse = false;
        fit.result.minos.clear();
        fit.result.minosValidity.clear();
      }
    }
  } catch (...) {
    data()->SetFit(false);
    data()->SetErrorAnalysis(false);
    throw;
  }
  data()->SetFit(false);
  data()->SetErrorAnalysis(false);
  // The errors as the program writes them into param.sav: Minuit's (after HESSE) for a MIGRAD row, the ones the
  // least-squares fitters set otherwise.
  const ROOT::Minuit2::MnUserParameters &end = (row->minos && fit.result.minimum && !fit.result.fellBackToMigrad &&
                                                fit.result.status != "stopped")
                                                   ? fit.result.minimum->UserParameters()
                                                   : params.GetMinuitParams();
  fit.result.errors.clear();
  for (unsigned int i = 0; i < end.Params().size(); i++) fit.result.errors.push_back(end.Error(i));
  // The end point evaluated once on the session's data, as the program's final evaluation.
  fit.objective = objective(fit.result.x);
  fit.log = log.str();
  return fit;
}
