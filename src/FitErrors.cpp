// The error analysis on a FitResult of a MIGRAD row: HESSE and MINOS (FitErrors.h).
#include "FitErrors.h"

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MinosError.h>
#include <Minuit2/MnHesse.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnUserCovariance.h>

#include "AZURECalc.h"
#include "CNuc.h"
#include "AZUREParams.h"
#include "WorkingCopy.h"

bool HasFreeParameterForHesse(const FitResult &result) {
  return result.minimum && result.minimum->UserState().VariableParameters() > 0;
}

bool RunHesse(AZURECalc &objective, FitResult &result) {
  if (!HasFreeParameterForHesse(result)) return false;
  ROOT::Minuit2::MnHesse hesse;
  try {
    hesse(objective, *result.minimum);
  } catch (...) {}
  if (!result.minimum->HasValidCovariance()) return false;
  const ROOT::Minuit2::MnUserCovariance &cov = result.minimum->UserCovariance();
  const unsigned int n = (unsigned int)result.freeIndices.size();
  if (cov.Nrow() == n) {
    result.covariance.assign((size_t)n * n, 0.0);
    for (unsigned int a = 0; a < n; a++)
      for (unsigned int b = 0; b < n; b++) result.covariance[(size_t)a * n + b] = cov(a, b);
    result.covarianceFromHesse = true;
  }
  return true;
}

void RunMinos(AZURECalc &objective, const AZUREParams &params, double up, FitResult &result, std::ostream &log) {
  objective.data()->SetErrorAnalysis(true);
  objective.SetErrorDef(up);
  ROOT::Minuit2::MnMinos minos(objective, *result.minimum);
  const ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  result.minos.clear();
  result.minosValidity.clear();
  for (int i = 0; i < (int)mp.Params().size(); i++) {
    log << "\tParameter " << i + 1 << "..." << std::endl;
    MinosValidity validity;
    if (!mp.Parameter(i).IsFixed()) {
      const ROOT::Minuit2::MinosError error = minos.Minos(i);
      result.minos.push_back(error());
      validity.lower = MinosSide{error.LowerValid(), error.AtLowerLimit(), error.AtLowerMaxFcn(), error.LowerNewMin()};
      validity.upper = MinosSide{error.UpperValid(), error.AtUpperLimit(), error.AtUpperMaxFcn(), error.UpperNewMin()};
    } else {
      result.minos.push_back(std::pair<double, double>(0., 0.));
    }
    result.minosValidity.push_back(validity);
  }
}
