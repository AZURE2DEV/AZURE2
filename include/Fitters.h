#ifndef FITTERS_H
#define FITTERS_H

// The table of the engine's fitters (docs/upgrade/fit.md section 3).  One row per fitter: its name, the program's
// switch, what it needs and yields, and its function.  The program's choice (FitterOfConfig), a pyazr session's
// fit_with_engine(method=) and pyazr.kinds()["fitters"] all read this table.  A new fitter is a row here, its
// function and its test (docs/source/developer/extending.rst, "A fitter").

#include <string>
#include <vector>

#include "FitControl.h"

class AZURECalc;
class AZUREParams;
class Config;

/// A fitter: minimises `objective` from `params` (start, limits, fixed flags); on return `params` holds the end point
/// with the errors the fitter sets.  Messages go to control.log, progress to control.progress.
using FitFunction = FitResult (*)(AZURECalc &objective, AZUREParams &params, const FitControl &control);

struct FitterRow {
  const char *name;           ///< the method's name (fit_with_engine(method=), kinds()["fitters"])
  const char *option;         ///< the program's switch that selects it; "" for the default
  unsigned int configBit;     ///< the Config::paramMask bit the switch sets; 0 for the default
  bool gradient;              ///< uses the analytic gradient of the objective
  bool residualRows;          ///< uses the residual rows and their analytic Jacobian
  bool bounds;                ///< honours the parameters' limits
  bool covariance;            ///< yields a covariance of the free parameters
  bool minos;                 ///< MINOS can follow it (a Minuit2 minimum)
  const char *progressLabel;  ///< the label of its iteration records ("" when it reports evaluations only)
  const char *doc;            ///< one line for the manual and kinds()
  FitFunction run;
};

/// The fitters of this build, in the program's precedence: the program runs the first whose switch is given.
const std::vector<FitterRow> &Fitters();
/// The row named `name`, nullptr for none.
const FitterRow *FitterNamed(const std::string &name);
/// The fitter the program's switches select (Config::paramMask).
const FitterRow &FitterOfConfig(const Config &configure);

/// One MIGRAD from `start` (FitMigrad.cpp), with or without the analytic gradient, and its report on control.log: what
/// the MIGRAD rows repeat (the refreeze loop) and the least-squares rows fall back to.  result.minimum holds the state.
namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT
FitResult MigradPass(AZURECalc &objective, const ROOT::Minuit2::MnUserParameters &start, bool analyticGradient,
                     const FitControl &control);

// The functions of the rows (FitMigrad.cpp, FitLeastSquares.cpp, FitGslTrustRegion.cpp, FitNLopt.cpp).
FitResult FitMigrad(AZURECalc &objective, AZUREParams &params, const FitControl &control);
FitResult FitMigradGradient(AZURECalc &objective, AZUREParams &params, const FitControl &control);
FitResult FitLevenbergMarquardt(AZURECalc &objective, AZUREParams &params, const FitControl &control);
FitResult FitGslTrustRegion(AZURECalc &objective, AZUREParams &params, const FitControl &control);
#ifdef USE_NLOPT
FitResult FitNLopt(AZURECalc &objective, AZUREParams &params, const FitControl &control);
#endif

#endif
