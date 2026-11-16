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
  bool gradient;              ///< uses the analytic gradient of the objective
  bool residualRows;          ///< uses the residual rows and their analytic Jacobian
  bool bounds;                ///< honours the parameters' limits
  bool covariance;            ///< yields a covariance of the free parameters
  bool minos;                 ///< MINOS can follow it (a Minuit2 minimum)
  bool callLimit;             ///< honours FitControl::maxCalls (--max-calls, max_calls=); otherwise a notice says so
  const char *fallback;       ///< the row that runs in its place for a model without an analytic Jacobian; "" none
  const char *progressLabel;  ///< the label of its iteration records ("" when it reports evaluations only)
  const char *doc;            ///< one line for the manual, kinds() and the program's usage text
  FitFunction run;
};

/// A fitter of the engine that this build leaves out (a build option): its switch is refused with how to get it.
struct AbsentFitterRow {
  const char *name;    ///< the method's name, as its row would have it
  const char *option;  ///< the program's switch
  const char *build;   ///< what a build needs to have it
};

/// The fitters of this build, in the program's precedence: the program runs the first whose switch is given.
const std::vector<FitterRow> &Fitters();
/// The fitters of the engine that this build leaves out.
const std::vector<AbsentFitterRow> &AbsentFitters();
/// The row named `name`, nullptr for none.
const FitterRow *FitterNamed(const std::string &name);
/// The fitter the program's switches chose: the row named Config::fitter, else the row without a switch.
const FitterRow &FitterOfConfig(const Config &configure);
/*!
 * Choose the row named `name` for the program's fit (Config::fitter) unless the row already chosen comes earlier in
 * Fitters(): after any sequence of choices the earliest row of them is chosen, whatever their order (the precedence
 * of the program's switches; docs/upgrade/fit.md section 5).  Throws std::logic_error for a name that is no row.
 */
void ChooseFitter(Config &configure, const std::string &name);
/// The row of the fitter that ran for `result` of `row`: the fallback row when it fell back, else `row`.
const FitterRow &FitterThatRan(const FitterRow &row, const FitResult &result);

/*!
 * `row`'s fitter on `objective` from `params` (FitterRow::run), then the notices of the settings of `control` the
 * fitter that ran ignores (a call limit where FitterRow::callLimit is false), each also on control.log.  The program
 * and a pyazr session run a fitter through this function.
 */
FitResult RunFitter(const FitterRow &row, AZURECalc &objective, AZUREParams &params, const FitControl &control);

/// The program's usage lines of the fitters' switches (one per row with a switch) and of --max-calls.
std::string FitterSwitchesHelp();
/*!
 * The program's switch `option` of a fitter: chooses its row (ChooseFitter) and returns true; false for no fitter
 * switch.  The switch of a fitter this build leaves out throws std::invalid_argument with what, the switches
 * allowed and how to get the fitter.
 */
bool ApplyFitterSwitch(Config &configure, const std::string &option);

/// One MIGRAD from `start` (FitMigrad.cpp), with or without the analytic gradient, and its report on control.log: what
/// the MIGRAD rows repeat (the refreeze loop) and the least-squares rows fall back to.  result.minimum holds the state.
namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT
FitResult MigradPass(AZURECalc &objective, const ROOT::Minuit2::MnUserParameters &start, bool analyticGradient,
                     const FitControl &control);
/// `params` take the end of a MIGRAD `result`: Minuit2's state (errors, limits; not for a stopped fit) with the values
/// of result.x, which at a call limit can be a point Minuit2 evaluated but did not accept.
void AdoptMigradEndPoint(AZUREParams &params, const FitResult &result);

// The functions of the rows (FitMigrad.cpp, FitLeastSquares.cpp, FitGslTrustRegion.cpp, FitNLopt.cpp).
FitResult FitMigrad(AZURECalc &objective, AZUREParams &params, const FitControl &control);
FitResult FitMigradGradient(AZURECalc &objective, AZUREParams &params, const FitControl &control);
FitResult FitLevenbergMarquardt(AZURECalc &objective, AZUREParams &params, const FitControl &control);
FitResult FitGslTrustRegion(AZURECalc &objective, AZUREParams &params, const FitControl &control);
#ifdef USE_NLOPT
FitResult FitNLopt(AZURECalc &objective, AZUREParams &params, const FitControl &control);
#endif

#endif
