#ifndef FITCONTROL_H
#define FITCONTROL_H

// What a fit is given and what it returns (docs/upgrade/fit.md section 3): the progress record, the controls a caller
// hands to a fitter (call limit, progress sink, checkpoint, stop request, log), and the FitResult record.  The same
// records serve the program's fit (src/ProgramFit.cpp) and a pyazr session's (api/src/AZUREAPIFit.cpp).

#include <atomic>
#include <functional>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "CovarianceBand.h"

namespace ROOT {
namespace Minuit2 {
class FunctionMinimum;
}
}  // namespace ROOT

/// One step of a fit, handed to the caller's sink.  An evaluation record comes from every evaluation of the objective
/// (AZURECalc::operator()); an iteration record from a fitter that counts its own steps (LM).
struct FitProgress {
  enum class Kind { kEvaluation, kIteration };
  Kind kind = Kind::kEvaluation;
  const char *fitter = "";  ///< the row's name (Fitters.h); "" for an evaluation outside a fit
  const char *label = "";   ///< an iteration record: the row's progressLabel ("LM iteration")
  const char *stage = "fit";  ///< what runs: "fit", or the error analysis after it, "hesse" or "minos"
  long evaluation = 0;      ///< an evaluation record: the run's evaluation number (WorkingCopy::Iterate)
  long iteration = 0;       ///< an iteration record: the fitter's iteration number
  long calls = 0;           ///< objective evaluations of this fit so far
  double objective = 0.0;   ///< the objective (chi^2 plus priors) at this record's point
  double best = std::numeric_limits<double>::quiet_NaN();     ///< lowest objective of this fit so far
  double measure = std::numeric_limits<double>::quiet_NaN();  ///< the fitter's own measure, NaN where it has none
  bool errorAnalysis = false;  ///< the evaluation belongs to MINOS (WorkingCopy::IsErrorAnalysis)
  const std::vector<double> *x = nullptr;  ///< the record's full parameter vector; valid during the call only
};

/// The sink a caller gives.  Evaluation records come from several threads at once: a sink must be thread-safe.
using FitProgressSink = std::function<void(const FitProgress &)>;
/// The program's intermediate output files at a full parameter vector (AZURECalc::WriteIterationOutput).
using FitCheckpoint = std::function<void(const std::vector<double> &)>;

/// What AZURECalc::operator() reports to the one who runs it: every evaluation, and the stop request it obeys.
/// Owned by the caller; AZURECalc holds a pointer (AZURECalc::Observe).
struct EvaluationObserver {
  /// (evaluation number, objective, full parameter vector, in MINOS): called after every evaluation; thread-safe.
  std::function<void(long, double, const std::vector<double> &, bool)> evaluated;
  /// The caller's stop request: when set, an evaluation inside a fit returns 0 at once (and the gradient 0).
  const std::atomic<bool> *stop = nullptr;
};

/// The controls of one fit (FitterRow::run).
struct FitControl {
  unsigned int maxCalls = 0;                 ///< the call limit of a fitter that honours one (FitterRow::callLimit);
                                             ///< 0: none
  std::string callLimitSource;               ///< appended to the limit in messages: " (--max-calls)" in the program
  int refreezePasses = 0;                    ///< MIGRAD reruns with a re-frozen relvector data covariance
  FitProgressSink progress;                  ///< every progress record; may be empty
  FitCheckpoint checkpoint;                  ///< the program's intermediate files; empty in a session
  const std::atomic<bool> *stop = nullptr;   ///< the caller's stop request; nullptr: none
  std::ostream *log = nullptr;               ///< the fitter's printed lines (the program: standard output)
};

/// One side of a MINOS interval as Minuit2's MinosError reports it (MnCross).
struct MinosSide {
  bool valid = false;       ///< the crossing of objective minimum + Up was found
  bool atLimit = false;     ///< the scan reached the parameter's limit before the crossing
  bool callLimit = false;   ///< the scan ran out of calls
  bool newMinimum = false;  ///< the scan found a point below the minimum: the fit was not at its minimum
};

/// The validity of one parameter's MINOS interval; both sides default (invalid) for a fixed parameter.
struct MinosValidity {
  MinosSide lower;
  MinosSide upper;
};

/// How a fit ended (FitterRow::run) and what its error analysis added (FitErrors.h).
struct FitResult {
  std::string fitter;                 ///< the row's name (the one asked; fellBackToMigrad says MIGRAD ran)
  std::vector<double> x;              ///< full parameter vector at the end (fit basis, Minuit order): the lowest
                                      ///< objective the fitter reached (at a call limit: the lowest it evaluated)
  double objective = std::numeric_limits<double>::quiet_NaN();  ///< the fitter's value of the objective at x
  std::string status;                 ///< the fitter's words: converged, call_limit, invalid, stopped, no_decrease,
                                      ///< iteration_limit, converged_<reason> (NLopt), failed_<reason>; FitOutcome
  bool converged = false;             ///< a minimum by the fitter's own criterion
  double edm = std::numeric_limits<double>::quiet_NaN();      ///< MIGRAD's estimated distance to the minimum
  double edmGoal = std::numeric_limits<double>::quiet_NaN();  ///< MIGRAD's EDM goal (0.002 x tolerance x Up)
  long calls = 0;                     ///< objective evaluations (MIGRAD, NLopt); LM: residual passes
  std::vector<int> freeIndices;       ///< full indices of the rows and columns of covariance
  std::vector<double> covariance;     ///< covariance over freeIndices, row-major; empty: none
  bool covarianceFromHesse = false;   ///< MIGRAD's matrix refined by HESSE (RunHesse)
  std::vector<std::pair<double, double>> minos;  ///< MINOS (lower, upper) per full index, (0, 0) for fixed; RunMinos
  std::vector<MinosValidity> minosValidity;      ///< per full index, beside minos; RunMinos
  std::vector<double> errors;         ///< per full index, the errors param.sav gets (Minuit's, or sqrt diag of LM's)
  std::vector<std::string> notices;   ///< what the caller should know of this fit, one sentence each: a fallback to
                                      ///< another fitter, a setting the fitter ignored, a best point kept at a limit
  bool fellBackToMigrad = false;      ///< LM: no analytic Jacobian, MIGRAD (numerical) ran instead
  std::shared_ptr<ROOT::Minuit2::FunctionMinimum> minimum;  ///< MIGRAD's state, for HESSE and MINOS
  BandCovariance leastSquaresBand;    ///< LM: the band block of (J^T J + priors)^-1; empty: none
};

/// The outcome of `result` in the words every fitter shares (fitStatus.out's status): converged (FitResult::converged),
/// call_limit, iteration_limit, stopped, invalid, or failed (any failed_<reason>).
std::string FitOutcome(const FitResult &result);

#endif
