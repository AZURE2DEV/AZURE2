#ifndef PROGRAMRUN_H
#define PROGRAMRUN_H

// One run of the program as named stages (docs/upgrade/fit.md section 3): AZUREMain::operator()() calls them in
// order on one ProgramRun, which holds what a stage hands to the next.

#include <memory>
#include <mutex>
#include <vector>

#include "AZURECalc.h"
#include "AZUREParams.h"
#include "CovarianceBand.h"
#include "FitControl.h"
#include "ParameterLimitsManager.h"
#include "WorkingData.h"

class CNuc;
class Config;
class PriorCovariance;
class WorkingCopy;

/// The program's sink of fit progress: today's lines on standard output, byte for byte, and the intermediate files
/// (AZURECalc::WriteIterationOutput) as the checkpoint.  The observer of the objective for the whole run, so that the
/// evaluations outside the fit (HESSE, the final evaluation, the band) print as the fit's do.
class ProgramProgress {
 public:
  ProgramProgress(AZURECalc &objective, std::ostream &out);
  ~ProgramProgress();
  ProgramProgress(const ProgramProgress &) = delete;
  ProgramProgress &operator=(const ProgramProgress &) = delete;
  /// The controls of the program's fit: this sink, the checkpoint, the call limit and the refreeze passes of `config`.
  FitControl Control(const Config &configure) const;

 private:
  void Print(const FitProgress &record);
  AZURECalc &objective_;
  std::ostream &out_;
  std::mutex mutex_;  ///< the stream (its width state and buffer) is shared by every evaluating thread
  FitControl runControl_;
  EvaluationObserver observer_;
};

struct ProgramRun {
  ProgramRun(const Config &configure, CNuc *compound, WorkingCopy *&data) :
    configure(configure), compound(compound), data(data) {}
  const Config &configure;
  CNuc *compound;  ///< AZUREMain's: written by the set-up (SetUpRun) only, read as constant after it
  WorkingCopy *&data;  ///< AZUREMain's: replaced when an adaptive extrapolation grid is refined
  // The run's levels (LevelState.h, LevelInputs.h): made by the compound's set-up, written by the evaluations on the
  // run's data (the calculation, the final evaluation), the band and the rate, and read by the outputs.
  LevelState state;
  LevelInputs inputs;
  FitterStates fitterStates;  ///< the fit's pooled level states (AZURECalc), before the objective that borrows them
  AZUREParams params;
  std::shared_ptr<const PriorCovariance> priorCov;
  bool refinesAdaptiveGrids = false;
  // The data runs (SetUpObjective):
  std::unique_ptr<ParameterLimitsManager> limits;
  azure2::Results results;  ///< the buffer of the evaluations on data (the calculation, the final evaluation)
  std::unique_ptr<AZURECalc> objective;
  std::vector<bool> fitFixedMask;  ///< the fixed flags the minimizers and the band honour
  // The band covariance of this run's fit (RunFitStage), saved and used by WriteRunOutputs:
  BandCovariance bandCov;
  bool haveBandCov = false;
};

/// Fill the compound nucleus and the data, initialize, make the parameter vector (and read a parameter file), refine
/// adaptive extrapolation grids, match the prior covariance, write param.par.  0, or -1 with the reason printed.
int SetUpRun(ProgramRun &run);
/// The data runs: the EC amplitude cache, the data's initialization, <parameterSettings>, the objective.  0 or -1.
int SetUpObjective(ProgramRun &run);
/// The fit (Config::PERFORM_FIT): the fitter of the program's switches, HESSE, MINOS, fitStatus.out,
/// covariance_matrix.out, param.sav; else the calculation's line.
void RunFitStage(ProgramRun &run, const ProgramProgress &progress);
/// The evaluation at the run's parameters and, with data, the chi^2 per segment and the total on standard output.
void FinalEvaluation(ProgramRun &run);
/// The band (covariance.dat, BuildBandData), the output files, scattering lengths and <modelData>.
void WriteRunOutputs(ProgramRun &run);
/// Name the parameters sitting at a limit of `params` (the usual reason a fit has no covariance).
void ReportParametersAtLimits(const ProgramRun &run, const AZUREParams &params);
/// The band covariance from a Minuit covariance (packed lower triangle over the free parameters).
void CaptureMinuitCovariance(ProgramRun &run, const std::vector<double> &covData, const AZUREParams &atLimits);

#endif
