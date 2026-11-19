#ifndef PROGRAMRUN_H
#define PROGRAMRUN_H

// One run of the program as named stages (docs/upgrade/fit.md section 3): AZUREMain::operator()() calls them in
// order on one ProgramRun, which holds what a stage hands to the next.

#include <memory>
#include <mutex>
#include <vector>

#include "AZURECalc.h"
#include "CompoundSystemSetUp.h"
#include "AZUREParams.h"
#include "CovarianceBand.h"
#include "FitControl.h"
#include "ParameterLimitsManager.h"
#include "WorkingData.h"

class CompoundSystem;
class Config;
struct ProgramObserver;
#include "InferenceTypes.h"
class WorkingCopy;

/// The program's sink of fit progress: today's lines on standard output, byte for byte, and the intermediate files
/// (AZURECalc::WriteIterationOutput) as the checkpoint.  The observer of the objective for the whole run, so that the
/// evaluations outside the fit (HESSE, the final evaluation, the band) print as the fit's do.
class ProgramProgress {
 public:
  /// `observer` (Config::observer, may be null): its stop request is the fit's and the objective's, its progress sink
  /// gets every record after the line is printed.
  ProgramProgress(AZURECalc &objective, std::ostream &out, const ProgramObserver *observer = nullptr);
  ~ProgramProgress();
  ProgramProgress(const ProgramProgress &) = delete;
  ProgramProgress &operator=(const ProgramProgress &) = delete;
  /// The controls of the program's fit: this sink, the checkpoint, the call limit and the refreeze passes of `config`.
  FitControl Control(const Config &configure) const;

 private:
  void Print(const FitProgress &record);
  AZURECalc &objective_;
  std::ostream &out_;
  const ProgramObserver *observer_ = nullptr;
  std::mutex mutex_;  ///< the stream (its width state and buffer) is shared by every evaluating thread
  FitControl runControl_;
  EvaluationObserver evaluationObserver_;
};

struct ProgramRun {
  ProgramRun(const Config &configure, WorkingCopy *&data) : configure(configure), data(data) {}
  const Config &configure;
  /// The run's compound system, made finished and constant by the set-up (SetUpRun, SetUpCompoundSystem); every
  /// evaluation, the band, the rate and the outputs borrow it.  Declared before what borrows it, so it outlives them.
  std::unique_ptr<const CompoundSystem> compound;
  WorkingCopy *&data;  ///< AZUREMain's: replaced when an adaptive extrapolation grid is refined
  // The run's levels (LevelState.h, LevelInputs.h): made by the compound's set-up, written by the evaluations on the
  // run's data (the calculation, the final evaluation), the band and the rate, and read by the outputs.
  LevelState state;
  LevelInputs inputs;
  /// What the fit's evaluations write (AZURECalc): before the objective that borrows it, and freed with it when the
  /// data run ends (its copies of the working data are not held through the final transformation).
  std::unique_ptr<FitterStates> fitterStates{new FitterStates};
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

/// When the program's data check (Config::CHECK_DATA) is written: at the data's set-up, or after an adaptive
/// extrapolation grid has been refined (SetUpRun writes it then, on the refined grid).
enum class DataCheckAt { kDataSetUp, kAfterGridRefinement };
/// The program's set-up of the compound system and of its data into `data`, the run's (SetUpRun) and the sampler's
/// (runMCMC) alike: "Filling Compound Nucleus...", then SetUpCompoundSystem with the program's data set-up (the
/// compound check; "Filling Data Structures..." and the data or the extrapolations, or the rate's pairs; the data
/// check), and on a refusal the program's words (PrintSetUpRefusal).
CompoundSetUp SetUpProgramCompound(const Config &configure, WorkingCopy &data, DataCheckAt dataCheck);
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
