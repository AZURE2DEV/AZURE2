// The program's run: the stages of ProgramRun.h in order (docs/upgrade/fit.md section 3).
#include "AZUREMain.h"

#include "AZURECalc.h"
#include "CompoundSystemChecks.h"
#include "LevelReport.h"
#include "Config.h"
#include "ECAmplitudeCache.h"
#include "GSLException.h"
#include "ProgramRun.h"
#include "ReactionRate.h"

namespace {

// Calculate Reaction Rate
//  This uses the adaptive integration routines of GSL.  As the energy stepsize is
//  unknown until integration, AZURE is called not in segment control mode but as
//  an energy dependent function.  As every data point must be reinitialized (new energy
//  dependent terms calculated) the routine is slow.  This should be a tradeoff
//  for good accuracy.
void RunReactionRate(ProgramRun &run) {
  const Config &configure = run.configure;
  configure.outStream << "Performing reaction rate calculation..." << std::endl;
  ReactionRate reactionRate(run.compound, run.state, run.inputs, run.params.GetMinuitParams().Params(), configure,
                            configure.rateParams.entrancePair, configure.rateParams.exitPair);
  if (configure.rateParams.useFile)
    reactionRate.CalculateFileRates();
  else
    reactionRate.CalculateRates();
  reactionRate.WriteRates();
  // The run's parameters are the rate's: the final transformation (parameters.out) reads them.
  run.state = reactionRate.State();
}

}  // namespace

int AZUREMain::operator()() {
  ProgramRun run(configure(), compound_, data_);  // the run holds the levels, the compound the set-up
  if (SetUpRun(run) != 0) return -1;
  if (!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
    if (SetUpObjective(run) != 0) return -1;
    {
      const ProgramProgress progress(*run.objective, configure().outStream);
      RunFitStage(run, progress);
      FinalEvaluation(run);
      WriteRunOutputs(run);
    }
    // The objective and the parameter settings end with the data run, before the transformation out.
    run.objective.reset();
    run.limits.reset();
  } else {
    RunReactionRate(run);
  }

  configure().outStream << "Performing final parameter transformation..." << std::endl;
  ObservedLevelTable observed;
  try {
    observed = compound()->TransformOut(configure(), run.state, run.inputs);
  } catch (GSLException e) {
    configure().outStream << e.what() << std::endl;
    configure().outStream << "Problem with output transformation.  Aborting."
                          << std::endl;
    return -1;
  }
  WriteParametersOut(*compound(), run.state, observed, configure());

  // Cleanup shared EC amplitude cache
  if (configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
    CleanupECAmplitudeCache(configure());
  }

  return 0;
}
