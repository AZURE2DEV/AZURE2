#ifndef SAMPLERRUN_H
#define SAMPLERRUN_H

// The run of the engine's ensemble sampler (AZURECalcMCMC, numcmc) from a set-up model: what the command `mcmc` does
// after its set-up, as one function for the program (runMCMC, src/AZURE2.cpp) and for a session
// (api/src/AZUREAPISampler.cpp, pyazr's Session.sample_with_engine; docs/upgrade/pyazr3.md section 12).

#include <functional>
#include <vector>

class AZUREParams;
class CompoundSystem;
class Config;
class WorkingCopy;
class LevelInputs;
class LevelState;
struct RunRequest;

/// How a sampler run ended.
enum class SamplerOutcome {
  kSampled,     ///< it ran: completed, stopped on request, or the chain of samples.mcmc was complete already
  kRefused,     ///< refused before the sampler was made: no free parameter, a bad row of <mcmc> <parameters>
  kNotSampled,  ///< the sampler refused its start or failed (AZURECalcMCMC::RunMCMCSampling returned false)
};

/// What the caller is told during a run.
struct SamplerHooks {
  /// After every step: the step and the steps of this run, and the log-probability, log-likelihood and log-prior of
  /// the step's best walker.  Returns true to end the run after this step (AZURECalcMCMC::RequestStop).  May be
  /// empty.  Called from the sampler's serial section.
  std::function<bool(int step, int steps, double logProbability, double logLikelihood, double logPrior)> afterStep;
};

/*!
 * Runs the sampler on a set-up model.  `params` is the parameter vector at the start values (FillMnParams, then the
 * user's parameter file where one is given); `state` holds those values (FillCompoundFromParams, and
 * CalcShiftFunctions in Brune's parametrization); `request` gives the command's answers: walkers, steps,
 * parameterSpread (percent), energySpread (keV), samplerThreads, reducedWidthAmplitudes (the sampling basis) and
 * overwriteSamples.
 *
 * In order: the start vector in the sampling basis (the free values of `params`, or the physical values of the
 * levels with the data parameters); the Gaussian priors of the project's <mcmc> <parameters> rows; the table of the
 * free parameters and the run's configuration on configure.outStream; samples.mcmc removed when overwriteSamples;
 * the sampler (AZURECalcMCMC on `data`, `compound` and copies of `state` and `inputs`), its priors, the run.
 * Writes samples.mcmc and walkers.mcmc, and every 10 steps the best walker's output files and param.mcmc, into
 * configure.outputdir.  `samples` receives the chain (every step of the file for a continued run).
 */
SamplerOutcome RunSampler(const Config &configure, WorkingCopy *data, const CompoundSystem *compound, LevelState &state,
                          LevelInputs &inputs, AZUREParams &params, const RunRequest &request,
                          const SamplerHooks &hooks, std::vector<std::vector<double>> &samples);

#endif
