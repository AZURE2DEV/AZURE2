#ifndef AZUREAPISAMPLER_H
#define AZUREAPISAMPLER_H

// The engine's ensemble sampler on a session (docs/upgrade/pyazr3.md section 12; pyazr's
// Session.sample_with_engine): the command `mcmc` after its set-up (RunSampler, include/SamplerRun.h) on the
// session's model, and the chain read back from samples.mcmc.

#include <string>
#include <vector>

#include "RunRequest.h"
#include "SamplerRun.h"

class AZUREAPI;

/// The answers of the command `mcmc` from a session's keywords, read by the command line's own reader
/// (ReadRunCommand): its ranges, its defaults and its refusals.  `basis` is "rwa" or "physical".  Throws
/// std::invalid_argument with the reader's refusal (an answer out of its range, a build without the sampler).
RunRequest SamplerRequest(int walkers, int steps, int samplerThreads, const std::string &basis, double parameterSpread,
                          double energySpread, bool overwriteSamples);

/// What SampleWithEngine is asked.
struct EngineSampleRequest {
  RunRequest options;        ///< SamplerRequest's
  std::vector<double> x0;    ///< the free start vector (fit basis); empty: the session's values
  SamplerHooks hooks;        ///< after every step; may ask for the run to end
};

/// What it gives: how the run ended, the files, and the chain of samples.mcmc (every step of the file, those of an
/// earlier run it continued included), row-major.
struct EngineSample {
  SamplerOutcome outcome = SamplerOutcome::kNotSampled;
  std::string samplesFile, walkersFile, outputDirectory;  ///< as opened (Config::OpenPath)
  int steps = 0, walkers = 0, parameters = 0;             ///< of the file
  std::vector<double> chain;           ///< [step][walker][parameter]
  std::vector<double> logProbability;  ///< [step][walker]
  std::vector<double> logLikelihood;
  std::vector<double> logPrior;
};

/*!
 * The sampler's run on the session's model: the parameter vector as the program makes it (FillMnParams) at the
 * session's values (x0 for the free ones), a copy of the session's levels filled from it as the command fills its
 * own (FillCompoundFromParams, CalcShiftFunctions in Brune's parametrization), then RunSampler, the function the
 * command calls.  Writes what the command's sampler writes, into the session's output directory.  The session's
 * parameter values and levels are not changed.  Throws std::invalid_argument for a wrong x0 and std::runtime_error
 * for a samples.mcmc that cannot be read back.
 */
EngineSample SampleWithEngine(AZUREAPI &api, const EngineSampleRequest &request);

/// samples.mcmc read back into `sample` (steps, walkers, parameters, the chain and the three log densities): the
/// header "step,walker,logprob,loglikelihood,logprior,param0,..." and one row per walker per step, as
/// AZURECalcMCMC::RecordStep writes them.  The one reader of the file (pyazr's sample_with_engine, the GUI's sampler
/// results through azure2::front::ReadSamplerChain).  Throws std::runtime_error for a file that cannot be read.
void ReadSampleFile(const std::string &path, EngineSample &sample);

#endif
