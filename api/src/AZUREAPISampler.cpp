// The engine's ensemble sampler on a session (AZUREAPISampler.h; docs/upgrade/pyazr3.md section 12).  The run is
// RunSampler (src/SamplerRun.cpp), the function the command `mcmc` calls: not a copy.
#include "AZUREAPISampler.h"

#include <algorithm>
#include <fstream>
#include <locale>
#include <sstream>
#include <stdexcept>

#include "AZUREAPI.h"
#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "LevelInputs.h"
#include "LevelState.h"
#include "WorkingCopy.h"
#include "azure2/inference/LevelParameterLayout.h"

namespace {

/// A number as the command line gives it, to the last digit.
std::string Written(double x) {
  std::ostringstream s;
  s.imbue(std::locale::classic());
  s.precision(17);
  s << x;
  return s.str();
}

/// The values of a parameter vector made from `inputs` as the sampler makes its own (FillMnParams).
AZUREParams VectorOfInputs(AZUREAPI &api, const LevelState &state, const LevelInputs &inputs) {
  AZUREParams seen;
  api.compound()->FillMnParams(seen, state, inputs, &api.configure());
  api.data()->FillMnParams(seen);
  return seen;
}

/*!
 * The sampler makes its own parameter vector from the levels' inputs and the data (AZURECalcMCMC::
 * UpdateParameterVectors: FillMnParams; the program's run then reads the user's parameter file again).  From that
 * vector it takes the fixed values and the start levels it reports; the free start is given to it.  A session's
 * values that are not the project's (a parameter file given when it was made, x0) are therefore put into `inputs`,
 * the copy handed to the sampler, along the layout of the level block; the vector the sampler will make is then
 * checked against `wanted`.  Throws std::invalid_argument where the values cannot be carried: a project with
 * gamma_decay records whose level values differ from the project's, and a fixed data parameter (a normalization,
 * an energy shift) that differs from the project's, which the sampler reads from the data.
 */
void CarryValuesToInputs(AZUREAPI &api, const ROOT::Minuit2::MnUserParameters &wanted, const LevelState &state,
                         LevelInputs &inputs) {
  const int levelBlock = api.compound()->ParameterLayout().Size();
  const auto firstDifference = [&](const AZUREParams &seen, bool levels) {
    const ROOT::Minuit2::MnUserParameters &got = seen.GetMinuitParams();
    for (int i = 0; i < (int)wanted.Params().size() && i < (int)got.Params().size(); i++) {
      if ((i < levelBlock) != levels) continue;
      // the free start of a data parameter is given to the sampler; its fixed value is read from the data
      if (!levels && !wanted.Parameter(i).IsFixed()) continue;
      if (got.Parameter(i).Value() != wanted.Parameter(i).Value()) return i;
    }
    return -1;
  };
  const AZUREParams asProject = VectorOfInputs(api, state, inputs);
  const int data = firstDifference(asProject, false);
  if (data >= 0)
    throw std::invalid_argument("the engine's sampler reads the fixed data parameters from the project, and the "
                                "session's " + wanted.GetName(data) + " (from its parameter file) is not the "
                                "project's value; allowed: fixed normalizations and energy shifts as the project "
                                "has them; write the value into the project (Project) and open a session of it.");
  if (firstDifference(asProject, true) < 0) return;  // the session's level values are the project's
  if (!inputs.GammaDecayNumbers().empty())
    throw std::invalid_argument("the engine's sampler starts a project with gamma_decay records from the project's "
                                "level values only; allowed: a session without a parameter file and without x0; "
                                "write the start into the project (Project) and open a session of it.");
  for (const azure2::inference::LevelParameterLayout::Level &slots : api.compound()->ParameterLayout().Levels()) {
    inputs.SetE(slots.group, slots.level, wanted.Parameter(slots.energy).Value());
    for (int ch = 1; ch <= (int)slots.amplitude.size(); ch++)
      inputs.SetGamma(slots.group, slots.level, ch, wanted.Parameter(slots.amplitude[ch - 1]).Value());
  }
  const AZUREParams carried = VectorOfInputs(api, state, inputs);
  const int level = firstDifference(carried, true);
  bool sameFixed = carried.GetMinuitParams().Params().size() == wanted.Params().size();
  for (int i = 0; sameFixed && i < (int)wanted.Params().size(); i++)
    sameFixed = carried.GetMinuitParams().Parameter(i).IsFixed() == wanted.Parameter(i).IsFixed();
  if (level >= 0 || !sameFixed)
    throw std::invalid_argument("the session's value of " + (level >= 0 ? wanted.GetName(level) : "a parameter") +
                                " cannot be given to the engine's sampler (a start that changes which parameters "
                                "are fixed, such as a free amplitude at 0); allowed: the project's values, or a "
                                "start with every free amplitude other than 0.");
}

}  // namespace

// samples.mcmc into `sample` (AZUREAPISampler.h): the header "step,walker,logprob,loglikelihood,logprior,param0,..." and one row per
// walker per step, as AZURECalcMCMC::RecordStep writes them (17 digits: the doubles the sampler held).
void ReadSampleFile(const std::string &path, EngineSample &sample) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error("the sampler's chain " + path + " cannot be opened.");
  std::string line;
  if (!std::getline(file, line) || line.compare(0, 42, "step,walker,logprob,loglikelihood,logprior") != 0)
    throw std::runtime_error(path + " does not start with the header of samples.mcmc.");
  const int columns = 1 + static_cast<int>(std::count(line.begin(), line.end(), ','));
  const int parameters = columns - 5;
  std::vector<int> steps, walkers;
  std::vector<double> values;  // per row: logprob, loglikelihood, logprior, the parameters
  int lineNumber = 1;
  while (std::getline(file, line)) {
    lineNumber++;
    if (line.empty()) continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream in(line);
    in.imbue(std::locale::classic());
    int step = 0, walker = 0;
    in >> step >> walker;
    for (int k = 0; k < 3 + parameters; k++) {
      double v = 0.0;
      in >> v;
      values.push_back(v);
    }
    std::string rest;
    if (!in || (in >> rest) || step < 0 || walker < 0)
      throw std::runtime_error(path + ", line " + std::to_string(lineNumber) + ": not a row of " +
                               std::to_string(columns) + " numbers.");
    steps.push_back(step);
    walkers.push_back(walker);
  }
  const int rows = static_cast<int>(steps.size());
  const int nWalkers = rows ? 1 + *std::max_element(walkers.begin(), walkers.end()) : 0;
  const int nSteps = rows ? 1 + *std::max_element(steps.begin(), steps.end()) : 0;
  if (rows != nWalkers * nSteps)
    throw std::runtime_error(path + " holds " + std::to_string(rows) + " rows, not " + std::to_string(nSteps) +
                             " steps of " + std::to_string(nWalkers) + " walkers.");
  sample.steps = nSteps;
  sample.walkers = nWalkers;
  sample.parameters = parameters;
  sample.chain.assign(static_cast<size_t>(rows) * parameters, 0.0);
  sample.logProbability.assign(rows, 0.0);
  sample.logLikelihood.assign(rows, 0.0);
  sample.logPrior.assign(rows, 0.0);
  std::vector<char> seen(rows, 0);
  const size_t width = 3 + parameters;
  for (int r = 0; r < rows; r++) {
    const size_t at = static_cast<size_t>(steps[r]) * nWalkers + walkers[r];
    if (seen[at])
      throw std::runtime_error(path + ": step " + std::to_string(steps[r]) + ", walker " + std::to_string(walkers[r]) +
                               " is written twice.");
    seen[at] = 1;
    const double *row = &values[r * width];
    sample.logProbability[at] = row[0];
    sample.logLikelihood[at] = row[1];
    sample.logPrior[at] = row[2];
    std::copy(row + 3, row + width, sample.chain.begin() + at * parameters);
  }
}

RunRequest SamplerRequest(int walkers, int steps, int samplerThreads, const std::string &basis, double parameterSpread,
                          double energySpread, bool overwriteSamples) {
  std::vector<std::string> arguments = {"mcmc",
                                        "--walkers=" + std::to_string(walkers),
                                        "--steps=" + std::to_string(steps),
                                        "--sampler-threads=" + std::to_string(samplerThreads),
                                        "--basis=" + basis,
                                        "--parameter-spread=" + Written(parameterSpread),
                                        "--energy-spread=" + Written(energySpread)};
  if (overwriteSamples) arguments.push_back("--overwrite-samples");
  arguments.push_back("project.azr");  // the reader asks for a project; the session has its own
  RunRequest request;
  std::vector<std::string> rest;
  ReadRunCommand(arguments, request, rest);
  return request;
}

EngineSample SampleWithEngine(AZUREAPI &api, const EngineSampleRequest &request) {
  const std::vector<bool> fixed = api.params_fixed();
  const vector_r values = api.params_all_rwa();
  size_t nFree = 0;
  for (bool f : fixed) nFree += f ? 0 : 1;
  if (!request.x0.empty() && request.x0.size() != nFree)
    throw std::invalid_argument("x0 has " + std::to_string(request.x0.size()) + " values; the session has " +
                                std::to_string(nFree) + " free parameters.");
  const Config &configure = api.configure();

  // The parameter vector as the program makes it (limits, fixed flags), at the session's values.
  AZUREParams params;
  api.compound()->FillMnParams(params, api.state(), api.inputs(), &configure);
  api.data()->FillMnParams(params);
  ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  for (size_t i = 0, k = 0; i < values.size() && i < mp.Params().size(); i++) {
    double value = values[i];
    if (!fixed[i] && !request.x0.empty()) value = request.x0[k];
    if (!fixed[i]) k++;
    mp.SetValue(i, value);
  }
  // The levels at those values, as the command fills its own before the sampler copies them: on copies, so that the
  // session's levels stay what they are.
  LevelState state(api.state());
  LevelInputs inputs(api.inputs());
  CarryValuesToInputs(api, mp, state, inputs);
  api.compound()->FillCompoundFromParams(mp.Params(), state);
  if (configure.paramMask & Config::USE_BRUNE_FORMALISM) api.compound()->CalcShiftFunctions(configure, state);

  EngineSample sample;
  sample.outputDirectory = configure.OpenPath(configure.outputdir);
  sample.samplesFile = configure.OpenPath(configure.outputdir + "samples.mcmc");
  sample.walkersFile = configure.OpenPath(configure.outputdir + "walkers.mcmc");
  std::vector<std::vector<double>> samples;
  sample.outcome = RunSampler(configure, api.data(), api.compound(), state, inputs, params, request.options,
                              request.hooks, samples);
  if (sample.outcome == SamplerOutcome::kSampled) ReadSampleFile(sample.samplesFile, sample);
  return sample;
}
