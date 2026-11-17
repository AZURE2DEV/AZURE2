// The command mcmc from its set-up on (SamplerRun.h): the start vector in the sampling basis, the priors of the
// project's <mcmc> section, the table of the free parameters, the sampler and its run.  Called by the program
// (runMCMC, src/AZURE2.cpp) and by a session (api/src/AZUREAPISampler.cpp): one definition of what a run is.
#include "SamplerRun.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <istream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>

#include "AZURECalcMCMC.h"
#include "AZUREParams.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "RunRequest.h"
#include "WorkingCopy.h"

namespace {

// The name the MCMC table prints for a parameter (the GUI's MCMC tab writes the same), from its record.
std::string PrettyParameterName(const ParameterRecord &r, const std::string &name) {
  const std::string seg = "Segment " + std::to_string(r.segmentKey);
  switch (r.kind) {
    case ParamKind::LevelEnergy: return "Level " + std::to_string(r.levelNumber) + " Energy (MeV)";
    case ParamKind::Gamma:
      return "Level " + std::to_string(r.levelNumber) + " Channel " + std::to_string(r.channel) + " Width (eV)";
    case ParamKind::Norm: return seg + " Normalization";
    case ParamKind::Thickness: return seg + " Sample Thickness (atoms/b)";
    case ParamKind::EnergyShift: return seg + " Energy Shift (keV)";
    case ParamKind::EnergyShiftSqrt: return seg + " sqrt(E) Energy-Shift Coefficient (MeV^1/2)";
    case ParamKind::EnergyShiftFn: return seg + " Energy-Shift Function a" + std::to_string(r.fnTerm);
    case ParamKind::Screening: return "Screening U_e (eV), effect " + std::to_string(r.effect);
    case ParamKind::Density: return "Target density (atoms/cm^2), effect " + std::to_string(r.effect);
    case ParamKind::Doppler: return "Doppler temperature (K), effect " + std::to_string(r.effect);
    case ParamKind::Abundance:
      return "Abundance " + (r.sampleName.empty() ? std::string() : "sample " + r.sampleName + " ") +
             "Z=" + std::to_string(r.nucleusZ) + " A=" + std::to_string(r.nucleusA);
    case ParamKind::IsospinShift:
      return "Isospin Coulomb Shift (MeV) Z=" + std::to_string(r.nucleusZ) + " A=" + std::to_string(r.nucleusA);
    case ParamKind::IsospinMix: return "Isospin-Forbidden Amplitude (MeV^1/2), record " + r.isospinRecord;
    case ParamKind::ThmCoherent: case ParamKind::Apparatus: case ParamKind::RadiativeDecay: break;  // their names say what they are
  }
  return name;
}

/*!
 * The parameter a row of the <mcmc> <parameters> block names, or -1.  The row gives the parameter's own name (whole,
 * not a prefix), or one of the GUI's short forms E<n> (energy of level n), G<n><c> (width of level n, one digit, in
 * channel c), N<k> (normalization of segment k), S<k> (constant energy shift of segment k), which are resolved through
 * the parameters' records.
 */
int FindMcmcPriorParameter(const std::string &token, const AZUREParams &params) {
  const int exact = params.IndexOfName(token);
  if (exact >= 0) return exact;
  if (token.size() < 2 || !std::isdigit((unsigned char)token[1])) return -1;
  ParamKind kind = ParamKind::LevelEnergy;
  int first = 0, second = -1;
  try {
    if (token[0] == 'E' || token[0] == 'N' || token[0] == 'S') {
      kind = token[0] == 'E' ? ParamKind::LevelEnergy : token[0] == 'N' ? ParamKind::Norm : ParamKind::EnergyShift;
      first = std::stoi(token.substr(1));
    } else if (token[0] == 'G' && token.size() > 2) {
      kind = ParamKind::Gamma;
      first = std::stoi(token.substr(1, 1));
      second = std::stoi(token.substr(2));
      if (second <= 0) return -1;
    } else {
      return -1;
    }
  } catch (...) {
    return -1;
  }
  if (first <= 0) return -1;
  for (int i = 0; i < (int)params.Records().size(); i++) {
    const ParameterRecord &r = params.Records()[i];
    if (r.kind != kind) continue;
    if ((kind == ParamKind::LevelEnergy && r.levelNumber == first) ||
        (kind == ParamKind::Gamma && r.levelNumber == first && r.channel == second) ||
        ((kind == ParamKind::Norm || kind == ParamKind::EnergyShift) && r.segmentKey == first))
      return i;
  }
  return -1;
}

}  // namespace

SamplerOutcome RunSampler(const Config &configure, WorkingCopy *data, const CompoundSystem *compound, LevelState &state,
                          LevelInputs &inputs, AZUREParams &params, const RunRequest &request,
                          const SamplerHooks &hooks, std::vector<std::vector<double>> &samples) {
  // What the sampler's walkers are evaluated on (idle levels and working copies).
  SamplerStates walkerPool;

  // Extract initial parameter values (only non-fixed parameters)
  std::vector<double> initialParams;
  if (request.reducedWidthAmplitudes) {
    // Use RWA parameters directly
    for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
      if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
        initialParams.push_back(params.GetMinuitParams().Parameter(i).Value());
      }
    }
  } else {
    // Use physical parameters
    vector_r physicalParams = compound->GetTransformParams(state, compound->TransformOut(configure, state, inputs));
    data->Data().Abundances().AppendTo(params.GetMinuitParams().Params(), physicalParams);

    // Extend with data parameters if needed
    for (size_t i = physicalParams.size(); i < params.GetMinuitParams().Params().size(); i++) {
      physicalParams.push_back(params.GetMinuitParams().Parameter(i).Value());
    }

    // Extract non-fixed parameters
    for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
      if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
        initialParams.push_back(physicalParams[i]);
      }
    }
  }

  if (initialParams.empty()) {
    configure.outStream << "Error: No free parameters found for MCMC sampling." << std::endl;
    configure.outStream << "Please ensure at least one parameter is not fixed." << std::endl;
    return SamplerOutcome::kRefused;
  }

  // Read priors from AZURE2 file if they exist
  std::map<int, std::tuple<double, double, bool>> priorMap;  // parameter index -> (priorMean, priorStd, usePrior)
  const std::unique_ptr<std::istream> priorSection = configure.SectionStream("mcmc");
  if (priorSection) {
    std::istream &priorFile = *priorSection;
    std::string line;
    bool inMCMCSection = false;
    bool inParametersSection = false;

    while (std::getline(priorFile, line)) {
      // Trim whitespace
      size_t start = line.find_first_not_of(" \t\r\n");
      if (start != std::string::npos) {
        line = line.substr(start);
        size_t end = line.find_last_not_of(" \t\r\n");
        if (end != std::string::npos) line = line.substr(0, end + 1);
      }

      if (line == "<mcmc>") {
        inMCMCSection = true;
        continue;
      } else if (line == "</mcmc>") {
        inMCMCSection = false;
        break;
      }

      if (inMCMCSection) {
        if (line == "<parameters>") {
          inParametersSection = true;
          continue;
        } else if (line == "</parameters>") {
          inParametersSection = false;
          continue;
        }

        if (inParametersSection && !line.empty()) {
          // Parse: "Parameter (raw_param_name) currentValue priorMean priorStd useGaussianPrior"
          std::istringstream iss(line);
          std::vector<std::string> parts;
          std::string part;
          while (iss >> part) parts.push_back(part);

          if (parts.size() >= 5) {
            // Last 4 parts are: currentValue, priorMean, priorStd, useGaussianPrior
            // Everything before that is the parameter name (possibly "Parameter (raw_name)")
            std::string fileParamName;
            for (size_t i = 0; i < parts.size() - 4; i++) {
              if (i > 0) fileParamName += " ";
              fileParamName += parts[i];
            }

            // Extract raw parameter name from "Parameter (raw_name)" format
            std::string rawParamName;
            size_t openParen = fileParamName.find('(');
            size_t closeParen = fileParamName.find(')');
            if (openParen != std::string::npos && closeParen != std::string::npos && closeParen > openParen) {
              // Extract the raw name from inside parentheses
              rawParamName = fileParamName.substr(openParen + 1, closeParen - openParen - 1);
            } else {
              // No parentheses - use the whole name
              rawParamName = fileParamName;
            }

            const int paramOfRow = FindMcmcPriorParameter(rawParamName, params);

            double priorMean = std::atof(parts[parts.size() - 3].c_str());
            double priorStd = std::atof(parts[parts.size() - 2].c_str());
            bool usePrior = (std::atoi(parts[parts.size() - 1].c_str()) == 1);

            if (paramOfRow < 0) {
              // A prior that is silently not applied is a wrong chain that looks right.
              configure.outStream << "ERROR: <mcmc> <parameters> row \"" << line << "\": \"" << rawParamName
                                  << "\" names no parameter of this model.  A row names a parameter by its whole "
                                  << "name as in param.sav (energy_1, width_1_2, segment_1_norm, ...), optionally "
                                  << "as \"label (name)\", or by the short forms E<n>, G<n><c>, N<k>, S<k>, "
                                  << "followed by: value, prior mean, prior sigma, use (0 or 1)." << std::endl;
              return SamplerOutcome::kRefused;
            }
            priorMap[paramOfRow] = std::make_tuple(priorMean, priorStd, usePrior);
          } else {
            configure.outStream << "ERROR: <mcmc> <parameters> row \"" << line << "\": expected a parameter name "
                                << "followed by: value, prior mean, prior sigma, use (0 or 1)." << std::endl;
            return SamplerOutcome::kRefused;
          }
        }
      }
    }
  }

  // Build prior vectors for MCMC calculator (matching GUI approach)
  std::vector<double> priorMeans;
  std::vector<double> priorStds;
  std::vector<bool> usePriors;

  // Print parameter information
  configure.outStream << std::endl;
  configure.outStream << "Free Parameters for MCMC:" << std::endl;
  configure.outStream << std::setw(40) << std::left << "Parameter Name"
                      << std::setw(15) << "Current Value"
                      << std::setw(15) << "Prior Mean"
                      << std::setw(15) << "Prior Std"
                      << std::setw(20) << "Use Prior" << std::endl;
  configure.outStream << std::string(105, '-') << std::endl;

  int paramIndex = 0;
  int priorsFound = 0;
  for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
    if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
      std::string paramName = PrettyParameterName(params.Record(i), params.GetMinuitParams().GetName(i));
      double currentValue = initialParams[paramIndex];

      // Check if prior exists for this parameter
      std::string priorMeanStr = "N/A";
      std::string priorStdStr = "N/A";
      std::string usePriorStr = "No (uniform)";

      double priorMean = 0.0;
      double priorStd = 0.0;
      bool usePrior = false;

      if (priorMap.find(i) != priorMap.end()) {
        priorMean = std::get<0>(priorMap[i]);
        priorStd = std::get<1>(priorMap[i]);
        usePrior = std::get<2>(priorMap[i]);

        // Only mark as found if usePrior is true AND std is positive
        if (usePrior && priorStd > 0.0) {
          std::ostringstream meanStream, stdStream;
          meanStream << std::scientific << std::setprecision(4) << priorMean;
          stdStream << std::scientific << std::setprecision(4) << priorStd;
          priorMeanStr = meanStream.str();
          priorStdStr = stdStream.str();
          usePriorStr = "Yes (Gaussian)";
          priorsFound++;
        } else {
          // Found in map but usePrior is false or std is invalid - reset to defaults
          priorMean = 0.0;
          priorStd = 0.0;
          usePrior = false;
        }
      }

      // Add to prior vectors - these must be added for ALL parameters to maintain index alignment
      priorMeans.push_back(priorMean);
      priorStds.push_back(priorStd);
      usePriors.push_back(usePrior);

      configure.outStream << std::setw(40) << std::left << paramName
                          << std::setw(15) << std::scientific << std::setprecision(6) << currentValue
                          << std::setw(15) << priorMeanStr
                          << std::setw(15) << priorStdStr
                          << std::setw(20) << usePriorStr << std::endl;
      paramIndex++;
    }
  }
  configure.outStream << std::endl;
  if (priorsFound > 0) {
    configure.outStream << "Found " << priorsFound << " Gaussian prior(s) in AZURE2 file." << std::endl;
    configure.outStream << "Priors will be applied during MCMC sampling." << std::endl;
  } else {
    configure.outStream << "No priors found in AZURE2 file. Using uniform priors for all parameters." << std::endl;
    configure.outStream << "To use Gaussian priors, use the GUI to configure them and save the .azr file." << std::endl;
  }
  configure.outStream << std::endl;

  // Handle overwriting existing samples.mcmc file if requested
  if (request.overwriteSamples) {
    std::string samplesFile = configure.outputdir + "samples.mcmc";
    std::ifstream checkFile(configure.OpenPath(samplesFile).c_str());
    if (checkFile.good()) {
      checkFile.close();
      if (std::remove(configure.OpenPath(samplesFile).c_str()) == 0) {
        configure.outStream << "Removed existing samples.mcmc file." << std::endl;
      } else {
        configure.outStream << "Warning: Could not remove existing samples.mcmc file." << std::endl;
      }
    }
  }

  configure.outStream << std::endl;
  configure.outStream << "MCMC Configuration:" << std::endl;
  configure.outStream << "  Number of walkers: " << request.walkers << std::endl;
  configure.outStream << "  Number of steps: " << request.steps << std::endl;
  configure.outStream << "  Initial parameter spread: " << request.parameterSpread << "%" << std::endl;
  configure.outStream << "  Level energy spread: " << request.energySpread << " keV" << std::endl;
  configure.outStream << "  Number of threads: " << request.samplerThreads << std::endl;
  configure.outStream << "  Parameter mode: " << (request.reducedWidthAmplitudes ? "RWA" : "Physical") << std::endl;
  configure.outStream << "  Number of free parameters: " << initialParams.size() << std::endl;
  configure.outStream << std::endl;

  // Create MCMC calculator
  AZURECalcMCMC mcmcCalc(data, compound, state, inputs, walkerPool, configure);

  // Set priors in MCMC calculator (same as GUI does in AZUREMCMCThread.cpp)
  if (priorsFound > 0) {
    mcmcCalc.SetPriors(priorMeans, priorStds, usePriors);
    configure.outStream << "Applied " << priorsFound << " Gaussian prior(s) to MCMC calculator." << std::endl;
    configure.outStream << std::endl;
  }

  // After every step (the sampler calls the progress callback only when an iteration callback is set): the caller's
  // hook, which may ask for the run to end after the step.
  mcmcCalc.SetIterationCallback([](int, int) {});
  mcmcCalc.SetProgressCallback([&hooks, &mcmcCalc](int step, int total, double logprob, double loglik,
                                                   double logprior) {
    if (hooks.afterStep && hooks.afterStep(step, total, logprob, loglik, logprior)) mcmcCalc.RequestStop();
  });

  // Run MCMC sampling
  const bool sampled = mcmcCalc.RunMCMCSampling(request.walkers, request.steps, initialParams,
                                                samples, request.parameterSpread, request.samplerThreads,
                                                request.reducedWidthAmplitudes, request.energySpread);
  return sampled ? SamplerOutcome::kSampled : SamplerOutcome::kNotSampled;
}
