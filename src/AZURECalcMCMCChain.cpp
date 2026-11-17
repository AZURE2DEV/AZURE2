// AZURECalcMCMC, the chain's files and the reports of a step: the accepted positions of every walker written
// once per step (RecordStep), the progress callbacks and the output files of the best walker (AfterStep,
// StartObservedLevels), and the samples of an earlier run read back (LoadExistingSamples, robust_stod).  One unit
// because these are everything the sampler writes or reads beside the walker positions, all from the sampler's
// serial section.
#include "Objective.h"
#include "AZURECalcMCMC.h"
#include "ModelPass.h"
#include "ParameterLabel.h"
#include "Config.h"
#include "CompoundSystem.h"
#include "LevelReport.h"
#include "WorkingCopy.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "PriorCovarianceFile.h"
#include "GSLException.h"
#include "ALevel.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include "AZURECalcMCMCRun.h"
#include <clocale>

// Records one completed MCMC step: every walker's accepted position, exactly
// once. Called from the sampler's serial section, so it needs no locking.
void AZURECalcMCMC::RecordStep(Run &run, int current_step, const std::vector<std::vector<double>> &positions,
                               const std::vector<double> &logps) const {
  // Absolute, 0-based step index, continuing an earlier run when resuming.
  const int absStep = run.stepOffset + current_step - 1;

  int bestWalker = -1;
  double bestLogProb = -std::numeric_limits<double>::infinity();

  for (size_t k = 0; k < positions.size(); k++) {
    // Splitting the sampled log-posterior back into likelihood and prior is
    // exact: the prior is analytic, so this reports the very values the
    // sampler used rather than a re-evaluation.
    const double logProbability = logps[k];
    const double logPrior = CalculateLogPrior(positions[k]);
    const double logLikelihood = logProbability - logPrior;

    if (logProbability > bestLogProb) {
      bestLogProb = logProbability;
      bestWalker = (int)k;
      run.lastLikelihood = logLikelihood;
      run.lastPrior = logPrior;
    }

    if (run.sampleFile && run.sampleFile->is_open()) {
      // Write with full double precision (17 digits) to preserve accuracy
      *run.sampleFile << absStep << "," << k << ","
                     << std::scientific << std::setprecision(17)
                     << logProbability << "," << logLikelihood << "," << logPrior;
      for (const double &param : positions[k]) {
        *run.sampleFile << "," << param;
      }
      *run.sampleFile << std::defaultfloat << "\n";
    }
  }

  if (run.sampleFile && run.sampleFile->is_open()) run.sampleFile->flush();  // Ensure data is written immediately

  // Best walker of the step drives the GUI readouts and the periodic output
  // files, which is more informative than whichever proposal happened to be
  // evaluated last.
  if (bestWalker >= 0) {
    run.lastParams = positions[bestWalker];
    run.lastLogprob = bestLogProb;
  }
  run.currentStep = absStep + 1;
}

// Custom robust string-to-double parser for scientific notation
double robust_stod(const std::string &str) {
  // Temporarily set numeric locale to "C"
  std::string old_locale = std::setlocale(LC_NUMERIC, nullptr);
  std::setlocale(LC_NUMERIC, "C");

  const char *cstr = str.c_str();
  char *end = nullptr;

  double value = std::strtod(cstr, &end);

  // Restore previous locale
  std::setlocale(LC_NUMERIC, old_locale.c_str());

  // Skip trailing spaces
  while (*end && std::isspace(static_cast<unsigned char>(*end))) ++end;

  if (end == cstr || *end != '\0') {
    throw std::invalid_argument("Invalid double format: '" + str + "'");
  }

  return value;
}

// After every step: the caller's progress callbacks, and every 10 steps the output files of the best walker.
ObservedLevelTable AZURECalcMCMC::StartObservedLevels() const {
  // As the set-up of the parameter vectors transforms the start vector (FillCompoundFromParams, CalcShiftFunctions,
  // TransformOut), on a copy.
  LevelState start(levels_.state);
  compound()->FillCompoundFromParams(all_rwa_, start);
  compound()->CalcShiftFunctions(configure(), start);
  return compound()->TransformOut(configure(), start, levels_.inputs);
}

void AZURECalcMCMC::AfterStep(const Run &run, int current_step, int total_steps) const {
  if (iterationCallback_) {
    // The values come from RecordStep(), which has just recorded this
    // step. An earlier version re-opened samples.mcmc and parsed its last line
    // on every iteration, which cost a file read per step and reported whichever
    // proposal was written last rather than an accepted state.
    int total_expected_samples = total_steps * run.nwalkers;

    iterationCallback_(run.currentStep, total_expected_samples);

    // For progress callback (time estimation), use the original numcmc step values
    // This keeps the time estimation logic working correctly
    if (progressCallback_) {
      progressCallback_(current_step, total_steps, run.lastLogprob, run.lastLikelihood, run.lastPrior);
    }
  }

  // Write output files every 10 iterations (moved from progress callback)
  if (current_step % 10 == 0) {
    try {
      // Get the most recent parameter set from the wrapper and reconstruct full parameter array
      if (!run.lastParams.empty()) {
        // Reconstruct full parameter array using the appropriate method based on parameter type
        vector_r fullParams;
        if (run.useRWA) {
          fullParams = ReconstructFullParameters(run.lastParams);
        } else {
          fullParams = ReconstructFullParametersPhysical(run.lastParams);
        }

        if (!fullParams.empty()) {
          const CompoundSystem *localCompound = compound();
          WalkerLevels snapshotLevels(levels_);  // this snapshot's copy of the sampler's levels
          WalkerLevels *walker = &snapshotLevels;
          WorkingCopy *localData = data()->Clone();
          azure2::Results localResults;  // this snapshot's buffer

          bool isValid = true;
          vector_r rwaFullParams;  // full RWA parameter set for param.mcmc

          if (run.useRWA) {
            // For RWA parameters, use the same approach as CalculateLogLikelihood
            localCompound->FillCompoundFromParams(fullParams, walker->state);
            localData->FillNormsFromParams(fullParams);
            localData->FillEnergyShiftsFromParams(fullParams, localResults, localCompound, &walker->state,
                                                  &walker->inputs, &configure());
            if (configure().paramMask & Config::USE_BRUNE_FORMALISM) {
              localCompound->CalcShiftFunctions(configure(), walker->state);
            }
            rwaFullParams = fullParams;  // already RWA
          } else {
            // For physical parameters, use the same approach as CalculateLogLikelihoodPhysical
            localCompound->FillCompoundFromParamsPhysical(fullParams, walker->state, walker->inputs);
            localData->FillAbundancesFromVector(fullParams);
            isValid = localCompound->TransformIn(configure(), walker->state, walker->inputs);

            if (isValid) {
              AZUREParams params;
              localCompound->FillMnParams(params, walker->state, walker->inputs, &configure());
              localData->FillMnParams(params);
              localCompound->FillCompoundFromParams(params.GetMinuitParams().Params(), walker->state);
              localData->FillNormsFromParams(fullParams);
              // After the compound: the sub-point grids are laid out from its fit values.
              localData->FillEnergyShiftsFromParams(fullParams, localResults, localCompound, &walker->state,
                                                  &walker->inputs, &configure());
              rwaFullParams = params.GetMinuitParams().Params();  // RWA params after TransformIn
              // FillMnParams read the clone's starting norms and shifts; the
              // sample's are in the tail of the physical vector, the same in
              // both bases from the first data parameter on, so param.mcmc
              // gets the sampled values rather than the start (2026-10-26).
              for (size_t i = (size_t)localData->Data().GetNormParamOffset(); i < rwaFullParams.size() && i < fullParams.size(); i++)
                rwaFullParams[i] = fullParams[i];

              if (configure().paramMask & Config::USE_BRUNE_FORMALISM) {
                localCompound->CalcShiftFunctions(configure(), walker->state);
              }
            }
          }

          if (isValid) {
            // Calculate cross sections for all data points
            for (const Segment &seg : localData->Data().GetSegments())
              if (!seg.HasComponents())
                for (Unit &point : seg.GetPoints(*localData))
                  point.Calculate(localData, localCompound, walker->state, configure(), localResults);
            // The model values (ModelPass.h): the points above are calculated; the loop skips the component
            // segments, so those of a segment that assembles its value are calculated by the model pass.
            for (int s = 1; s <= localData->Data().NumSegments(); s++) {
              const Segment *seg = localData->Data().GetSegment(s);
              azure2::EvaluateModelSegment(*localData, *seg, localCompound, walker->state, configure(), localResults,
                                           seg->NeedsAssembly() ? azure2::LeafValues::Calculate : azure2::LeafValues::Ready);
            }

            localData->WriteOutputFiles(configure(), localResults, true);

            WriteParametersOut(*localCompound, walker->state, run.useRWA ? StartObservedLevels() : localCompound->TransformOut(configure(), walker->state, walker->inputs),
                               configure());

            // Write param.mcmc snapshot: all RWA parameters (including fixed), analogous to param.sav
            if (!rwaFullParams.empty() && rwaFullParams.size() == paramNames_.size()) {
              std::string paramMCMCFile = configure().outputdir + "param.mcmc";
              std::ofstream paramOut(configure().OpenPath(paramMCMCFile));
              if (paramOut) {
                paramOut.precision(7);
                // Tagged like param.sav (the amplitudes are this run's basis);
                // untagged before 2026-11-04, see AZURECalc::WriteParameters.
                paramOut << std::setw(20) << "#parametrization" << std::scientific << std::setw(20)
                         << (double)AZUREParams::BasisForConfig(configure())
                         << std::scientific << std::setw(20) << 0.0 << std::endl;
                for (size_t i = 0; i < rwaFullParams.size(); i++) {
                  paramOut << std::setw(20) << paramNames_[i]
                           << std::scientific << std::setw(20) << rwaFullParams[i]
                           << std::scientific << std::setw(20) << 0.0 << std::endl;
                }
              }
            }
          }

          delete localData;
        }
      }
    } catch (const std::exception &e) {
      // Ignore errors in output file writing to not interrupt sampling
    }
  }
}

void AZURECalcMCMC::LoadExistingSamples(const std::string &filename, std::vector<std::vector<double>> &samples) const {
  samples.clear();

  std::ifstream file(configure().OpenPath(filename));
  if (!file.good()) {
    configure().outStream << "Warning: Could not open samples file for loading: " << filename << "\n";
    return;
  }

  std::string line;
  bool headerSkipped = false;

  while (std::getline(file, line)) {
    // Skip header line
    if (!headerSkipped && line.find("step,walker,") == 0) {
      headerSkipped = true;
      continue;
    }

    if (!headerSkipped) continue;

    // Parse CSV line: step,walker,logprob,param0,param1,...
    std::vector<double> sample;
    std::stringstream ss(line);
    std::string cell;

    int colIndex = 0;
    while (std::getline(ss, cell, ',')) {
      if (colIndex >= 3) {  // Skip step, walker, logprob columns
        try {
          double value = robust_stod(cell);
          sample.push_back(value);
        } catch (const std::exception &) {
          // Skip malformed lines
          break;
        }
      }
      colIndex++;
    }

    if (!sample.empty()) {
      samples.push_back(sample);
    }
  }

  configure().outStream << "Loaded " << samples.size() << " existing samples from " << filename << "\n";
}
