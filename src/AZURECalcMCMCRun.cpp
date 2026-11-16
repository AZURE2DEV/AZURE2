// AZURECalcMCMC, the run of the chain: RunMCMCSampling (the start positions, a resumed chain, the ensemble
// sampler nu::Mcmc with its team of threads, the final walker positions), the walker-position files that let a
// resumed run continue the same ensemble, the check of the start vector against Park's J > 0
// (ParkNonPositiveLevels), and the caller's callbacks.  One unit because this is the one place the sampler is
// driven; the state of a run is AZURECalcMCMC::Run (include/AZURECalcMCMCRun.h), owned by the call.
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
#ifdef _OPENMP
#include <omp.h>
#endif

/* Walker-position persistence.
 *
 * The chain file records where the walkers have been; these two functions
 * record where they *are*, so that resuming continues the same ensemble. Without
 * them a resumed run re-scatters the walkers around the starting parameters and
 * splices a fresh burn-in into the middle of the chain file, which quietly
 * corrupts any posterior summary taken over the whole file.
 */
static bool load_walker_state(const std::string &filename, int nwalkers, int ndim,
                              std::vector<std::vector<double>> &positions) {
  std::ifstream file(filename);
  if (!file.is_open()) return false;

  std::vector<std::vector<double>> loaded;
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string token;
    std::vector<double> row;
    while (std::getline(ss, token, ',')) {
      try {
        row.push_back(robust_stod(token));
      } catch (...) {
        return false;
      }
    }
    if ((int)row.size() != ndim) return false;
    loaded.push_back(row);
  }

  // A state saved for a different walker count or model cannot be reused.
  if ((int)loaded.size() != nwalkers) return false;

  positions = loaded;
  return true;
}

static bool save_walker_state(const std::string &filename,
                              const std::vector<std::vector<double>> &positions) {
  std::ofstream file(filename);
  if (!file.is_open()) return false;

  for (size_t k = 0; k < positions.size(); k++) {
    for (size_t j = 0; j < positions[k].size(); j++) {
      // Full precision: a rounded restart position is a different model.
      file << std::scientific << std::setprecision(17) << positions[k][j];
      if (j + 1 != positions[k].size()) file << ",";
    }
    file << "\n";
  }
  return true;
}

void AZURECalcMCMC::SetProgressCallback(ProgressCallback callback) { progressCallback_ = std::move(callback); }

void AZURECalcMCMC::SetIterationCallback(IterationCallback callback) { iterationCallback_ = std::move(callback); }

std::vector<std::string> AZURECalcMCMC::ParkNonPositiveLevels(const std::vector<double> &p, bool useRWA,
                                                             bool &filled) const {
  std::vector<std::string> labels;
  filled = false;
  if (!parametersInitialized_) return labels;
  const CompoundSystem *lc = compound();
  WalkerLevels *walker = GetPooledLevels();
  try {
    if (useRWA) {
      lc->FillCompoundFromParams(ReconstructFullParameters(p), walker->state);
      filled = true;
    } else {
      // As CalculateLogLikelihoodPhysical: observed widths -> Park's amplitudes.
      lc->FillCompoundFromParamsPhysical(ReconstructFullParametersPhysical(p), walker->state, walker->inputs);
      if (lc->TransformIn(configure(), walker->state, walker->inputs)) {
        AZUREParams params;
        lc->FillMnParams(params, walker->state, walker->inputs, &configure());
        lc->FillCompoundFromParams(params.GetMinuitParams().Params(), walker->state);
        filled = true;
      }
    }
    if (filled) {
      lc->CalcShiftFunctions(configure(), walker->state);
      for (int j = 1; j <= lc->NumJGroups(); j++) {
        const JGroup *jg = lc->GetJGroup(j);
        if (!jg->IsInRMatrix()) continue;
        for (int la = 1; la <= jg->NumLevels(); la++) {
          const ALevel *level = jg->GetLevel(la);
          if (!level->IsInRMatrix() || walker->state.ParkNorm(j, la) > 0.0) continue;
          std::ostringstream text;
          text << AZURELabel::Level(jg, walker->inputs.E(j, la), j, la) << ": J = " << walker->state.ParkNorm(j, la);
          labels.push_back(text.str());
        }
      }
    }
  } catch (...) {
    filled = false;
    labels.clear();
  }
  ReturnPooledLevels(walker);
  return labels;
}

bool AZURECalcMCMC::RunMCMCSampling(int nwalkers, int nsteps, const std::vector<double> &initialParams,
                                    std::vector<std::vector<double>> &samples, double chainSpreadPercent, int nthreads, bool useRWA,
                                    double energySpreadKeV) {
  try {
    int ndim = initialParams.size();

    if (!PrepareRun(initialParams, useRWA)) {
      configure().outStream << "MCMC sampling not started." << std::endl;
      return false;
    }

    if (ndim == 0) {
      configure().outStream << "Error: No parameters provided for MCMC sampling\n";
      return false;
    }

    // Set up CSV file for saving samples with enhanced format
    std::string samplesFile = configure().outputdir + "samples.mcmc";
    std::ofstream csvFile;
    bool resuming = false;
    int startStep = 0;

    // Check if samples file exists for resuming
    std::ifstream existingFile(configure().OpenPath(samplesFile));
    if (existingFile.good()) {
      configure().outStream << "Found existing samples file: " << samplesFile << "\n";
      configure().outStream << "Counting existing samples for resume capability...\n";

      std::string line;
      int lineCount = 0;
      bool hasHeader = false;
      while (std::getline(existingFile, line)) {
        if (!hasHeader && line.find("step,walker,logprob,loglikelihood,logprior") == 0) {
          hasHeader = true;
          continue;
        }
        if (hasHeader) lineCount++;
      }
      existingFile.close();

      if (lineCount > 0) {
        int samplesPerStep = nwalkers;
        startStep = lineCount / samplesPerStep;

        configure().outStream << "Found " << lineCount << " existing samples (" << startStep << " steps)\n";

        if (startStep >= nsteps) {
          configure().outStream << "Sampling already completed! Loading existing samples.\n";
          LoadExistingSamples(samplesFile, samples);
          return true;
        }

        resuming = true;
        configure().outStream << "Resuming from step " << startStep << "\n";
        csvFile.open(configure().OpenPath(samplesFile), std::ios::app);
      }
    }

    if (!resuming) {
      // Create new CSV file with enhanced header
      csvFile.open(configure().OpenPath(samplesFile), std::ios::out);
      csvFile << "step,walker,logprob,loglikelihood,logprior";
      for (int i = 0; i < ndim; i++) {
        csvFile << ",param" << i;
      }
      csvFile << "\n";
      configure().outStream << "Created new samples file: " << samplesFile << "\n";
    }

    configure().outStream << "Initializing numcmc with " << ndim << " parameters, "
                          << nwalkers << " walkers, " << (nsteps - startStep) << " remaining steps\n";

    // Initialize walkers in a ball around the initial parameters
    std::vector<std::vector<double>> initial_positions;
    std::default_random_engine generator;
    double spreadFraction = chainSpreadPercent / 100.0;

    // Per-parameter initial spread. A single percentage cannot serve every
    // parameter kind: 1% of a 5 MeV level energy is 50 keV, which scatters the
    // walkers across unrelated resonance structures and the ensemble never
    // contracts. Level energies therefore get an absolute spread in keV, and
    // norms/shifts are held inside their own prior width.
    const double energySpreadMeV = std::max(0.0, energySpreadKeV) / 1000.0;
    std::vector<double> paramSpread(ndim, 0.0);
    int nEnergyCapped = 0, nPriorCapped = 0;

    for (int j = 0; j < ndim; j++) {
      double relative = fabs(initialParams[j]) * spreadFraction;
      double spread = relative;

      const bool known = j < (int)freeKinds_.size();

      if (known && freeKinds_[j] == ParamKind::LevelEnergy) {
        spread = energySpreadMeV;
        nEnergyCapped++;
      } else if (known && KindRow(freeKinds_[j]).ownPrior) {
        // Starting outside the prior only wastes steps walking back into it.
        if (j < (int)priorStds_.size() && j < (int)usePriors_.size() &&
            usePriors_[j] && priorStds_[j] > 0.0 && priorStds_[j] < spread) {
          spread = priorStds_[j];
          nPriorCapped++;
        }
      }

      if (spread <= 0.0) spread = 1e-6;
      paramSpread[j] = spread;
    }

    // Resuming continues the ensemble where it stopped, when that state is
    // available and matches this model; otherwise fall back to a fresh ball.
    const std::string walkerStateFile = configure().outputdir + "walkers.mcmc";
    if (resuming && load_walker_state(configure().OpenPath(walkerStateFile), nwalkers, ndim, initial_positions)) {
      configure().outStream << "Resumed walker positions from " << walkerStateFile << "\n";
    } else {
      if (resuming) {
        configure().outStream << "Warning: no usable walker state in " << walkerStateFile
                              << "; the walkers restart from a fresh spread, so the "
                              << "steps that follow are a new burn-in rather than a "
                              << "continuation of the earlier chain.\n";
      }
      if (nEnergyCapped > 0) {
        configure().outStream << "Initial spread: " << nEnergyCapped
                              << " level-energy parameter" << (nEnergyCapped == 1 ? "" : "s")
                              << " limited to " << energySpreadKeV << " keV";
        if (nPriorCapped > 0) {
          configure().outStream << ", " << nPriorCapped
                                << " norm/shift parameter" << (nPriorCapped == 1 ? "" : "s")
                                << " limited to their prior width";
        }
        configure().outStream << "\n";
      }

      initial_positions = DrawInitialEnsemble(nwalkers, initialParams, paramSpread, generator, useRWA);
    }

    // Check the start before any step (see RunMCMCSampling in the header).
    auto report = [&](const std::string &text) {
      configure().outStream << text;
      configure().outStream.flush();
    };
    // The initial ensemble (evaluated first, so that the transformation
    // warnings of its walkers come before the report).
    int nFinite = 0;
    {
      std::vector<double> startLogP(nwalkers, -std::numeric_limits<double>::infinity());
#pragma omp parallel for num_threads(std::max(1, nthreads)) schedule(dynamic)
      for (int k = 0; k < nwalkers; k++) {
        const std::vector<double> &x = initial_positions[k];
        const double lp = CalculateLogPrior(x);
        if (std::isfinite(lp)) startLogP[k] = lp + (useRWA ? LogLikelihood(x) : LogLikelihoodPhysical(x));
      }
      for (int k = 0; k < nwalkers; k++)
        if (std::isfinite(startLogP[k])) nFinite++;
    }
    // Park: name the levels with J <= 0 at the starting point.
    if (configure().paramMask & Config::USE_PARK_FORMALISM) {
      bool filled = true;
      std::vector<std::string> badLevels = ParkNonPositiveLevels(initialParams, useRWA, filled);
      std::ostringstream msg;
      if (!badLevels.empty()) {
        msg << "WARNING: the MCMC starting point is outside Park's parameter space: "
            << badLevels.size() << " level" << (badLevels.size() == 1 ? " has" : "s have")
            << " J = 1 - sum gamma^2 dS/dE <= 0\n";
        for (size_t i = 0; i < badLevels.size(); ++i) msg << "    " << badLevels[i] << "\n";
        msg << "  Their widths exceed what the channel radii allow; the posterior is zero "
            << "there (J < 0 is rejected), so walkers near this point are rejected.\n"
            << "  Start from a fit under --use-park (its penalty keeps J > 0), or from a "
            << "Brune parameter file (converted on read), or reduce these widths.\n";
      } else if (!filled) {
        msg << "WARNING: the MCMC starting point could not be transformed to "
            << "Park's amplitudes; its posterior probability is zero.\n";
      }
      if (!msg.str().empty()) report(msg.str());
    }
    // Every walker of the initial ensemble outside the support: refuse.
    {
      if (nFinite == 0) {
        std::ostringstream msg;
        msg << "ERROR: none of the " << nwalkers << " walkers of the initial "
            << "ensemble has a finite posterior probability";
        if (configure().paramMask & Config::USE_PARK_FORMALISM)
          msg << " (under Park, a level with J < 0 is outside the prior)";
        msg << ".\n  The sampler cannot move from there: every walker would stay "
            << "at its start with acceptance 0.  MCMC not started.\n";
        report(msg.str());
        csvFile.close();
        if (!resuming) std::remove(configure().OpenPath(samplesFile).c_str());  // only its header was written
        return false;
      }
    }

    // This run's state, owned here; the callbacks below reach it and this sampler through their captures.
    Run run;
    run.sampleFile = &csvFile;
    run.nwalkers = nwalkers;
    run.currentStep = startStep;
    run.stepOffset = startStep;  // Continue the chain numbering when resuming
    run.useRWA = useRWA;
    stopRequested_ = false;  // a stop requested before this run does not stop it

    // Create numcmc sampler
    nu::Mcmc mcmc_sampler(nwalkers, ndim, initial_positions);
    // Deterministic given the starting step, so a run is reproducible while a
    // resumed run still draws a fresh stretch of the random stream.
    mcmc_sampler.set_seed(0x9E3779B97F4A7C15ull ^ (std::uint64_t)(startStep + 1));

    configure().outStream << "Starting MCMC sampling with enhanced progress tracking";
    if (nthreads > 1) {
      configure().outStream << " using " << nthreads << " parallel threads";
    }
    configure().outStream << "...\n";
    configure().outStream.flush();

    int remainingSteps = nsteps - startStep;

    configure().outStream << "Using " << (useRWA ? "reduced width amplitudes (RWA)" : "physical parameters") << " for MCMC fitting\n";

    const nu::Mcmc::LogProbability logProbability = [this, useRWA](std::vector<double> &x) {
      return LogProbability(x, useRWA);
    };
    const nu::Mcmc::StopCallback shouldStop = [this]() { return stopRequested_.load(); };
    const nu::Mcmc::IterationCallback afterStep = [this, &run](int step, int total) { AfterStep(run, step, total); };
    const nu::Mcmc::SampleCallback recordStep = [this, &run](int step, const std::vector<std::vector<double>> &positions,
                                                             const std::vector<double> &logps) {
      RecordStep(run, step, positions, logps);
    };
    int result;
    if (nthreads > 1) {
      result = mcmc_sampler.run_parallel_with_callback(logProbability, remainingSteps, nthreads, nullptr, shouldStop,
                                                       afterStep, recordStep);
    } else {
      result = mcmc_sampler.run_with_callback(logProbability, remainingSteps, nullptr, shouldStop, afterStep,
                                              recordStep);
    }

    // Save the ensemble so a later run resumes it rather than re-scattering.
    // Also written when the user stopped early: that is precisely the case where
    // continuing from the same positions matters.
    {
      std::vector<nu::Walker> finalWalkers = mcmc_sampler.get_walkers();
      std::vector<std::vector<double>> finalPositions;
      for (size_t k = 0; k < finalWalkers.size(); k++) {
        finalPositions.push_back(finalWalkers[k].getPos());
      }
      if (!save_walker_state(configure().OpenPath(walkerStateFile), finalPositions)) {
        configure().outStream << "Warning: could not write walker state to "
                              << walkerStateFile << "; a later resume will restart "
                              << "the walkers instead of continuing them.\n";
      }
    }

    if (mcmc_sampler.num_invalid_initial() > 0) {
      configure().outStream << "Warning: " << mcmc_sampler.num_invalid_initial()
                            << " of " << nwalkers << " walkers started at a position "
                            << "with zero posterior probability and can never move. "
                            << "Reduce the initial spread or check the priors.\n";
    }
    configure().outStream << "Acceptance fraction: "
                          << std::fixed << std::setprecision(3)
                          << mcmc_sampler.acceptance_fraction()
                          << std::defaultfloat << "\n";
    if (mcmc_sampler.acceptance_fraction() < 0.05) {
      configure().outStream << "  Warning: a very low acceptance fraction means the "
                            << "chain is barely moving; the samples are not yet a "
                            << "usable posterior.\n";
    }
    configure().outStream.flush();

    csvFile.close();

    if (result == 1) {
      configure().outStream << "MCMC sampling stopped early by user request.\n";
    } else if (result != 0) {
      configure().outStream << "MCMC sampling failed with error code: " << result << "\n";
      return false;
    }

    // Get samples from the sampler
    std::vector<std::vector<double>> newSamples = mcmc_sampler.get_chain();

    // Return all samples
    if (resuming) {
      std::vector<std::vector<double>> existingSamples;
      LoadExistingSamples(samplesFile, existingSamples);
      samples = existingSamples;
    } else {
      samples = newSamples;
    }

    configure().outStream << "MCMC sampling completed successfully. Generated "
                          << samples.size() << " total samples.\n";
    configure().outStream << "Samples saved to: " << samplesFile << "\n";
    return true;

  } catch (const std::exception &e) {
    configure().outStream << "Error during MCMC sampling: " << e.what() << "\n";
  } catch (...) {
    configure().outStream << "Unknown error during MCMC sampling\n";
  }
  return false;
}
