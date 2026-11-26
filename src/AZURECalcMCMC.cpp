// AZURECalcMCMC, the posterior of one walker: the log-likelihood of a parameter vector on a pooled level state and
// working copy (CalculateLogLikelihood, CalculateLogLikelihoodPhysical, LogLikelihoodOfFilled), the map between
// the sampler's varying parameters and the full vector in both bases (UpdateParameterVectors,
// ReconstructFullParameters*), and the log-probability the ensemble sampler asks for (LogProbability*).  One
// unit because this is all a proposal evaluates, from several threads at once and without a side effect.  The
// other jobs of the class are files of their own: AZURECalcMCMCPriors.cpp (the priors), AZURECalcMCMCRun.cpp
// (the run of the chain), AZURECalcMCMCChain.cpp (the chain's files and the step's reports),
// AZURECalcMCMCPool.cpp (the idle states of the caller's SamplerStates).
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

double AZURECalcMCMC::CalculateLogLikelihood(const vector_r &p) const {
  // Outside the support (the fit's limits of any parameter): zero probability, nothing evaluated.
  if (!InsideFitLimits(p)) return -std::numeric_limits<double>::infinity();
  const CompoundSystem *localCompound = compound();
  WalkerLevels *walker = nullptr;
  WorkingCopy *localData = NULL;
  WorkingData *working = nullptr;

  try {
    // An evaluation of its own (made when no idle one is left: one per walker evaluated at the same time).
    walker = GetPooledLevels();
    working = GetPooledData();
    localData = working->data.get();

    // Objects are already cloned, no need for deep copying

    // Fill Compound Nucleus From Parameters
    AZUREParams params;
    localCompound->FillCompoundFromParams(p, walker->state);
    localData->FillNormsFromParams(p);
    localData->FillEnergyShiftsFromParams(p, working->results, localCompound, &walker->state, &walker->inputs, &configure());
    if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure(), walker->state, CompoundSystem::ShiftOrders::kFirstDerivativesOnly);
  } catch (GSLException &e) {
    // Clean up and return bad likelihood for GSL errors
    if (walker) ReturnPooledLevels(walker);
    if (working) ReturnPooledData(working);
    return -std::numeric_limits<double>::infinity();
  } catch (...) {
    // Clean up and return bad likelihood for any other errors
    if (walker) ReturnPooledLevels(walker);
    if (working) ReturnPooledData(working);
    return -std::numeric_limits<double>::infinity();
  }

  return LogLikelihoodOfFilled(walker, working, p);
}

double AZURECalcMCMC::LogLikelihoodOfFilled(WalkerLevels *walker, WorkingData *working, const vector_r &p) const {
  const CompoundSystem *localCompound = compound();
  azure2::Results &results = working->results;
  ObjectiveContext context{localCompound, walker->state, walker->inputs, working->data.get(), results, p, configure()};
  // A support of the posterior (Park's J > 0, docs/upgrade/objective.md D2): -infinity outside.
  if (Objective::Engine().ExcludesForSampler(context)) {
    ReturnPooledLevels(walker);
    ReturnPooledData(working);
    return -std::numeric_limits<double>::infinity();
  }
  // The model value of every point, then the fit's objective as the sampler takes it (Objective::SamplerChi2,
  // objective.md section 6): every term at its value, the data with THM norms profiled as the fit does (D1), except
  // the priors this sampler normalises (CalculateLogPrior).  A term added to the objective arrives here unedited.
  // The points first in one pass over the compound energies, as a fit's evaluation (AZURECalc::PrecalculatePoints): on
  // this thread alone when the walkers are evaluated at once, every record it writes being this evaluation's own.
  std::vector<const Segment *> segments;
  for (int i = 1; i <= working->data->Data().NumSegments(); i++)
    if (const Segment *segment = working->data->Data().GetSegment(i)) segments.push_back(segment);
  const bool pre = working->data->CalculatePointsInParallel(segments, localCompound, walker->state, configure(), results) >= 0;
  azure2::EvaluateModel(*working->data, localCompound, walker->state, configure(), results,
                        pre ? azure2::LeafValues::Ready : azure2::LeafValues::Calculate);
  const double chiSquared = Objective::Engine().SamplerChi2(context);

  ReturnPooledLevels(walker);
  ReturnPooledData(working);
  return -0.5 * chiSquared;
}

double AZURECalcMCMC::CalculateLogLikelihoodPhysical(const vector_r &params_) const {
  const CompoundSystem *localCompound = compound();
  WalkerLevels *walker = nullptr;
  WorkingCopy *localData = NULL;
  WorkingData *working = nullptr;
  AZUREParams params;

  try {
    // An evaluation of its own (made when no idle one is left: one per walker evaluated at the same time).
    walker = GetPooledLevels();
    working = GetPooledData();
    localData = working->data.get();

    // Objects are already cloned, no need for deep copying

    localCompound->FillCompoundFromParamsPhysical(params_, walker->state, walker->inputs);
    localData->FillAbundancesFromVector(params_);
    bool isValid = localCompound->TransformIn(configure(), walker->state, walker->inputs);

    if (!isValid) {
      ReturnPooledLevels(walker);
      ReturnPooledData(working);
      return -std::numeric_limits<double>::infinity();
    }
  } catch (GSLException &e) {
    // Clean up and return bad likelihood for GSL errors
    if (walker) ReturnPooledLevels(walker);
    if (working) ReturnPooledData(working);
    return -std::numeric_limits<double>::infinity();
  } catch (...) {
    // Clean up and return bad likelihood for any other errors
    if (walker) ReturnPooledLevels(walker);
    if (working) ReturnPooledData(working);
    return -std::numeric_limits<double>::infinity();
  }

  // The fit's vector at this point (amplitudes, the data parameters from the tail of the observed vector), so that the
  // support and a term reading the vector read the fit's: the limits act on the converted amplitudes and energies.
  const vector_r fitVector = FitVectorOfObserved(*walker, localData, params_);
  if (!InsideFitLimits(fitVector)) {
    ReturnPooledLevels(walker);
    ReturnPooledData(working);
    return -std::numeric_limits<double>::infinity();
  }
  localCompound->FillMnParams(params, walker->state, walker->inputs, &configure());
  localData->FillMnParams(params);
  localCompound->FillCompoundFromParams(params.GetMinuitParams().Params(), walker->state);
  // The tail of the physical vector (normalizations, energy shifts, screening)
  // is the same in both bases.  The norms were never copied out of it, so the
  // pooled WorkingCopy kept its starting norms: the sampler moved the norm
  // parameters under their priors while the likelihood never saw them (norm
  // posterior = prior, R-matrix posterior at frozen norms; 2026-10-26).
  localData->FillNormsFromParams(params_);
  // After the compound: the sub-point grids are laid out from its fit values.
  localData->FillEnergyShiftsFromParams(params_, working->results, localCompound, &walker->state, &walker->inputs, &configure());

  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure(), walker->state, CompoundSystem::ShiftOrders::kFirstDerivativesOnly);
  return LogLikelihoodOfFilled(walker, working, fitVector);
}

void AZURECalcMCMC::UpdateParameterVectors(const vector_r &physicalParams) {
  if (parametersInitialized_) return;

  all_indexes.clear();
  all_physical_.clear();
  all_rwa_.clear();
  rwa_.clear();
  physical_.clear();
  fixed_.clear();
  paramNames_.clear();
  paramKinds_.clear();

  AZUREParams params;
  compound()->FillMnParams(params, levels_.state, levels_.inputs, &configure());
  data()->FillMnParams(params);
  if (configure().paramMask & Config::USE_PREVIOUS_PARAMETERS) {
    params.ReadUserParameters(configure());
    params.ReconcileBasis(compound(), configure());
    params.ApplyFileWidthMask(compound(), levels_.state, configure());
  }
  // The support: the fit's limits, as AZUREMain sets them before a fit (the kinds' own and <parameterSettings>).
  support_ = PosteriorSupportOfFit(configure(), compound(), levels_.state, levels_.inputs, data(), params);

  compound()->FillCompoundFromParams(params.GetMinuitParams().Params(), levels_.state);
  compound()->CalcShiftFunctions(configure(), levels_.state);
  const ObservedLevelTable observed = compound()->TransformOut(configure(), levels_.state, levels_.inputs);

  for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
    all_rwa_.push_back(params.GetMinuitParams().Parameter(i).Value());
    fixed_.push_back(params.GetMinuitParams().Parameter(i).IsFixed());
    paramNames_.push_back(params.GetMinuitParams().GetName(i));
    paramKinds_.push_back(params.Kind(i));
    // If not fixed, add to rwa
    const ROOT::Minuit2::MinuitParameter &par = params.GetMinuitParams().Parameter(i);
    if (!par.IsFixed()) rwa_.push_back(par.Value());
  }

  all_physical_ = compound()->GetTransformParams(levels_.state, observed);
  data()->Data().Abundances().AppendTo(all_rwa_, all_physical_);
  // Add missing parameters
  for (size_t i = all_physical_.size(); i < all_rwa_.size(); i++) {
    all_physical_.push_back(all_rwa_[i]);
  }

  // Save the not fixed physical parameters
  int k = 0;
  for (size_t i = 0; i < fixed_.size(); ++i) {
    if (!fixed_[i]) {
      physical_.push_back(all_physical_[i]);
      k++;
    }
  }

  parametersInitialized_ = true;
}

vector_r AZURECalcMCMC::ReconstructFullParametersPhysical(const std::vector<double> &varyingParams) const {
  if (!parametersInitialized_) {
    return vector_r();  // Return empty if not initialized
  }

  // Same logic as LogLikelihoodPhysical
  int k = 0;
  vector_r fullParams = all_physical_;  // Start with all parameters
  for (int i = 0; i < all_physical_.size(); ++i) {
    if (!fixed_[i]) {
      if (k < varyingParams.size()) {
        fullParams[i] = varyingParams[k];  // Only substitute non-fixed parameters
        ++k;
      }
    }
  }

  return fullParams;
}

vector_r AZURECalcMCMC::ReconstructFullParameters(const std::vector<double> &varyingParams) const {
  if (!parametersInitialized_) {
    return vector_r();  // Return empty if not initialized
  }

  // Same logic but for RWA parameters
  int k = 0;
  vector_r fullParams = all_rwa_;  // Start with all RWA parameters
  for (int i = 0; i < all_rwa_.size(); ++i) {
    if (!fixed_[i]) {
      if (k < varyingParams.size()) {
        fullParams[i] = varyingParams[k];  // Only substitute non-fixed parameters
        ++k;
      }
    }
  }

  return fullParams;
}

double AZURECalcMCMC::LogLikelihoodPhysical(const std::vector<double> &physicalParams) const {
  int k = 0;
  vector_r params_ = all_physical_;  // Start with all parameters (like AZUREAPI::UpdateSegments)
  for (int i = 0; i < all_physical_.size(); ++i) {
    if (!fixed_[i]) {
      if (k < physicalParams.size()) {
        params_[i] = physicalParams[k];  // Only substitute non-fixed parameters
        ++k;
      }
    }
  }

  // Calculate log-likelihood from data only (no priors)
  return CalculateLogLikelihoodPhysical(params_);
}

double AZURECalcMCMC::LogLikelihood(const std::vector<double> &rwaParams) const {
  int k = 0;
  vector_r params_ = all_rwa_;  // Start with all RWA parameters
  for (int i = 0; i < all_rwa_.size(); ++i) {
    if (!fixed_[i]) {
      if (k < rwaParams.size()) {
        params_[i] = rwaParams[k];  // Only substitute non-fixed parameters
        ++k;
      }
    }
  }

  // Calculate log-likelihood from data only (no priors) - use standard RWA calculation
  return CalculateLogLikelihood(params_);
}

double AZURECalcMCMC::LogProbabilityPhysical(const std::vector<double> &physicalParams) const {
  // Calculate log-likelihood from data
  double logLikelihood = LogLikelihoodPhysical(physicalParams);

  // Add log-prior contribution for Bayesian inference
  double logPrior = CalculateLogPrior(physicalParams);

  return logLikelihood + logPrior;
}

/* The log-probability must stay free of side effects.
 *
 * It is called for every *proposal*, several times per walker per step and
 * from several threads at once.  Recording the chain here -- as an earlier
 * version did -- writes rejected proposals into the samples file with the same
 * weight as accepted ones and omits the repeats that a rejection contributes,
 * so the file is a log of evaluations rather than a Markov chain and posterior
 * summaries taken from it are wrong.  The chain is written once per step by
 * RecordStep() instead.
 */
double AZURECalcMCMC::LogProbability(const std::vector<double> &params, bool useRWA) const {
  // Priors are analytic and cheap; evaluating them first lets a point that
  // the prior already excludes skip the chi-squared entirely.
  const double logPrior = CalculateLogPrior(params);
  if (!std::isfinite(logPrior)) return -std::numeric_limits<double>::infinity();
  const double logLikelihood = useRWA ? LogLikelihood(params) : LogLikelihoodPhysical(params);
  // The posterior, not the likelihood: returning the likelihood alone made
  // the sampler ignore every prior the user set.
  return logLikelihood + logPrior;
}
