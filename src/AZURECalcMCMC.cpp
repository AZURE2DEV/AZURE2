#include "AZURECalcMCMC.h"
#include "ParameterLabel.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "GSLException.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Global function pointer to access the MCMC calculator instance
static const AZURECalcMCMC* g_mcmc_calc = nullptr;
static std::ofstream* g_sample_file = nullptr;
static int g_sample_count = 0;
static int g_nwalkers = 0;
static int g_current_step = 0;
static volatile int g_samples_written_to_file = 0; // Atomic counter of samples actually written to file
static volatile bool g_stop_requested = false;
// Step index this run started from, so a resumed run keeps numbering the chain
// where the previous one stopped.
static int g_step_offset = 0;
// Best-walker values of the most recent completed step, for progress callbacks
static double g_last_likelihood = 0.0;
static double g_last_prior = 0.0;
static double g_last_logprob = 0.0;
static std::vector<double> g_last_params;
static std::vector<std::string> g_param_names;
#ifdef _OPENMP
#include <omp.h>
#endif

double AZURECalcMCMC::CalculateLogLikelihood(const vector_r& p) const {

  CNuc* localCompound = NULL;
  EData* localData = NULL;

  try {
    // Initialize pools on first use
    if (!pools_initialized_) {
      InitializePools();
    }
    
    // Get objects from pool
    localCompound = GetPooledCNuc();
    localData = GetPooledEData();
    
    // Objects are already cloned, no need for deep copying

    //Fill Compound Nucleus From Parameters
    AZUREParams params;
    localCompound->FillCompoundFromParams(p);
    localData->FillNormsFromParams(p);
    localData->FillEnergyShiftsFromParams(p,localData,localCompound,&configure());
    if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());
    
    // Sub-segments are now integrated into ESegment, no separate initialization needed
  } catch (GSLException& e) {
    // Clean up and return bad likelihood for GSL errors
    if(localCompound) ReturnPooledCNuc(localCompound);
    if(localData) ReturnPooledEData(localData);
    return -std::numeric_limits<double>::infinity();
  } catch (...) {
    // Clean up and return bad likelihood for any other errors
    if(localCompound) ReturnPooledCNuc(localCompound);
    if(localData) ReturnPooledEData(localData);
    return -std::numeric_limits<double>::infinity();
  }

  bool isFit=true;
  
  // Process segments with components - use new integrated calculation method
  double chiSquared=0.0;
  for(int i = 1; i <= localData->NumSegments(); i++) {
    ESegment* segment = localData->GetSegment(i);
    if(segment) {
      // Recalculate points using the new combined calculation method
      for(int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        double theoreticalValue = segment->CalculateTheoreticalCrossSection(pointIdx, localCompound, configure(), localData);
        EPoint* point = segment->GetPoint(pointIdx + 1);
        if(point) {
          point->SetFitCrossSection(theoreticalValue);
        }
      }
      
      // Recalculate chi-squared for this segment with components
      double segmentChiSquared = 0.0;
      for(int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        EPoint* point = segment->GetPoint(pointIdx + 1);
        if(point) {
          double residual = point->GetFitCrossSection() - point->GetCMCrossSection() * segment->GetNorm();
          double error = point->GetCMCrossSectionError() * segment->GetNorm();
          if(error != 0.0) {
            segmentChiSquared += (residual * residual) / (error * error);
          }
        }
      }

      segment->SetSegmentChiSquared(segmentChiSquared);
      chiSquared += segmentChiSquared;
    }
  }

  // Return objects to pool instead of deleting
  ReturnPooledCNuc(localCompound);
  ReturnPooledEData(localData);

  // Convert chi-squared to log-likelihood: ln(L) = -0.5 * chi^2
  return -0.5 * chiSquared;
  
}

double AZURECalcMCMC::CalculateLogLikelihoodPhysical(const vector_r& params_) const {

  CNuc * localCompound = NULL;
  EData *localData = NULL;
  AZUREParams params;
  
  try {
    // Initialize pools on first use
    if (!pools_initialized_) {
      InitializePools();
    }
    
    // Get objects from pool
    localCompound = GetPooledCNuc();
    localData = GetPooledEData();
    
    // Objects are already cloned, no need for deep copying

    localCompound->FillCompoundFromParamsPhysical(params_);
    bool isValid = localCompound->TransformIn( configure( ) );
    
    if(!isValid) {
      ReturnPooledCNuc(localCompound);
      ReturnPooledEData(localData);
      return -std::numeric_limits<double>::infinity();
    }
  } catch (GSLException& e) {
    // Clean up and return bad likelihood for GSL errors
    if(localCompound) ReturnPooledCNuc(localCompound);
    if(localData) ReturnPooledEData(localData);
    return -std::numeric_limits<double>::infinity();
  } catch (...) {
    // Clean up and return bad likelihood for any other errors
    if(localCompound) ReturnPooledCNuc(localCompound);
    if(localData) ReturnPooledEData(localData);
    return -std::numeric_limits<double>::infinity();
  }

  localCompound->FillMnParams(params.GetMinuitParams(), &configure());
  localData->FillMnParams(params.GetMinuitParams());
  localData->FillEnergyShiftsFromParams(params_,localData,localCompound,&configure());
  localCompound->FillCompoundFromParams(params.GetMinuitParams( ).Params( ));

  //Fill Compound Nucleus From Minuit Parameters
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());
  
  // Sub-segments are now integrated into ESegment, no separate initialization needed
  
  bool isFit = true; // For MCMC, we are always fitting
  
  // Process segments with components - use new integrated calculation method
  double chiSquared=0.0;
  for(int i = 1; i <= localData->NumSegments(); i++) {
    ESegment* segment = localData->GetSegment(i);
    if(segment) {
      // Recalculate points using the new combined calculation method
      for(int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        double theoreticalValue = segment->CalculateTheoreticalCrossSection(pointIdx, localCompound, configure(), localData);
        EPoint* point = segment->GetPoint(pointIdx + 1);
        if(point) {
          point->SetFitCrossSection(theoreticalValue);
        }
      }
      
      // Recalculate chi-squared for this segment with components
      double segmentChiSquared = 0.0;
      for(int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        EPoint* point = segment->GetPoint(pointIdx + 1);
        if(point) {
          double residual = point->GetFitCrossSection() - point->GetCMCrossSection() * segment->GetNorm();
          double error = point->GetCMCrossSectionError() * segment->GetNorm();
          if(error != 0.0) {
            segmentChiSquared += (residual * residual) / (error * error);
          }
        }
      }

      segment->SetSegmentChiSquared(segmentChiSquared);
      chiSquared += segmentChiSquared;
    }
  }

  // Return objects to pool instead of deleting
  ReturnPooledCNuc(localCompound);
  ReturnPooledEData(localData);

  // Convert chi-squared to log-likelihood: ln(L) = -0.5 * chi^2
  return -0.5 * chiSquared;
}

void AZURECalcMCMC::UpdateParameterVectors(const vector_r& physicalParams) const {

  if (parametersInitialized_) return;

  all_indexes.clear();
  all_physical_.clear();
  all_rwa_.clear();
  rwa_.clear();
  physical_.clear();
  fixed_.clear();
  g_param_names.clear();
  
  AZUREParams params;
  compound()->FillMnParams(params.GetMinuitParams(), &configure());
  data()->FillMnParams(params.GetMinuitParams());
  if(configure().paramMask & Config::USE_PREVIOUS_PARAMETERS) {
    params.ReadUserParameters(configure());
  }
  
  compound()->FillCompoundFromParams(params.GetMinuitParams().Params());
  compound()->CalcShiftFunctions(configure());
  compound()->TransformOut(configure());
  
  for(int i = 0; i < params.GetMinuitParams().Params().size(); i++){
    all_rwa_.push_back(params.GetMinuitParams().Parameter(i).Value());
    fixed_.push_back(params.GetMinuitParams().Parameter(i).IsFixed());
    g_param_names.push_back(params.GetMinuitParams().GetName(i));
    // If not fixed, add to rwa
    if (!params.GetMinuitParams().Parameter(i).IsFixed()) {
      rwa_.push_back(params.GetMinuitParams().Parameter(i).Value());
    }
  }
  
  all_physical_ = compound()->GetTransformParams(configure());
  // Add missing parameters
  for (size_t i = all_physical_.size(); i < all_rwa_.size(); i++) {
    all_physical_.push_back(all_rwa_[i]);
  }

  // Save the not fixed physical parameters
  int k =0;
  for (size_t i = 0; i < fixed_.size(); ++i) {
    if (!fixed_[i]) {
      physical_.push_back(all_physical_[i]);
      k++;
    }
  }

  const_cast<AZURECalcMCMC*>(this)->parametersInitialized_ = true;
}

vector_r AZURECalcMCMC::ReconstructFullParametersPhysical(const std::vector<double>& varyingParams) const {
  if (!parametersInitialized_) {
    return vector_r(); // Return empty if not initialized
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

vector_r AZURECalcMCMC::ReconstructFullParameters(const std::vector<double>& varyingParams) const {
  if (!parametersInitialized_) {
    return vector_r(); // Return empty if not initialized
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

double AZURECalcMCMC::LogLikelihoodPhysical(const std::vector<double>& physicalParams) const {
  
  int k = 0;
  vector_r params_ = all_physical_;  // Start with all parameters (like AZUREAPI::UpdateSegments)
  for( int i = 0; i < all_physical_.size( ); ++i ){
    if( !fixed_[i] ){
      if(k < physicalParams.size()) {
        params_[i] = physicalParams[k];  // Only substitute non-fixed parameters
        ++k;
      }
    }
  }

  // Calculate log-likelihood from data only (no priors)
  return CalculateLogLikelihoodPhysical(params_);
}

double AZURECalcMCMC::LogLikelihood(const std::vector<double>& rwaParams) const {
  
  int k = 0;
  vector_r params_ = all_rwa_;  // Start with all RWA parameters
  for( int i = 0; i < all_rwa_.size( ); ++i ){
    if( !fixed_[i] ){
      if(k < rwaParams.size()) {
        params_[i] = rwaParams[k];  // Only substitute non-fixed parameters
        ++k;
      }
    }
  }

  // Calculate log-likelihood from data only (no priors) - use standard RWA calculation
  return CalculateLogLikelihood(params_);
}

double AZURECalcMCMC::LogProbabilityPhysical(const std::vector<double>& physicalParams) const {
  
  // Calculate log-likelihood from data
  double logLikelihood = LogLikelihoodPhysical(physicalParams);
  
  // Add log-prior contribution for Bayesian inference
  double logPrior = CalculateLogPrior(physicalParams);
  
  return logLikelihood + logPrior;
}

void AZURECalcMCMC::SetPriors(const std::vector<double>& priorMeans,
                              const std::vector<double>& priorStds,
                              const std::vector<bool>& usePriors) {
  priorMeans_ = priorMeans;
  priorStds_ = priorStds;
  usePriors_ = usePriors;
}

void AZURECalcMCMC::BuildAutoPriors() const {

  if(!parametersInitialized_) {
    configure().outStream << "Warning: BuildAutoPriors() called before parameters were "
                          << "initialized; automatic priors skipped.\n";
    return;
  }

  // Walk the parameters in exactly the order CNuc::FillMnParams and then
  // EData::FillMnParams add them, which is the order fixed_ is indexed in:
  //   per J-group, per level: one energy, then one width per channel
  //   one norm per segment with IsVaryNorm()
  //   one energy shift per segment (all segments)
  // AZUREAPI::GetParameterInfo() walks the same sequence.
  std::vector<int> allKinds;
  std::vector<double> allAutoMean;
  std::vector<double> allAutoStd;   // <= 0 means "no automatic prior available"

  CNuc* nuc = compound();
  for(int j = 1; j <= nuc->NumJGroups(); ++j) {
    JGroup* jgroup = nuc->GetJGroup(j);
    for(int la = 1; la <= jgroup->NumLevels(); ++la) {
      allKinds.push_back(PARAM_ENERGY);
      allAutoMean.push_back(0.0);
      allAutoStd.push_back(-1.0);
      for(int ch = 1; ch <= jgroup->NumChannels(); ++ch) {
        allKinds.push_back(PARAM_WIDTH);
        allAutoMean.push_back(0.0);
        allAutoStd.push_back(-1.0);
      }
    }
  }

  std::vector<ESegment>& segments = data()->GetSegments();

  // Normalizations. GetNormError() is a percentage of the nominal norm; this
  // reproduces the penalty AZURECalc::operator() adds during a Minuit fit.
  for(size_t s = 0; s < segments.size(); ++s) {
    if(segments[s].IsVaryNorm()) {
      double nominal = segments[s].GetNominalNorm();
      allKinds.push_back(PARAM_NORM);
      allAutoMean.push_back(nominal);
      allAutoStd.push_back(nominal / 100.0 * segments[s].GetNormError());
    }
  }

  // Energy shifts, one per segment whether it varies or not.
  for(size_t s = 0; s < segments.size(); ++s) {
    allKinds.push_back(PARAM_SHIFT);
    allAutoMean.push_back(segments[s].GetNominalEnergyShift());
    allAutoStd.push_back(segments[s].GetEnergyShiftError());
  }

  if(allKinds.size() != fixed_.size()) {
    configure().outStream << "Warning: parameter classification produced "
                          << allKinds.size() << " entries but there are "
                          << fixed_.size() << " parameters; automatic priors "
                          << "skipped to avoid mis-assigning them.\n";
    freeKinds_.clear();
    return;
  }

  // Compress to the varying parameters, which is how the priors and the
  // sampled vector are indexed.
  freeKinds_.clear();
  std::vector<double> freeAutoMean;
  std::vector<double> freeAutoStd;
  std::vector<int> freeAllIndex;   // varying index -> index among all parameters
  for(size_t i = 0; i < fixed_.size(); ++i) {
    if(!fixed_[i]) {
      freeKinds_.push_back(allKinds[i]);
      freeAutoMean.push_back(allAutoMean[i]);
      freeAutoStd.push_back(allAutoStd[i]);
      freeAllIndex.push_back((int)i);
    }
  }

  const size_t nfree = freeKinds_.size();

  // Any priors the caller set cover the varying parameters; grow them to the
  // full length so the loop below can index them unconditionally.
  priorMeans_.resize(nfree, 0.0);
  priorStds_.resize(nfree, 0.0);
  usePriors_.resize(nfree, false);

  int nAutoNorm = 0, nAutoShift = 0, nNoError = 0, nUserNeeded = 0, nUserSet = 0;
  std::vector<std::string> noErrorLabels;
  std::vector<std::string> noPriorLabels;

  for(size_t i = 0; i < nfree; ++i) {
    const int kind = freeKinds_[i];

    if(kind == PARAM_NORM || kind == PARAM_SHIFT) {
      // Derived from the data, so it always wins over whatever the GUI holds.
      if(freeAutoStd[i] > 0.0) {
        priorMeans_[i] = freeAutoMean[i];
        priorStds_[i] = freeAutoStd[i];
        usePriors_[i] = true;
        if(kind == PARAM_NORM) nAutoNorm++; else nAutoShift++;
      } else {
        // No experimental error quoted: nothing to build a prior from, so the
        // parameter stays uniform rather than getting an invented width.
        usePriors_[i] = false;
        nNoError++;
        noErrorLabels.push_back(AZURELabel::Parameter(compound(), data(), freeAllIndex[i]));
      }
    } else {
      // Level energies and widths: the user's business.
      nUserNeeded++;
      if(usePriors_[i] && priorStds_[i] > 0.0) nUserSet++;
      else noPriorLabels.push_back(AZURELabel::Parameter(compound(), data(), freeAllIndex[i]));
    }
  }

  configure().outStream << "Automatic priors: " << nAutoNorm << " normalization"
                        << (nAutoNorm == 1 ? "" : "s") << ", " << nAutoShift
                        << " energy shift" << (nAutoShift == 1 ? "" : "s")
                        << " taken from the quoted experimental errors\n";
  if(nNoError > 0) {
    configure().outStream << "  " << nNoError << " normalization/energy-shift "
                          << "parameter" << (nNoError == 1 ? " has" : "s have")
                          << " no quoted error and " << (nNoError == 1 ? "is" : "are")
                          << " sampled with a uniform prior:\n";
    for(size_t i = 0; i < noErrorLabels.size(); ++i)
      configure().outStream << "      " << noErrorLabels[i] << "\n";
  }
  configure().outStream << "  " << nUserSet << " of " << nUserNeeded
                        << " level-energy/width parameters have a user-defined prior\n";
  if(nUserSet < nUserNeeded) {
    configure().outStream << "  Warning: " << (nUserNeeded - nUserSet)
                          << " level-energy/width parameter"
                          << ((nUserNeeded - nUserSet) == 1 ? " is" : "s are")
                          << " sampled with an unbounded uniform prior. Define "
                          << "priors for them in the MCMC tab if the chain wanders:\n";
    // Cap the list so a large model cannot bury the rest of the log.
    const size_t maxListed = 10;
    for(size_t i = 0; i < noPriorLabels.size() && i < maxListed; ++i)
      configure().outStream << "      " << noPriorLabels[i] << "\n";
    if(noPriorLabels.size() > maxListed)
      configure().outStream << "      ...and " << (noPriorLabels.size() - maxListed)
                            << " more\n";
  }
  configure().outStream.flush();
}

double AZURECalcMCMC::CalculateLogPrior(const std::vector<double>& params) const {
  if(priorMeans_.empty() || priorStds_.empty() || usePriors_.empty()) {
    return 0.0; // No priors specified - uniform priors
  }
  
  double logPrior = 0.0;
  
  for(size_t i = 0; i < params.size() && i < usePriors_.size(); i++) {
    if(usePriors_[i] && i < priorMeans_.size() && i < priorStds_.size()) {
      double mean = priorMeans_[i];
      double std = priorStds_[i];
      
      if(std > 0.0) {
        // Gaussian prior: log(prior) = -0.5 * ((x - mean) / std)^2 - log(std * sqrt(2*pi))
        double deviation = (params[i] - mean) / std;
        logPrior += -0.5 * deviation * deviation - log(std * sqrt(2.0 * M_PI));
      }
    }
  }
  
  return logPrior;
}

/* The log-probability wrappers below must stay free of side effects.
 *
 * They are called for every *proposal*, several times per walker per step and
 * from several threads at once.  Recording the chain here -- as an earlier
 * version did -- writes rejected proposals into the samples file with the same
 * weight as accepted ones and omits the repeats that a rejection contributes,
 * so the file is a log of evaluations rather than a Markov chain and posterior
 * summaries taken from it are wrong.  The chain is written once per step by
 * mcmc_sample_callback() instead.
 */

// C-style function wrapper for numcmc library - uses RWA parameters
double mcmc_log_probability_wrapper_rwa(std::vector<double>& params) {
  if(g_mcmc_calc) {
    // Priors are analytic and cheap; evaluating them first lets a point that
    // the prior already excludes skip the chi-squared entirely.
    double logPrior = g_mcmc_calc->CalculateLogPrior(params);
    if(!std::isfinite(logPrior)) return -std::numeric_limits<double>::infinity();

    double logLikelihood = g_mcmc_calc->LogLikelihood(params);

    // The posterior, not the likelihood: returning the likelihood alone made
    // the sampler ignore every prior the user set.
    return logLikelihood + logPrior;
  }
  return -std::numeric_limits<double>::infinity();
}

// C-style function wrapper for numcmc library - uses physical parameters
double mcmc_log_probability_wrapper(std::vector<double>& params) {
  if(g_mcmc_calc) {
    double logPrior = g_mcmc_calc->CalculateLogPrior(params);
    if(!std::isfinite(logPrior)) return -std::numeric_limits<double>::infinity();

    double logLikelihood = g_mcmc_calc->LogLikelihoodPhysical(params);

    return logLikelihood + logPrior;
  }
  return -std::numeric_limits<double>::infinity();
}

double robust_stod(const std::string& str);

/* Walker-position persistence.
 *
 * The chain file records where the walkers have been; these two functions
 * record where they *are*, so that resuming continues the same ensemble. Without
 * them a resumed run re-scatters the walkers around the starting parameters and
 * splices a fresh burn-in into the middle of the chain file, which quietly
 * corrupts any posterior summary taken over the whole file.
 */
static bool load_walker_state(const std::string& filename, int nwalkers, int ndim,
                              std::vector<std::vector<double>>& positions) {
  std::ifstream file(filename);
  if(!file.is_open()) return false;

  std::vector<std::vector<double>> loaded;
  std::string line;
  while(std::getline(file, line)) {
    if(line.empty()) continue;
    std::stringstream ss(line);
    std::string token;
    std::vector<double> row;
    while(std::getline(ss, token, ',')) {
      try {
        row.push_back(robust_stod(token));
      } catch(...) {
        return false;
      }
    }
    if((int)row.size() != ndim) return false;
    loaded.push_back(row);
  }

  // A state saved for a different walker count or model cannot be reused.
  if((int)loaded.size() != nwalkers) return false;

  positions = loaded;
  return true;
}

static bool save_walker_state(const std::string& filename,
                              const std::vector<std::vector<double>>& positions) {
  std::ofstream file(filename);
  if(!file.is_open()) return false;

  for(size_t k = 0; k < positions.size(); k++) {
    for(size_t j = 0; j < positions[k].size(); j++) {
      // Full precision: a rounded restart position is a different model.
      file << std::scientific << std::setprecision(17) << positions[k][j];
      if(j + 1 != positions[k].size()) file << ",";
    }
    file << "\n";
  }
  return true;
}

// Records one completed MCMC step: every walker's accepted position, exactly
// once. Called from the sampler's serial section, so it needs no locking.
void mcmc_sample_callback(int current_step,
                          const std::vector<std::vector<double>>& positions,
                          const std::vector<double>& logps) {
  if(!g_mcmc_calc) return;

  // Absolute, 0-based step index, continuing an earlier run when resuming.
  const int absStep = g_step_offset + current_step - 1;

  int bestWalker = -1;
  double bestLogProb = -std::numeric_limits<double>::infinity();

  for(size_t k = 0; k < positions.size(); k++) {
    // Splitting the sampled log-posterior back into likelihood and prior is
    // exact: the prior is analytic, so this reports the very values the
    // sampler used rather than a re-evaluation.
    const double logProbability = logps[k];
    const double logPrior = g_mcmc_calc->CalculateLogPrior(positions[k]);
    const double logLikelihood = logProbability - logPrior;

    if(logProbability > bestLogProb) {
      bestLogProb = logProbability;
      bestWalker = (int)k;
      g_last_likelihood = logLikelihood;
      g_last_prior = logPrior;
    }

    if(g_sample_file && g_sample_file->is_open()) {
      // Write with full double precision (17 digits) to preserve accuracy
      *g_sample_file << absStep << "," << k << ","
                     << std::scientific << std::setprecision(17)
                     << logProbability << "," << logLikelihood << "," << logPrior;
      for(const double& param : positions[k]) {
        *g_sample_file << "," << param;
      }
      *g_sample_file << std::defaultfloat << "\n";
    }
  }

  if(g_sample_file && g_sample_file->is_open()) {
    g_sample_file->flush(); // Ensure data is written immediately
    g_samples_written_to_file += (int)positions.size();
    g_sample_count += (int)positions.size();
  }

  // Best walker of the step drives the GUI readouts and the periodic output
  // files, which is more informative than whichever proposal happened to be
  // evaluated last.
  if(bestWalker >= 0) {
    g_last_params = positions[bestWalker];
    g_last_logprob = bestLogProb;
  }
  g_current_step = absStep + 1;
}

// Global pointer to the GUI progress handler
static void (*g_gui_progress_callback)(int, int, double, double, double) = nullptr;
static void (*g_gui_iteration_callback)(int, int) = nullptr;
static void (*g_gui_results_callback)(int, int, const std::vector<std::vector<double>>&) = nullptr;
static const AZURECalcMCMC* g_mcmc_for_callbacks = nullptr;
static nu::Mcmc* g_mcmc_sampler = nullptr;
static bool g_using_rwa_parameters = false;

// Custom robust string-to-double parser for scientific notation
double robust_stod(const std::string& str) {
    // Temporarily set numeric locale to "C"
    std::string old_locale = std::setlocale(LC_NUMERIC, nullptr);
    std::setlocale(LC_NUMERIC, "C");

    const char* cstr = str.c_str();
    char* end = nullptr;

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

// Iteration callback for GUI updates (called every iteration)
void mcmc_iteration_callback(int current_step, int total_steps) {
  if (g_gui_iteration_callback) {
    // The values come from mcmc_sample_callback(), which has just recorded this
    // step. An earlier version re-opened samples.mcmc and parsed its last line
    // on every iteration, which cost a file read per step and reported whichever
    // proposal was written last rather than an accepted state.
    int total_expected_samples = total_steps * g_nwalkers;

    g_gui_iteration_callback(g_current_step, total_expected_samples);

    // For progress callback (time estimation), use the original numcmc step values
    // This keeps the time estimation logic working correctly
    if (g_gui_progress_callback) {
      g_gui_progress_callback(current_step, total_steps,
                              g_last_logprob, g_last_likelihood, g_last_prior);
    }
  }
  
  // Write output files every 10 iterations (moved from progress callback)
  if (g_mcmc_for_callbacks && (current_step % 10 == 0)) {
    try {
      // Get the most recent parameter set from the wrapper and reconstruct full parameter array
      if (!g_last_params.empty()) {
        
        // Reconstruct full parameter array using the appropriate method based on parameter type
        vector_r fullParams;
        if (g_using_rwa_parameters) {
          fullParams = g_mcmc_for_callbacks->ReconstructFullParameters(g_last_params);
        } else {
          fullParams = g_mcmc_for_callbacks->ReconstructFullParametersPhysical(g_last_params);
        }
        
        if (!fullParams.empty()) {
          CNuc* localCompound = g_mcmc_for_callbacks->compound()->Clone();
          EData* localData = g_mcmc_for_callbacks->data()->Clone();

          bool isValid = true;
          vector_r rwaFullParams; // full RWA parameter set for param.mcmc

          if (g_using_rwa_parameters) {
            // For RWA parameters, use the same approach as CalculateLogLikelihood
            localCompound->FillCompoundFromParams(fullParams);
            localData->FillNormsFromParams(fullParams);
            localData->FillEnergyShiftsFromParams(fullParams, localData, localCompound, &g_mcmc_for_callbacks->configure());
            if (g_mcmc_for_callbacks->configure().paramMask & Config::USE_BRUNE_FORMALISM) {
              localCompound->CalcShiftFunctions(g_mcmc_for_callbacks->configure());
            }
            rwaFullParams = fullParams; // already RWA
          } else {
            // For physical parameters, use the same approach as CalculateLogLikelihoodPhysical
            localCompound->FillCompoundFromParamsPhysical(fullParams);
            isValid = localCompound->TransformIn(g_mcmc_for_callbacks->configure());

            if (isValid) {
              AZUREParams params;
              localCompound->FillMnParams(params.GetMinuitParams(), &g_mcmc_for_callbacks->configure());
              localData->FillMnParams(params.GetMinuitParams());
              localData->FillEnergyShiftsFromParams(fullParams, localData, localCompound, &g_mcmc_for_callbacks->configure());
              localCompound->FillCompoundFromParams(params.GetMinuitParams().Params());
              rwaFullParams = params.GetMinuitParams().Params(); // RWA params after TransformIn

              if (g_mcmc_for_callbacks->configure().paramMask & Config::USE_BRUNE_FORMALISM) {
                localCompound->CalcShiftFunctions(g_mcmc_for_callbacks->configure());
              }
            }
          }

          if (isValid) {
            // Calculate cross sections for all data points
            for(EDataIterator data=localData->begin(); data!=localData->end(); data++) {
              if(!data.point()->IsMapped()) {
                data.point()->Calculate(localCompound, g_mcmc_for_callbacks->configure());
              }
            }

            localData->WriteOutputFiles(g_mcmc_for_callbacks->configure(), true);

            if (!g_using_rwa_parameters) {
              localCompound->TransformOut(g_mcmc_for_callbacks->configure());
            }
            localCompound->PrintTransformParams(g_mcmc_for_callbacks->configure());

            // Write param.mcmc snapshot: all RWA parameters (including fixed), analogous to param.sav
            if (!rwaFullParams.empty() && rwaFullParams.size() == g_param_names.size()) {
              std::string paramMCMCFile = g_mcmc_for_callbacks->configure().outputdir + "param.mcmc";
              std::ofstream paramOut(paramMCMCFile);
              if (paramOut) {
                paramOut.precision(7);
                for (size_t i = 0; i < rwaFullParams.size(); i++) {
                  paramOut << std::setw(20) << g_param_names[i]
                           << std::scientific << std::setw(20) << rwaFullParams[i]
                           << std::scientific << std::setw(20) << 0.0 << std::endl;
                }
              }
            }
          }

          delete localCompound;
          delete localData;
        }
      }
    } catch (const std::exception& e) {
      // Ignore errors in output file writing to not interrupt sampling
    }
  }
}

// Progress callback for GUI updates
void mcmc_progress_callback(int current_step, int total_steps, double log_prob) {
  g_current_step = current_step;
  
  // No longer call g_gui_progress_callback here - let the iteration callback handle all GUI updates
  // This prevents double updates and flickering values
  
  // Call results callback every 10 steps for results tab updates
  if (g_gui_results_callback && g_mcmc_sampler && (current_step % 10 == 0 || current_step == total_steps)) {
    try {
      auto current_samples = g_mcmc_sampler->get_chain();
      g_gui_results_callback(current_step, total_steps, current_samples);
    } catch (const std::exception& e) {
      // Ignore errors in results callback to not interrupt sampling
    }
  }
  
  // Output file writing moved to iteration callback to ensure it runs
}

// Stop callback for early termination
bool mcmc_should_stop() {
  return g_stop_requested;
}

// Functions to control the stop flag
void mcmc_request_stop() {
  g_stop_requested = true;
}

void mcmc_clear_stop() {
  g_stop_requested = false;
}

// Static methods for external access
void AZURECalcMCMC::RequestStop() {
  mcmc_request_stop();
}

void AZURECalcMCMC::ClearStop() {
  mcmc_clear_stop();
}

void AZURECalcMCMC::SetGUIProgressCallback(void (*callback)(int, int, double, double, double)) {
  g_gui_progress_callback = callback;
}

void AZURECalcMCMC::SetGUIIterationCallback(void (*callback)(int, int)) {
  g_gui_iteration_callback = callback;
}

void AZURECalcMCMC::SetGUIResultsCallback(void (*callback)(int, int, const std::vector<std::vector<double>>&)) {
  g_gui_results_callback = callback;
}

void AZURECalcMCMC::RunMCMCSampling(int nwalkers, int nsteps, const std::vector<double>& initialParams,
                                   std::vector<std::vector<double>>& samples, double chainSpreadPercent, int nthreads, bool useRWA,
                                   double energySpreadKeV) const {
  try {
    int ndim = initialParams.size();

    // Update parameter vectors if needed
    const_cast<AZURECalcMCMC*>(this)->UpdateParameterVectors(initialParams);

    // Classify the varying parameters and derive the normalization and
    // energy-shift priors from the data. Must follow UpdateParameterVectors()
    // (it needs fixed_) and precede any log-probability evaluation.
    BuildAutoPriors();

    if(ndim == 0) {
      configure().outStream << "Error: No parameters provided for MCMC sampling\n";
      return;
    }

    // Set up CSV file for saving samples with enhanced format
    std::string samplesFile = configure().outputdir + "samples.mcmc";
    std::ofstream csvFile;
    bool resuming = false;
    int startStep = 0;
    
    // Check if samples file exists for resuming
    std::ifstream existingFile(samplesFile);
    if(existingFile.good()) {
      configure().outStream << "Found existing samples file: " << samplesFile << "\n";
      configure().outStream << "Counting existing samples for resume capability...\n";
      
      std::string line;
      int lineCount = 0;
      bool hasHeader = false;
      while(std::getline(existingFile, line)) {
        if(!hasHeader && line.find("step,walker,logprob,loglikelihood,logprior") == 0) {
          hasHeader = true;
          continue;
        }
        if(hasHeader) lineCount++;
      }
      existingFile.close();
      
      if(lineCount > 0) {
        int samplesPerStep = nwalkers;
        startStep = lineCount / samplesPerStep;
        
        configure().outStream << "Found " << lineCount << " existing samples (" << startStep << " steps)\n";
        
        if(startStep >= nsteps) {
          configure().outStream << "Sampling already completed! Loading existing samples.\n";
          LoadExistingSamples(samplesFile, samples);
          return;
        }
        
        resuming = true;
        configure().outStream << "Resuming from step " << startStep << "\n";
        csvFile.open(samplesFile, std::ios::app);
      }
    }
    
    if(!resuming) {
      // Create new CSV file with enhanced header
      csvFile.open(samplesFile, std::ios::out);
      csvFile << "step,walker,logprob,loglikelihood,logprior";
      for(int i = 0; i < ndim; i++) {
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

    for(int j = 0; j < ndim; j++) {
      double relative = fabs(initialParams[j]) * spreadFraction;
      double spread = relative;

      const int kind = (j < (int)freeKinds_.size()) ? freeKinds_[j] : -1;

      if(kind == PARAM_ENERGY) {
        spread = energySpreadMeV;
        nEnergyCapped++;
      } else if(kind == PARAM_NORM || kind == PARAM_SHIFT) {
        // Starting outside the prior only wastes steps walking back into it.
        if(j < (int)priorStds_.size() && j < (int)usePriors_.size() &&
           usePriors_[j] && priorStds_[j] > 0.0 && priorStds_[j] < spread) {
          spread = priorStds_[j];
          nPriorCapped++;
        }
      }

      if(spread <= 0.0) spread = 1e-6;
      paramSpread[j] = spread;
    }

    // Resuming continues the ensemble where it stopped, when that state is
    // available and matches this model; otherwise fall back to a fresh ball.
    const std::string walkerStateFile = configure().outputdir + "walkers.mcmc";
    if(resuming && load_walker_state(walkerStateFile, nwalkers, ndim, initial_positions)) {
      configure().outStream << "Resumed walker positions from " << walkerStateFile << "\n";
    } else {
      if(resuming) {
        configure().outStream << "Warning: no usable walker state in " << walkerStateFile
                              << "; the walkers restart from a fresh spread, so the "
                              << "steps that follow are a new burn-in rather than a "
                              << "continuation of the earlier chain.\n";
      }
      if(nEnergyCapped > 0) {
        configure().outStream << "Initial spread: " << nEnergyCapped
                              << " level-energy parameter" << (nEnergyCapped == 1 ? "" : "s")
                              << " limited to " << energySpreadKeV << " keV";
        if(nPriorCapped > 0) {
          configure().outStream << ", " << nPriorCapped
                                << " norm/shift parameter" << (nPriorCapped == 1 ? "" : "s")
                                << " limited to their prior width";
        }
        configure().outStream << "\n";
      }

      for(int i = 0; i < nwalkers; i++) {
        std::vector<double> pos = initialParams;
        for(int j = 0; j < ndim; j++) {
          std::normal_distribution<double> init_distribution(0.0, paramSpread[j]);
          pos[j] += init_distribution(generator);
        }
        initial_positions.push_back(pos);
      }
    }
    
    // Set up global variables for sampling
    g_mcmc_calc = this;
    g_sample_file = &csvFile;
    g_mcmc_for_callbacks = this;
    g_sample_count = resuming ? startStep * nwalkers : 0;
    g_nwalkers = nwalkers;
    g_current_step = startStep;
    g_samples_written_to_file = resuming ? startStep * nwalkers : 0; // Initialize sample counter
    g_stop_requested = false; // Clear stop flag
    g_last_likelihood = 0.0;
    g_last_prior = 0.0;
    g_last_logprob = -std::numeric_limits<double>::infinity();
    g_step_offset = startStep; // Continue the chain numbering when resuming
    g_using_rwa_parameters = useRWA; // Set parameter type flag

    // Create numcmc sampler
    nu::Mcmc mcmc_sampler(nwalkers, ndim, initial_positions);
    // Deterministic given the starting step, so a run is reproducible while a
    // resumed run still draws a fresh stretch of the random stream.
    mcmc_sampler.set_seed(0x9E3779B97F4A7C15ull ^ (std::uint64_t)(startStep + 1));
    g_mcmc_sampler = &mcmc_sampler;

    configure().outStream << "Starting MCMC sampling with enhanced progress tracking";
    if (nthreads > 1) {
        configure().outStream << " using " << nthreads << " parallel threads";
    }
    configure().outStream << "...\n";
    configure().outStream.flush();
    
    int remainingSteps = nsteps - startStep;
    
    // Use enhanced run method with callbacks (parallel or serial based on nthreads)
    // Choose appropriate wrapper function based on parameter type
    auto wrapper_func = useRWA ? mcmc_log_probability_wrapper_rwa : mcmc_log_probability_wrapper;
    
    configure().outStream << "Using " << (useRWA ? "reduced width amplitudes (RWA)" : "physical parameters") << " for MCMC fitting\n";
    
    int result;
    if (nthreads > 1) {
        result = mcmc_sampler.run_parallel_with_callback(
            wrapper_func,
            remainingSteps,
            nthreads,
            mcmc_progress_callback,
            mcmc_should_stop,
            mcmc_iteration_callback,
            mcmc_sample_callback
        );
    } else {
        result = mcmc_sampler.run_with_callback(
            wrapper_func,
            remainingSteps,
            mcmc_progress_callback,
            mcmc_should_stop,
            mcmc_iteration_callback,
            mcmc_sample_callback
        );
    }

    // Save the ensemble so a later run resumes it rather than re-scattering.
    // Also written when the user stopped early: that is precisely the case where
    // continuing from the same positions matters.
    {
      std::vector<nu::Walker> finalWalkers = mcmc_sampler.get_walkers();
      std::vector<std::vector<double>> finalPositions;
      for(size_t k = 0; k < finalWalkers.size(); k++) {
        finalPositions.push_back(finalWalkers[k].getPos());
      }
      if(!save_walker_state(walkerStateFile, finalPositions)) {
        configure().outStream << "Warning: could not write walker state to "
                              << walkerStateFile << "; a later resume will restart "
                              << "the walkers instead of continuing them.\n";
      }
    }

    if(mcmc_sampler.num_invalid_initial() > 0) {
      configure().outStream << "Warning: " << mcmc_sampler.num_invalid_initial()
                            << " of " << nwalkers << " walkers started at a position "
                            << "with zero posterior probability and can never move. "
                            << "Reduce the initial spread or check the priors.\n";
    }
    configure().outStream << "Acceptance fraction: "
                          << std::fixed << std::setprecision(3)
                          << mcmc_sampler.acceptance_fraction()
                          << std::defaultfloat << "\n";
    if(mcmc_sampler.acceptance_fraction() < 0.05) {
      configure().outStream << "  Warning: a very low acceptance fraction means the "
                            << "chain is barely moving; the samples are not yet a "
                            << "usable posterior.\n";
    }
    configure().outStream.flush();

    // Clean up global pointers
    g_mcmc_calc = nullptr;
    g_sample_file = nullptr;
    g_mcmc_for_callbacks = nullptr;
    g_mcmc_sampler = nullptr;
    g_samples_written_to_file = 0; // Reset sample counter
    g_using_rwa_parameters = false; // Reset parameter type flag
    
    csvFile.close();
    
    if(result == 1) {
      configure().outStream << "MCMC sampling stopped early by user request.\n";
    } else if(result != 0) {
      configure().outStream << "MCMC sampling failed with error code: " << result << "\n";
      return;
    }
    
    // Get samples from the sampler
    std::vector<std::vector<double>> newSamples = mcmc_sampler.get_chain();
    
    // Return all samples
    if(resuming) {
      std::vector<std::vector<double>> existingSamples;
      LoadExistingSamples(samplesFile, existingSamples);
      samples = existingSamples;
    } else {
      samples = newSamples;
    }
    
    configure().outStream << "MCMC sampling completed successfully. Generated " 
                         << samples.size() << " total samples.\n";
    configure().outStream << "Samples saved to: " << samplesFile << "\n";
    
  } catch(const std::exception& e) {
    configure().outStream << "Error during MCMC sampling: " << e.what() << "\n";
    g_mcmc_calc = nullptr;
    g_sample_file = nullptr;
    g_using_rwa_parameters = false;
  } catch(...) {
    configure().outStream << "Unknown error during MCMC sampling\n";
    g_mcmc_calc = nullptr;
    g_sample_file = nullptr;
    g_using_rwa_parameters = false;
  }
}

void AZURECalcMCMC::LoadExistingSamples(const std::string& filename, std::vector<std::vector<double>>& samples) const {
  samples.clear();
  
  std::ifstream file(filename);
  if(!file.good()) {
    configure().outStream << "Warning: Could not open samples file for loading: " << filename << "\n";
    return;
  }
  
  std::string line;
  bool headerSkipped = false;
  
  while(std::getline(file, line)) {
    // Skip header line
    if(!headerSkipped && line.find("step,walker,") == 0) {
      headerSkipped = true;
      continue;
    }
    
    if(!headerSkipped) continue;
    
    // Parse CSV line: step,walker,logprob,param0,param1,...
    std::vector<double> sample;
    std::stringstream ss(line);
    std::string cell;
    
    int colIndex = 0;
    while(std::getline(ss, cell, ',')) {
      if(colIndex >= 3) { // Skip step, walker, logprob columns
        try {
          double value = robust_stod(cell);
          sample.push_back(value);
        } catch(const std::exception&) {
          // Skip malformed lines
          break;
        }
      }
      colIndex++;
    }
    
    if(!sample.empty()) {
      samples.push_back(sample);
    }
  }
  
  configure().outStream << "Loaded " << samples.size() << " existing samples from " << filename << "\n";
}

/*!
 * Initialize object pools with pre-allocated CNuc and EData objects
 */
void AZURECalcMCMC::InitializePools() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);
  if (pools_initialized_) return;
  
  // Calculate pool size based on available hardware threads
  const int pool_size = std::max(4, static_cast<int>(std::thread::hardware_concurrency() * 2));
  
  // Pre-allocate CNuc objects by cloning once
  for (int i = 0; i < pool_size; ++i) {
    cnuc_pool_.push(std::unique_ptr<CNuc>(compound_->Clone()));
  }
  
  // Pre-allocate EData objects by cloning once
  for (int i = 0; i < pool_size; ++i) {
    edata_pool_.push(std::unique_ptr<EData>(data_->Clone()));
  }
  
  pools_initialized_ = true;
}

/*!
 * Get a CNuc object from the pool, creating new if pool is empty
 */
CNuc* AZURECalcMCMC::GetPooledCNuc() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);
  
  if (!cnuc_pool_.empty()) {
    auto obj = std::move(cnuc_pool_.top());
    cnuc_pool_.pop();
    return obj.release();
  }
  
  // Fallback: create new if pool is empty (shouldn't happen often)
  return compound_->Clone();
}

/*!
 * Get an EData object from the pool, creating new if pool is empty
 */
EData* AZURECalcMCMC::GetPooledEData() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);
  
  if (!edata_pool_.empty()) {
    auto obj = std::move(edata_pool_.top());
    edata_pool_.pop();
    return obj.release();
  }
  
  // Fallback: create new if pool is empty (shouldn't happen often)
  return data_->Clone();
}

/*!
 * Return a CNuc object to the pool for reuse
 */
void AZURECalcMCMC::ReturnPooledCNuc(CNuc* obj) const {
  if (!obj) return;
  
  std::lock_guard<std::mutex> lock(pool_mutex_);
  cnuc_pool_.push(std::unique_ptr<CNuc>(obj));
}

/*!
 * Return an EData object to the pool for reuse  
 */
void AZURECalcMCMC::ReturnPooledEData(EData* obj) const {
  if (!obj) return;
  
  std::lock_guard<std::mutex> lock(pool_mutex_);
  edata_pool_.push(std::unique_ptr<EData>(obj));
}
