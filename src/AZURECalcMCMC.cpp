#include "AZURECalcMCMC.h"
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
// Store the last calculated likelihood and prior for progress callbacks
static double g_last_likelihood = 0.0;
static double g_last_prior = 0.0;
static std::vector<double> g_last_params;
// Mutex for thread-safe file writing
#ifdef _OPENMP
#include <omp.h>
static omp_lock_t g_file_lock;
#endif

double AZURECalcMCMC::CalculateLogLikelihood(const vector_r& p) const {

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  
  try {
    localCompound = compound()->Clone();
    localData = data()->Clone();

    //Fill Compound Nucleus From Parameters
    AZUREParams params;
    localCompound->FillCompoundFromParams(p);
    //localData->FillNormsFromParams(p);
    localData->FillEnergyShiftsFromParams(p,localData,localCompound,&configure());
    if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());
  } catch (GSLException& e) {
    // Clean up and return bad likelihood for GSL errors
    if(localCompound) delete localCompound;
    if(localData) delete localData;
    return -std::numeric_limits<double>::infinity();
  } catch (...) {
    // Clean up and return bad likelihood for any other errors
    if(localCompound) delete localCompound;
    if(localData) delete localData;
    return -std::numeric_limits<double>::infinity();
  }

  bool isFit=true;
  
  // Calculate chi-squared (same logic as original but without parameter filling)
  double chiSquared=0.0;
  double segmentChiSquared=0.0;
  ESegmentIterator firstSumIterator = localData->GetSegments().end();
  ESegmentIterator lastSumIterator = localData->GetSegments().end();
  
  try {
    for(EDataIterator data=localData->begin();data!=localData->end();data++) {
      if(data.segment()->GetPoints().begin()==data.point()) {
        segmentChiSquared=0.0;
        if(data.segment()->IsTotalCapture()) {
          firstSumIterator=data.segment();
          lastSumIterator=data.segment()+data.segment()->IsTotalCapture()-1;
        } 
      }
      if(!data.point()->IsMapped()) {
        try {
          data.point()->Calculate(localCompound,configure());
        } catch (GSLException& e) {
          // Skip this point if GSL calculation fails
          continue;
        }
      }
    if(firstSumIterator!=localData->GetSegments().end()&&
       data.segment()!=lastSumIterator) continue;
    double fitCrossSection=data.point()->GetFitCrossSection();
    ESegmentIterator thisSegment = data.segment();
    if(data.segment()==lastSumIterator) {
      int pointIndex=data.point()-data.segment()->GetPoints().begin()+1;
      for(ESegmentIterator it=firstSumIterator;it<data.segment();it++) 
        fitCrossSection+=it->GetPoint(pointIndex)->GetFitCrossSection();
      thisSegment = firstSumIterator;
    }
    double dataNorm=thisSegment->GetNorm();
    double CrossSection=data.point()->GetCMCrossSection()*dataNorm;
    double CrossSectionError=data.point()->GetCMCrossSectionError()*dataNorm;
    double chi=(fitCrossSection-CrossSection)/CrossSectionError;
    double pointChiSquared=pow(chi,2.0);
    segmentChiSquared+=pointChiSquared;
    if(data.segment()->GetPoints().end()-1==data.point()) {
      if(!isFit) thisSegment->SetSegmentChiSquared(segmentChiSquared);
      if(data.segment()==lastSumIterator) {
        firstSumIterator=localData->GetSegments().end();
        lastSumIterator=localData->GetSegments().end();
      }
      chiSquared+=segmentChiSquared;
    }

    }
  } catch (GSLException& e) {
    // Clean up and return bad likelihood for GSL errors during calculation
    delete localCompound;
    delete localData;
    return -std::numeric_limits<double>::infinity();
  } catch (...) {
    // Clean up and return bad likelihood for any other errors during calculation
    delete localCompound;
    delete localData;
    return -std::numeric_limits<double>::infinity();
  }

  std::cout << "Chi-squared: " << chiSquared << std::endl;

  delete localCompound;
  delete localData;

  // Convert chi-squared to log-likelihood: ln(L) = -0.5 * chi^2
  return -0.5 * chiSquared;
  
}

double AZURECalcMCMC::CalculateLogLikelihoodPhysical(const vector_r& params_) const {

  CNuc * localCompound = NULL;
  EData *localData = NULL;
  AZUREParams params;
  
  try {
    localCompound = compound()->Clone();
    localData = data()->Clone();

    localCompound->FillCompoundFromParamsPhysical(params_);
    bool isValid = localCompound->TransformIn( configure( ) );
    
    if(!isValid) {
      delete localCompound;
      delete localData;
      return -std::numeric_limits<double>::infinity();
    }
  } catch (GSLException& e) {
    // Clean up and return bad likelihood for GSL errors
    if(localCompound) delete localCompound;
    if(localData) delete localData;
    return -std::numeric_limits<double>::infinity();
  } catch (...) {
    // Clean up and return bad likelihood for any other errors
    if(localCompound) delete localCompound;
    if(localData) delete localData;
    return -std::numeric_limits<double>::infinity();
  }

  localCompound->FillMnParams(params.GetMinuitParams());
  localData->FillMnParams(params.GetMinuitParams());
  localData->FillEnergyShiftsFromParams(params_,localData,localCompound,&configure());
  localCompound->FillCompoundFromParams(params.GetMinuitParams( ).Params( ));

  //Fill Compound Nucleus From Minuit Parameters
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());
  
  bool isFit = true; // For MCMC, we are always fitting
  
  // Calculate chi-squared (same logic as original but without parameter filling)
  double chiSquared=0.0;
  double segmentChiSquared=0.0;
  ESegmentIterator firstSumIterator = localData->GetSegments().end();
  ESegmentIterator lastSumIterator = localData->GetSegments().end();
  
  try {
    for(EDataIterator data=localData->begin();data!=localData->end();data++) {
      if(data.segment()->GetPoints().begin()==data.point()) {
        segmentChiSquared=0.0;
        if(data.segment()->IsTotalCapture()) {
          firstSumIterator=data.segment();
          lastSumIterator=data.segment()+data.segment()->IsTotalCapture()-1;
        } 
      }
      if(!data.point()->IsMapped()) {
        try {
          data.point()->Calculate(localCompound,configure());
        } catch (GSLException& e) {
          // Skip this point if GSL calculation fails
          continue;
        }
      }
    if(firstSumIterator!=localData->GetSegments().end()&&
       data.segment()!=lastSumIterator) continue;
    double fitCrossSection=data.point()->GetFitCrossSection();
    ESegmentIterator thisSegment = data.segment();
    if(data.segment()==lastSumIterator) {
      int pointIndex=data.point()-data.segment()->GetPoints().begin()+1;
      for(ESegmentIterator it=firstSumIterator;it<data.segment();it++) 
        fitCrossSection+=it->GetPoint(pointIndex)->GetFitCrossSection();
      thisSegment = firstSumIterator;
    }
    double dataNorm=thisSegment->GetNorm();
    double CrossSection=data.point()->GetCMCrossSection()*dataNorm;
    double CrossSectionError=data.point()->GetCMCrossSectionError()*dataNorm;
    double chi=(fitCrossSection-CrossSection)/CrossSectionError;
    double pointChiSquared=pow(chi,2.0);
    segmentChiSquared+=pointChiSquared;
    if(data.segment()->GetPoints().end()-1==data.point()) {
      if(!isFit) thisSegment->SetSegmentChiSquared(segmentChiSquared);
      if(data.segment()==lastSumIterator) {
        firstSumIterator=localData->GetSegments().end();
        lastSumIterator=localData->GetSegments().end();
      }
      chiSquared+=segmentChiSquared;
    }

    }
  } catch (GSLException& e) {
    // Clean up and return bad likelihood for GSL errors during calculation
    delete localCompound;
    delete localData;
    return -std::numeric_limits<double>::infinity();
  } catch (...) {
    // Clean up and return bad likelihood for any other errors during calculation
    delete localCompound;
    delete localData;
    return -std::numeric_limits<double>::infinity();
  }

  delete localCompound;
  delete localData;

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
  
  AZUREParams params;
  compound()->FillMnParams(params.GetMinuitParams());
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

// C-style function wrapper for numcmc library - uses RWA parameters
double mcmc_log_probability_wrapper_rwa(std::vector<double>& params) {
  if(g_mcmc_calc) {
    // Use the new LogLikelihood function for RWA parameters
    double logLikelihood = g_mcmc_calc->LogLikelihood(params);
    double logPrior = g_mcmc_calc->CalculateLogPrior(params);
    double logProbability = logLikelihood + logPrior;
    
    // Store for GUI callbacks - these are the EXACT values used
    g_last_likelihood = logLikelihood;
    g_last_prior = logPrior;
    g_last_params = params;
    
    // Save sample to file if file is open (thread-safe)
    if(g_sample_file && g_sample_file->is_open()) {
#ifdef _OPENMP
      omp_set_lock(&g_file_lock);
#endif
      int walker_id = g_sample_count % g_nwalkers;
      int step = g_sample_count / g_nwalkers;
      
      *g_sample_file << step << "," << walker_id << "," << logProbability << "," << logLikelihood << "," << logPrior;
      for(const double& param : params) {
        *g_sample_file << "," << param;
      }
      *g_sample_file << "\n";
      g_sample_file->flush(); // Ensure data is written immediately
      g_sample_count++;
      
      // Increment the atomic counter AFTER the sample is actually written to file
      g_samples_written_to_file++;
#ifdef _OPENMP
      omp_unset_lock(&g_file_lock);
#endif
    }
    
    return logLikelihood;  // Return only likelihood for MCMC sampler
  }
  return -std::numeric_limits<double>::infinity();
}

// C-style function wrapper for numcmc library - uses physical parameters
double mcmc_log_probability_wrapper(std::vector<double>& params) {
  if(g_mcmc_calc) {

    double logLikelihood = g_mcmc_calc->LogLikelihoodPhysical(params);
    double logPrior = g_mcmc_calc->CalculateLogPrior(params);
    double logProbability = logLikelihood + logPrior;
    
    // Store for GUI callbacks - these are the EXACT values used
    g_last_likelihood = logLikelihood;
    g_last_prior = logPrior;
    g_last_params = params;
    
    // Save sample to file if file is open (thread-safe)
    if(g_sample_file && g_sample_file->is_open()) {
#ifdef _OPENMP
      omp_set_lock(&g_file_lock);
#endif
      int walker_id = g_sample_count % g_nwalkers;
      int step = g_sample_count / g_nwalkers;
      
      *g_sample_file << step << "," << walker_id << "," << logProbability << "," << logLikelihood << "," << logPrior;
      for(const double& param : params) {
        *g_sample_file << "," << param;
      }
      *g_sample_file << "\n";
      g_sample_file->flush(); // Ensure data is written immediately
      g_sample_count++;
      
      // Increment the atomic counter AFTER the sample is actually written to file
      g_samples_written_to_file++;
#ifdef _OPENMP
      omp_unset_lock(&g_file_lock);
#endif
    }
    
    return logLikelihood;  // Return only likelihood for MCMC sampler
  }
  return -std::numeric_limits<double>::infinity();
}

// Global pointer to the GUI progress handler
static void (*g_gui_progress_callback)(int, int, double, double, double) = nullptr;
static void (*g_gui_iteration_callback)(int, int) = nullptr;
static void (*g_gui_results_callback)(int, int, const std::vector<std::vector<double>>&) = nullptr;
static const AZURECalcMCMC* g_mcmc_for_callbacks = nullptr;
static nu::Mcmc* g_mcmc_sampler = nullptr;
static bool g_using_rwa_parameters = false;

// Helper function to read last line of file and extract values (thread-safe)
struct LastSampleInfo {
  int step;
  int walker;
  double logprob;
  double loglikelihood;
  double logprior;
  bool valid;
};

LastSampleInfo read_last_sample_from_file(const std::string& filename) {
  LastSampleInfo info = {0, 0, 0.0, 0.0, 0.0, false};
  
#ifdef _OPENMP
  omp_set_lock(&g_file_lock);
#endif
  
  std::ifstream file(filename);
  if (!file.is_open()) {
#ifdef _OPENMP
    omp_unset_lock(&g_file_lock);
#endif
    return info;
  }
  
  // Seek to end and read backwards to find last line
  file.seekg(-1, std::ios::end);
  if (file.tellg() <= 0) {
    file.close();
#ifdef _OPENMP
    omp_unset_lock(&g_file_lock);
#endif
    return info;
  }
  
  // Skip any trailing newlines
  char ch;
  while (file.tellg() > 0) {
    file.get(ch);
    if (ch != '\n' && ch != '\r') {
      file.seekg(-1, std::ios::cur);
      break;
    }
    file.seekg(-2, std::ios::cur);
  }
  
  // Read backwards to find start of last line
  while (file.tellg() > 0) {
    file.seekg(-1, std::ios::cur);
    file.get(ch);
    if (ch == '\n') {
      break;
    }
    file.seekg(-1, std::ios::cur);
  }
  
  // Read the last line
  std::string lastLine;
  std::getline(file, lastLine);
  file.close();
  
#ifdef _OPENMP
  omp_unset_lock(&g_file_lock);
#endif
  
  // Parse CSV: step,walker,logprob,loglikelihood,logprior,param1,param2,...
  if (!lastLine.empty() && lastLine.find("step,walker") == std::string::npos) { // Skip header
    std::stringstream ss(lastLine);
    std::string token;
    int field = 0;
    
    while (std::getline(ss, token, ',') && field < 5) {
      try {
        switch (field) {
          case 0: info.step = std::stoi(token); break;
          case 1: info.walker = std::stoi(token); break;
          case 2: info.logprob = std::stod(token); break;
          case 3: info.loglikelihood = std::stod(token); break;
          case 4: info.logprior = std::stod(token); break;
        }
      } catch (...) {
        return info; // Invalid format
      }
      field++;
    }
    
    if (field >= 5) {
      info.valid = true;
    }
  }
  
  return info;
}

// Iteration callback for GUI updates (called every iteration)
void mcmc_iteration_callback(int current_step, int total_steps) {
  if (g_gui_iteration_callback) {
    // Read the last line from the file to get current step and probability values
    std::string samplesFile = "";
    if (g_mcmc_for_callbacks) {
      samplesFile = g_mcmc_for_callbacks->configure().outputdir + "samples.mcmc";
    }
    
    LastSampleInfo lastSample = {0, 0, 0.0, 0.0, 0.0, false};
    if (!samplesFile.empty()) {
      lastSample = read_last_sample_from_file(samplesFile);
    }
    
    // Use the step from the file for iteration display
    int current_iteration = lastSample.valid ? lastSample.step : current_step;
    int total_expected_samples = total_steps * g_nwalkers;
    
    g_gui_iteration_callback(current_iteration, total_expected_samples);
    
    // For progress callback (time estimation), use the original numcmc step values
    // This keeps the time estimation logic working correctly
    if (g_gui_progress_callback && lastSample.valid) {
      g_gui_progress_callback(current_step, total_steps, 
                              lastSample.logprob, lastSample.loglikelihood, lastSample.logprior);
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
          
          if (g_using_rwa_parameters) {
            // For RWA parameters, use the same approach as CalculateLogLikelihood
            localCompound->FillCompoundFromParams(fullParams);
            localData->FillNormsFromParams(fullParams);
            localData->FillEnergyShiftsFromParams(fullParams, localData, localCompound, &g_mcmc_for_callbacks->configure());
            if (g_mcmc_for_callbacks->configure().paramMask & Config::USE_BRUNE_FORMALISM) {
              localCompound->CalcShiftFunctions(g_mcmc_for_callbacks->configure());
            }
          } else {
            // For physical parameters, use the same approach as CalculateLogLikelihoodPhysical
            localCompound->FillCompoundFromParamsPhysical(fullParams);
            isValid = localCompound->TransformIn(g_mcmc_for_callbacks->configure());
            
            if (isValid) {
              AZUREParams params;
              localCompound->FillMnParams(params.GetMinuitParams());
              localData->FillMnParams(params.GetMinuitParams());
              localData->FillEnergyShiftsFromParams(fullParams, localData, localCompound, &g_mcmc_for_callbacks->configure());
              localCompound->FillCompoundFromParams(params.GetMinuitParams().Params());
              
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

bool AZURECalcMCMC::Initialize( ){

  std::string file;
  if( configure().paramMask & Config::CALCULATE_WITH_DATA ) file = configure().outputdir + "intEC.dat";
  else file = configure().outputdir + "intEC.extrap";

  // Initialize EC Integral caching system

  std::string cacheFile;
  if (configure().paramMask & Config::CALCULATE_WITH_DATA) {
    cacheFile = configure().outputdir + "intEC_cache.dat";
  } else {
    cacheFile = configure().outputdir + "intEC_cache.extrap";
  }

  // FIXME: It crashes on Linux (but fine on Mac)
  //if( compound_ != nullptr ) delete compound_;
  //if( data_ != nullptr ) delete data_;

  delete compound_;
  delete data_;
  data_ = new EData( );
  compound_ = new CNuc( );

  //configure().outStream << "Filling Compound Nucleus..." << std::endl;
  if(compound()->Fill(configure())==-1) {
    //configure().outStream << "Could not fill compound nucleus from file." << std::endl;
    return -1;
  } else if(compound()->NumPairs()==0 || compound()->NumJGroups()==0) {
    //configure().outStream << "No nuclear data exists. Calculation not possible." << std::endl; 
    return -1;
  } 
  if((configure().screenCheckMask|configure().fileCheckMask) & 
     Config::CHECK_COMPOUND_NUCLEUS) compound()->PrintNuc(configure());

  if(!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
    //Fill the data object from the segments and data file
    //  Compound object is passed to the function for pair key verification and
    //  center of mass conversions, s-factor conversions, etc.
    //configure().outStream << "Filling Data Structures..." << std::endl;
    if(configure().paramMask & Config::CALCULATE_WITH_DATA) {
      if(data()->Fill(configure(),compound())==-1) {
	//configure().outStream << "Could not fill data object from file." << std::endl;
	return -1;
      } else if(data()->NumSegments()==0) {
	//configure().outStream << "There is no data provided." << std::endl;
	return -1;
      }
    } else {
      if(data()->MakePoints(configure(),compound())==-1) {
	//configure().outStream << "Could not fill data object from file." << std::endl;
	return -1;
      } else if(data()->NumSegments()==0) {
	//configure().outStream << "Extrapolation segments produce no data." << std::endl;
	return -1;
      }
    } 
    if((configure().fileCheckMask|configure().screenCheckMask) & Config::CHECK_DATA)
      data()->PrintData(configure());
  } else {
    if(!compound()->IsPairKey(configure().rateParams.entrancePair)||!compound()->IsPairKey(configure().rateParams.exitPair)) {
      //configure().outStream << "Reaction rate pairs do not exist in compound nucleus." << std::endl;
      return -1;
    } else {
      compound()->GetPair(compound()->GetPairNumFromKey(configure().rateParams.entrancePair))->SetEntrance();
    }
  }

  //Initialize compound nucleus object
  try {
    compound()->Initialize(configure());
  } catch (GSLException e) {
    //configure().outStream << e.what() << std::endl;
    //configure().outStream << "Calculation was aborted." << std::endl;
    return -1;
  }

  // Create new parameters for minuit, fill them from compound nucleus object and data file.
  AZUREParams params;
  compound()->FillMnParams(params.GetMinuitParams());
  data()->FillMnParams(params.GetMinuitParams());
  if(!(configure().paramMask & Config::USE_PREVIOUS_PARAMETERS)) {
    //configure().outStream << "Creating New param.par File..." << std::endl;
    params.WriteUserParameters(configure(),false);
  } else {
    //configure().outStream << "Reading User Parameter File..." << std::endl;
    params.ReadUserParameters(configure());
  }

  if(data()->Initialize(compound(),configure())==-1) return -1;

  return 0;
  
}

void AZURECalcMCMC::RunMCMCSampling(int nwalkers, int nsteps, const std::vector<double>& initialParams, 
                                   std::vector<std::vector<double>>& samples, double chainSpreadPercent, int nthreads, bool useRWA) const {
  try {
    int ndim = initialParams.size();

    // Update parameter vectors if needed
    const_cast<AZURECalcMCMC*>(this)->UpdateParameterVectors(initialParams);

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

    std::cout << "Spread Fraction: " << spreadFraction << std::endl;

    for(int i = 0; i < nwalkers; i++) {
      std::vector<double> pos = initialParams;
      for(int j = 0; j < ndim; j++) {
        double paramValue = initialParams[j];
        double spread = fabs(paramValue) * spreadFraction;
        if(spread == 0.0) spread = 0.01;
        std::normal_distribution<double> init_distribution(0.0, spread);
        pos[j] += init_distribution(generator);
      }
      initial_positions.push_back(pos);
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
    g_using_rwa_parameters = useRWA; // Set parameter type flag
    
    // Initialize OpenMP lock for thread-safe file writing
#ifdef _OPENMP
    omp_init_lock(&g_file_lock);
#endif
    
    // Create numcmc sampler
    nu::Mcmc mcmc_sampler(nwalkers, ndim, initial_positions);
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
            mcmc_iteration_callback
        );
    } else {
        result = mcmc_sampler.run_with_callback(
            wrapper_func, 
            remainingSteps,
            mcmc_progress_callback,
            mcmc_should_stop,
            mcmc_iteration_callback
        );
    }
    
    // Clean up global pointers and OpenMP lock
#ifdef _OPENMP
    omp_destroy_lock(&g_file_lock);
#endif
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
          double value = std::stod(cell);
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
