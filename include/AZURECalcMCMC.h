#ifndef AZURECALCMCMC_H
#define AZURECALCMCMC_H

#include "numcmc/mcmc.h"
#include "Constants.h"
#include <vector>

class Config;
class EData;
class CNuc;

///A function class to perform MCMC Bayesian calculation of parameters

/*!
 * The AZURECalcMCMC function class calculates the log-likelihood based on a 
 * parameter set for all available data, and returns a log-probability value.
 * This function class is what the MCMC sampler calls repeatedly during the
 * Bayesian sampling process to explore the parameter space.
 */

class AZURECalcMCMC {
 public:
  /*!
   * The AZURECalcMCMC object is created with reference to an EData and CNuc object.
   * The runtime configurations are also passed through a Config structure.
   */
  AZURECalcMCMC(EData* data, CNuc* compound, const Config& configure) 
    : configure_(configure), parametersInitialized_(false) {
    data_=data;
    compound_=compound;
  };
  
  ~AZURECalcMCMC() {};
  
  /*!
   * Log-likelihood function for MCMC sampler with physical parameters (no priors).
   * Handles physical parameters like AZUREAPI::UpdateSegments.
   */
  double LogLikelihoodPhysical(const std::vector<double>& physicalParams) const;
  
  /*!
   * Log-likelihood function for MCMC sampler with RWA parameters (no priors).
   * Handles RWA parameters directly.
   */
  double LogLikelihood(const std::vector<double>& rwaParams) const;
  
  /*!
   * Log-probability function for MCMC sampler with physical parameters.
   * Handles physical parameters like AZUREAPI::UpdateSegments.
   */
  double LogProbabilityPhysical(const std::vector<double>& physicalParams) const;
  
  /*!
   * Run MCMC sampling with specified parameters.
   */
  void RunMCMCSampling(int nwalkers, int nsteps, const std::vector<double>& initialParams, 
                       std::vector<std::vector<double>>& samples, double chainSpreadPercent = 1.0, int nthreads = 1, bool useRWA = false) const;

  /*!
   * Set prior information for parameters.
   */
  void SetPriors(const std::vector<double>& priorMeans, 
                 const std::vector<double>& priorStds,
                 const std::vector<bool>& usePriors);
  
  /*!
   * Calculate log-prior contribution for given parameters.
   */
  double CalculateLogPrior(const std::vector<double>& params) const;

  /*!
   * Load existing samples from CSV file for resume capability.
   */
  void LoadExistingSamples(const std::string& filename, std::vector<std::vector<double>>& samples) const;

  bool Initialize( );
  
  /*!
   * Request stop for currently running MCMC
   */
  static void RequestStop();
  
  /*!
   * Clear stop flag for new MCMC run
   */
  static void ClearStop();
  
  /*!
   * Set GUI progress callback function
   */
  static void SetGUIProgressCallback(void (*callback)(int, int, double, double, double));
  
  /*!
   * Set GUI iteration callback function (called every iteration)
   */
  static void SetGUIIterationCallback(void (*callback)(int, int));
  
  /*!
   * Set GUI results callback function (called every 100 iterations for results update)
   */
  static void SetGUIResultsCallback(void (*callback)(int, int, const std::vector<std::vector<double>>&));
  
  /*!
   * Returns a reference to the Config structure.
   */
  const Config &configure() const {return configure_;};
  /*!
   * Returns a pointer to the EData object.
   */
  EData *data() const {return data_;};
  /*!
   * Returns a pointer to the CNuc object.
   */
  CNuc *compound() const {return compound_;};
 
  /*!
   * See Minuit2 documentation for an explanation of this function.
   */
  void SetErrorDef(double def) {theErrorDef=def;};
  
  /*!
   * Calculate log-likelihood from chi-squared
   */
  double CalculateLogLikelihood(const vector_r& p) const;
  
  /*!
   * Calculate log-likelihood from physical parameters (like AZUREAPI::UpdateSegments)
   */
  double CalculateLogLikelihoodPhysical(const vector_r& physicalParams) const;
  
  /*!
   * Update parameter vectors for parameter transformation handling
   */
  void UpdateParameterVectors(const vector_r& initialParams) const;
  
  /*!
   * Reconstruct full parameter array from varying physical parameters (for output file writing)
   */
  vector_r ReconstructFullParametersPhysical(const std::vector<double>& varyingParams) const;
  
  /*!
   * Reconstruct full parameter array from varying RWA parameters (for output file writing)
   */
  vector_r ReconstructFullParameters(const std::vector<double>& varyingParams) const;


 private:
  const Config &configure_;
  EData *data_;
  CNuc *compound_;
  double theErrorDef;
  
  // Parameter vectors for physical parameter handling (mutable for lazy initialization)
  mutable vector_r all_physical_;
  mutable vector_r physical_;
  mutable vector_r rwa_;
  mutable vector_r all_rwa_;
  mutable vector_r all_indexes;
  mutable std::vector<bool> fixed_;
  mutable bool parametersInitialized_;
  
  // Prior information for Bayesian analysis
  mutable std::vector<double> priorMeans_;
  mutable std::vector<double> priorStds_;
  mutable std::vector<bool> usePriors_;
};

#endif
