#ifndef AZURECALCMCMC_H
#define AZURECALCMCMC_H

#ifdef USE_MCMC
#include "numcmc/mcmc.h"
#endif
#include "Minuit2/FCNBase.h"
#include "Constants.h"
#include <vector>

class Config;
class EData;
class CNuc;
class ParameterLimitsManager;

///A function class to perform MCMC Bayesian calculation of parameters

/*!
 * The AZURECalcMCMC function class calculates the log-likelihood based on a 
 * parameter set for all available data, and returns a log-probability value.
 * This function class is what the MCMC sampler calls repeatedly during the
 * Bayesian sampling process to explore the parameter space.
 */

class AZURECalcMCMC : public ROOT::Minuit2::FCNBase {
 public:
  /*!
   * The AZURECalcMCMC object is created with reference to an EData and CNuc object.
   * The runtime configurations are also passed through a Config structure.
   */
  AZURECalcMCMC(EData* data, CNuc* compound, const Config& configure, ParameterLimitsManager* limitsManager = nullptr) 
    : configure_(configure), parametersInitialized_(false) {
    data_=data;
    compound_=compound;
    limitsManager_=limitsManager;
  };
  
  ~AZURECalcMCMC() {};
  /*!
   * See Minuit2 documentation for an explanation of this function.
   */
  virtual double Up() const {return theErrorDef;};
  /*!
   * Overloaded operator to make the class instance callable as a function. 
   * A parameter array is passed as the dependent variable. The function
   * returns the log-likelihood value for MCMC sampling.
   */
  virtual double operator()(const vector_r&) const;
  
#ifdef USE_MCMC
  
  /*!
   * Log-likelihood function for MCMC sampler with physical parameters (no priors).
   * Handles physical parameters like AZUREAPI::UpdateSegments.
   */
  double LogLikelihoodPhysical(const std::vector<double>& physicalParams) const;
  
  /*!
   * Log-probability function for MCMC sampler with physical parameters.
   * Handles physical parameters like AZUREAPI::UpdateSegments.
   */
  double LogProbabilityPhysical(const std::vector<double>& physicalParams) const;
  
  /*!
   * Run MCMC sampling with specified parameters.
   */
  void RunMCMCSampling(int nwalkers, int nsteps, const std::vector<double>& initialParams, 
                       std::vector<std::vector<double>>& samples, double chainSpreadPercent = 1.0, int nthreads = 1) const;

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
#endif
  
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
   * Calculate nuisance parameter chi-squared contribution
   */
  double CalculateNuisanceChiSquared(const vector_r& p) const;
  
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
   * Reconstruct full parameter array from varying parameters (for output file writing)
   */
  vector_r ReconstructFullParameters(const std::vector<double>& varyingParams) const;


 private:
  const Config &configure_;
  EData *data_;
  CNuc *compound_;
  ParameterLimitsManager *limitsManager_;
  double theErrorDef;
  
  // Parameter vectors for physical parameter handling (mutable for lazy initialization)
  mutable vector_r all_physical_;
  mutable vector_r physical_;
  mutable vector_r rwa_;
  mutable vector_r all_rwa_;
  mutable vector_r all_indexes;
  mutable std::vector<bool> fixed_;
  mutable bool parametersInitialized_;
  
#ifdef USE_MCMC
  // Prior information for Bayesian analysis
  mutable std::vector<double> priorMeans_;
  mutable std::vector<double> priorStds_;
  mutable std::vector<bool> usePriors_;
#endif
};

#endif
