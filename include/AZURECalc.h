#ifndef AZURECALC_H
#define AZURECALC_H

#include "Minuit2/FCNBase.h"
#include "Minuit2/FCNGradientBase.h"
#include "Constants.h"
#include "AZUREParams.h"
#include <vector>
#include <memory>
#include <mutex>
#include <stack>

class Config;
class EData;
class CNuc;
class ParameterLimitsManager;

///A function class to perform the calculation of the chi-squared value

/*!
 * The AZURECalc function class calculates the cross section based on a 
 * parameter set for all available data, and returns a chi-squared value.
 * This function class is what Minuit calls repeatedly during the fitting
 * process to perform the minimization.
 */

class AZURECalc : public ROOT::Minuit2::FCNGradientBase {
 public:
  /*!
   * The AZURECalc object is created with reference to an EData and CNuc object.
   *. The runtime configurations are also passed through a Config structure.
   */
  AZURECalc(EData* data,CNuc* compound, const Config& configure, ParameterLimitsManager* limitsManager = nullptr) : configure_(configure), pools_initialized_(false) {
    data_=data;
    compound_=compound;
    limitsManager_=limitsManager;
  };
  
  ~AZURECalc() {};
  /*!
   * See Minuit2 documentation for an explanation of this function.
   */
  virtual double Up() const {return theErrorDef;};
  /*!
   * Overloaded operator to make the class instance callable as a function. 
   * A Minuit parameter array is passed as the dependent variable.  The function
   * returns the total chi-squared value.
   */
  virtual double operator()(const vector_r&) const;

  /*!
   * Analytic gradient of the total chi-squared with respect to the (external)
   * Minuit parameters, for Minuit2's FCNGradientBase interface.  The energy /
   * reduced-width block is computed by the shared reverse-mode adjoint
   * (AccumulateEGammaGradient -- the same engine AZUREAPI uses); norms, energy
   * shifts and any analytically-unsupported configuration fall back to central
   * finite differences of the chi-squared.
   */
  std::vector<double> Gradient(const std::vector<double>&) const override;
  /*!
   * Skip Minuit's start-up numerical check of the analytic gradient (trusted,
   * separately validated).
   */
  bool CheckGradient() const override { return false; }
  /*!
   * Side-effect-free chi-squared evaluation (no iteration counter / file output
   * / object pools), used by the finite-difference part of Gradient().
   */
  double Chi2Value(const vector_r& p, bool thmOnly = false) const;

  /*!
   * Standardized residuals and analytic residual Jacobian at the full parameter
   * vector `full`.  `residuals` is N_res; `jac` is row-major N_res x nFree;
   * `packedToFull[a]` is the full-vector index of Jacobian column a (the a-th
   * non-fixed parameter).  Returns false if the analytic path is unsupported.
   */
  bool ResidualJacobian(const vector_r& full, vector_r& residuals,
                        vector_r& jac, std::vector<int>& packedToFull) const;

  /*!
   * Levenberg-Marquardt / Gauss-Newton minimization of the total chi-squared
   * using the analytic residual Jacobian.  Converges in far fewer iterations
   * than MIGRAD for least-squares problems.  Updates `params` in place with the
   * best-fit values and (J^T J)^{-1} parameter errors, respecting parameter
   * limits.  Returns the final chi-squared, or a negative value if the analytic
   * Jacobian is unsupported (caller should fall back to MIGRAD).
   */
  double RunLevenbergMarquardt(AZUREParams& params, int maxIter = 200) const;

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
   * Add the analytic gradient of the nuisance-parameter penalty to `grad`,
   * mirroring CalculateNuisanceChiSquared exactly (so it is consistent with the
   * objective).  Used when the energy/gamma block is taken analytically (which
   * covers only the data chi-squared).
   */
  void AddNuisanceGradient(const vector_r& p, std::vector<double>& grad) const;
  
  /*!
   * Object pool management methods
   */
  void InitializePools() const;
  CNuc* GetPooledCNuc() const;
  EData* GetPooledEData() const;
  void ReturnPooledCNuc(CNuc* obj) const;
  void ReturnPooledEData(EData* obj) const;

  /*!
   * Write parameters to file
   */
  void WriteParameters(AZUREParams& params, const Config& configure) const;

 private:
  const Config &configure_;
  EData *data_;
  CNuc *compound_;
  ParameterLimitsManager *limitsManager_;
  double theErrorDef;
  
  // Object pools for memory reuse - store cloned objects directly
  mutable std::stack<std::unique_ptr<CNuc>> cnuc_pool_;
  mutable std::stack<std::unique_ptr<EData>> edata_pool_;
  mutable std::mutex pool_mutex_;
  mutable bool pools_initialized_;
};

#endif
