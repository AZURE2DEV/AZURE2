#ifndef AZURECALC_H
#define AZURECALC_H

#include "Minuit2/FCNBase.h"
#include "Minuit2/FCNGradientBase.h"
#include "Constants.h"
#include "AZUREParams.h"
#include "WorkingData.h"
#include <vector>
#include <memory>
#include <atomic>
#include <mutex>
#include <stack>

class Config;
class WorkingCopy;
class CNuc;
class ParameterLimitsManager;

/// A function class to perform the calculation of the chi-squared value

/*!
 * The AZURECalc function class calculates the cross section based on a
 * parameter set for all available data, and returns a chi-squared value.
 * This function class is what Minuit calls repeatedly during the fitting
 * process to perform the minimization.
 */

class AZURECalc : public ROOT::Minuit2::FCNGradientBase {
 public:
  /*!
   * The AZURECalc object is created with reference to an WorkingCopy and CNuc object.
   *. The runtime configurations are also passed through a Config structure.
   */
  /*!
   * `results` is the buffer of the evaluations on `data` itself (a calculation, the final evaluation after a fit),
   * owned by the caller, which then reads it (the output files).  A fit's evaluations run on pooled working copies,
   * each with a buffer of its own (WorkingData).
   */
  AZURECalc(WorkingCopy *data, azure2::Results &results, CNuc *compound, const Config &configure,
            ParameterLimitsManager *limitsManager = nullptr) :
    configure_(configure),
    results_(results) {
    data_ = data;
    compound_ = compound;
    limitsManager_ = limitsManager;
  };

  ~AZURECalc() {};
  /*!
   * See Minuit2 documentation for an explanation of this function.
   */
  virtual double Up() const override { return theErrorDef; };
  /*!
   * Overloaded operator to make the class instance callable as a function.
   * A Minuit parameter array is passed as the dependent variable.  The function
   * returns the total chi-squared value.
   */
  virtual double operator()(const vector_r &) const override;
  /*!
   * Count every evaluation of operator() in `counter` (owned by the caller, which outlives the counting; nullptr
   * stops it).  The fit reports MIGRAD's calls from it: Minuit2's own MnFcn::fNumCall is a plain int that its
   * OpenMP numerical gradient increments from several threads at once, so under threads it can come out a few low.
   */
  void CountCallsInto(std::atomic<long> *counter) { callCounter_ = counter; }

  /*!
   * Analytic gradient of the total chi-squared with respect to the (external)
   * Minuit parameters, for Minuit2's FCNGradientBase interface.  The energy /
   * reduced-width block is computed by the shared reverse-mode adjoint
   * (AccumulateEGammaGradient -- the same engine AZUREAPI uses); norms, energy
   * shifts and any analytically-unsupported configuration fall back to central
   * finite differences of the chi-squared.
   */
  std::vector<double> Gradient(const std::vector<double> &) const override;
  /// How GradientBy takes the gradient: Gradient() takes the defaults.
  struct GradientRoute {
    /// The adjoint gradient even when G2 is on (AZURE_GRAD_CHECK's comparison with the Jacobian path).
    bool adjointOnly = false;
    /// AZURE_GRAD_CHECK on the Jacobian path: that gradient, compared with the same finite differences as the adjoint
    /// one; the check's exit is then the caller's, after its own comparison.
    const std::vector<double> *jacobianPathGradient = nullptr;
  };
  std::vector<double> GradientBy(const std::vector<double> &p, const GradientRoute &route) const;
  /*!
   * Skip Minuit's start-up numerical check of the analytic gradient (trusted,
   * separately validated).
   */
  bool CheckGradient() const override { return false; }
  /*!
   * Diagonal second derivatives for Minuit2 (Gauss-Newton: 2 * sum_i (dr_i/dp)^2
   * from the analytic residual Jacobian, plus the norm and energy-shift penalty
   * curvatures).  Without it Minuit2 seeds MIGRAD with a full numerical
   * gradient + G2 pass and NegativeG2 line searches, which costs as much as the
   * numerical fit it is meant to replace.  Enabled by EnableG2() only after a
   * successful probe (an unsupported configuration keeps the numerical seed).
   */
  std::vector<double> G2(const std::vector<double> &) const override;
  bool HasG2() const override { return g2Enabled_; }
  /// Probe the analytic Jacobian at `p`; enables G2 for Minuit2 when it works.
  bool EnableG2(const std::vector<double> &p);
  /// AZURE_DUMP_NONFINITE diagnostic (see AZURECalc.cpp).
  void DumpNonFinite(const vector_r &p, double chi, CNuc *lc, WorkingCopy *ld, const azure2::Results &results) const;
  /// True when the limits manager declares at least one nuisance parameter.
  bool HasNuisanceParameters() const;
  /*!
   * Calculate every non-mapped point of every segment and of its leaves in
   * one OpenMP pass (a point's Calculate also fills its mapped copies), so the
   * per-segment loops only combine.  Returns false -- and does nothing -- when
   * already inside a parallel region (Minuit's parallel numerical gradient) or
   * when AZURE_SERIAL_CHI2 is set; callers then use the per-point path.
   */
  bool PrecalculatePoints(CNuc *lc, WorkingCopy *ld, azure2::Results &results, const Config &config,
                          bool thmOnly = false) const;
  /// chi^2 of `p` evaluated on the given working copies, writing `results` (optionally per segment).
  double Chi2On(CNuc *lc, WorkingCopy *ld, azure2::Results &results, const vector_r &p, std::vector<double> *segChis = nullptr,
                bool thmOnly = false) const;
  /// The same under another Config than the run's (the pool check's serial reference).
  double Chi2OnWith(const Config &config, CNuc *lc, WorkingCopy *ld, azure2::Results &results, const vector_r &p,
                    std::vector<double> *segChis, bool thmOnly) const;
  /*!
   * Side-effect-free chi-squared evaluation (no iteration counter / file output
   * / object pools), used by the finite-difference part of Gradient().
   */
  double Chi2Value(const vector_r &p, bool thmOnly = false) const;

  /*!
   * Write the intermediate output files -- parameters, calculated segments and
   * the transformed physical parameters -- for the parameter vector \p p.
   *
   * Called periodically during a fit so a long run leaves usable output behind
   * before it finishes, and so a fit that is stopped or crashes is not a total
   * loss. MIGRAD reaches this through operator(); the analytic-Jacobian
   * minimizers call it directly, since they evaluate through the
   * side-effect-free Chi2Value and would otherwise write nothing until the end.
   */
  void WriteIterationOutput(const vector_r &p) const;

  //! Chi-squared evaluations between intermediate writes (MIGRAD, NLopt, LM).
  static const int kOutputInterval = 100;

  /*! Solver iterations between intermediate writes for the analytic-Jacobian
   *  minimizers (LM, GSL trust-region).
   *
   *  Their iterations are whole Gauss-Newton steps, not single function
   *  evaluations, and they converge in tens of them -- a converged LM fit runs
   *  11 to 30 iterations, well under a hundred evaluations. Counting
   *  evaluations here would mean the snapshot never fired at all. */
  static const int kIterationOutputInterval = 10;

  /*!
   * Standardized residuals and analytic residual Jacobian at the full parameter
   * vector `full`.  `residuals` is N_res; `jac` is row-major N_res x nFree;
   * `packedToFull[a]` is the full-vector index of Jacobian column a (the a-th
   * non-fixed parameter).  Returns false if the analytic path is unsupported.
   */
  bool ResidualJacobian(const vector_r &full, vector_r &residuals,
                        vector_r &jac, std::vector<int> &packedToFull) const;

  /*!
   * The fixed flags of the parameters being minimized, in Minuit order: those of
   * the AZUREParams handed to the minimizer, i.e. the .azr's, then an external
   * param.par's ("fixed", and the release/fix of a zero value), then the
   * <parameterSettings>.  MIGRAD reads them off the parameters themselves; the
   * least-squares solvers (ResidualJacobian's columns), the gradient's
   * finite-difference skip and the band covariance take them from here.  Set by
   * RunLevenbergMarquardt / RunGSLNonlinear from their parameters and by the
   * caller before MIGRAD.  Unset (or of the wrong size) falls back to the flags
   * a fresh FillMnParams gives, the .azr's alone.
   */
  void SetFixedMask(const std::vector<bool> &fixed) const { fixedMask_ = fixed; }
  std::vector<bool> FixedMask() const;

  /*!
   * Standardized residuals only (no Jacobian) at `full`, in the same row order
   * as ResidualJacobian.  A plain forward pass, used for the trust-region trial
   * evaluations of the GSL solver.  Returns false if the evaluation fails.
   */
  bool ResidualsOnly(const vector_r &full, vector_r &residuals) const;

  /*!
   * Levenberg-Marquardt / Gauss-Newton minimization of the total chi-squared
   * using the analytic residual Jacobian.  Converges in far fewer iterations
   * than MIGRAD for least-squares problems.  Updates `params` in place with the
   * best-fit values and (J^T J)^{-1} parameter errors, respecting parameter
   * limits.  Returns the final chi-squared, or a negative value if the analytic
   * Jacobian is unsupported (caller should fall back to MIGRAD).
   */
  double RunLevenbergMarquardt(AZUREParams &params, int maxIter = 200,
                               struct BandCovariance *bandCovOut = nullptr) const;

  /*!
   * Nonlinear least-squares minimization via GSL's trust-region solver
   * (gsl_multifit_nlinear, geodesic-accelerated).  Uses the same analytic
   * residual Jacobian as RunLevenbergMarquardt but factorizes J directly (QR).
   * Gaussian priors enter as extra residual rows; parameter limits are enforced
   * by projection.  Updates `params` in place with the best-fit values and
   * (J^T J + priors)^{-1} errors.  Returns the final chi-squared, or a negative
   * value if the analytic Jacobian is unsupported (caller falls back to MIGRAD).
   */
  double RunGSLNonlinear(AZUREParams &params, int maxIter = 200,
                         struct BandCovariance *bandCovOut = nullptr) const;

  /*!
   * Returns a reference to the Config structure.
   */
  const Config &configure() const { return configure_; };
  /*!
   * Returns a pointer to the WorkingCopy object.
   */
  WorkingCopy *data() const { return data_; };
  /*!
   * Returns a pointer to the CNuc object.
   */
  CNuc *compound() const { return compound_; };

  /*!
   * See Minuit2 documentation for an explanation of this function.
   */
  void SetErrorDef(double def) override { theErrorDef = def; };

  /*!
   * Calculate nuisance parameter chi-squared contribution
   */
  double CalculateNuisanceChiSquared(const vector_r &p) const;
  /// The nuisance prior (nominal, sigma) of full index f, false for none or one the prior covariance replaces (D3).
  bool NuisancePriorOf(int f, double &nominal, double &sigma) const;
  /// Is full index f a parameter of the prior covariance (<config> priorCovariance=)?
  bool CoveredByPriorCovariance(int f) const;

  /*!
   * Add the analytic gradient of the nuisance-parameter penalty to `grad`,
   * mirroring CalculateNuisanceChiSquared exactly (so it is consistent with the
   * objective).  Used when the energy/gamma block is taken analytically (which
   * covers only the data chi-squared).
   */
  void AddNuisanceGradient(const vector_r &p, std::vector<double> &grad) const;

  /*!
   * The copies an evaluation of a fit runs on: an idle one, or a new clone of compound() and data() when every copy
   * made so far is taken (made when needed, so a fit holds one per evaluation running at the same time: one with the
   * analytic gradient, one per thread of Minuit2's parallel numerical gradient).  Returned after the evaluation.
   */
  CNuc *GetPooledCNuc() const;
  WorkingData *GetPooledData() const;
  void ReturnPooledCNuc(CNuc *obj) const;
  void ReturnPooledData(WorkingData *obj) const;

  /*!
   * Write parameters to file
   */
  void WriteParameters(AZUREParams &params, const Config &configure) const;

 private:
  std::atomic<long> *callCounter_ = nullptr;  ///< CountCallsInto: not owned
  /*!
   * Shared setup for the least-squares solvers.  Runs one ResidualJacobian at
   * `full` (returned in `r`/`J`/`p2f`) to fix the packed column order, and fills
   * the free-parameter start values `x`, projection limits `lo`/`hi`, and the
   * Gaussian priors `pen_nom` / `pen_inv2` (zero where a param has no prior).
   * Returns the number of free parameters, or -1 if the Jacobian is unsupported.
   */
  int PrepareFreeParams(const vector_r &full,
                        const ROOT::Minuit2::MnUserParameters &mp,
                        vector_r &r, vector_r &J, std::vector<int> &p2f,
                        vector_r &x, std::vector<double> &lo,
                        std::vector<double> &hi, std::vector<double> &pen_nom,
                        std::vector<double> &pen_inv2) const;

  /*!
   * Shared covariance / error finalization for the least-squares solvers.  At
   * the best-fit `full`, builds H = J^T J (+ priors), inverts it (ridge fallback
   * if singular), writes the errors into `mp`, reports poorly-constrained
   * directions, and, if `bandCovOut` is set, exports the R-matrix sub-block
   * covariance for the uncertainty band.
   */
  void FinalizeLeastSquaresCovariance(const vector_r &full,
                                      const std::vector<int> &p2f,
                                      const std::vector<double> &pen_inv2,
                                      ROOT::Minuit2::MnUserParameters &mp,
                                      struct BandCovariance *bandCovOut) const;

  const Config &configure_;
  WorkingCopy *data_;
  azure2::Results &results_;  ///< the caller's buffer for the evaluations on data_
  CNuc *compound_;
  ParameterLimitsManager *limitsManager_;
  double theErrorDef;
  mutable std::vector<bool> fixedMask_;  ///< see FixedMask()

  // The idle copies of GetPooledCNuc / GetPooledData (no value of an evaluation survives in them that the next one
  // does not overwrite: tests/upgrade/history_test.cpp, AZURE_TRACE_FCN=fresh).
  mutable std::stack<std::unique_ptr<CNuc>> cnuc_pool_;
  mutable std::stack<std::unique_ptr<WorkingData>> edata_pool_;
  mutable std::mutex pool_mutex_;
  bool g2Enabled_ = false;
  // Gauss-Newton cache: the residual Jacobian taken in Gradient() also gives G2;
  // Minuit2 asks for G2 at the same point right after the gradient.
  mutable std::vector<double> g2CacheP_, g2Cache_;
  mutable std::vector<double> lastGradP_, lastGrad_;  ///< parameters of the last Gradient() call (AZURE_DUMP_NONFINITE)
  mutable std::mutex lastGradMutex_;
  /// Gradient and Gauss-Newton G2 from one residual Jacobian; false if unsupported.
  bool GradientFromJacobian(const std::vector<double> &p, std::vector<double> &grad) const;

  // Working copies for WriteIterationOutput, cloned on first use and reused for
  // every later snapshot.  Cloning the whole data set afresh every 100
  // evaluations (several GB for a large target-effect model) and freeing it
  // again fragments the per-thread malloc arenas: the freed copy is not reused
  // by the next snapshot, so resident memory climbed by roughly one data-set
  // copy per snapshot until long fits were OOM-killed.
  mutable std::unique_ptr<CNuc> output_compound_;
  mutable std::unique_ptr<WorkingCopy> output_data_;
  mutable azure2::Results output_results_;
  mutable std::mutex output_mutex_;
};

#endif
