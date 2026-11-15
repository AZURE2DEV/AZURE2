#ifndef AZURECALC_H
#define AZURECALC_H

#include "Minuit2/FCNBase.h"
#include "Minuit2/FCNGradientBase.h"
#include "Constants.h"
#include "AZUREParams.h"
#include "WorkingData.h"
#include "FitControl.h"
#include "LevelInputs.h"
#include "LevelState.h"
#include "StatePool.h"
#include <vector>
#include <memory>
#include <atomic>
#include <mutex>
#include <stack>

class Config;
class WorkingCopy;
class CompoundSystem;
class ParameterLimitsManager;

/// What a fitter's evaluations run on and leave behind, owned by whoever runs the fitter (the program's run, a
/// session's fit, a sweep of an adaptive grid, a test) and lent to it as its azure2::Results are: the evaluations
/// are const methods (Minuit2's interface) and write nothing that is the fitter's own.  One per fitter: the copies
/// are copies of that fitter's level state and working copy.
struct FitterStates {
  FitterStates();
  ~FitterStates();  // where WorkingCopy is complete (AZURECalc.cpp)

  /// The idle level states and working copies (each with the buffer it writes) of the evaluations that run at the
  /// same time.  Each is a copy of the caller's, made when no idle one is left; no value of an evaluation survives
  /// in a working copy that the next one does not overwrite (tests/upgrade/history_test.cpp, AZURE_TRACE_FCN=fresh).
  StatePool<LevelState> idle;
  StatePool<WorkingData> idleData;

  /// What the intermediate output files are written from (WriteIterationOutput, under outputMutex): copies made at
  /// the first write and used again by every later one.  Cloning the whole data set afresh every 100 evaluations
  /// (several GB for a large target-effect model) and freeing it again fragments the per-thread malloc arenas: the
  /// freed copy is not reused by the next snapshot, so resident memory climbed by roughly one data-set copy per
  /// snapshot until long fits were OOM-killed.
  std::unique_ptr<LevelState> output;
  std::unique_ptr<WorkingCopy> outputData;
  azure2::Results outputResults;
  std::mutex outputMutex;

  /// The last gradient: the point it was taken at and its value (AZURE_DUMP_NONFINITE reads them; both under
  /// `mutex`), and the Gauss-Newton G2 the same residual Jacobian gave, with its point (Minuit2 asks for G2 at the
  /// same point right after the gradient).
  struct LastGradient {
    std::vector<double> at, value;
    std::mutex mutex;
    std::vector<double> g2At, g2;
  };
  LastGradient gradient;
};

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
   * The AZURECalc object is created with reference to an WorkingCopy and CompoundSystem object.
   *. The runtime configurations are also passed through a Config structure.
   */
  /*!
   * `results` is the buffer of the evaluations on `data` itself (a calculation, the final evaluation after a fit),
   * owned by the caller, which then reads it (the output files); `state` is the caller's level state, which those
   * evaluations write as they write `results`, and `inputs` the caller's level inputs (read only).  A fit's
   * evaluations run on pooled working copies, each with a buffer of its own (WorkingData), and on pooled copies of
   * `state`, all held in the caller's `states`.  The compound is shared and constant; the caller outlives this
   * object.
   */
  AZURECalc(WorkingCopy *data, azure2::Results &results, const CompoundSystem *compound, LevelState &state,
            const LevelInputs &inputs, FitterStates &states, const Config &configure,
            ParameterLimitsManager *limitsManager = nullptr) :
    configure_(configure),
    results_(results),
    state_(state),
    inputs_(inputs),
    states_(states) {
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
   * Hand every evaluation to `observer` (owned by the caller, which outlives the observing; nullptr stops it) and
   * obey its stop request during a fit (FitControl.h).  The program's observer prints the progress lines and writes
   * the intermediate files (src/ProgramFit.cpp); a fit puts its own in place while it runs (Fitters.h).
   */
  void Observe(const EvaluationObserver *observer) { observer_ = observer; }
  const EvaluationObserver *Observer() const { return observer_; }
  /// The observer's stop request is set.
  bool StopRequested() const {
    return observer_ && observer_->stop && observer_->stop->load(std::memory_order_relaxed);
  }

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
  void DumpNonFinite(const vector_r &p, double chi, const CompoundSystem *lc, WorkingCopy *ld, const azure2::Results &results) const;
  /// True when the limits manager declares at least one nuisance parameter.
  bool HasNuisanceParameters() const;
  /*!
   * Calculate every non-mapped point of every segment and of its leaves in
   * one OpenMP pass (a point's Calculate also fills its mapped copies), so the
   * per-segment loops only combine.  Returns false -- and does nothing -- when
   * already inside a parallel region (Minuit's parallel numerical gradient) or
   * when AZURE_SERIAL_CHI2 is set; callers then use the per-point path.
   */
  bool PrecalculatePoints(const LevelState &state, WorkingCopy *ld, azure2::Results &results, const Config &config,
                          bool thmOnly = false) const;
  /// chi^2 of `p` evaluated on the given working copies, writing `results` (optionally per segment).
  double Chi2On(LevelState &state, WorkingCopy *ld, azure2::Results &results, const vector_r &p, std::vector<double> *segChis = nullptr,
                bool thmOnly = false) const;
  /// The same under another Config than the run's (the pool check's serial reference).
  double Chi2OnWith(const Config &config, LevelState &state, WorkingCopy *ld, azure2::Results &results, const vector_r &p,
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
   * loss. It is the program's FitControl::checkpoint (src/ProgramFit.cpp):
   * its observer calls it every kOutputInterval evaluations (MIGRAD, NLopt), and
   * LM and GSL-LM call the checkpoint every kIterationOutputInterval iterations,
   * since they evaluate through the side-effect-free Chi2Value.  A session's fit
   * has no checkpoint and writes nothing.
   */
  void WriteIterationOutput(const vector_r &p) const;

  //! Evaluations of the run between intermediate writes (MIGRAD, NLopt, the evaluations outside a fit).
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
  void SetFixedMask(const std::vector<bool> &fixed) { fixedMask_ = fixed; }
  std::vector<bool> FixedMask() const;

  /*!
   * Standardized residuals only (no Jacobian) at `full`, in the same row order
   * as ResidualJacobian.  A plain forward pass, used for the trust-region trial
   * evaluations of the GSL solver.  Returns false if the evaluation fails.
   */
  bool ResidualsOnly(const vector_r &full, vector_r &residuals) const;

  /*!
   * Returns a reference to the Config structure.
   */
  const Config &configure() const { return configure_; };
  /*!
   * Returns a pointer to the WorkingCopy object.
   */
  WorkingCopy *data() const { return data_; };
  /*!
   * Returns a pointer to the CompoundSystem object.
   */
  const CompoundSystem *compound() const { return compound_; };
  /// The caller's level inputs (the labels of the parameters read them).
  const LevelInputs &inputs() const { return inputs_; }

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
   * The copies an evaluation of a fit runs on: an idle one, or a new copy of the caller's level state and of data()
   * when every copy made so far is taken (made when needed, so a fit holds one per evaluation running at the same
   * time: one with the analytic gradient, one per thread of Minuit2's parallel numerical gradient).  Returned after
   * the evaluation.
   */
  std::unique_ptr<LevelState> GetPooledState() const;
  WorkingData *GetPooledData() const;
  void ReturnPooledState(std::unique_ptr<LevelState> state) const;
  void ReturnPooledData(WorkingData *obj) const;

  /*!
   * Write parameters to file
   */
  void WriteParameters(AZUREParams &params, const Config &configure) const;

 private:
  std::atomic<long> *callCounter_ = nullptr;  ///< CountCallsInto: not owned
  const EvaluationObserver *observer_ = nullptr;  ///< Observe: not owned
  const Config &configure_;
  WorkingCopy *data_;
  azure2::Results &results_;  ///< the caller's buffer for the evaluations on data_
  LevelState &state_;          ///< the caller's level state: those evaluations' state, the seed of the pooled ones
  const LevelInputs &inputs_;  ///< the caller's level inputs
  FitterStates &states_;       ///< the caller's: the pooled copies, the output copies, the last gradient
  const CompoundSystem *compound_;
  ParameterLimitsManager *limitsManager_;
  double theErrorDef;
  std::vector<bool> fixedMask_;  ///< see FixedMask()
  bool g2Enabled_ = false;
  /// Gradient and Gauss-Newton G2 from one residual Jacobian (the G2 into the caller's FitterStates::gradient);
  /// false if unsupported.
  bool GradientFromJacobian(const std::vector<double> &p, std::vector<double> &grad) const;
};

#endif
