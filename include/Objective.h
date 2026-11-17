#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <functional>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "AZUREGrad.h"
#include "DiagonalPriors.h"
#include "CompoundSystem.h"
#include "Config.h"
#include "Constants.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "TargetEffect.h"
#include "azure2/inference/TermList.h"

/*!
 * The fit objective as one ordered list of terms (docs/upgrade/objective.md, Phase 2).  Every path that needs the
 * chi-squared of a parameter vector -- the Minuit FCN (AZURECalc::operator()), the finite-difference and LM cost
 * (AZURECalc::Chi2On) -- sums this list, in this order, after the model values of every point are in place.
 *
 * The rule of the list -- the sum in order, the sum with a breakdown, the sampler's sum and support, the terms of the
 * model alone -- is the inference/ library's (azure2/inference/TermList.h: Term<Context>, TermList<Context>), over
 * the context below.  The terms themselves read the engine's objects and are engine code (src/Objective.cpp;
 * docs/upgrade/inference.md section 2).  A term states in its class its value, its gradient, its least-squares rows,
 * its total and its field of chiSquared.out (ObjectiveTerm below; docs/upgrade/objective_arrays.md section 8); the
 * walks over them (Objective::AddGradient, AppendModelRows, WriteTotals) name no term.
 *
 * The order of the value is the one AZURECalc::operator() always used, so its chi-squared keeps its bits:
 *   T1-T4  data chi^2 per segment, with each segment's normalization prior and energy-shift priors
 *          (THM experiments and profiled norms as before);
 *   T2     the data covariance correction;
 *   T5     nuisance priors (<parameterSettings> useAsNuisance);
 *   T6     Park's J > 0 wall;
 *   T8     electron-screening, target-density, Doppler-temperature, sample-thickness and apparatus-parameter priors
 *          (five terms; each the subtotal of its priors as AppendTargetPriors, DiagonalPriors.h, gathers them);
 *   T7     abundance priors and sum-to-one;
 *   T9     the model data (ModelData.h): scattering lengths and the level quantities of <modelData>;
 *   T10    the prior parameter covariance;
 *   T11    the priors of the gamma_decay records and the wall on a negative remainder;
 *   T12    the priors of the isospin shifts and mixes and the soft ties of the isospin levels.
 * The gradient, the rows and the writer have orders of their own (stated in the Objective constructor).
 */
/// Scale of Park's J > 0 wall (T6): a level at J = -kParkNormScale costs one unit of chi-squared, so the wall is steep
/// on the scale of J but smooth.
constexpr double kParkNormScale = 1.0e-3;

/*!
 * T6's gather (inference/'s ParkWall.h holds the arithmetic): one ParkLevel per level in the R matrix (of a J group in
 * the R matrix), in group and level order, with its Park norm J from `state`.  With `channels`, a level on the wall
 * (OnParkWall) also carries its particle channels ('P' channels of a particle pair) with gamma_c, S'_c and S''_c from
 * `state` (the gradient and the rows; the value reads J only, and S'' is not filled for it).
 */
std::vector<ParkLevel> ParkLevelsOf(const CompoundSystem *compound, const LevelState &state, bool channels);

struct ObjectiveContext {
  const CompoundSystem *compound;
  /// The levels at the evaluated vector (the model data's finite differences write and restore them) and the
  /// levels' input values of the evaluation's owner.
  LevelState &state;
  const LevelInputs &inputs;
  WorkingCopy *data;
  /// The buffer of the evaluation: the model value of every point of `data` (azure2::EvaluateModel) in, what the terms
  /// store out (Results::objective: profiled norms, THM profiles, segment chi^2, totals; the writer reads them).
  azure2::Results &results;
  const vector_r &p;
  const Config &config;
  /// Per-segment data chi^2 out (index = segment number, the THM experiment's at its last segment), or null.
  std::vector<double> *segmentChi2 = nullptr;
  /// Only the THM segments' data chi^2 (T1-T4 over THM segments; no other term).
  bool thmOnly = false;
  /// Add the normalization and energy-shift priors (T3, T4) in the segment loop.  The sampler takes them as
  /// normalised log-priors (AZURECalcMCMC::CalculateLogPrior) and evaluates the rest (Objective::SamplerChi2).
  bool segmentPriors = true;
  /// Optional breakdown (index = segment number): the normalization prior and the three energy-shift prior parts
  /// of each segment, as added by the segment loop (pyazr's penalties(), the writer).
  struct SegmentPriors {
    std::vector<double> norm, shift, shiftSqrt, shiftFn;
  };
  SegmentPriors *segmentPriorsOut = nullptr;
  /// Optional breakdown: each term's contribution, in Objective::Terms() order.
  std::vector<double> *termsOut = nullptr;
  /// T5, the nuisance chi^2 (AZURECalc::CalculateNuisanceChiSquared), or empty when there is no limits manager.
  std::function<double(const vector_r &)> nuisance;
};

/// What a term's AddGradient did: nothing (no gradient of its own here), added its part, or could not (the caller
/// then takes finite differences of the whole objective).
enum class TermGradient { kNone, kAdded, kUnavailable };

/// Which terms a walk visits: all, those of the model alone (ModelOnly: pyazr's calculate_chi2), or the others.
enum class TermSet { kAll, kModelOnly, kOthers };

/*!
 * Where a term's least-squares rows go (ObjectiveTerm::AppendRows): each row is a residual r and, with a Jacobian,
 * its row over nCols packed columns of pmap: a GradAccum of dr/d(E, gamma, abundance, ...) scattered through pmap,
 * plus (full index, dr/dx) pairs added as they are.  Without pmap or jac only the residuals are appended.
 */
struct ModelRows {
  const ParamIndexMap *pmap;
  int nCols;
  vector_r &residuals;
  vector_r *jac;
  void Append(double r, const GradAccum *acc, const std::vector<std::pair<int, double>> *direct);
};

/*!
 * One term of the engine's objective, in one place (docs/upgrade/objective_arrays.md section 8).  From
 * inference/'s Term: Name, Add (the value, added to the running sum; a term with a total stores it with StoreTotal),
 * ModelOnly, AddForSampler, ExcludesForSampler (TermList.h).  Here, each with a default that does nothing:
 *   AddGradient  d chi^2/d theta into accum (Objective::AddGradient walks it);
 *   AppendRows   its least-squares rows (Objective::AppendModelRows);
 *   WriteTotal   its field(s) of chiSquared.out (Objective::WriteTotals), TotalLabel the label of the field that
 *                holds the whole term ("" when none does).
 */
class ObjectiveTerm : public azure2::inference::Term<ObjectiveContext> {
 public:
  /*!
   * The term's gradient at c (the levels of c.state, the working copy's values) into accum: level energies and
   * widths, abundances, screening, Doppler, thickness, and any full index through accum.direct.  shiftDeriv: the
   * Brune shift-function derivative table, or null.
   */
  virtual TermGradient AddGradient(const ObjectiveContext &, const ParamIndexMap &, const vector_matrix_r *,
                                   GradAccum &) const {
    return TermGradient::kNone;
  }
  /// Its rows (sums of squares only); the number of rows must not depend on the parameters.
  virtual void AppendRows(const ObjectiveContext &, ModelRows &) const {}
  /// Its field(s) of chiSquared.out, each " <label>: <value>", when the term's condition holds.
  virtual void WriteTotal(const WorkingCopy &, const Config &, const azure2::Results &, std::ostream &) const {}
  virtual std::string TotalLabel() const { return ""; }
  /// The term's place in Objective::Terms(), its slot in Results::objective.totals.
  int Slot() const { return slot_; }

 protected:
  void StoreTotal(const ObjectiveContext &c, double value) const { c.results.objective.StoreTotal(slot_, value); }
  /// The usual field: " <TotalLabel()>: <the stored total>".
  void WriteStoredTotal(const azure2::Results &results, std::ostream &out) const;

 private:
  friend class Objective;
  int slot_ = -1;  // set once by the Objective's constructor
};

class Objective {
 public:
  Objective();
  Objective(const Objective &) = delete;
  Objective &operator=(const Objective &) = delete;
  /// The one objective of the engine (it holds no state).
  static const Objective &Engine();
  /// The chi-squared of c.p (model values already computed).
  double Chi2(const ObjectiveContext &c) const;
  /// The data part alone, T1 (with the THM profiling) and T2, without the segment priors: the API's data
  /// chi-squared (AZUREAPI::EvaluateFilledChi2, pyazr's calculate_chi2_rwa).
  double DataChi2(ObjectiveContext c) const;
  /// The ensemble sampler's -2 ln L: every term's AddForSampler, in the table's order (segment priors off).
  double SamplerChi2(ObjectiveContext c) const;
  /// Whether a term puts c's point outside the sampler's support (ExcludesForSampler).
  bool ExcludesForSampler(const ObjectiveContext &c) const;
  /// The terms of the model alone (T6 Park, T7 abundances, T9 model data, T11), in the engine's order: what the API
  /// adds to its data chi^2 (AZUREAPI::FilledModelPenalties, the analytic-gradient path, the writer's totals).
  double ModelChi2(const ObjectiveContext &c) const;
  /*!
   * T9, the model data (ModelData.h): its chi-squared alone (the MCMC likelihood; 0 without model data), and its
   * residual rows in row order, without and with their Jacobian rows over nCols packed columns of pmap (false where
   * the Jacobian is not available).  These rows sit with the data rows, before the prior covariance's.
   */
  double ModelDataChi2(const ObjectiveContext &c) const;
  void ModelDataResiduals(const ObjectiveContext &c, vector_r &residuals) const;
  bool AppendModelDataRows(const ObjectiveContext &c, const ParamIndexMap &pmap, const vector_matrix_r *shiftDeriv,
                           int nCols, vector_r &residuals, vector_r &jac) const;
  /*!
   * The gradient walk: every term of `which`'s AddGradient into accum, in the gradient's order.  Returns false as soon
   * as a term's gradient is unavailable (the caller differentiates numerically); `added`, if given, is set to whether
   * any term added something.
   */
  bool AddGradient(const ObjectiveContext &c, const ParamIndexMap &pmap, const vector_matrix_r *shiftDeriv,
                   GradAccum &accum, TermSet which = TermSet::kAll, bool *added = nullptr) const;
  /*!
   * The rows walk: every term's AppendRows, in the rows' order, appended to `residuals` (and, with `jac` and `pmap`,
   * their Jacobian rows over nCols packed columns).  Today Park's J > 0 wall (one row per R-matrix level, Park
   * formalism only; the shift functions and their derivatives must be current), the gamma_decay records' priors and
   * wall, the abundance sums (objective.md D2, D7).  The diagonal priors' rows are the caller's (DiagonalPriors.h).
   */
  void AppendModelRows(const ObjectiveContext &c, const ParamIndexMap *pmap, int nCols, vector_r &residuals,
                       vector_r *jac) const;
  /// The writer's walk: every term's WriteTotal, in the order of the fields of chiSquared.out.
  void WriteTotals(const WorkingCopy &data, const Config &config, const azure2::Results &results,
                   std::ostream &out) const;
  /// T10's total stored for chiSquared.out by a path whose sum leaves T10 out (the API's calculate_chi2).
  void RecordPriorCovariance(const ObjectiveContext &c) const;
  const std::vector<std::unique_ptr<ObjectiveTerm> > &Terms() const { return terms_.Terms(); }

 private:
  /// Appends a term to the value's list and gives it its slot.
  template <class T, class... A>
  const ObjectiveTerm *Append(A &&...args);
  azure2::inference::TermList<ObjectiveContext, ObjectiveTerm> terms_;
  /// The orders of the walks (the constructor: today's order of the terms that take part, then the rest).
  std::vector<const ObjectiveTerm *> gradientWalk_, rowWalk_, reportWalk_;
  const ObjectiveTerm *priorCovariance_ = nullptr;
};

/*!
 * The prior covariance (<config> priorCovariance=) replaces the diagonal prior of each of its parameters
 * (docs/upgrade/objective.md D3): marks those priors on the segments, target effects and abundances (nuisance priors
 * are skipped by index, AZURECalc::NuisancePriorOf).  Returns the names of the parameters whose own prior it
 * replaced, for a note.  Call once, after the covariance is set (WorkingCopy::SetPriorCovariance).
 */
std::vector<std::string> ReplaceDiagonalPriors(const CompoundSystem *compound, WorkingCopy *data, const PriorCovariance &pc);

#endif
