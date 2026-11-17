#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <functional>
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
 * docs/upgrade/inference.md section 2).
 *
 * The order is the one AZURECalc::operator() always used, so its chi-squared keeps its bits:
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
 */
/// Scale of Park's J > 0 wall (T6): a level at J = -kParkNormScale costs one unit of chi-squared, so the wall is steep
/// on the scale of J but smooth.
constexpr double kParkNormScale = 1.0e-3;

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

/// One term of the engine's objective: Name, Add, ModelOnly, AddForSampler, ExcludesForSampler (TermList.h).
using ObjectiveTerm = azure2::inference::Term<ObjectiveContext>;

class Objective {
 public:
  Objective();
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
  /// The terms of the model alone (T6 Park, T7 abundances, T9 model data), in the engine's order: what the API adds
  /// to its data chi^2 (AZUREAPI::FilledModelPenalties, the analytic-gradient path, the writer's totals).
  double ModelChi2(const ObjectiveContext &c) const;
  /*!
   * T9, the model data (ModelData.h), for the paths that need more than its value, at the parameters c.compound
   * holds: its chi-squared alone (the MCMC likelihood; 0 without model data), its residual rows in row order, its
   * gradient d chi2/d(E, gamma) into accum, and its rows with their Jacobian rows over nCols packed columns of pmap.
   * The last two return false where the gradient is not available here (the caller differentiates numerically).
   */
  double ModelDataChi2(const ObjectiveContext &c) const;
  void ModelDataResiduals(const ObjectiveContext &c, vector_r &residuals) const;
  bool ModelDataGradient(const ObjectiveContext &c, const ParamIndexMap &pmap, const vector_matrix_r *shiftDeriv,
                         GradAccum &accum) const;
  bool AppendModelDataRows(const ObjectiveContext &c, const ParamIndexMap &pmap, const vector_matrix_r *shiftDeriv,
                           int nCols, vector_r &residuals, vector_r &jac) const;
  const std::vector<std::unique_ptr<ObjectiveTerm> > &Terms() const { return terms_.Terms(); }

 private:
  azure2::inference::TermList<ObjectiveContext> terms_;
};

/*!
 * The least-squares rows of the model terms that are sums of squares but not single-parameter priors
 * (docs/upgrade/objective.md D2, D7), appended to `residuals` (and, with `jac` and `pmap`, their Jacobian rows over
 * nCols packed columns): Park's J > 0 wall, one row per R-matrix level, J/s while J < 0 and 0 otherwise (Park
 * formalism only; the shift functions and their derivatives must be current); the rows of the gamma_decay records'
 * priors and wall (GammaDecaySet::PenaltyRows); the abundance sum-to-one term, one row per sample with a width,
 * ((sum w) - c)/s.  The number of rows does not depend on the parameters.  The diagonal priors' rows are the caller's,
 * from ForEachDiagonalPrior (DiagonalPriors.h).
 */
void AppendModelRows(const CompoundSystem *compound, const LevelState &state, const WorkingCopy &data, const Config &config, const ParamIndexMap *pmap, int nCols,
                     vector_r &residuals, vector_r *jac);

/*!
 * The prior covariance (<config> priorCovariance=) replaces the diagonal prior of each of its parameters
 * (docs/upgrade/objective.md D3): marks those priors on the segments, target effects and abundances (nuisance priors
 * are skipped by index, AZURECalc::NuisancePriorOf).  Returns the names of the parameters whose own prior it
 * replaced, for a note.  Call once, after the covariance is set (WorkingCopy::SetPriorCovariance).
 */
std::vector<std::string> ReplaceDiagonalPriors(const CompoundSystem *compound, WorkingCopy *data, const PriorCovariance &pc);

#endif
