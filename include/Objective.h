#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "AZUREGrad.h"
#include "CNuc.h"
#include "Config.h"
#include "Constants.h"
#include "EData.h"
#include "ESegment.h"
#include "TargetEffect.h"

/*!
 * The fit objective as one ordered list of terms (docs/upgrade/objective.md, Phase 2).  Every path that needs the
 * chi-squared of a parameter vector -- the Minuit FCN (AZURECalc::operator()), the finite-difference and LM cost
 * (AZURECalc::Chi2On) -- sums this list, in this order, after the model values of every point are in place.
 *
 * The order is the one AZURECalc::operator() always used, so its chi-squared keeps its bits:
 *   T1-T4  data chi^2 per segment, with each segment's normalization prior and energy-shift priors
 *          (THM experiments and profiled norms as before);
 *   T2     the data covariance correction;
 *   T5     nuisance priors (<parameterSettings> useAsNuisance);
 *   T6     Park's J > 0 wall;
 *   T8     electron-screening, target-density, Doppler-temperature and sample-thickness priors (four terms);
 *   T7     abundance priors and sum-to-one;
 *   T9     measured scattering lengths;
 *   T10    the prior parameter covariance.
 */
struct ObjectiveContext {
  CNuc *compound;
  EData *data;
  const vector_r &p;
  const Config &config;
  /// Per-segment data chi^2 out (index = segment number, the THM experiment's at its last segment), or null.
  std::vector<double> *segmentChi2 = nullptr;
  /// Only the THM segments' data chi^2 (T1-T4 over THM segments; no other term).
  bool thmOnly = false;
  /// Store the segment chi^2 and the Park, abundance, scattering-length and prior-covariance totals into `data`
  /// (the writer reads them): the evaluation the fit reports, not a probe.
  bool storeTotals = false;
  /// T5, the nuisance chi^2 (AZURECalc::CalculateNuisanceChiSquared), or empty when there is no limits manager.
  std::function<double(const vector_r &)> nuisance;
};

class ObjectiveTerm {
 public:
  virtual ~ObjectiveTerm() = default;
  virtual std::string Name() const = 0;
  /// Adds this term to chi2 (some terms add in several steps, in the order the engine always used).
  virtual void Add(const ObjectiveContext &c, double &chi2) const = 0;
};

class Objective {
 public:
  Objective();
  /// The one objective of the engine (it holds no state).
  static const Objective &Engine();
  /// The chi-squared of c.p (model values already computed).
  double Chi2(const ObjectiveContext &c) const;
  const std::vector<std::unique_ptr<ObjectiveTerm> > &Terms() const { return terms_; }

 private:
  std::vector<std::unique_ptr<ObjectiveTerm> > terms_;
};


/*
 * The diagonal priors (x - nominal)^2 / sigma^2 of single parameters, as (full parameter index, nominal, sigma), for
 * the code that needs their derivatives or their least-squares rows rather than their value: d/dx = 2 (x - nominal) /
 * sigma^2, curvature 2 / sigma^2.  Exactly the terms the objective adds (T3, T4, T8 and T7's per-nucleus priors), so
 * the gradients, G2 and the least-squares rows stay consistent with chi^2.
 */

// The shift penalties as (full parameter index, nominal, sigma), for the code
// that needs a penalty's derivatives rather than its value: d/dx = 2 (x - nominal)
// / sigma^2, curvature 2 / sigma^2.  Exactly the terms ESegment::AddShiftPenalty
// adds, so the gradient, G2 and the least-squares priors stay consistent with
// chi^2.  The nominal is 0 for every shift term since 2026-10-25 (the file's value
// is only the start).
template <class F>
inline void ForEachShiftPenalty(EData *data, const ParamIndexMap &pmap, F f) {
  for (int s = 1; s <= data->NumSegments(); s++) {
    ESegment *seg = data->GetSegment(s);
    if (!seg) continue;
    int idx;
    if (seg->IsVaryEnergyShift() && seg->GetEnergyShiftError() != 0.0 && !seg->IsShiftPriorReplaced() &&
        (idx = pmap.EnergyShiftIndex(s)) >= 0)
      f(idx, seg->GetNominalEnergyShift(), seg->GetEnergyShiftError());  // 0 unless a prior_centre row
    if (seg->IsVaryEnergyShiftSqrt() && seg->GetEnergyShiftSqrtError() != 0.0 && !seg->IsShiftSqrtPriorReplaced() &&
        (idx = pmap.EnergyShiftSqrtIndex(s)) >= 0)
      f(idx, 0.0, seg->GetEnergyShiftSqrtError());
    if (seg->OwnsShiftFn() && seg->HasShiftFn())
      for (int j = 0; j < seg->NumShiftFnParams(); j++)
        if (seg->IsVaryShiftFn(j) && seg->GetShiftFnErrors()[j] != 0.0 && !seg->IsShiftFnPriorReplaced(j) &&
            (idx = pmap.EnergyShiftFnIndex(s, j)) >= 0)
          f(idx, 0.0, seg->GetShiftFnErrors()[j]);
  }
}

// The electron-screening priors ((U_e - nominal)/sigma)^2, U_e in eV: f(full index, nominal, sigma) for each
// screening parameter with sigma > 0 (EData::ScreeningPenalty is their sum).
template <class F>
inline void ForEachScreeningPrior(EData *data, const ParamIndexMap &pmap, F f) {
  for (int k = 1; k <= data->NumScreeningParams(); k++) {
    const TargetEffect *te = data->GetScreeningTargetEffect(k);
    if (te->GetScreeningSigma() <= 0.0 || !te->IsScreeningVaried() || te->IsScreeningPriorReplaced()) continue;
    const int idx = pmap.ScreeningIndex(k);
    if (idx >= 0) f(idx, te->GetScreeningNominal(), te->GetScreeningSigma());
  }
}

// The fitted-density priors ((N - N0)/sigma)^2, N in atoms/cm^2: f(full index, nominal, sigma) for each density
// parameter with sigma > 0 (EData::DensityPenalty is their sum).
template <class F>
inline void ForEachDensityPrior(EData *data, const ParamIndexMap &pmap, F f) {
  for (int k = 1; k <= data->NumDensityParams(); k++) {
    const TargetEffect *te = data->GetDensityTargetEffect(k);
    if (te->GetDensitySigma() <= 0.0 || !te->IsDensityVaried() || te->IsDensityPriorReplaced()) continue;
    const int idx = pmap.DensityIndex(k);
    if (idx >= 0) f(idx, te->GetDensityNominal(), te->GetDensitySigma());
  }
}

// The fitted-Doppler-temperature priors ((T - T0)/sigma)^2, T in K: f(full index, nominal, sigma) for each Doppler
// parameter with sigma > 0 (EData::DopplerPenalty is their sum).
template <class F>
inline void ForEachDopplerPrior(EData *data, const ParamIndexMap &pmap, F f) {
  for (int k = 1; k <= data->NumDopplerParams(); k++) {
    const TargetEffect *te = data->GetDopplerTargetEffect(k);
    if (te->GetDopplerSigma() <= 0.0 || !te->IsDopplerVaried() || te->IsDopplerPriorReplaced()) continue;
    const int idx = pmap.DopplerIndex(k);
    if (idx >= 0) f(idx, te->GetDopplerNominal(), te->GetDopplerSigma());
  }
}

// The fitted-sample-thickness priors ((n - n0)/sigma)^2, n in atoms/barn: f(full index, nominal, sigma) for each
// thickness parameter with sigma > 0 that owns its value (EData::ThicknessPenalty is their sum; a thickness-group
// follower's parameter is fixed and carries no prior).
template <class F>
inline void ForEachThicknessPrior(EData *data, const ParamIndexMap &pmap, F f) {
  for (int k = 1; k <= data->NumThicknessParams(); k++) {
    const ESegment *seg = data->GetSegment(data->GetThicknessSegment(k));
    if (seg->GetThicknessSigma() <= 0.0 || !seg->IsThicknessVaried() || seg->IsThicknessPriorReplaced()) continue;
    const int idx = pmap.ThicknessIndex(k);
    if (idx >= 0) f(idx, seg->GetNominalThickness(), seg->GetThicknessSigma());
  }
}

// Every diagonal prior: normalizations (T3), energy shifts (T4), screening, density, Doppler and thickness (T8),
// abundances (T7, per nucleus; the sum-to-one term couples them and is not diagonal).
template <class F>
inline void ForEachDiagonalPrior(EData *data, CNuc *compound, const ParamIndexMap &pmap, F f) {
  for (int s = 1; s <= data->NumSegments(); s++) {
    ESegment *seg = data->GetSegment(s);
    int idx;
    if (!seg || !seg->IsVaryNorm() || seg->IsNormPriorReplaced() || (idx = pmap.NormIndex(s)) < 0) continue;
    const double n0 = seg->GetNominalNorm(), nerr = n0 / 100.0 * seg->GetNormError();
    if (nerr != 0.0) f(idx, n0, nerr);
  }
  ForEachShiftPenalty(data, pmap, f);
  ForEachScreeningPrior(data, pmap, f);
  ForEachDensityPrior(data, pmap, f);
  ForEachDopplerPrior(data, pmap, f);
  ForEachThicknessPrior(data, pmap, f);
  if (compound->HasAbundanceBlock())
    for (int k = 1; k <= compound->NumAbundanceParams(); k++) {
      const int idx = pmap.AbundanceIndex(k);
      if (idx < 0) continue;
      int sample, nucleus;
      compound->AbundanceOf(k, sample, nucleus);
      const SampleAbundance &e = compound->GetSample(sample)->entries[nucleus - 1];
      if (e.sigma > 0.0 && !e.priorReplaced) f(idx, e.nominal, e.sigma);
    }
}

/*!
 * The prior covariance (<config> priorCovariance=) replaces the diagonal prior of each of its parameters
 * (docs/upgrade/objective.md D3): marks those priors on the segments, target effects and abundances (nuisance priors
 * are skipped by index, AZURECalc::NuisancePriorOf).  Returns the names of the parameters whose own prior it
 * replaced, for a note.  Call once, after the covariance is set (EData::SetPriorCovariance).
 */
class PriorCovariance;
std::vector<std::string> ReplaceDiagonalPriors(CNuc *compound, EData *data, const PriorCovariance &pc);

#endif
