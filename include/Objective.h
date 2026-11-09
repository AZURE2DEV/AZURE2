#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "Constants.h"

class CNuc;
class Config;
class EData;

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

#endif
