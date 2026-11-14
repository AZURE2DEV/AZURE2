#ifndef POSTERIORSUPPORT_H
#define POSTERIORSUPPORT_H

#include <vector>

#include "Constants.h"

class AZUREParams;
class CompoundSystem;
class LevelInputs;
class LevelState;
class Config;
class WorkingCopy;

/*!
 * The support of the posterior in the parameters (docs/upgrade/objective.md section 6, "The support"): the fit vectors
 * whose free components lie inside the limits the fit gives them, of every parameter kind.  The one rule of the
 * ensemble sampler (AZURECalcMCMC, both sampling bases: it tests the fit vector it evaluates the objective at) and of
 * any other log-probability of the objective.
 */
struct PosteriorSupport {
  /// The limits of one free parameter: its index in the full vector, lower and upper (+-infinity where it has none).
  struct Limit {
    int index;
    double lower, upper;
  };
  std::vector<Limit> limits;

  /// The fit's limits of the free parameters of `params` (filled by CompoundSystem and WorkingCopy::FillMnParams, after a
  /// parameter file if any): the kinds' own, then those of <parameterSettings>, which this applies to `params` as
  /// AZUREMain does before a fit.
  static PosteriorSupport OfFit(const Config &configure, const CompoundSystem *compound, const LevelState &state,
                                const LevelInputs &inputs, WorkingCopy *data, AZUREParams &params);
  /// The limits `params` holds for its free parameters, as they stand (OfFit after it has applied the settings; a
  /// caller that has applied them itself, the API's session).
  static PosteriorSupport OfParameters(const AZUREParams &params);
  /// Whether the full fit vector `fitVector` lies inside every limit (a vector too short for a limit is outside).
  bool Contains(const vector_r &fitVector) const;
};

#endif  // POSTERIORSUPPORT_H
