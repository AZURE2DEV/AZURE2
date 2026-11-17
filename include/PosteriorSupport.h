#ifndef POSTERIORSUPPORT_H
#define POSTERIORSUPPORT_H

#include "InferenceTypes.h"

class AZUREParams;
class CompoundSystem;
class LevelInputs;
class LevelState;
class Config;
class WorkingCopy;

/*!
 * The makers of the posterior's support (azure2/inference/PosteriorSupport.h: the limits and the rule Contains) from
 * the fit's parameters, which are Minuit2's and the engine's (docs/upgrade/objective.md section 6, "The support").
 */

/// The fit's limits of the free parameters of `params` (filled by CompoundSystem and WorkingCopy::FillMnParams, after a
/// parameter file if any): the kinds' own, then those of <parameterSettings>, which this applies to `params` as
/// AZUREMain does before a fit.
PosteriorSupport PosteriorSupportOfFit(const Config &configure, const CompoundSystem *compound, const LevelState &state,
                                       const LevelInputs &inputs, WorkingCopy *data, AZUREParams &params);
/// The limits `params` holds for its free parameters, as they stand (PosteriorSupportOfFit after it has applied the
/// settings; a caller that has applied them itself, the API's session).
PosteriorSupport PosteriorSupportOfParameters(const AZUREParams &params);

#endif  // POSTERIORSUPPORT_H
