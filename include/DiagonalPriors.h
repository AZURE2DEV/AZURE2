#ifndef DIAGONALPRIORS_H
#define DIAGONALPRIORS_H

#include "InferenceTypes.h"

class WorkingCopy;

/*!
 * The engine's gathers of the objective's diagonal priors (docs/upgrade/objective_arrays.md section 2.1): walks of the
 * working copy that hold the rules of eligibility and fill the list of inference/'s DiagonalPriors.h, whose arithmetic
 * the terms, the gradients and the chiSquared.out writer call.  Nothing here forms a chi^2.
 */

/*!
 * T8: the priors of one kind of target or apparatus parameter (`source` one of kScreening, kDensity, kDoppler,
 * kThickness, kApparatus), in the order of its parameters k = 1..N, appended to `out`.  A parameter has a prior when
 * it is varied, its prior is not replaced by the prior covariance and its width is positive (screening, density and
 * Doppler: also the effect's own flag, and a width that is not <= 0).  The value is the working copy's, as the
 * hand-over left it (clamped into the layout range for density, Doppler and apparatus; at least 0 for a thickness);
 * the index is the parameter's in `pmap`, or -1 without one.
 */
void AppendTargetPriors(const WorkingCopy &data, PriorSource source, const ParamIndexMap *pmap, DiagonalPriors &out);

#endif  // DIAGONALPRIORS_H
