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

/*!
 * T3: the normalization prior of segment `s` (1-based), when it has one: the norm varied and not profiled (a THM
 * segment's free norm is profiled, and a segment of a THM experiment pays none), not replaced by the prior covariance,
 * and sigma_n = n0 sigma% / 100 != 0 (a norm-group follower's sigma is 0).  Value: the working copy's norm.
 */
void AppendNormPrior(const WorkingCopy &data, int s, const ParamIndexMap *pmap, DiagonalPriors &out);

/*!
 * T4: the energy-shift priors of segment `s`, in the order the objective adds them: the constant shift (MeV), the
 * sqrt(E) coefficient (MeV^1/2), then each parameter of the shift function the segment owns (a function shared by a
 * group is penalised once, on its owner); each when varied, its sigma != 0 and not replaced.  Values: the working
 * copy's.
 */
void AppendShiftPriors(const WorkingCopy &data, int s, const ParamIndexMap *pmap, DiagonalPriors &out);

/// T3 and T4 of segment `s`: the norm, then the shifts (the order SegmentTerms adds them in).
void AppendSegmentPriors(const WorkingCopy &data, int s, const ParamIndexMap *pmap, DiagonalPriors &out);

/// T7's per-nucleus priors, in the order of the abundance parameters: sigma > 0 and not replaced.  (The sum-to-one
/// term couples the abundances of a sample and is not diagonal.)
void AppendAbundancePriors(const WorkingCopy &data, const ParamIndexMap *pmap, DiagonalPriors &out);

/*!
 * Every diagonal prior of the objective with its index in `pmap`, in the order the least-squares rows of the API
 * expect (AZUREAPIObjective.cpp): the norms of all segments, the shifts of all segments, screening, density, Doppler,
 * thickness, apparatus, abundances.  An entry whose parameter `pmap` does not hold (index -1) is left out.
 */
DiagonalPriors DiagonalPriorsOf(const WorkingCopy &data, const ParamIndexMap &pmap);

/// f(full index, centre, sigma) for each entry of DiagonalPriorsOf, in its order: the derivatives, the curvature and
/// the least-squares rows of the diagonal priors (each path keeps its own form of the expression).
template <class F>
inline void ForEachDiagonalPrior(const WorkingCopy *data, const ParamIndexMap &pmap, F f) {
  for (const DiagonalPrior &e : DiagonalPriorsOf(*data, pmap)) f(e.index, e.centre, e.sigma);
}

#endif  // DIAGONALPRIORS_H
