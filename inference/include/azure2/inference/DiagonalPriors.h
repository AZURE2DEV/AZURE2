#ifndef AZURE2_INFERENCE_DIAGONALPRIORS_H
#define AZURE2_INFERENCE_DIAGONALPRIORS_H

#include <vector>

#include "azure2/inference/GaussianPrior.h"

namespace azure2 {
namespace inference {

/*!
 * The diagonal priors of the objective as one list (docs/upgrade/objective_arrays.md section 2.1): each entry one
 * Gaussian prior ((x - x0)/sigma)^2 on one parameter (GaussianPrior.h), every number in the unit of the parameter's
 * kind (ParamKindTable(): a normalization dimensionless, an energy shift in MeV, U_e in eV, an areal density in
 * atoms/cm^2, a temperature in K, a thickness in atoms/b, an apparatus parameter in its kernel's engine unit).
 *
 * Which parameter carries a prior, and the value it is formed on, is the engine's gather (include/DiagonalPriors.h):
 * this header holds no rule of eligibility, only the arithmetic over the list.  The value of the objective is formed
 * on `value` (the working copy's, clamped into the parameter's layout range where the hand-over clamps); a gradient
 * and a curvature on the component p[index] of the vector (objective_arrays.md R3).  The sums run in list order, so
 * the bits of a total are those of the loop the list replaced (R1).
 */
enum class PriorSource {
  kNorm,       ///< T3, a segment's normalization (group: the segment, 1-based)
  kShift,      ///< T4, a segment's constant energy shift (group: the segment)
  kShiftSqrt,  ///< T4, a segment's sqrt(E) coefficient (group: the segment)
  kShiftFn,    ///< T4, a parameter of the shift function a segment owns (group: the segment)
  kScreening,  ///< T8, an electron-screening potential
  kDensity,    ///< T8, a fitted areal density
  kDoppler,    ///< T8, a fitted Doppler temperature
  kThickness,  ///< T8, a fitted sample thickness
  kApparatus,  ///< T8, an apparatus parameter (group: its row of ApparatusRows())
  kAbundance,  ///< T7, the prior of one nucleus in one sample
  kNuisance    ///< T5, a <parameterSettings> nuisance prior
};

/// Any group, in PriorsChi2.
constexpr int kAnyGroup = -1;

struct DiagonalPrior {
  int index = -1;       ///< the full-vector index (ParamIndexMap), -1 where the gather had no map
  double value = 0.0;   ///< the value the objective's chi^2 reads
  double centre = 0.0;  ///< x0
  double sigma = 0.0;   ///< the width; the gather enters a prior only where the engine's rule gives it one
  PriorSource source = PriorSource::kNorm;
  int group = 0;        ///< the segment (T3, T4), the apparatus row (T8 apparatus), else 0
};
using DiagonalPriors = std::vector<DiagonalPrior>;

/// chi2 += ((value - x0)/sigma)^2 for each entry, in list order (a term that adds each prior to its running sum).
inline void AddPriorsChi2(const DiagonalPriors &list, double &chi2) {
  for (const DiagonalPrior &e : list) chi2 += GaussianPriorChi2(e.value, e.centre, e.sigma);
}

/// The sum from 0, in list order, of the entries of `source` (and of `group`, unless kAnyGroup).
inline double PriorsChi2(const DiagonalPriors &list, PriorSource source, int group = kAnyGroup) {
  double chi2 = 0.0;
  for (const DiagonalPrior &e : list)
    if (e.source == source && (group == kAnyGroup || e.group == group)) chi2 += GaussianPriorChi2(e.value, e.centre, e.sigma);
  return chi2;
}

/// The sum from 0, in list order, of every entry.
inline double PriorsChi2(const DiagonalPriors &list) {
  double chi2 = 0.0;
  AddPriorsChi2(list, chi2);
  return chi2;
}

/// grad[index] += 2 (p[index] - x0)/sigma^2 for each entry whose index lies in both p and grad.
inline void AddPriorsGradient(const DiagonalPriors &list, const std::vector<double> &p, std::vector<double> &grad) {
  for (const DiagonalPrior &e : list)
    if (e.index >= 0 && e.index < (int)p.size() && e.index < (int)grad.size())
      grad[e.index] += GaussianPriorGradient(p[e.index], e.centre, e.sigma);
}

/// g2[index] += 2/sigma^2 for each entry whose index lies in g2.
inline void AddPriorsCurvature(const DiagonalPriors &list, std::vector<double> &g2) {
  for (const DiagonalPrior &e : list)
    if (e.index >= 0 && e.index < (int)g2.size()) g2[e.index] += GaussianPriorCurvature(e.sigma);
}

}  // namespace inference
}  // namespace azure2

#endif  // AZURE2_INFERENCE_DIAGONALPRIORS_H
