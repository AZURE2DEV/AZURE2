// The apparatus fold of a point: the kernels of its data set's effect (azure2/experiment/Kernels.h) applied to the
// values of its sub-points.
#ifndef POINTFOLD_H
#define POINTFOLD_H

#include <utility>
#include <vector>

#include "azure2/experiment/Kernels.h"

class Config;
class Unit;
class WorkingCopy;
namespace azure2 {
struct Results;
}

/// The folded value of `point` from its sub-points' cross-section slots, written into its own slot.
void FoldPoint(Unit &point, const WorkingCopy &data, const Config &, azure2::Results &results);
/// The folded value of given sub-point values (0-based, one per sub-point), nothing written (the transform of a
/// transmission or yield segment, a Coulomb-reference leaf).  kAlreadyApplied: the caller broadened them.
double FoldValues(const Unit &point, const WorkingCopy &data, const Config &, const std::vector<double> &values,
                  azure2::experiment::Broadening = azure2::experiment::Broadening::kApply);
/// The exact weights c (0-based, one per sub-point) of the fold, value = sum_i c[i] values[i]: the kernels evaluated
/// on a linear form.  False if the fold is not a pure linear form for this point.
bool FoldWeights(const Unit &point, const WorkingCopy &data, const Config &, std::vector<double> &c,
                 azure2::experiment::Broadening = azure2::experiment::Broadening::kApply);
/*!
 * The fold weights shared between a gradient call's forward folds and its cotangents.  The weights of a point depend
 * on its nodes (energies, stopping), its kernels' parameters and its own energy, never on the sub-point values, so
 * within one gradient call -- the forward pass, then the cotangents, at one set of parameters, with no hand-over, no
 * energy shift and no mesh re-placement in between -- they are computed once: while a KeepFoldWeights is alive on
 * `results`, FoldKept computes a point's weights at its first fold, keeps them, and returns sum_i c[i] values[i]
 * (equal to FoldValues to rounding); KeptWeights hands the same weights to the gradient.  The scope is the whole
 * invalidation rule: it is opened after the hand-over of the call's parameters and closed before anything can move a
 * node, a kernel parameter (density, thickness, width, temperature) or a point (the energy-shift columns), and the
 * weights are dropped with it.  One vector of NumSubPoints doubles per folded point.  Not opened under
 * --fold-weights-by-differences.
 */
class KeepFoldWeights {
 public:
  KeepFoldWeights(azure2::Results &results, const WorkingCopy &data, const Config &);
  ~KeepFoldWeights();
  KeepFoldWeights(const KeepFoldWeights &) = delete;
  KeepFoldWeights &operator=(const KeepFoldWeights &) = delete;

 private:
  azure2::Results &results_;
};
/// FoldValues (broadening applied) of `point`, through the kept weights while a KeepFoldWeights is alive; the fold
/// itself otherwise, and whenever the product is not finite (a failed point is then found exactly as without them).
double FoldKept(const Unit &point, const WorkingCopy &data, const Config &, const std::vector<double> &values,
                azure2::Results &results);
/// The weights FoldKept keeps for `point`; null if none are kept.
const std::vector<double> *KeptWeights(const Unit &point, const azure2::Results &results);
/// d/dT of the fold for the sub-point values (0-based), T the effect's Doppler temperature (doppler_<n>, the analytic
/// gradient); 0 without a Doppler step.
double FoldTemperatureDerivative(const Unit &point, const WorkingCopy &data, const Config &,
                                 const std::vector<double> &values);
/// The Doppler weights at c.m. energy eCM on the point's sub-points (1-based), and dw/dT if asked: the linear map the
/// transform of a transmission or yield applies before its own transform.  Empty without a Doppler step.
void FoldBroadeningWeights(const Unit &point, const WorkingCopy &data, const Config &, double eCM,
                           std::vector<std::pair<int, double>> &w, std::vector<std::pair<int, double>> *dwdT);

#endif
