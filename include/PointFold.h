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
/// d/dT of the fold for the sub-point values (0-based), T the effect's Doppler temperature (doppler_<n>, the analytic
/// gradient); 0 without a Doppler step.
double FoldTemperatureDerivative(const Unit &point, const WorkingCopy &data, const Config &,
                                 const std::vector<double> &values);
/// The Doppler weights at c.m. energy eCM on the point's sub-points (1-based), and dw/dT if asked: the linear map the
/// transform of a transmission or yield applies before its own transform.  Empty without a Doppler step.
void FoldBroadeningWeights(const Unit &point, const WorkingCopy &data, const Config &, double eCM,
                           std::vector<std::pair<int, double>> &w, std::vector<std::pair<int, double>> *dwdT);

#endif
