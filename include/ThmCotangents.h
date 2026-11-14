#ifndef AZURE2_THMCOTANGENTS_H
#define AZURE2_THMCOTANGENTS_H

#include <memory>
#include <vector>

#include "GradEnergyTasks.h"

class CNuc;
class LevelState;
class Config;
struct ThmAdjointContext;
class WorkingCopy;
namespace azure2 {
struct Results;
}

/*!
 * THM data in the scalar gradient pass (GradEnergyTasks.h; the design is the comment of ComputeTHMRows, AZUREGrad.h).
 *
 * The chi^2 of THM data compares s m_i + b(E_i) with d_i, the scale s (and an experiment's background coefficients) a
 * closed-form least-squares profile of the models (Segment::ProfileNormChiSquared, ThmProfile.h), or s = 1/n for a
 * norm that is fixed or fitted.  A profiled linear parameter is at the minimum of chi^2 over itself, so d chi^2/d s = 0
 * there and the total derivative of chi^2 in a level parameter is its partial derivative at fixed s and b (the
 * envelope theorem): the cotangent of point i's model is
 *
 *   d chi^2/d m_i = 2 r_i s f_i / e_i,   r_i = (s f_i m_i + b(E_i) - d_i)/e_i,
 *
 * f_i the factor that takes the model to the frame of the datum (AppendInDataFrame).  From there a THM point is a
 * datum like any other: its fold's weights take the cotangent to the sites it reads, the cotangents of every reader
 * of a kHoes site are summed, and the HOES adjoint (THMMatrixFunc::THMPointAdjoint) is taken once per site, in the
 * site's energy task (ScalarTask).
 */
namespace azure2 {

/// Whether the scalar pass takes the THM data (ThmModelCotangents, the HOES adjoint per site): every THM data set of
/// `data` has no components, the experiments carry no coherent-background parameter, the formalism is not RMC, and
/// AZURE2_THM_FD is not set.  Otherwise the THM data stay outside it (ComputeTHMRows, central differences).
bool ThmThroughEnergyTasks(const WorkingCopy &data, const Config &config);

/// The context of the HOES adjoint for the scalar pass (BuildThmAdjointContext at the compound as filled, shiftDeriv
/// its Brune table or null) where ThmThroughEnergyTasks holds and `data` has THM data; null otherwise.
std::unique_ptr<ThmAdjointContext> ThmGradientContext(CNuc *compound, const LevelState &state, const WorkingCopy &data, const Config &config,
                                                      const vector_matrix_r *shiftDeriv);

/// The cotangents d chi^2/d m (above) of the points of `points` that belong to THM data sets, into fitBar (the same
/// index), from their models `models` (model frame); the others are left alone.  The THM points' model values are
/// set in `results` and each experiment's profile is taken again there (WorkingCopy::ProfileThmGroup).
void ThmModelCotangents(WorkingCopy &data, Results &results, const std::vector<GradPoint> &points,
                        const std::vector<double> &models, std::vector<double> &fitBar);

}  // namespace azure2

#endif
