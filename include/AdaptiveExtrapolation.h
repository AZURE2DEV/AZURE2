#ifndef ADAPTIVEEXTRAPOLATION_H
#define ADAPTIVEEXTRAPOLATION_H

#include <memory>
#include <vector>

#include "Constants.h"

class CNuc;
class Config;
class LevelInputs;
class LevelState;
class WorkingCopy;

/*!
 * The extrapolation run's adaptive energy grids (docs/upgrade/project_format3.md section 3.8).  An adaptive grid is a
 * function of the parameters, so it is made in the run, at the run's parameters, before the extrapolations are
 * computed: the extrapolation data set the run then computes and writes is made on the refined grids, as on `list`
 * grids of the same energies.
 *
 * For each active extrapolation with an adaptive grid (in file order): its start grid (Build: the uniform grid of
 * start_step) with, as seeds, E_lambda + k Gamma_lambda for k = 0, +-1/2, +-1, +-2, +-4 of every level of the R-matrix
 * of the entrance pair's compound nucleus (the observed energy E_lambda as c.m. energy of the entrance pair and the
 * total width Gamma_lambda of the open channels, CNuc::ObservedWidthsAt; the formal energy alone if the observed
 * values cannot be computed; no seeds for a beta-delayed extrapolation) that fall inside its range, is refined by
 * experiment::RefineGrid.  Each sweep's batch of energies is computed by the point pass of the run (AZURECalc on an
 * extrapolation data set made by WorkingCopy::MakePoints from the project's specs, the batch as this extrapolation's grid
 * and every other extrapolation reduced to its first point, so that the parameter vector's layout is the run's), in
 * a scratch directory for the external-capture integrals and without messages.  The value refined is the point
 * pass's model value at each angle of the extrapolation (the union grid: an energy is added at every angle when any
 * angle needs it): a c.m. cross section, differential cross section, phase shift or analyzing power, before the
 * output's frame and S-factor conversions.
 */

/// Whether the project has an active extrapolation with an adaptive grid (a format-3 project only).
bool HasAdaptiveGrids(const Config &configure);

/// A pyazr session evaluates at many parameter vectors and an adaptive grid is made at one, by the extrapolation run:
/// true, after an error message saying so, when the project has an adaptive grid (AZUREAPI's set-up refuses it).
bool SessionRefusesAdaptiveGrids(const Config &configure);

/// The extrapolation data set on the refined grids, made by WorkingCopy::MakePoints and laid out on the parameter vector as
/// the run's data set, for the parameters `params` (the run's vector); nullptr after an error message.  The
/// compound nucleus is initialized and its pairs marked as the run's MakePoints marked them.
std::unique_ptr<WorkingCopy> MakeRefinedExtrapolations(const Config &configure, const CNuc *compound, LevelState &state, const LevelInputs &inputs, const vector_r &params);

#endif
