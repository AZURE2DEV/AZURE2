#ifndef AZURE2_EXPERIMENT_LAYERCHECKS_H
#define AZURE2_EXPERIMENT_LAYERCHECKS_H

#include <string>

#include "azure2/experiment/Diagnostics.h"

namespace azure2 {
namespace experiment {

/*!
 * The refusals of a yield on a target layer (quantity yield, target kind layer; docs/upgrade/energy_grid.md section 5)
 * that only the engine's set-up can find: it evaluates the target's stopping on each point's sub-point grid and finds
 * where the beam leaves the layer (include/LayerExit.h).  The validator cannot: it neither evaluates the formula nor
 * reads the data energies.  The engine prints the diagnostic and refuses the project; the messages say what to write
 * instead.  All energies are lab MeV, the stopping MeV cm^2 per active nucleus.
 */
struct LayerPoint {
  std::string target;    ///< the target record's name
  std::string dataFile;  ///< the data set's file (or "" for an extrapolation), to name the point
  double beamLab = 0.0;  ///< the beam energy E_0 of the point
};

/// The stopping is not positive and finite at the lab energy \p atLab inside the layer (\p value as evaluated).
Diagnostic LayerStoppingUndefined(const LayerPoint &p, double atLab, double value);
/// The beam does not leave the layer above \p lowestLab, the lowest energy the model is evaluated at (1 keV c.m.,
/// TargetEffect::minIntegrationEnergy): the areal density \p density (atoms/cm^2) is more than the beam can cross.
Diagnostic LayerExitBelowModel(const LayerPoint &p, double density, double lowestLab);

}  // namespace experiment
}  // namespace azure2

#endif
