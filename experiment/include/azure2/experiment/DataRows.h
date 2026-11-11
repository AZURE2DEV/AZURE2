#ifndef AZURE2_EXPERIMENT_DATAROWS_H
#define AZURE2_EXPERIMENT_DATAROWS_H

// One row of a data file into the engine's units and frame (SetupSpecs.h, ColumnConversion): the one place where the
// unit of a column of a data file is applied.  The engine calls it on every row it reads, before the point is made.

#include <vector>

#include "azure2/experiment/SetupSpecs.h"

namespace azure2 {
namespace experiment {

/*!
 * The row (energy, angle, value, error, extra columns) as written -> in MeV, deg and the engine unit of the quantity:
 * every factor through Units.h Rescale (a power of ten divides, so a value exact in the written unit gives the
 * correctly rounded engine value); an E_cm column is then taken to the lab energy of the beam on the heavy particle at
 * rest (Frames.h BeamLabFromCM, the inverse of the engine's map).  The extra columns are energies (window_lo,
 * window_hi of a beam-profile kernel), converted with the energy factor.  The identity conversion leaves every bit.
 */
void RowToEngine(const ColumnConversion &c, double &energy, double &angle, double &value, double &error,
                 std::vector<double> &extra);

}  // namespace experiment
}  // namespace azure2

#endif
