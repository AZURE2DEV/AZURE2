// A row of a data file into the engine's units and frame (DataRows.h).
#include "azure2/experiment/DataRows.h"

#include "azure2/experiment/Frames.h"
#include "azure2/experiment/Units.h"

namespace azure2 {
namespace experiment {

void RowToEngine(const ColumnConversion &c, double &energy, double &angle, double &value, double &error,
                 std::vector<double> &extra) {
  if (c.Identity()) return;
  energy = Rescale(energy, c.energy);
  if (c.energyCM) {
    Partition entrance;
    entrance.m1 = c.lightMass;
    entrance.m2 = c.heavyMass;
    energy = BeamLabFromCM(entrance, energy, c.relativistic);
  }
  for (double &x : extra) x = Rescale(x, c.energy);
  angle = Rescale(angle, c.angle);
  value = Rescale(value, c.value);
  error = Rescale(error, c.value);
}

}  // namespace experiment
}  // namespace azure2
