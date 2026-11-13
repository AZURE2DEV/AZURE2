#ifndef AZURE2_EXPERIMENT_GRIDS_H
#define AZURE2_EXPERIMENT_GRIDS_H

#include <string>
#include <vector>

#include "azure2/experiment/LevelTables.h"
#include "azure2/experiment/Values.h"

namespace azure2 {
namespace experiment {

/*!
 * The values of an extrapolation's grid (docs/upgrade/project_format3.md section 3.8), in the engine's units (MeV,
 * deg) and, for energies, the lab frame of the entrance pair: the one implementation of the grid forms, used by
 * Schema::Build and by pyazr (pyazr.grid).
 *
 *   range lo hi <unit> <frame> step s        x_k = lo + k s, k = 0 ... floor((hi - lo)/s + 1e-9) (decision F10)
 *   range lo hi <unit> <frame> points n      x_k = lo + k (hi - lo)/(n - 1), k = 0 ... n - 1
 *   range lo hi <unit> <frame> points n log  x_k = lo (hi/lo)^(k/(n - 1)); lo > 0
 *   list x1 x2 ... <unit> <frame>            as written
 *   file "path" <unit> <frame>               one value per line ("#" starts a comment; blank lines skipped), the
 *                                            path relative to the project's directory as a data file's
 *   adaptive lo hi <unit> <frame> ...        the start grid: x_k = lo + k start_step while x_k < hi, then hi; the
 *                                            refinement is the extrapolation run's (AdaptiveGrid.h)
 *
 * Each value is taken to the engine's unit by the unit table (Units.h), in the frame it is written in; an energy
 * grid in the c.m. frame is then taken to lab energies with the entrance pair's masses (BeamLabFromCM, Frames.h), the
 * inverse of the map the engine applies to every point, so the engine computes the c.m. energies asked for (to
 * rounding).  An angle grid is kept in its frame: the observable's frame follows from it (Build).
 */
struct GridSetting {
  bool energies = true;               ///< an energy grid (else angles)
  const PairInfo *entrance = nullptr;  ///< the entrance pair: its masses take a c.m. energy grid to lab
  bool relativistic = false;          ///< the project's kinematics (<config> relativistic)
  std::string projectDir;             ///< where a grid file is read from ("" = the working directory)
};

struct GridValues {
  bool ok = false;
  std::vector<double> values;  ///< MeV (lab) or deg, in order
  double lo = 0.0, hi = 0.0;   ///< the grid's bounds, as the segment reports them (range: lo and hi as written)
  std::string file;            ///< a file grid: the path as written (opened relative to GridSettings::projectDir)
  const char *code = nullptr;  ///< when not ok: the diagnostic code (Diagnostics.h)
  std::string error;           ///< when not ok: why
};

/// The values of `grid` (see above); an adaptive grid gives its start grid.
GridValues ReadGridValues(const Grid &grid, const GridSetting &setting);

}  // namespace experiment
}  // namespace azure2

#endif
