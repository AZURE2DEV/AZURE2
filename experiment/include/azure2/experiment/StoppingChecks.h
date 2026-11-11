#ifndef AZURE2_EXPERIMENT_STOPPINGCHECKS_H
#define AZURE2_EXPERIMENT_STOPPINGCHECKS_H

#include <string>

#include "azure2/experiment/Diagnostics.h"

namespace azure2 {
namespace experiment {

/*!
 * The refusals of a target's stopping that only the engine's set-up can find (docs/upgrade/energy_grid.md section 6):
 * it builds the material of `stopping evaluation` or `stopping table` with the stopping/ library, which the validator
 * does not link (experiment/ uses physics/ and io/'s records only), and evaluates it at the energies of each point's
 * sub-point grid.  The engine prints the diagnostic and refuses the project.  Energies are lab MeV.
 */

/// The library cannot give the stopping at all: SRIM's data files are neither compiled in nor in AZURE2_SRIM_DATA_DIR,
/// or the table file cannot be read; \p why is the library's own message.
Diagnostic StoppingDataMissing(const std::string &target, const std::string &why);
/// The stopping of \p target for \p projectile (a label, "1H") is asked for at lab energy \p atLab outside its range
/// (SRIM's 0 < E/M <= 1e7 keV/u, or the ends of a table); \p why is the library's message.
Diagnostic StoppingOutOfRange(const std::string &target, const std::string &projectile, double atLab, const std::string &why);

}  // namespace experiment
}  // namespace azure2

#endif
