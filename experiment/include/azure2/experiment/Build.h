#ifndef AZURE2_EXPERIMENT_BUILD_H
#define AZURE2_EXPERIMENT_BUILD_H

#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/SetupSpecs.h"

namespace azure2 {
namespace experiment {

/*!
 * What Schema::Build needs besides the validated records, and what it makes (docs/upgrade/project_format3.md
 * step 3b): the pairs of the project's <levels> section, against which reaction labels are resolved (section 3.5),
 * and the engine's set-up specs (SetupSpecs.h) with the diagnostics of the checks the validator leaves to it.
 */

/// One particle pair of <levels> (the columns of its rows that name the pair; io::LevelRow).
struct PairRow {
  int key = 0;    ///< the pair key (column 6)
  int pType = 0;  ///< 0 particles, 10 photon + nucleus, 20 beta-delayed
  int z1 = 0, z2 = 0;
  double m1 = 0.0, m2 = 0.0;  ///< u; the mass number is the nearest integer
  double j1 = 0.0, j2 = 0.0;  ///< spins of the light and the heavy particle
  double e1 = 0.0, e2 = 0.0;  ///< MeV: excitation of the light and of the heavy particle
};

struct BuildOptions {
  /// The directory data files are read from, as the engine opens them ("" = the working directory, the engine's).
  std::string projectDir;
  /// Read every active data set's file (it must exist and hold the columns its roles name).  Off only in tests of
  /// projects whose files are not at hand.
  bool readDataFiles = true;
  /// The project's relativistic flag (<config>): THM data sets are refused with it.
  bool relativistic = false;
};

/// The specs of a format-3 project.  When diagnostics hold an error, the specs are incomplete and not set up.
struct BuiltProject {
  DataSpecs data;
  ExtrapolationSpecs extrapolations;
  AbundanceSpecs abundances;
  ThmMembershipSpec thm;
  /// Diagnostics of the build checks, ordered by line.
  std::vector<Diagnostic> diagnostics;
  /// The checks that need a calculation and are made by the engine's set-up, which this build does not do (the
  /// quantities' refusals marked "build" that Build does not make here); 0 when every check was made.
  int checksLeftToEngine = 0;
  bool Ok() const;
};

/// The pair a reaction label names (section 3.5): the light and heavy particle and their states.  Used by Build and
/// by `--describe` to print the pairs a reaction resolved to.
std::string PairLabel(const PairRow &pair);

}  // namespace experiment
}  // namespace azure2

#endif
