#ifndef AZURE2_EXPERIMENT_BUILD_H
#define AZURE2_EXPERIMENT_BUILD_H

#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/LevelTables.h"
#include "azure2/experiment/SetupSpecs.h"

namespace azure2 {
namespace experiment {

/*!
 * What Schema::Build needs besides the validated records, and what it makes (docs/upgrade/project_format3.md
 * step 3b): the pairs and levels of the project's <levels> section (LevelTables.h), against which reactions and
 * model data are resolved (section 3.5), and the engine's set-up specs (SetupSpecs.h) with the diagnostics of the checks the validator leaves to it.
 */

/// Whether the run the project is built for uses its data (BuildOptions::dataFiles).
enum class DataFiles {
  Required,  ///< a fit or calculation with data: an unreadable data file of an active data set is E_DATA_FILE
  NotUsed    ///< an extrapolation (no data): an unreadable data file is W_DATA_FILE_NOT_READ, and nothing of its
             ///< contents is checked
};

struct BuildOptions {
  /// The directory data files are read from, as the engine opens them ("" = the working directory, the engine's).
  std::string projectDir;
  /// Read every active data set's file (it must exist and hold the columns its roles name).  Off only in tests of
  /// projects whose files are not at hand.
  bool readDataFiles = true;
  /// Whether this run uses the data: --validate, --describe and a data run require readable files; the engine's
  /// extrapolation run (and pyazr's) does not.
  DataFiles dataFiles = DataFiles::Required;
  /// The project's relativistic flag (<config>): THM data sets are refused with it.
  bool relativistic = false;
};

/// The specs of a format-3 project.  When diagnostics hold an error, the specs are incomplete and not set up.
struct BuiltProject {
  DataSpecs data;
  ExtrapolationSpecs extrapolations;
  AbundanceSpecs abundances;
  ThmMembershipSpec thm;
  std::vector<ModelDatum> modelData;  ///< the `datum` records of <modelData>, in record order (BuildModelData)
  /// Diagnostics of the build checks, ordered by line.
  std::vector<Diagnostic> diagnostics;
  /// The checks that need a calculation and are made by the engine's set-up, which this build does not do (the
  /// quantities' refusals marked "build" that Build does not make here); 0 when every check was made.
  int checksLeftToEngine = 0;
  bool Ok() const;
};

}  // namespace experiment
}  // namespace azure2

#endif
