#ifndef AZURE2_EXPERIMENT_BUILD_H
#define AZURE2_EXPERIMENT_BUILD_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/ChannelRadii.h"
#include "azure2/experiment/Isospin.h"
#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/GammaDecay.h"
#include "azure2/experiment/HybridPotentialSettings.h"
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
  /// Also make the checks whose refusal the engine's set-up makes itself, by the same function (the final states of
  /// external capture, FinalStates.h): on for --validate, --describe and pyazr; off in the engine's own build, so a
  /// run refuses them once, at set-up, with the set-up's message.
  bool setUpChecks = true;
  /// Build every data set and extrapolation as if it were active (`active no` ignored): the names its parameters get
  /// when it is switched on (pyazr's Project.rename).  Never for a run.
  bool allRecordsActive = false;
};

/*!
 * Where the value of a fit parameter is written in the project (ParameterValues.h): the `occurrence`-th field `field`
 * of record `record` (0 for a field given once), and within a structured field its part `part` ("" when the field's
 * value is the nuisance itself).  E.g. a data set's own norm {"s01", "norm", 0, ""}, a shared normalization
 * {"D", "value", 0, ""}, a data set's own offset {"s01", "shift", k, "value"} where the k-th shift field is the offset,
 * the j-th parameter of an energy scale's function {"G", "param", j, "value"}.
 */
struct ParameterSource {
  std::string record;
  std::string field;
  int occurrence = 0;
  std::string part;
  /// The engine's unit of the value where it is not the base unit of the nuisance's dimension ("eV" for a screening
  /// potential, "atoms/b" for a sample thickness); "" = the base unit (ToBase), or the number as written without one.
  std::string engineUnit;
};

/// The specs of a format-3 project.  When diagnostics hold an error, the specs are incomplete and not set up.
struct BuiltProject {
  DataSpecs data;
  ExtrapolationSpecs extrapolations;
  AbundanceSpecs abundances;
  ThmMembershipSpec thm;
  std::vector<ModelDatum> modelData;  ///< the `datum` records of <modelData>, in record order (BuildModelData)
  std::vector<GammaDecaySpec> gammaDecays;  ///< the `gamma_decay` records of <gammaDecays> (BuildGammaDecays)
  /// The `potential` records of <potentials>: one setting per pair key, no default (a pair without a record has no
  /// potential).
  HybridPotentialSettings potentials;
  /// The `channel_radius` records of <channelRadii>, in record order (BuildChannelRadii): the channels that take a
  /// radius of their own.  Which channel takes which is the set-up's (ChannelRadiusFor), which knows the channels.
  std::vector<ChannelRadiusEntry> channelRadii;
  /// The records of <isospin> (BuildIsospin): the isospin levels across compound nuclei, resolved against the pairs
  /// and level IDs of <levels>.  Which level lines, channels and parameters they tie is the set-up's (IsospinSet).
  IsospinSpec isospin;
  /// The record of each pair's potential: pair key -> (record name, line).
  std::map<int, std::pair<std::string, int>> potentialRecords;
  /// The defaults of the table the potential records used (matching radius, step), one line each for the log.
  std::vector<std::string> potentialNotes;
  /// The record field each data-set parameter's value is written in (norm, shift offset and sqrt, the shift function's
  /// parameters, a target isotope's abundance), by parameter name; a group's parameter only for its first active
  /// member, whose parameter it is.
  std::map<std::string, ParameterSource> parameterSources;
  /// Diagnostics of the build checks, ordered by line.
  std::vector<Diagnostic> diagnostics;
  /// The checks that need a calculation and are made by the engine's set-up, which this build does not do (the
  /// quantities' refusals marked "build" that Build does not make here); 0 when every check was made.
  int checksLeftToEngine = 0;
  bool Ok() const;
};

/*!
 * Every fit-parameter name a built project gives that is made from a record's name, in record order, each once: the
 * data sets' (norm, shift offset and sqrt, shift-function parameters, sample thickness, a group follower's
 * <record>.<field>[<data set>]), the extrapolations', the isospin records' (<isospin_shift>.shift,
 * <isospin_mix>.amplitude), the isotope abundances' (AbundanceParameterName) and the target
 * effects' (<target>.screening, .areal_density, .doppler: the names src/EDataParameters.cpp EffectParameterName gives).
 * pyazr's Project.rename maps the names before and after a rename by position.
 */
std::vector<std::string> ParameterNamesOf(const BuiltProject &built);

}  // namespace experiment
}  // namespace azure2

#endif
