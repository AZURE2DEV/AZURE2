#ifndef AZURE2_EXPERIMENT_PROJECTREADER_H
#define AZURE2_EXPERIMENT_PROJECTREADER_H

#include <memory>
#include <string>
#include <vector>

#include "azure2/experiment/Schema.h"

namespace azure2 {
namespace io {
class Document;
class LevelRow;
struct LevelSection;
}
namespace experiment {

/*!
 * The one reading path of a format-3 project (docs/upgrade/project_format3.md section 7.3, step 3b): its record
 * sections read with io::ReadRecords, validated with the table (Schema::Validate), the pairs of its <levels> section
 * (io::ReadLevelSection), and the specs built from both (Schema::Build).  Library azure2_experiment_reader (it links
 * io/; the table library azure2_experiment does not).  Used by the engine (Config::ReadConfigFile), by
 * `AZURE2 --validate` / `--describe`, by pyazr.validate and by the converter on its own output.
 */
struct Format3Project {
  std::vector<Diagnostic> diagnostics;  ///< the reader's, the validator's and the build's, each with the file name
  ValidationResult validated;
  LevelTables levels;  ///< the pairs and levels of <levels> (LevelTablesOfRows)
  /// The <levels> section as read once (null when the file has none): the engine fills its compound nucleus from
  /// these rows (CNuc::Fill), the tables above are made from them.
  std::shared_ptr<const io::LevelSection> levelSection;
  std::string projectDir;  ///< the directory the data files were read from ("" = the working directory)
  BuiltProject built;  ///< empty when the validation found an error (Build needs a valid project)
  bool levelsRead = false;  ///< <levels> was read without a malformed line
  /// The checks the table leaves to the build that were not made: all of them when the validation failed or no
  /// <levels> could be read, else Build's own count (0).
  int buildChecksSkipped = 0;
  bool Ok() const;  ///< no diagnostic of severity error
};

struct Format3ReadOptions {
  std::string fileName;    ///< set on every diagnostic
  std::string projectDir;  ///< data files are read relative to it
  /// Run Schema::Build (needs <levels> and the data files); false: the reader and the validator only.
  bool build = true;
  bool readDataFiles = true;
  DataFiles dataFiles = DataFiles::Required;  ///< BuildOptions::dataFiles
  bool setUpChecks = true;                     ///< BuildOptions::setUpChecks
  bool allRecordsActive = false;               ///< BuildOptions::allRecordsActive
};

Format3Project ReadFormat3Project(const io::Document &doc, const Format3ReadOptions &options);
Format3Project ReadFormat3Text(const std::string &text, const Format3ReadOptions &options);

/// The rows of the <levels> section of a document.  `ok` false when the section is absent or holds a malformed line.
std::vector<io::LevelRow> LevelRowsOf(const io::Document &doc, bool &ok);

/*!
 * The pairs and levels of <levels> rows as the plain tables of LevelTables.h, made once and read by Schema::Build,
 * BuildModelData and the engine's model data.  The indices of a level are the ones CNuc::Fill gives: J^pi groups
 * numbered by their first active row (a group is J, parity and compound nucleus), channels of a group by their first
 * active row (pair, l, s), levels of a group by io::LevelGrouper.  Pairs: one per key, from its first row, active or
 * not, in row order.
 */
LevelTables LevelTablesOfRows(const std::vector<io::LevelRow> &rows);

/// The pair of one row (its pair columns).
PairInfo PairInfoOfRow(const io::LevelRow &row);

/// The project's relativistic flag (<config> key relativistic=1|true).
bool RelativisticFlag(const io::Document &doc);

}  // namespace experiment
}  // namespace azure2

#endif
