#ifndef AZURE2_EXPERIMENT_PARAMETERFILENAMES_H
#define AZURE2_EXPERIMENT_PARAMETERFILENAMES_H

#include <string>
#include <vector>

#include "azure2/experiment/Diagnostics.h"

namespace azure2 {
namespace experiment {

/*!
 * The names of a parameter file (a param.sav or param.par given as the start of a run: one "name value error" line per
 * parameter) against the parameters of the run it is read into.  A file is applied by name, so a name the run does not
 * have is skipped and a parameter the file does not name keeps the project's value; both were silent, and a file
 * written for other names (the old segment_<k>_norm of a project since converted to format 3) changed a chi^2 without
 * a word.  The match says which names, so that the run can say so.
 */
struct ParameterFileNames {
  std::vector<std::string> inFile;     ///< the file's parameter names, in file order (not "#parametrization")
  std::vector<std::string> ofRun;      ///< the parameters of this run, in parameter order
  /// The run's parameters a file may leave out without a word: data-set and extrapolation nuisances the project
  /// fixes (a fit's param.sav has no extrapolation shifts, an extrapolation's param.par no data-set norms).
  std::vector<std::string> optionalOfRun;
  /// Names the project gives to the parameters of another run (a format-3 project's data-set parameters in an
  /// extrapolation run, its extrapolations' in a data run: ParameterNamesOf): not applied here, and not unknown.
  std::vector<std::string> ofOtherRuns;
};

/*!
 * The diagnostics of the match, at most one, with the file named in `file`:
 *  - no name of the file is a parameter of the run (or the file names none): E_PARAMETER_FILE_NO_MATCH, an error,
 *    listing the file's names; the run is refused, since nothing of the file would be applied;
 *  - some names match: W_PARAMETER_FILE_NAMES, listing the file's names the project does not have and the run's
 *    parameters (not optional ones) the file does not name;
 *  - every name matches: none.
 */
std::vector<Diagnostic> MatchParameterFileNames(const std::string &file, const ParameterFileNames &names);

/// A row of <parameterSettings> stored under a parameter's own name (the one-token form), with its line in the file.
struct ParameterSettingRow {
  std::string name;
  int line = 0;  ///< 1-based, in the project file
};

/*!
 * The <parameterSettings> rows that name no parameter of the project (`known`: the run's parameters and those the
 * project gives other runs, ParameterNamesOf): W_PARAMETER_SETTING_UNKNOWN for each, with its line and name.  Such a
 * row is ignored; before 2026-11-11 it was ignored without a word (a row naming an old parameter name, as
 * screening_1 after the conversion to format 3, tested nothing: tests/nuisance_kinds).  A warning, not an error: the
 * GUI writes a row for every parameter it shows and leaves rows of parameters since removed.
 */
std::vector<Diagnostic> UnknownParameterSettings(const std::string &file, const std::vector<ParameterSettingRow> &rows,
                                                 const std::vector<std::string> &known);

}  // namespace experiment
}  // namespace azure2

#endif
