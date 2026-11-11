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

}  // namespace experiment
}  // namespace azure2

#endif
