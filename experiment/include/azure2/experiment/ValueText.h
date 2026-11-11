#ifndef AZURE2_EXPERIMENT_VALUETEXT_H
#define AZURE2_EXPERIMENT_VALUETEXT_H

#include <string>

#include "azure2/experiment/Build.h"
#include "azure2/experiment/Schema.h"
#include "azure2/experiment/Values.h"

namespace azure2 {
namespace experiment {

/*!
 * The text a writer puts into a format-3 project for a value (docs/upgrade/project_format3.md sections 3.3, 4.1):
 * the one writer of a nuisance's text, used by the converter and by pyazr (pyazr.Nuisance), next to its one parser
 * (ValueParsers.cpp); and the edit that writes a fitted parameter's value into the nuisance it comes from (pyazr's
 * save_fit of a format-3 project).
 */

/// The words of a nuisance (section 3.3), each as it is to be written.
struct NuisanceWords {
  std::string value;  ///< the number as written; empty for a profiled nuisance
  std::string unit;   ///< a unit token, or "" (a dimensionless nuisance)
  std::string mode;   ///< the mode's word (NuisanceModeName: fixed, fit, profiled), written as given
  std::string priorCentre;  ///< with priorWidth: "prior <centre> <width>"; both empty: no prior
  std::string priorWidth;   ///< in the nuisance's unit, or a percentage of the centre ("5%")
  std::string gridMax;      ///< "" or the largest value the sub-point grids must span
};

/// "<value> [<unit>] <mode> [prior <centre> <width>[%]] [grid_max <value>]", or the mode alone without a value
/// ("profiled").  The text is not checked: the validator reads it (pyazr validates every record it writes).
std::string NuisanceText(const NuisanceWords &words);

/// The shortest decimal text that reads back (strtod) as `x`.
std::string ShortestNumberText(double x);

/// An edit of one token of a project's text: the token `old` at (line, col) (1-based, the file's) becomes `text`.
struct ValueEdit {
  int line = 0;
  int col = 0;
  std::string old;
  std::string text;
};

/*!
 * The edit that writes `engineValue` (in the engine's unit: MeV for an energy shift, MeV^1/2 for a sqrt term) as the
 * value of the nuisance `source` names (Build.h, ParameterSource): in the unit the nuisance is written in, as the
 * shortest decimal (ShortestNumberText), its mode and prior unchanged.  False with `why` when the source names no
 * nuisance with a value (an unknown record or field, a profiled norm).
 */
bool ParameterValueEdit(const ValidationResult &validated, const ParameterSource &source, double engineValue,
                        ValueEdit &edit, std::string &why);

}  // namespace experiment
}  // namespace azure2

#endif
