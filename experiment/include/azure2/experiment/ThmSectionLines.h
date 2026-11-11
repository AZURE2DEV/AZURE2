#ifndef AZURE2_EXPERIMENT_THMSECTIONLINES_H
#define AZURE2_EXPERIMENT_THMSECTIONLINES_H

#include <string>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * One line of a project's <thm> section as far as the names go (docs/upgrade/project_format3.md section 3.9): whether
 * it is an experiment line `experiment[<name>] key=value ...` and which experiment it names, else the key of a global
 * line `key=value`.  The one reading of an experiment's name: the reading path (ProjectReader: which names a `thm`
 * field may use, the old member lists) and the engine's <thm> reader and writer (src/ConfigThm.cpp,
 * src/ThmBlockEdits.cpp) use it.  The values of the keys are the engine's THM grammar (src/ThmExperiment.cpp).
 */
struct ThmSectionLine {
  std::string text;                ///< the line before '#', without the blanks at either end
  bool experiment = false;         ///< text starts with "experiment["
  std::string name;                ///< the experiment's name: the text between '[' and the first ']' ("" without ']')
  std::vector<std::string> keys;   ///< an experiment line: the key (text before '=') of each blank-separated token
                                   ///< after ']'; a global line: its one key (text before '=', trimmed), if it has '='
};

ThmSectionLine ReadThmSectionLine(const std::string &line);

/// The key of an old member or weight list that format 3 does not have, with the record field that replaces it, or
/// "" for a line without one: `segments=` on an experiment line (a data set names its experiment: `thm <name>`), a
/// global `weight[<k>]=` / `weightTest[<k>]=` (the data set's or extrapolation's `thm_weight <file>`).
std::string ThmOldListKey(const ThmSectionLine &line, std::string &field);

}  // namespace experiment
}  // namespace azure2

#endif
