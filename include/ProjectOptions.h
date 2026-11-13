#ifndef PROJECT_OPTIONS_H
#define PROJECT_OPTIONS_H

#include <string>
#include <vector>

class Config;
namespace azure2 {
namespace io {
struct ConfigSection;
}
}  // namespace azure2

/*!
 * The options of the calculation that are a property of the project (docs/upgrade/project_format3.md section 3.12):
 * the option rows of io::ConfigKeys() (name = <config> key = pyazr keyword, flags, default), where the engine keeps
 * each one in Config (src/ProjectOptions.cpp), and the rule that combines them: a flag given on the command line or a
 * keyword given to a session wins over the project's key, the key over the default.  Config::optionSources says where
 * each value came from.
 */
enum class OptionSource { kDefault, kProject, kCommandLine, kKeyword };

/// "default", "project", "command line", "keyword".
const char *OptionSourceName(OptionSource source);

/// The value of option `name` in `config`; std::invalid_argument for a name that is not an option with a storage.
bool OptionIsOn(const Config &config, const std::string &name);
/// Where the value of option `name` came from (kDefault when nothing set it).
OptionSource OptionSourceOf(const Config &config, const std::string &name);

/// An option given by a flag or a keyword: sets it, then the options it implies (io::ConfigOptionImplications) that no
/// flag or keyword stated itself, each marked as given by `source`.  std::invalid_argument for a name that is not an option.
void GiveOption(Config &config, const std::string &name, bool on, OptionSource source);

/// A flag of an option row (its flag or its negation): GiveOption(..., kCommandLine); false when `arg` is neither.
bool ApplyOptionFlag(Config &config, const std::string &arg);

/*!
 * The project's option keys, for the options not given by a flag or keyword (marked kProject), then the
 * contradictions (io::ConfigOptionContradictions), whatever the source of the two values.  False when there is one:
 * each is printed as an ERROR line (E_CONFIG_OPTIONS_CONTRADICT: both options with their sources, what is allowed, and
 * the fix in the words of the front end that gave values: a flag, a keyword, or a key) and kept in
 * Config::refusedOptions.  `tagLine`: the 1-based line of the <config> tag in the file.  A key
 * whose value is not one the engine reads is skipped (ReadFormat3Project refused it before).
 */
bool ApplyProjectOptions(Config &config, const azure2::io::ConfigSection &section, int tagLine);

/// The options of the calculation for the project text `projectText` as a run with these command-line arguments would
/// have them: each flag of an option row in `args` (other arguments ignored), then the project's keys.  False (ERROR
/// lines printed) when a flag and a key contradict each other.
bool OptionsOfRun(const std::string &projectText, const std::vector<std::string> &args, Config &config);
/// A description of the project (experiment::DescribeText or DescribeJson) with the options of `config` after it: the
/// JSON member "options" (OptionsJson), in the text the lines "options of the calculation (value, source):" and
/// OptionsText, for every project (an empty text, which describes nothing, stays empty).  AZURE2 --describe and
/// pyazr.describe give the same text.
std::string WithOptions(const std::string &described, bool json, const Config &config);

/// "  use_brune = 1 (default)" per option, one line each, in the table's order.
std::string OptionsText(const Config &config);
/// {"use_brune": {"value": true, "source": "default"}, ...} in the table's order.
std::string OptionsJson(const Config &config);

#endif
