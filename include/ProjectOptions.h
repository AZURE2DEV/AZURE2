#ifndef PROJECT_OPTIONS_H
#define PROJECT_OPTIONS_H

#include <string>

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

/// An option given by a flag or a keyword: sets it, then the options it implies (io::ConfigOptionImplications), each
/// marked as given by `source`.  std::invalid_argument for a name that is not an option.
void GiveOption(Config &config, const std::string &name, bool on, OptionSource source);

/// A flag of an option row (its flag or its negation): GiveOption(..., kCommandLine); false when `arg` is neither.
bool ApplyOptionFlag(Config &config, const std::string &arg);

/*!
 * The project's option keys, for the options not given by a flag or keyword (marked kProject), then the
 * contradictions (io::ConfigOptionContradictions) in which at least one value came from the project.  False when
 * there is one: each is printed as an ERROR line (E_CONFIG_OPTIONS_CONTRADICT, naming both options and their
 * sources) and kept in Config::refusedOptions.  `tagLine`: the 1-based line of the <config> tag in the file.  A key
 * whose value is not one the engine reads is skipped (ReadFormat3Project refused it before).
 */
bool ApplyProjectOptions(Config &config, const azure2::io::ConfigSection &section, int tagLine);

/// "  use_brune = 1 (default)" per option, one line each, in the table's order.
std::string OptionsText(const Config &config);
/// {"use_brune": {"value": true, "source": "default"}, ...} in the table's order.
std::string OptionsJson(const Config &config);

#endif
