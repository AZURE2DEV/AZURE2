#ifndef AZURE2_IO_CONFIGSECTION_H
#define AZURE2_IO_CONFIGSECTION_H

#include <string>
#include <vector>

namespace azure2 {
namespace io {

/*!
 * The <config> section: its one reader (Config::ReadConfigFile reads the project through it, pyazr through the
 * binding) and the line of an edited entry.  The grammar, as the engine has always read it:
 *
 *   - after the <config> tag, eleven positional entries: `amatrix` (the first token of its line: true = A-matrix
 *     calculation), `output_dir` and `checks_dir` (the whole line before its last '#', trimmed; may be empty), then
 *     the eight checks `check_compound_nucleus`, `check_boundary_conditions`, `check_data`,
 *     `check_energy_dependence`, `check_legendre`, `check_coulomb_amplitudes`, `check_pathways`,
 *     `check_angular_distributions` (the first token: none | screen | file).  An entry read as a token skips blank
 *     lines before it (the stream's >>); the rest of its line is ignored.
 *   - then `key=value` lines up to </config>: the text before the first '#', split at the first '=', key and value
 *     trimmed; a line without '=' is skipped; the last line of a key wins.  The keys the engine reads are
 *     ConfigKeys(); an unknown key is skipped by the engine, refused by an edit.
 */
struct ConfigEntry {
  std::string name;   ///< a positional name above, or the key
  std::string value;  ///< as the engine reads it (trimmed; a token entry's first token)
  int line = 0;       ///< index into the lines given to ReadConfigSection (0 = the <config> tag)
  bool positional = false;
};

struct ConfigSection {
  std::vector<ConfigEntry> entries;  ///< the positional entries, then every key=value line in file order
  bool closed = false;               ///< the eleven entries were read and the </config> line was reached
  /// The value of the last entry named `name`, or nullptr.
  const ConfigEntry *Find(const std::string &name) const;
};

/// One of the eleven positional entries: how it is read and what a user is told about it.  The words (`title`,
/// `group`, `rare`) are for the forms that set it (pyazr.kinds()["config"], the web page's run settings); the engine
/// reads none of them.
struct ConfigPositional {
  const char *name;
  bool token;          ///< read as its line's first token (one of `values`); else the whole line (a directory)
  const char *values;  ///< a token entry's values, separated by blanks, in the order a form offers them
  const char *title;   ///< the old setup program's words for it
  const char *group;   ///< the group of settings it belongs to (ConfigKey::group)
  const char *doc;     ///< what it sets
  bool rare;           ///< set rarely (ConfigKey::rare)
};
/// The positional entries, in file order.
const std::vector<ConfigPositional> &ConfigPositionals();
/// The names of the eleven positional entries, in file order.
const std::vector<std::string> &ConfigPositionalNames();

/// One key of the key=value part the engine reads.
struct ConfigKey {
  const char *name;
  /// bool (0 | 1 | true | false), int (>= minimum), number (>= 0), path (one line, no '#'), version (not editable)
  enum Kind { Bool, Int, Number, Path, Version } kind;
  int minimum;      ///< Int: the smallest value taken
  const char *doc;  ///< what it sets
  /*!
   * An option of the calculation (docs/upgrade/project_format3.md section 3.12): a bool key that is also a flag of the
   * program and a keyword of a pyazr session, under the same name.  A flag or keyword given wins over the key, the key
   * over `defaultOn`.  `flag` sets it on, `negation` off ("" where the program has none); `effect` says what it
   * changes: "physics", "convention" or "numerics".
   */
  bool option = false;
  const char *flag = "";
  const char *negation = "";
  bool defaultOn = false;
  const char *effect = "";
  // What a form that sets the key tells a user (pyazr.kinds()["config"], the web page's run settings); the engine
  // reads none of these.
  const char *title = "";   ///< the old setup program's words for it where it had the setting, else a short name
  const char *group = "";   ///< the group of settings a user finds it in ("Runtime options", "Physics", ...)
  const char *absent = "";  ///< what the engine does when the key is absent, in words (an option: its default)
  bool rare = false;        ///< set rarely: a form shows it one click away ("more"), unless the project sets it
};
/// The keys, in the order the engine applies them: the options in the order a session applies its keywords (use_brune
/// before use_park, so that use_park=True with use_brune=False gives Park, as before the keys).
const std::vector<ConfigKey> &ConfigKeys();
/// The key named `name` (an option or not), or nullptr when <config> has no such key.
const ConfigKey *ConfigKeyNamed(const std::string &name);
/// The option named `name`, or nullptr when it is not an option.
const ConfigKey *ConfigOptionNamed(const std::string &name);
/// An option key's value: 1 | true is on, 0 | false off; false when the text is neither.
bool ConfigOptionValue(const std::string &text, bool &on);

/// A flag or keyword that also sets another option (a project's key implies nothing): `name` = `value` gives `implied`
/// = `impliedValue`.
struct ConfigOptionImplication {
  const char *name;
  bool value;
  const char *implied;
  bool impliedValue;
};
const std::vector<ConfigOptionImplication> &ConfigOptionImplications();

/// Two option values that cannot hold together (E_CONFIG_OPTIONS_CONTRADICT): `a` = `aValue` with `b` = `bValue`.  Refused
/// whatever the source of the two values (a flag, a keyword, a key or the default); an option one flag or keyword
/// implies counts as given, so that --use-park alone (Brune implied on) is no contradiction.  `why` says why in words.
struct ConfigOptionContradiction {
  const char *a;
  bool aValue;
  const char *b;
  bool bValue;
  const char *why;
};
const std::vector<ConfigOptionContradiction> &ConfigOptionContradictions();

/// How the front end that meets a contradiction spells an option's value in the refusal's fix: the program's flag
/// ("--no-brune"), a session's keyword ("use_brune=False"), the project's key ("use_brune=0").
enum class OptionSpelling { kFlag, kKeyword, kKey };

/// A key <config> had and has no more.  A file that still holds the line is refused with either value
/// (E_CONFIG_KEY_REMOVED, experiment/'s reader), and so are the flag the key had as an option of the calculation
/// (`flag`, "" when it had none) and its keyword of a pyazr session: the key says nothing any more, and going on
/// without it would hide that from the user who wrote it.  `why` says in a sentence what replaced it.
struct RemovedConfigKey {
  const char *name;
  const char *flag;
  const char *why;
};
const std::vector<RemovedConfigKey> &RemovedConfigKeys();
/// The removed key named `name`, or nullptr.
const RemovedConfigKey *RemovedConfigKeyNamed(const std::string &name);
/// The removed key whose flag was `flag` ("--gsl-coul"), or nullptr.
const RemovedConfigKey *RemovedConfigKeyOfFlag(const std::string &flag);
/// The one refusal of a removed key, in the spelling of the front end that met it: "<name> is no longer a key of
/// <config>: delete the line.  <why>" (kKey), "<flag> is no longer an option of AZURE2: leave it out.  <why>" (kFlag),
/// "<name> is no longer a keyword of pyazr.Session: leave it out.  <why>" (kKeyword).
std::string RemovedConfigKeyText(const RemovedConfigKey &key, OptionSpelling spelling);

/// The refusal's text in three parts, one line: what ("use_rmc=1 (command line) with use_brune=1 (the default) cannot
/// hold together: <why>"), "allowed: ..." and "fix: ..." (b's other value, else a's, in `spelling`).  `aWords` and
/// `bWords` are the two values with their sources as the caller names them; `aGiven` and `bGiven` say whether the
/// front end of `spelling` gave the value (a flag to replace rather than one to add).
std::string ConfigOptionContradictionText(const ConfigOptionContradiction &c, const std::string &aWords,
                                          const std::string &bWords, OptionSpelling spelling, bool aGiven, bool bGiven);

/// Reads the section from its lines as Config::SectionStream gives them: the <config> tag line first, then the body,
/// then the </config> line (absent when the section is not closed).
ConfigSection ReadConfigSection(const std::vector<std::string> &lines);

class Document;
/// The lines of a document's <config> section as Config::SectionStream gives them (the tags written "<config>",
/// "</config>"), and in `open` the index of its tag line in the document; empty when there is no such section.
std::vector<std::string> ConfigLinesOf(const Document &doc, std::size_t &open);

/// The edit of one entry: replace line `line` (index as above) by `text`, or remove it (`remove`), or insert `text`
/// before line `line` (`insert`).
struct ConfigLineEdit {
  int line = 0;
  std::string text;
  bool remove = false;
  bool insert = false;
};

/*!
 * The edits that set entry `name` to `value` in the section of `lines` (as for ReadConfigSection): a positional entry
 * rewrites its line, the value where the old one began and the rest of the line (its comment) kept in its column; a
 * key rewrites the line that wins (`key=value`, its comment kept) and removes the earlier lines of the key, or adds
 * `key=value` before </config>; ConfigChange::Remove (keys only) removes every line of the key.  False with `why` when the name is
 * neither a positional entry nor a key of ConfigKeys(), or the value is not one the engine reads (a word of the
 * entry's list, a number, one line without '#').
 */
enum class ConfigChange { Set, Remove };
bool ConfigEdits(const std::vector<std::string> &lines, const std::string &name, const std::string &value,
                 ConfigChange change, std::vector<ConfigLineEdit> &edits, std::string &why);

}  // namespace io
}  // namespace azure2

#endif
