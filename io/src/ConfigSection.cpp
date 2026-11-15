// The <config> section: its reader and the lines of an edited entry (ConfigSection.h).
#include "azure2/io/ConfigSection.h"

#include <algorithm>
#include <cstdlib>

#include "azure2/io/Document.h"

namespace azure2 {
namespace io {

namespace {

std::string Trim(const std::string &s, const char *blanks) {
  const std::size_t a = s.find_first_not_of(blanks);
  if (a == std::string::npos) return std::string();
  return s.substr(a, s.find_last_not_of(blanks) - a + 1);
}

bool BlankLine(const std::string &s) { return s.find_first_not_of(" \t\r\n\f\v") == std::string::npos; }

// The first token of a line and where it lies (the stream's >> on that line).
bool FirstToken(const std::string &line, std::size_t &a, std::size_t &b) {
  a = line.find_first_not_of(" \t\r\n\f\v");
  if (a == std::string::npos) return false;
  b = line.find_first_of(" \t\r\n\f\v", a);
  if (b == std::string::npos) b = line.size();
  return true;
}

}  // namespace

// The positional entries and how each is read: a token (>> then the rest of the line skipped) or a whole line.
const std::vector<ConfigPositional> &ConfigPositionals() {
  static const char *const kChecks = "none screen file";
  static const char *const kCheckDoc = "none: not written; screen: printed while the program runs; file: written into "
                                       "the checks directory";
  static const std::vector<ConfigPositional> rows = {
      {"amatrix", true, "true false", "Formalism: A-Matrix", "Runtime options",
       "true: the A-matrix calculation (needed to ignore external widths); false: the R-matrix calculation", false},
      {"output_dir", false, "", "Output directory", "Directories",
       "where the program writes its output files, relative to the project (empty: the project's directory)", false},
      {"checks_dir", false, "", "Checks directory", "Directories",
       "where the checks set to 'file' are written, relative to the project", true},
      {"check_compound_nucleus", true, kChecks, "Compound nucleus", "Checks", kCheckDoc, true},
      {"check_boundary_conditions", true, kChecks, "Boundary conditions", "Checks", kCheckDoc, true},
      {"check_data", true, kChecks, "Data", "Checks", kCheckDoc, true},
      {"check_energy_dependence", true, kChecks, "L-Matrix, phases, penetrabilities", "Checks", kCheckDoc, true},
      {"check_legendre", true, kChecks, "Legendre polynomials", "Checks", kCheckDoc, true},
      {"check_coulomb_amplitudes", true, kChecks, "Coulomb amplitudes", "Checks", kCheckDoc, true},
      {"check_pathways", true, kChecks, "Reaction pathways", "Checks", kCheckDoc, true},
      {"check_angular_distributions", true, kChecks, "Angular distributions", "Checks", kCheckDoc, true}};
  return rows;
}

namespace {

// "a b c" as "a, b or c".
std::string ValuesInWords(const std::string &values) {
  std::vector<std::string> words;
  std::size_t a = 0;
  while ((a = values.find_first_not_of(' ', a)) != std::string::npos) {
    const std::size_t b = std::min(values.find(' ', a), values.size());
    words.push_back(values.substr(a, b - a));
    a = b;
  }
  std::string out;
  for (std::size_t i = 0; i < words.size(); i++)
    out += (i == 0 ? "" : (i + 1 == words.size() ? " or " : ", ")) + words[i];
  return out;
}

// Whether `value` is one of a token entry's values.
bool PositionalValueTaken(const ConfigPositional &p, const std::string &value) {
  const std::string values = std::string(" ") + p.values + " ";
  return !value.empty() && value.find(' ') == std::string::npos && values.find(" " + value + " ") != std::string::npos;
}

// key=value of a line after the positional ones: false when it has no '=' before its comment.
bool KeyValue(const std::string &line, std::string &key, std::string &value) {
  std::string entry = line;
  const std::size_t hash = entry.find('#');
  if (hash != std::string::npos) entry = entry.substr(0, hash);
  const std::size_t eq = entry.find('=');
  if (eq == std::string::npos) return false;
  key = Trim(entry.substr(0, eq), " \t\r");
  value = Trim(entry.substr(eq + 1), " \t\r");
  return true;
}

// What follows a value that changed length by -`shrink` characters: when a comment follows after blanks, the comment
// keeps its column (blanks added, or removed down to one); otherwise unchanged.
std::string Refit(const std::string &rest, long shrink) {
  const std::size_t hash = rest.find('#');
  if (hash == std::string::npos || rest.find_first_not_of(" \t") != hash || shrink == 0) return rest;
  const long blanks = (long)hash + shrink;
  return std::string((std::size_t)std::max(1L, blanks), ' ') + rest.substr(hash);
}

const ConfigKey *KeyNamed(const std::string &name) {
  for (const ConfigKey &k : ConfigKeys())
    if (name == k.name) return &k;
  return nullptr;
}

bool OneLine(const std::string &value, std::string &why) {
  if (value.find_first_of("#\n\r") != std::string::npos) {
    why = "'" + value + "': a <config> value is one line without '#' (the line is cut there)";
    return false;
  }
  if (Trim(value, " \t") != value) {
    why = "'" + value + "': no blank before or after a <config> value (the engine trims them)";
    return false;
  }
  return true;
}

bool ValueTaken(const ConfigKey &k, const std::string &value, std::string &why) {
  if (!OneLine(value, why)) return false;
  char *end = nullptr;
  switch (k.kind) {
    case ConfigKey::Bool:
      if (value == "0" || value == "1" || value == "true" || value == "false") return true;
      why = std::string(k.name) + "=" + value + ": the value is 0, 1, true or false";
      return false;
    case ConfigKey::Int: {
      const long v = std::strtol(value.c_str(), &end, 10);
      if (!value.empty() && *end == '\0' && v >= k.minimum) return true;
      why = std::string(k.name) + "=" + value + ": the value is a whole number >= " + std::to_string(k.minimum);
      return false;
    }
    case ConfigKey::Number: {
      const double v = std::strtod(value.c_str(), &end);
      if (!value.empty() && *end == '\0' && v >= 0.0) return true;
      why = std::string(k.name) + "=" + value + ": the value is a number >= 0";
      return false;
    }
    case ConfigKey::Path:
      if (!value.empty()) return true;
      why = std::string(k.name) + ": the path is empty (remove the key instead)";
      return false;
    case ConfigKey::Version:
      why = std::string(k.name) + " is the file's format, written by the converter; it is not edited";
      return false;
  }
  return false;
}

}  // namespace

const std::vector<std::string> &ConfigPositionalNames() {
  static const std::vector<std::string> names = [] {
    std::vector<std::string> n;
    for (const ConfigPositional &p : ConfigPositionals()) n.push_back(p.name);
    return n;
  }();
  return names;
}

namespace {

ConfigKey Key(const char *name, ConfigKey::Kind kind, int minimum, const char *doc) {
  ConfigKey k;
  k.name = name;
  k.kind = kind;
  k.minimum = minimum;
  k.doc = doc;
  return k;
}

// An option of the calculation: its key (= the pyazr keyword), the program's flags, its default and what it changes.
ConfigKey Option(const char *name, const char *flag, const char *negation, bool defaultOn, const char *effect,
                 const char *doc) {
  ConfigKey k = Key(name, ConfigKey::Bool, 0, doc);
  k.option = true;
  k.flag = flag;
  k.negation = negation;
  k.defaultOn = defaultOn;
  k.effect = effect;
  return k;
}

// What a form tells a user about a key: its title, group, what its absence means, and whether it is set rarely.
ConfigKey Told(ConfigKey k, const char *title, const char *group, const char *absent, bool rare) {
  k.title = title;
  k.group = group;
  k.absent = absent;
  k.rare = rare;
  return k;
}
constexpr bool kRare = true;
constexpr bool kOften = false;

}  // namespace

// The groups, in the order a form shows them: the old setup program's "Runtime Options" first.
const std::vector<ConfigKey> &ConfigKeys() {
  static const std::vector<ConfigKey> keys = {
      Told(Key("includeAllJPi", ConfigKey::Bool, 0, "every J^pi up to maxOrbitalMomentum gets a hard-sphere channel"),
           "Hard-sphere scattering in every J^pi", "Physics", "0: only the J^pi of the levels", kOften),
      Told(Key("maxOrbitalMomentum", ConfigKey::Int, 0, "the largest l of the hard-sphere channels"),
           "Maximum orbital momentum", "Physics", "the largest l of any particle channel of the levels", kOften),
      Told(Key("abundanceSumSigma", ConfigKey::Number, 0,
               "the width of the penalty on the sum of the abundances (0: none)"),
           "Abundance-sum penalty width", "Fit", "0: no penalty", kRare),
      Told(Key("relativistic", ConfigKey::Bool, 0, "relativistic two-body kinematics"), "Relativistic kinematics",
           "Physics", "0: non-relativistic kinematics", kOften),
      Told(Key("relativisticReferencePair", ConfigKey::Int, 1, "the pair whose masses define the relativistic frames"),
           "Reference pair of the relativistic frames", "Physics", "the first pair", kRare),
      Told(Key("priorCovariance", ConfigKey::Path, 0, "the prior parameter covariance file, relative to the project"),
           "Prior covariance file", "Fit", "no prior covariance", kOften),
      Told(Key("format", ConfigKey::Version, 0, "the project format (3)"), "Project format", "File",
           "the converter writes it", kRare),
      // The options of the calculation (project_format3.md section 3.12); src/ProjectOptions.cpp says where the engine
      // keeps each one.
      Told(Option("use_brune", "--use-brune", "--no-brune", true, "convention",
                  "Brune's alternative level matrix (on-shell boundary conditions); off: Lane-Thomas constant "
                  "boundary conditions"),
           "Use Brune formalism", "Runtime options", "1", kOften),
      Told(Option("use_park", "--use-park", "", false, "convention",
                  "the fit's amplitudes are the observed reduced width amplitudes (T.-S. Park, PRC 104, 064612); "
                  "needs use_brune"),
           "Use Park's observed amplitudes", "Runtime options", "0", kOften),
      Told(Option("transform", "", "--no-transform", true, "convention",
                  "the input parameters are transformed to the run's basis at start"),
           "Perform parameter transformations", "Runtime options", "1", kOften),
      Told(Option("ignore_externals", "--ignore-externals", "--no-ignore-externals", true, "physics",
                  "the external capture of a photon channel whose input width is exactly zero is ignored"),
           "Ignore external width if internal width is zeroed", "Runtime options", "1", kOften),
      Told(Option("use_long_wavelength", "", "--no-long-wavelength", true, "physics",
                  "the long-wavelength approximation for EL capture"),
           "Long-wavelength approximation for capture", "Physics", "1", kRare),
      Told(Option("ec_convection_current", "--ec-convection-current", "--no-ec-convection-current", true, "physics",
                  "with use_long_wavelength off: external capture with the convection-current part of the EL operator "
                  "too"),
           "Convection current in external capture", "Physics", "1", kRare),
      Told(Option("use_rmc", "--use-rmc", "", false, "physics",
                  "the Reich-Moore approximation for capture (neutron capture only); needs use_brune off"),
           "Use RMC capture formalism (neutron capture only)", "Runtime options", "0", kOften),
      Told(Option("use_gsl_coul", "--gsl-coul", "", false, "numerics", "GSL's Coulomb functions instead of AZURE2's own"),
           "Use GSL Coulomb functions", "Runtime options", "0", kOften),
      Told(Option("use_coulomb_table", "--coulomb-table", "--no-coulomb-table", true, "numerics",
                  "P, S and the hard-sphere phase of charged channels from the universal Coulomb table, exact outside "
                  "it"),
           "Coulomb functions from the table", "Numerics", "1", kRare),
      Told(Option("use_ec_table", "--ec-table", "--no-ec-table", true, "numerics",
                  "external-capture integrals from per-term energy tables, interpolated"),
           "External-capture integrals from tables", "Numerics", "1", kRare),
      Told(Option("replace_grids", "--grid-replacement", "--no-grid-replacement", true, "numerics",
                  "the adaptive sub-point grids are laid out again as levels and energy shifts move"),
           "Adaptive grids follow the levels", "Numerics", "1", kRare)};
  return keys;
}

const ConfigKey *ConfigKeyNamed(const std::string &name) { return KeyNamed(name); }

const ConfigKey *ConfigOptionNamed(const std::string &name) {
  const ConfigKey *k = KeyNamed(name);
  return k && k->option ? k : nullptr;
}

bool ConfigOptionValue(const std::string &text, bool &on) {
  if (text == "1" || text == "true") return on = true, true;
  if (text == "0" || text == "false") return on = false, true;
  return false;
}

const std::vector<ConfigOptionImplication> &ConfigOptionImplications() {
  static const std::vector<ConfigOptionImplication> rows = {{"use_park", true, "use_brune", true},
                                                            {"use_brune", false, "use_park", false}};
  return rows;
}

const std::vector<ConfigOptionContradiction> &ConfigOptionContradictions() {
  static const std::vector<ConfigOptionContradiction> rows = {
      {"use_park", true, "use_brune", false, "Park's observed amplitudes are defined with Brune's shift functions"},
      {"use_rmc", true, "use_brune", true,
       "the Reich-Moore capture is computed with constant boundary conditions only (Brune is on by default)"},
      {"ec_convection_current", false, "use_long_wavelength", true,
       "the convection current is a part of the exact EL operator only; with the long-wavelength approximation the "
       "key changes nothing"}};
  return rows;
}

namespace {

// `name` = `on` as `spelling` writes it; "" for a flag that does not exist (an option without a flag for that value).
std::string Spelled(const char *name, bool on, OptionSpelling spelling) {
  if (spelling == OptionSpelling::kKeyword) return std::string(name) + "=" + (on ? "True" : "False");
  if (spelling == OptionSpelling::kKey) return std::string(name) + "=" + (on ? "1" : "0");
  const ConfigKey *k = ConfigOptionNamed(name);
  return k ? std::string(on ? k->flag : k->negation) : std::string();
}

// One way out of the contradiction: `name` to `to` (it is at `!to` now), in `spelling`; `given` says whether the call
// gave the present value (a flag to replace, a keyword to change) or not (one to add).
std::string WayOut(const char *name, bool to, OptionSpelling spelling, bool given) {
  const std::string want = Spelled(name, to, spelling), now = Spelled(name, !to, spelling);
  if (spelling == OptionSpelling::kFlag) {
    if (want.empty())  // no flag for that value: drop the flag given, else the project's key says it
      return given && !now.empty() ? "leave out " + now : WayOut(name, to, OptionSpelling::kKey, false);
    return given && !now.empty() ? want + " in place of " + now : "add " + want;
  }
  if (given) return want + " in place of " + now;
  return spelling == OptionSpelling::kKeyword ? "pass " + want : "write " + want + " in <config>";
}

}  // namespace

std::string ConfigOptionContradictionText(const ConfigOptionContradiction &c, const std::string &aWords,
                                          const std::string &bWords, OptionSpelling spelling, bool aGiven, bool bGiven) {
  auto value = [](const char *name, bool on) { return std::string(name) + "=" + (on ? "1" : "0"); };
  return aWords + " with " + bWords + " cannot hold together: " + c.why + "; allowed: " + value(c.a, c.aValue) +
         " with " + value(c.b, !c.bValue) + ", or " + value(c.a, !c.aValue) + "; fix: " +
         WayOut(c.b, !c.bValue, spelling, bGiven) + " (or " + WayOut(c.a, !c.aValue, spelling, aGiven) + ")";
}

const ConfigEntry *ConfigSection::Find(const std::string &name) const {
  const ConfigEntry *found = nullptr;
  for (const ConfigEntry &e : entries)
    if (e.name == name) found = &e;
  return found;
}

ConfigSection ReadConfigSection(const std::vector<std::string> &lines) {
  ConfigSection c;
  std::size_t next = 1;  // line 0 is the <config> tag
  for (const ConfigPositional &p : ConfigPositionals()) {
    ConfigEntry e;
    e.name = p.name;
    e.positional = true;
    if (p.token) {
      while (next < lines.size() && BlankLine(lines[next])) next++;
      std::size_t a = 0, b = 0;
      if (next >= lines.size() || !FirstToken(lines[next], a, b)) return c;  // the stream ran out: not closed
      e.value = lines[next].substr(a, b - a);
    } else {
      if (next >= lines.size()) return c;
      std::string temp = lines[next];
      const std::size_t hash = temp.find_last_of('#');
      if (hash != std::string::npos) temp = temp.substr(0, hash);
      e.value = Trim(temp, " \n\t\r");
    }
    e.line = (int)next++;
    c.entries.push_back(e);
  }
  for (; next < lines.size(); next++) {
    if (lines[next] == "</config>") {
      c.closed = true;
      break;
    }
    ConfigEntry e;
    if (!KeyValue(lines[next], e.name, e.value)) continue;
    e.line = (int)next;
    c.entries.push_back(e);
  }
  return c;
}

std::vector<std::string> ConfigLinesOf(const Document &doc, std::size_t &open) {
  std::vector<std::string> lines;
  const Document::Section *s = doc.Find("config");
  if (!s) return lines;
  open = s->open;
  const std::size_t end = s->close == Document::Section::npos ? doc.NumLines() : s->close + 1;
  for (std::size_t i = s->open; i < end; i++) {
    bool closing = false;
    const std::string tag = TagOf(doc.Line(i), closing);
    lines.push_back((i == s->open || (i == s->close && !tag.empty())) ? (closing ? "</" : "<") + tag + ">" : doc.Line(i));
  }
  return lines;
}

bool ConfigEdits(const std::vector<std::string> &lines, const std::string &name, const std::string &value,
                 ConfigChange change, std::vector<ConfigLineEdit> &edits, std::string &why) {
  edits.clear();
  const ConfigSection c = ReadConfigSection(lines);
  if (!c.closed) {
    why = "the <config> section is not complete (eleven entries, then </config>)";
    return false;
  }
  int k = 0;
  for (const ConfigPositional &p : ConfigPositionals()) {
    if (name != p.name) {
      k++;
      continue;
    }
    if (change == ConfigChange::Remove) {
      why = name + " is a positional entry of <config>: it is set, not removed";
      return false;
    }
    if (!OneLine(value, why)) return false;
    if (p.token && !PositionalValueTaken(p, value)) {
      why = name + "=" + value + ": the value is " + ValuesInWords(p.values);
      return false;
    }
    const ConfigEntry &e = c.entries[k];
    const std::string &old = lines[e.line];
    // The span [a, b) of the value as read; a token entry's value is its first token, a directory's the text before
    // the last '#', trimmed (an empty value stands at the start of the line).
    std::size_t a = 0, b = 0;
    if (p.token) {
      FirstToken(old, a, b);
    } else {
      const std::size_t hash = old.find_last_of('#');
      const std::string before = hash == std::string::npos ? old : old.substr(0, hash);
      const std::size_t first = before.find_first_not_of(" \n\t\r");
      if (first != std::string::npos) {
        a = first;
        b = before.find_last_not_of(" \n\t\r") + 1;
      }
    }
    const std::string text = old.substr(0, a) + value + Refit(old.substr(b), (long)(b - a) - (long)value.size());
    edits.push_back({e.line, text, false, false});
    return true;
  }
  const ConfigKey *key = KeyNamed(name);
  if (!key) {
    std::string list;
    for (const std::string &n : ConfigPositionalNames()) list += (list.empty() ? "" : ", ") + n;
    for (const ConfigKey &x : ConfigKeys()) list += std::string(", ") + x.name;
    why = "'" + name + "' is not an entry of <config> (the entries: " + list + ")";
    return false;
  }
  if (change == ConfigChange::Set && !ValueTaken(*key, value, why)) return false;
  if (change == ConfigChange::Remove && key->kind == ConfigKey::Version) {
    why = name + " is the file's format, written by the converter; it is not removed";
    return false;
  }
  std::vector<int> own;
  for (const ConfigEntry &e : c.entries)
    if (!e.positional && e.name == name) own.push_back(e.line);
  int close = 0;
  for (std::size_t i = 1; i < lines.size(); i++)
    if (lines[i] == "</config>") close = (int)i;
  if (change == ConfigChange::Remove) {
    for (int line : own) edits.push_back({line, "", true, false});
    return true;
  }
  if (own.empty()) {
    edits.push_back({close, name + "=" + value, false, true});
    return true;
  }
  for (std::size_t i = 0; i + 1 < own.size(); i++) edits.push_back({own[i], "", true, false});
  // The value's span in the winning line: after '=', trimmed, before the comment.
  const std::string &old = lines[own.back()];
  const std::size_t eq = old.find('='), end = std::min(old.find('#'), old.size());
  std::size_t a = old.find_first_not_of(" \t\r", eq + 1), b = a;
  if (a == std::string::npos || a > end) {
    a = b = eq + 1;
  } else {
    b = old.find_last_not_of(" \t\r", end - 1) + 1;
  }
  edits.push_back({own.back(), old.substr(0, a) + value + Refit(old.substr(b), (long)(b - a) - (long)value.size()), false,
                   false});
  return true;
}

}  // namespace io
}  // namespace azure2
