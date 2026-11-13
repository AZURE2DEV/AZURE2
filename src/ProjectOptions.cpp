// The options of the calculation that are a property of the project: where the engine keeps each, and how flags,
// keywords and the project's keys combine (include/ProjectOptions.h; docs/upgrade/project_format3.md section 3.12).
#include "ProjectOptions.h"

#include <sstream>
#include <stdexcept>
#include <vector>

#include "Config.h"
#include "azure2/experiment/Diagnostics.h"
#include "azure2/io/ConfigSection.h"
#include "azure2/io/Document.h"

namespace {

// Where the engine keeps one option: a bit of Config::paramMask, or a bool member of Config.  One row per option row
// of io::ConfigKeys() (tests/upgrade/project_options_test.cpp fails for a row without one).
struct OptionStorage {
  const char *name;
  unsigned int maskBit;
  bool Config::*member;
};
const std::vector<OptionStorage> &OptionStorages() {
  static const std::vector<OptionStorage> rows = {
      {"use_brune", Config::USE_BRUNE_FORMALISM, nullptr},
      {"use_park", Config::USE_PARK_FORMALISM, nullptr},
      {"transform", Config::TRANSFORM_PARAMETERS, nullptr},
      {"ignore_externals", Config::IGNORE_ZERO_WIDTHS, nullptr},
      {"use_long_wavelength", Config::USE_LONGWAVELENGTH_APPROX, nullptr},
      {"ec_convection_current", 0, &Config::ecConvectionCurrent},
      {"use_rmc", Config::USE_RMC_FORMALISM, nullptr},
      {"use_gsl_coul", Config::USE_GSL_COULOMB_FUNC, nullptr},
      {"use_coulomb_table", 0, &Config::useCoulombTable},
      {"use_ec_table", 0, &Config::useECTable},
      {"replace_grids", 0, &Config::replaceSubPointGrids}};
  return rows;
}

const OptionStorage &StorageOf(const std::string &name) {
  for (const OptionStorage &s : OptionStorages())
    if (name == s.name && azure2::io::ConfigOptionNamed(name)) return s;
  throw std::invalid_argument("'" + name + "' is not an option of the calculation with a place in Config");
}

void Store(Config &config, const std::string &name, bool on) {
  const OptionStorage &s = StorageOf(name);
  if (s.member)
    config.*s.member = on;
  else if (on)
    config.paramMask |= s.maskBit;
  else
    config.paramMask &= ~s.maskBit;
}

// "use_brune=1 (command line)": an option's value and its source, for a message.
std::string OptionWords(const Config &config, const char *name) {
  const OptionSource source = OptionSourceOf(config, name);
  return std::string(name) + "=" + (OptionIsOn(config, name) ? "1" : "0") + " (" +
         (source == OptionSource::kProject   ? "this file"
          : source == OptionSource::kDefault ? "the default"
                                             : OptionSourceName(source)) +
         ")";
}

}  // namespace

const char *OptionSourceName(OptionSource source) {
  switch (source) {
    case OptionSource::kDefault:
      return "default";
    case OptionSource::kProject:
      return "project";
    case OptionSource::kCommandLine:
      return "command line";
    case OptionSource::kKeyword:
      return "keyword";
  }
  return "default";
}

bool OptionIsOn(const Config &config, const std::string &name) {
  const OptionStorage &s = StorageOf(name);
  return s.member ? config.*s.member : (config.paramMask & s.maskBit) != 0;
}

OptionSource OptionSourceOf(const Config &config, const std::string &name) {
  const auto it = config.optionSources.find(name);
  return it == config.optionSources.end() ? OptionSource::kDefault : it->second;
}

void GiveOption(Config &config, const std::string &name, bool on, OptionSource source) {
  Store(config, name, on);
  config.optionSources[name] = source;
  config.optionsStated.insert(name);
  // An implied value fills an option the call did not state itself; a stated one stays, so that a contradiction
  // between two stated values is seen (ApplyProjectOptions) whatever their order.
  for (const azure2::io::ConfigOptionImplication &i : azure2::io::ConfigOptionImplications())
    if (name == i.name && on == i.value && !config.optionsStated.count(i.implied)) {
      Store(config, i.implied, i.impliedValue);
      config.optionSources[i.implied] = source;
    }
}

bool ApplyOptionFlag(Config &config, const std::string &arg) {
  for (const azure2::io::ConfigKey &k : azure2::io::ConfigKeys()) {
    if (!k.option) continue;
    if (*k.flag && arg == k.flag) {
      GiveOption(config, k.name, true, OptionSource::kCommandLine);
      return true;
    }
    if (*k.negation && arg == k.negation) {
      GiveOption(config, k.name, false, OptionSource::kCommandLine);
      return true;
    }
  }
  return false;
}

bool ApplyProjectOptions(Config &config, const azure2::io::ConfigSection &section, int tagLine) {
  config.refusedOptions.clear();
  for (const azure2::io::ConfigKey &k : azure2::io::ConfigKeys()) {
    if (!k.option) continue;
    const OptionSource source = OptionSourceOf(config, k.name);
    if (source == OptionSource::kCommandLine || source == OptionSource::kKeyword) continue;
    const azure2::io::ConfigEntry *e = section.Find(k.name);
    bool on = false;
    if (!e || e->positional || !azure2::io::ConfigOptionValue(e->value, on)) continue;
    Store(config, k.name, on);
    config.optionSources[k.name] = OptionSource::kProject;
  }
  // The front end that gave values spells the fix: a session's keywords, else the program's flags, else the keys.
  azure2::io::OptionSpelling spelling = azure2::io::OptionSpelling::kKey;
  for (const auto &given : config.optionSources)
    if (given.second == OptionSource::kKeyword)
      spelling = azure2::io::OptionSpelling::kKeyword;
    else if (given.second == OptionSource::kCommandLine && spelling == azure2::io::OptionSpelling::kKey)
      spelling = azure2::io::OptionSpelling::kFlag;
  auto givenIn = [&config, spelling](const char *name) {
    const OptionSource source = OptionSourceOf(config, name);
    return spelling == azure2::io::OptionSpelling::kKey ? source == OptionSource::kProject
                                                        : source == OptionSource::kCommandLine || source == OptionSource::kKeyword;
  };
  // Every contradicting pair is refused, whatever the source of its two values (a pair of defaults never contradicts).
  for (const azure2::io::ConfigOptionContradiction &c : azure2::io::ConfigOptionContradictions()) {
    if (OptionIsOn(config, c.a) != c.aValue || OptionIsOn(config, c.b) != c.bValue) continue;
    const char *key = OptionSourceOf(config, c.a) == OptionSource::kProject   ? c.a
                      : OptionSourceOf(config, c.b) == OptionSource::kProject ? c.b
                                                                              : c.a;
    const azure2::io::ConfigEntry *e = OptionSourceOf(config, key) == OptionSource::kProject ? section.Find(key) : nullptr;
    azure2::experiment::Diagnostic d;
    d.code = azure2::experiment::codes::kConfigOptionsContradict;
    d.file = config.configfile;
    d.line = e ? tagLine + e->line : 0;
    d.section = "config";
    d.field = key;
    d.message = azure2::io::ConfigOptionContradictionText(c, OptionWords(config, c.a), OptionWords(config, c.b), spelling,
                                                          givenIn(c.a), givenIn(c.b));
    config.refusedOptions += d.ToText() + "\n";
    config.outStream << "ERROR: " << d.ToText() << std::endl;
  }
  return config.refusedOptions.empty();
}

bool OptionsOfRun(const std::string &projectText, const std::vector<std::string> &args, Config &config) {
  for (const std::string &arg : args) ApplyOptionFlag(config, arg);
  const azure2::io::Document doc = azure2::io::Document::Parse(projectText);
  std::size_t open = 0;
  return ApplyProjectOptions(config, azure2::io::ReadConfigSection(azure2::io::ConfigLinesOf(doc, open)), (int)open + 1);
}

std::string WithOptions(const std::string &described, bool json, const Config &config) {
  std::string out = described;
  const std::size_t close = out.rfind('}');
  if (json && close != std::string::npos)
    out.insert(close, ",\n \"options\": " + OptionsJson(config) + "\n");
  else if (!json && !described.empty())  // an empty text describes nothing (a record pyazr asks for that has no entry)
    out += "options of the calculation (value, source):\n" + OptionsText(config);
  return out;
}

std::string OptionsText(const Config &config) {
  std::ostringstream out;
  for (const azure2::io::ConfigKey &k : azure2::io::ConfigKeys())
    if (k.option)
      out << "  " << k.name << " = " << (OptionIsOn(config, k.name) ? 1 : 0) << " ("
          << OptionSourceName(OptionSourceOf(config, k.name)) << ")\n";
  return out.str();
}

std::string OptionsJson(const Config &config) {
  std::ostringstream out;
  out << "{";
  bool first = true;
  for (const azure2::io::ConfigKey &k : azure2::io::ConfigKeys()) {
    if (!k.option) continue;
    out << (first ? "" : ", ") << "\"" << k.name << "\": {\"value\": " << (OptionIsOn(config, k.name) ? "true" : "false")
        << ", \"source\": \"" << OptionSourceName(OptionSourceOf(config, k.name)) << "\"}";
    first = false;
  }
  out << "}";
  return out.str();
}
