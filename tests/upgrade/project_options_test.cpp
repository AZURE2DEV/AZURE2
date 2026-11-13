// The options of the calculation (include/ProjectOptions.h; docs/upgrade/project_format3.md section 3.12): every option
// row of io::ConfigKeys() has its place in Config and the default of Config::Reset; flags, keywords and the project's
// keys combine as the rule says (given > key > default); a contradiction with a key is refused, naming both options.
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Config.h"
#include "ProjectOptions.h"
#include "azure2/io/ConfigSection.h"

namespace {

int failures = 0;
long checks = 0;
void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    failures++;
    std::printf("  FAIL  %s\n", what.c_str());
  }
}

// A <config> section with the eleven positional entries and these key lines.
azure2::io::ConfigSection Section(const std::vector<std::string> &keys) {
  std::vector<std::string> lines = {"<config>", "true", "output/", "checks/"};
  for (int i = 0; i < 8; i++) lines.push_back("none");
  for (const std::string &k : keys) lines.push_back(k);
  lines.push_back("format=3");
  lines.push_back("</config>");
  return azure2::io::ReadConfigSection(lines);
}

}  // namespace

int main() {
  std::ostringstream sink;
  // Every option row: a place in Config, the default of Config::Reset, source "default".
  {
    Config fresh(sink);
    int rows = 0;
    for (const azure2::io::ConfigKey &k : azure2::io::ConfigKeys()) {
      if (!k.option) continue;
      rows++;
      bool on = false, stored = true;
      try {
        on = OptionIsOn(fresh, k.name);
      } catch (const std::invalid_argument &) {
        stored = false;
      }
      check(std::string(k.name) + ": has a place in Config (src/ProjectOptions.cpp)", stored);
      check(std::string(k.name) + ": its row's default is Config::Reset's", stored && on == k.defaultOn);
      check(std::string(k.name) + ": source default", OptionSourceOf(fresh, k.name) == OptionSource::kDefault);
      // set and read back through the storage, both values
      for (const bool value : {true, false}) {
        Config c(sink);
        GiveOption(c, k.name, value, OptionSource::kKeyword);
        check(std::string(k.name) + "=" + (value ? "1" : "0") + ": read back", OptionIsOn(c, k.name) == value);
      }
    }
    check("eleven options", rows == 11);
    bool refused = false;
    try {
      OptionIsOn(fresh, "relativistic");
    } catch (const std::invalid_argument &) {
      refused = true;
    }
    check("a key that is not an option has no storage", refused);
  }
  // The flags: each sets its option for this run, marked "command line"; the implications as the flags had them.
  for (const azure2::io::ConfigKey &k : azure2::io::ConfigKeys()) {
    if (!k.option) continue;
    for (const bool on : {true, false}) {
      const std::string flag = on ? k.flag : k.negation;
      if (flag.empty()) continue;
      Config c(sink);
      check(flag + ": a flag of the table", ApplyOptionFlag(c, flag));
      check(flag + ": " + k.name + (on ? " on" : " off") + ", from the command line",
            OptionIsOn(c, k.name) == on && OptionSourceOf(c, k.name) == OptionSource::kCommandLine);
    }
  }
  {
    Config c(sink);
    check("--use-api is no option flag", !ApplyOptionFlag(c, "--use-api") && !ApplyOptionFlag(c, "--coulomb-table=x"));
  }
  auto both = [&sink](const std::vector<std::string> &flags, bool brune, bool park) {
    Config c(sink);
    for (const std::string &f : flags) ApplyOptionFlag(c, f);
    return OptionIsOn(c, "use_brune") == brune && OptionIsOn(c, "use_park") == park &&
           (c.paramMask & Config::USE_BRUNE_FORMALISM) == (brune ? Config::USE_BRUNE_FORMALISM : 0u);
  };
  check("--use-park: Park and Brune", both({"--use-park"}, true, true));
  check("--no-brune: neither", both({"--no-brune"}, false, false));
  check("--no-brune --use-park: both (the order decides, as before)", both({"--no-brune", "--use-park"}, true, true));
  check("--use-park --no-brune: neither", both({"--use-park", "--no-brune"}, false, false));
  {
    Config c(sink);
    ApplyOptionFlag(c, "--use-park");
    check("--use-park gives use_brune too (command line)", OptionSourceOf(c, "use_brune") == OptionSource::kCommandLine);
  }
  // The project's keys: a key over the default, a flag or keyword over the key.
  {
    Config c(sink);
    check("keys applied", ApplyProjectOptions(c, Section({"use_brune=0", "transform=false", "use_coulomb_table=0"}), 1));
    check("use_brune=0 from the project", !OptionIsOn(c, "use_brune") && OptionSourceOf(c, "use_brune") == OptionSource::kProject);
    check("transform=false from the project", !OptionIsOn(c, "transform") && OptionSourceOf(c, "transform") == OptionSource::kProject);
    check("use_coulomb_table=0 from the project", !c.useCoulombTable && OptionSourceOf(c, "use_coulomb_table") == OptionSource::kProject);
    check("an option without a key keeps its default", OptionIsOn(c, "ignore_externals") &&
                                                         OptionSourceOf(c, "ignore_externals") == OptionSource::kDefault);
  }
  {
    Config c(sink);
    ApplyOptionFlag(c, "--use-brune");
    GiveOption(c, "transform", true, OptionSource::kKeyword);
    check("flag and keyword over the key", ApplyProjectOptions(c, Section({"use_brune=0", "transform=0"}), 1) &&
                                               OptionIsOn(c, "use_brune") && OptionIsOn(c, "transform") &&
                                               OptionSourceOf(c, "use_brune") == OptionSource::kCommandLine &&
                                               OptionSourceOf(c, "transform") == OptionSource::kKeyword);
  }
  {
    Config c(sink), fresh(sink);
    check("no keys: nothing changes", ApplyProjectOptions(c, Section({"relativistic=1"}), 1) && c.paramMask == fresh.paramMask &&
                                          c.useCoulombTable == fresh.useCoulombTable && c.useECTable == fresh.useECTable &&
                                          c.ecConvectionCurrent == fresh.ecConvectionCurrent &&
                                          c.replaceSubPointGrids == fresh.replaceSubPointGrids && c.optionSources.empty());
  }
  // Contradictions with a key: refused, both named with their sources; without a key: as before (no check).
  {
    std::ostringstream out;
    Config c(out);
    c.configfile = "p.azr";
    ApplyOptionFlag(c, "--no-brune");
    const bool ok = ApplyProjectOptions(c, Section({"use_park=1"}), 1);
    // --no-brune also gives use_park off: the key use_park=1 is over-ruled, no contradiction
    check("--no-brune over use_park=1: Lane-Thomas, no refusal", ok && !OptionIsOn(c, "use_park"));
  }
  {
    std::ostringstream out;
    Config c(out);
    c.configfile = "p.azr";
    GiveOption(c, "use_brune", false, OptionSource::kKeyword);
    ApplyOptionFlag(c, "--use-rmc");
    const bool ok = ApplyProjectOptions(c, Section({"use_brune=1"}), 1);
    check("use_brune given off, the key use_brune=1 over-ruled: --use-rmc valid", ok);
  }
  {
    std::ostringstream out;
    Config c(out);
    c.configfile = "p.azr";
    ApplyOptionFlag(c, "--use-rmc");
    const bool ok = ApplyProjectOptions(c, Section({"use_brune=1"}), 1);
    check("--use-rmc with the key use_brune=1: refused, E_CONFIG_OPTIONS_CONTRADICT naming both and their sources",
          !ok && out.str().find("E_CONFIG_OPTIONS_CONTRADICT") != std::string::npos &&
              out.str().find("use_rmc=1 (command line)") != std::string::npos &&
              out.str().find("use_brune=1 (this file)") != std::string::npos && out.str().find("p.azr:13") != std::string::npos &&
              c.refusedOptions == out.str().substr(7));
  }
  {
    std::ostringstream out;
    Config c(out);
    ApplyOptionFlag(c, "--no-ec-convection-current");
    check("--no-ec-convection-current with the key use_long_wavelength=1: refused",
          !ApplyProjectOptions(c, Section({"use_long_wavelength=1"}), 1));
  }
  {
    std::ostringstream out;
    Config c(out);
    ApplyOptionFlag(c, "--use-rmc");
    check("--use-rmc without a key: no refusal here (the program's rule drops Brune later, as before)",
          ApplyProjectOptions(c, Section({}), 1) && out.str().empty());
  }
  {
    Config c(sink);
    ApplyOptionFlag(c, "--no-transform");
    const std::string text = OptionsText(c), json = OptionsJson(c);
    check("OptionsText: one line per option with its source",
          text.find("  transform = 0 (command line)\n") != std::string::npos &&
              text.find("  use_brune = 1 (default)\n") != std::string::npos);
    check("OptionsJson", json.find("\"transform\": {\"value\": false, \"source\": \"command line\"}") != std::string::npos);
  }
  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
