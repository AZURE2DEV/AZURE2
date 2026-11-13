// The options of the calculation (include/ProjectOptions.h; docs/upgrade/project_format3.md section 3.12): every option
// row of io::ConfigKeys() has its place in Config and the default of Config::Reset; flags, keywords and the project's
// keys combine as the rule says (given > key > default); a contradiction with a key is refused, naming both options.
#include <cstdio>
#include <map>
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
  // Two stated values stay as stated, whatever their order (an implication fills only what no flag stated): the pair
  // is refused below, not decided by the order.
  check("--no-brune --use-park: Park on, Brune off, as stated", both({"--no-brune", "--use-park"}, false, true));
  check("--use-park --no-brune: the same", both({"--use-park", "--no-brune"}, false, true));
  check("--no-brune --use-brune: the last stated wins", both({"--no-brune", "--use-brune"}, true, false));
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
  // Contradictions without a key: refused as well (one rule), the fix in the front end's words.
  auto refusal = [](const std::vector<std::string> &flags, const std::map<std::string, bool> &keywords) {
    std::ostringstream out;
    Config c(out);
    for (const auto &k : keywords) GiveOption(c, k.first, k.second, OptionSource::kKeyword);
    for (const std::string &f : flags) ApplyOptionFlag(c, f);
    return ApplyProjectOptions(c, Section({}), 1) ? std::string() : out.str();
  };
  const char *kRmc = "E_CONFIG_OPTIONS_CONTRADICT: use_rmc=1 (command line) with use_brune=1 (the default) cannot "
                     "hold together: the Reich-Moore capture is computed with constant boundary conditions only (Brune is "
                     "on by default); allowed: use_rmc=1 with use_brune=0, or use_rmc=0; fix: add --no-brune (or leave out "
                     "--use-rmc)\n";
  const std::string rmc = refusal({"--use-rmc"}, {}), rmcWant = kRmc;
  check("--use-rmc alone: refused, what, allowed and fix (--no-brune)",
        rmc.size() > rmcWant.size() && rmc.compare(rmc.size() - rmcWant.size(), std::string::npos, rmcWant) == 0);
  check("--use-rmc --no-brune: valid", refusal({"--use-rmc", "--no-brune"}, {}).empty());
  const std::string rmcKeyword = refusal({}, {{"use_rmc", true}});
  check("use_rmc=True alone (keyword): refused, fix use_brune=False",
        rmcKeyword.find("use_rmc=1 (keyword) with use_brune=1 (the default)") != std::string::npos &&
            rmcKeyword.find("; fix: pass use_brune=False (or use_rmc=False in place of use_rmc=True)") != std::string::npos);
  check("use_rmc=True, use_brune=False: valid", refusal({}, {{"use_rmc", true}, {"use_brune", false}}).empty());
  for (const std::vector<std::string> &order : {std::vector<std::string>{"--use-park", "--no-brune"},
                                                std::vector<std::string>{"--no-brune", "--use-park"}}) {
    const std::string text = refusal(order, {});
    check(order[0] + " " + order[1] + ": refused, either order",
          text.find("use_park=1 (command line) with use_brune=0 (command line)") != std::string::npos &&
              text.find("; fix: --use-brune in place of --no-brune (or leave out --use-park)") != std::string::npos);
  }
  check("use_park=True, use_brune=False (keywords): refused",
        refusal({}, {{"use_park", true}, {"use_brune", false}}).find("E_CONFIG_OPTIONS_CONTRADICT") != std::string::npos);
  check("--use-park alone (Brune implied on): valid", refusal({"--use-park"}, {}).empty());
  check("--no-brune alone (Park implied off): valid", refusal({"--no-brune"}, {}).empty());
  check("use_park=True alone, use_brune=False alone (keywords): valid",
        refusal({}, {{"use_park", true}}).empty() && refusal({}, {{"use_brune", false}}).empty());
  check("--no-ec-convection-current alone: refused, fix --no-long-wavelength",
        refusal({"--no-ec-convection-current"}, {}).find("; fix: add --no-long-wavelength") != std::string::npos);
  check("--no-long-wavelength --no-ec-convection-current: valid",
        refusal({"--no-long-wavelength", "--no-ec-convection-current"}, {}).empty());
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
