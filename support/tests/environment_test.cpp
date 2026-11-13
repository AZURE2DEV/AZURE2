// support_environment: the table of environment switches and the one reader (include/azure2/support/Environment.h).
// ctest runs it with AZURE_PROFILE=1 and AZURE_GRAD_CHECK_H=2.5e-7 set, so that the process's record is seen read.
#include "azure2/support/Environment.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <map>
#include <set>
#include <string>

using namespace azure2::support;

namespace {
int failures = 0;
void check(const std::string &name, bool ok) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", name.c_str());
  if (!ok) failures++;
}
}  // namespace

int main() {
  // The table: one row per switch, unique names, every name an AZURE variable, every column filled.
  std::set<std::string> names;
  bool filled = true;
  for (const SwitchRow &row : SwitchTable()) {
    names.insert(row.name);
    filled = filled && row.unset && *row.unset && row.what && *row.what && row.documented && *row.documented;
  }
  check("one row per switch, unique names", (int)names.size() == kSwitchCount);
  check("every row has its value when unset, what it does and its page", filled);
  check("RowOf finds the row", std::string(RowOf(Switch::kProfile).name) == "AZURE_PROFILE" &&
                                   std::string(RowOf(Switch::kTestFailExit).name) == "AZURE_TEST_FAIL_EXIT");
  check("kinds by name", std::string(KindName(SwitchKind::kNumbers)) == "numbers" &&
                             std::string(KindName(SwitchKind::kTest)) == "test");

  // A record made from given values: set at all, set and on, text, numbers, the home directory.
  const std::map<std::string, std::string> given = {{"AZURE_PROFILE", "0"},
                                                    {"AZURE_LEVEL_SPACE", ""},
                                                    {"AZURE_GRAD_CHECK", "exit"},
                                                    {"AZURE_GRAD_CHECK_H", "1e-5"},
                                                    {"AZURE2_EC_TABLE_NODES", "80"},
                                                    {"AZURE_DUMP_NONFINITE", ""},
                                                    {"HOME", ""},
                                                    {"USERPROFILE", "C:/Users/x"}};
  const Environment e = Environment::Read([&](const char *name) -> const char * {
    auto it = given.find(name);
    return it == given.end() ? nullptr : it->second.c_str();
  });
  check("'0' is set but not on", e.Set(Switch::kProfile) && !e.On(Switch::kProfile));
  check("'' is set but not on", e.Set(Switch::kLevelSpace) && !e.On(Switch::kLevelSpace));
  check("a word is on, with its text", e.On(Switch::kGradCheck) && e.Text(Switch::kGradCheck) == "exit");
  check("an unset switch: not set, empty text, the fallback number",
        !e.Set(Switch::kPoolCheck) && e.Text(Switch::kPoolCheck).empty() && e.Number(Switch::kGradCheckRows, 60) == 60);
  check("numbers by atof and atoi", e.Number(Switch::kGradCheckStep, 1e-6) == 1e-5 &&
                                        e.Integer(Switch::kEcTableNodes, 0) == 80);
  check("a set but empty number is 0, as atof(\"\")", e.Number(Switch::kDumpNonFinite, 1e10) == 0.0);
  check("HOME empty: USERPROFILE", e.Home() == "C:/Users/x");
  check("the set switches in the table's order",
        e.SetSwitches() == std::vector<Switch>{Switch::kEcTableNodes, Switch::kLevelSpace, Switch::kGradCheck,
                                               Switch::kGradCheckStep, Switch::kProfile, Switch::kDumpNonFinite});
  check("the text of the set switches for --describe: the other kinds (numbers have their own block)",
        SetSwitchesText(e).rfind("Environment switches set:\n  AZURE_LEVEL_SPACE", 0) == std::string::npos &&
            SetSwitchesText(e).rfind("Environment switches set:\n  AZURE_GRAD_CHECK=exit (prints)\n", 0) == 0);
  check("the block of the switches that move numbers: every one, value and source",
        NumberSwitchesText(e).find("  --ec-table-nodes=N (AZURE2_EC_TABLE_NODES) = 80 (environment)\n") !=
                std::string::npos &&
            NumberSwitchesText(e).find("  --refine-quadrature=K (AZURE_REFINE_QUADRATURE) = 1 (default)\n") !=
                std::string::npos);
  const Environment none = Environment::Read([](const char *) -> const char * { return nullptr; });
  check("nothing set: no switch, home '.', no text", none.SetSwitches().empty() && none.Home() == "." &&
                                                         SetSwitchesText(none).empty());

  // The options of the switches that move numbers (docs/upgrade/evaluation.md section 9): every row of kind numbers
  // has one, no other row; the parser takes a flag whole and a value row with its value; a value row without one is an
  // error; an argument that is no switch option is not one.
  bool numbersHaveOptions = true;
  for (const SwitchRow &row : SwitchTable())
    numbersHaveOptions = numbersHaveOptions && ((row.kind == SwitchKind::kNumbers) == (*row.option != '\0'));
  check("an option for each switch of kind numbers, none for the others", numbersHaveOptions);
  std::string value;
  const SwitchRow *row = ParseSwitchOption("--level-space", value);
  check("--level-space: the flag, value 1", row && row->id == Switch::kLevelSpace && value == "1");
  row = ParseSwitchOption("--refine-quadrature=4", value);
  check("--refine-quadrature=4: the value 4", row && row->id == Switch::kRefineQuadrature && value == "4");
  check("--level-spacex, --use-brune, --ec-table: no switch option",
        !ParseSwitchOption("--level-spacex", value) && !ParseSwitchOption("--use-brune", value) &&
            !ParseSwitchOption("--ec-table", value));
  check("--srim-data-dir and --coulomb-table bare: no switch option (an unknown option, a flag of another row)",
        !ParseSwitchOption("--srim-data-dir", value) && !ParseSwitchOption("--coulomb-table", value));
  for (const char *bare : {"--srim-data-dir=", "--coulomb-table="}) {
    bool refused = false;
    try {
      ParseSwitchOption(bare, value);
    } catch (const std::invalid_argument &) {
      refused = true;
    }
    check(std::string(bare) + ": refused, a value is needed", refused);
  }
  // Precedence: the command line over the options file over the variable over the default.
  const std::string optionsFile = "environment_test_options.txt";
  {
    std::ofstream f(optionsFile);
    f << "--refine-quadrature=3\n  --level-space\n--ec-table-nodes = 70\n--no-brune\nnot an option\n";
  }
  const std::map<std::string, std::string> variables = {
      {"AZURE_REFINE_QUADRATURE", "2"}, {"AZURE2_EC_TABLE_NODES", "60"}, {"AZURE_GRAD_DENSE", "1"},
      {"AZURE_OPTIONS_FILE", optionsFile}, {"OMP_NUM_THREADS", "3"}};
  auto lookup = [&](const char *name) -> const char * {
    auto it = variables.find(name);
    return it == variables.end() ? nullptr : it->second.c_str();
  };
  const Environment layered = Environment::Read(lookup, {"--refine-quadrature=5", "project.azr", "--use-brune"});
  check("the option over the file and the variable",
        layered.Text(Switch::kRefineQuadrature) == "5" && layered.Source(Switch::kRefineQuadrature) == SwitchSource::kOption);
  check("the file over the variable (its spaces removed)", layered.Text(Switch::kEcTableNodes) == "70" &&
                                                              layered.Source(Switch::kEcTableNodes) ==
                                                                  SwitchSource::kOptionsFile &&
                                                              layered.On(Switch::kLevelSpace));
  check("the variable over the default",
        layered.On(Switch::kGradDense) && layered.Source(Switch::kGradDense) == SwitchSource::kEnvironment &&
            layered.Source(Switch::kGradAdjoint) == SwitchSource::kDefault && !layered.Set(Switch::kGradAdjoint));
  check("the file's lines starting with --, in order, for the program's other options",
        layered.OptionsFileLines() == std::vector<std::string>{"--refine-quadrature=3", "--level-space",
                                                               "--ec-table-nodes=70", "--no-brune"});
  check("OMP_NUM_THREADS kept", layered.OmpNumThreads() == "3");
  const Environment byOption = Environment::Read(lookup, {"--options-file=no_such_file.txt"});
  check("--options-file= over AZURE_OPTIONS_FILE: the named file (unreadable here, said)",
        byOption.OptionsFileUnreadable() && byOption.Source(Switch::kOptionsFile) == SwitchSource::kOption &&
            byOption.Text(Switch::kRefineQuadrature) == "2");
  // The checked value: the same range and refusal whatever gave it, named after its source.
  auto refusal = [&](const std::vector<std::string> &arguments, const char *variable) {
    const Environment r = Environment::Read(
        [&](const char *name) -> const char * {
          return std::string(name) == "AZURE_REFINE_QUADRATURE" ? variable : nullptr;
        },
        arguments);
    try {
      CheckSwitchValues(r);
    } catch (const std::invalid_argument &x) {
      return std::string(x.what());
    }
    return std::string();
  };
  check("AZURE_REFINE_QUADRATURE=65 refused under its name",
        refusal({}, "65") == "AZURE_REFINE_QUADRATURE needs a whole number from 1 to 64, not '65'");
  check("--refine-quadrature=65 refused alike, under the option's name",
        refusal({"--refine-quadrature=65"}, nullptr) == "--refine-quadrature needs a whole number from 1 to 64, not '65'");
  check("--refine-quadrature=x and 0 refused", !refusal({"--refine-quadrature=x"}, nullptr).empty() &&
                                                   !refusal({"--refine-quadrature=0"}, nullptr).empty());
  check("64 and an unset one accepted", refusal({"--refine-quadrature=64"}, nullptr).empty() && refusal({}, nullptr).empty());
  std::remove(optionsFile.c_str());

  // The process's record: read once, so a variable set afterwards changes nothing.
  const Environment &p = ProcessEnvironment();
  check("the process's record sees AZURE_PROFILE=1 (ctest)", p.On(Switch::kProfile));
  check("... and AZURE_GRAD_CHECK_H=2.5e-7", p.Number(Switch::kGradCheckStep, 1e-6) == 2.5e-7);
#ifndef _WIN32
  setenv("AZURE_POOL_CHECK", "1", 1);
  check("a variable set after the first read is not seen", !ProcessEnvironment().Set(Switch::kPoolCheck) &&
                                                                &ProcessEnvironment() == &p);
#endif
  std::printf("%d failure(s)\n", failures);
  return failures ? 1 : 0;
}
