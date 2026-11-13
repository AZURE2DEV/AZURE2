// support_environment: the table of environment switches and the one reader (include/azure2/support/Environment.h).
// ctest runs it with AZURE_PROFILE=1 and AZURE_GRAD_CHECK_H=2.5e-7 set, so that the process's record is seen read.
#include "azure2/support/Environment.h"

#include <cstdio>
#include <cstdlib>
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
  const Environment none = Environment::Read([](const char *) -> const char * { return nullptr; });
  check("nothing set: no switch, home '.'", none.SetSwitches().empty() && none.Home() == ".");

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
