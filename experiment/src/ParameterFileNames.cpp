// A parameter file's names against the run's parameters (ParameterFileNames.h).
#include "azure2/experiment/ParameterFileNames.h"

#include <set>

namespace azure2 {
namespace experiment {

namespace {

// "a, b, c" (at most `shown` names, then "and <n> more").
std::string NameList(const std::vector<std::string> &names, std::size_t shown = 12) {
  std::string out;
  for (std::size_t i = 0; i < names.size() && i < shown; i++) out += (i ? ", " : "") + names[i];
  if (names.size() > shown) out += " and " + std::to_string(names.size() - shown) + " more";
  return out;
}

}  // namespace

std::vector<Diagnostic> MatchParameterFileNames(const std::string &file, const ParameterFileNames &names) {
  const std::set<std::string> run(names.ofRun.begin(), names.ofRun.end());
  const std::set<std::string> elsewhere(names.ofOtherRuns.begin(), names.ofOtherRuns.end());
  const std::set<std::string> optional(names.optionalOfRun.begin(), names.optionalOfRun.end());
  const std::set<std::string> inFile(names.inFile.begin(), names.inFile.end());
  std::vector<std::string> unknown, missing;
  int matched = 0;
  for (const std::string &n : names.inFile) {
    if (run.count(n))
      matched++;
    else if (!elsewhere.count(n))
      unknown.push_back(n);
  }
  for (const std::string &n : names.ofRun)
    if (!inFile.count(n) && !optional.count(n)) missing.push_back(n);
  Diagnostic d;
  d.file = file;
  if (matched == 0) {
    d.code = codes::kParameterFileNoMatch;
    d.severity = Severity::Error;
    d.message = names.inFile.empty()
                    ? "the parameter file names no parameter"
                    : "no name of the parameter file is a parameter of this project (" + NameList(names.inFile) +
                          "); nothing of it would be applied";
    return {d};
  }
  if (unknown.empty() && missing.empty()) return {};
  d.code = codes::kParameterFileNames;
  d.severity = Severity::Warning;
  d.message = std::to_string(matched) + " of its " + std::to_string(names.inFile.size()) + " names applied";
  if (!unknown.empty()) d.message += "; names this project has no parameter for (skipped): " + NameList(unknown);
  if (!missing.empty()) d.message += "; parameters it does not name (the project's values kept): " + NameList(missing);
  return {d};
}

}  // namespace experiment
}  // namespace azure2
