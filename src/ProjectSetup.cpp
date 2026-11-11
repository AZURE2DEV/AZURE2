// The set-up input of a project (ProjectSetup.h): the format-3 project as Config::ReadConfigFile read it once (records
// -> Validate -> <levels> tables -> Schema::Build, Config::project).  The engine's set-up (WorkingCopy::Fill,
// MakePoints, CNuc::SetAbundances) consumes its specs.
#include "ProjectSetup.h"

#include "CNuc.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "azure2/experiment/ProjectReader.h"

namespace {

// The data files of the project that could not be read (Config::ReadFormat3Project builds it for a run without data):
// a data run refuses them as E_DATA_FILE and refuses the project (Config::ReportFormat3Refused, the line of a project
// refused when it is read), a run without data warns.  True when there were any.
bool ReportUnreadDataFiles(const Config &configure, bool dataRun) {
  bool any = false;
  for (azure2::experiment::Diagnostic d : configure.project->diagnostics) {
    if (d.code != azure2::experiment::codes::kDataFileNotRead) continue;
    any = true;
    if (dataRun) {
      d.code = azure2::experiment::codes::kDataFile;
      d.severity = azure2::experiment::Severity::Error;
      d.message = d.recordKind + " " + d.record + ", field file: cannot read " + d.found;
    }
    configure.outStream << (dataRun ? "ERROR: " : "WARNING: ") << d.ToText() << std::endl;
  }
  if (any && dataRun) configure.ReportFormat3Refused();
  return any;
}

}  // namespace

int WorkingCopy::Fill(const Config &configure, CNuc *theCNuc) {
  if (ReportUnreadDataFiles(configure, true)) return kFillRefused;
  // Every refusal of the set-up prints its ERROR line before it returns: that line ends the run (WorkingCopy.h).
  return this->Fill(configure, theCNuc, configure.project->built.data) == 0 ? 0 : kFillRefused;
}

int WorkingCopy::MakePoints(const Config &configure, CNuc *theCNuc) {
  ReportUnreadDataFiles(configure, false);
  return this->MakePoints(configure, theCNuc, configure.project->built.extrapolations) == 0 ? 0 : kFillRefused;
}

int CNuc::ReadAbundances(const Config &configure) {
  return this->SetAbundances(configure, configure.project->built.abundances);
}

bool ProjectLevelSection(const Config &configure, azure2::io::LevelSection &out) {
  if (!configure.project->levelSection) return false;
  out = *configure.project->levelSection;
  return true;
}

std::set<int> ProjectExitPairKeys(const Config &configure) {
  std::set<int> keys;
  for (const auto &item : configure.project->built.data.datasets.items)
    if (item.spec.active) keys.insert(item.spec.reaction.exitKey);
  for (const auto &item : configure.project->built.extrapolations.extrapolations.items)
    if (item.spec.active) keys.insert(item.spec.reaction.exitKey);
  if (configure.paramMask & Config::CALCULATE_REACTION_RATE) keys.insert(configure.rateParams.exitPair);
  return keys;
}

std::vector<std::pair<int, int>> ProjectReactionPairs(const Config &configure, bool data) {
  std::vector<std::pair<int, int>> pairs;
  if (data) {
    for (const auto &item : configure.project->built.data.datasets.items)
      if (item.spec.active) pairs.push_back({item.spec.reaction.entranceKey, item.spec.reaction.exitKey});
  } else {
    for (const auto &item : configure.project->built.extrapolations.extrapolations.items)
      if (item.spec.active) pairs.push_back({item.spec.reaction.entranceKey, item.spec.reaction.exitKey});
  }
  return pairs;
}

bool ProjectHasThmData(const Config &configure) {
  for (const auto &item : configure.project->built.data.datasets.items)
    if (item.spec.active && item.spec.reaction.thm) return true;
  return false;
}
