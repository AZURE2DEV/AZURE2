// The set-up input of a project (ProjectSetup.h): the format-3 project as Config::ReadConfigFile read it once (records
// -> Validate -> <levels> tables -> Schema::Build, Config::project).  The engine's set-up (WorkingCopy::Fill,
// MakePoints, WorkingCopy::SetAbundanceSamples) consumes its specs.
#include "ProjectSetup.h"

#include "CompoundSystemSetUp.h"

#include "CompoundSystem.h"
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

int WorkingCopy::SetAbundanceSamples(std::ostream &messages, const Config &configure, const CompoundSystem &compound) {
  std::vector<std::pair<int, int> > nucleiZA;
  for (int k = 1; k <= compound.NumNuclei(); k++) nucleiZA.push_back({compound.GetNucleus(k)->Z, compound.GetNucleus(k)->A});
  const int read = Draft().abundances_.Set(messages, configure.project->built.abundances, nucleiZA);
  abundances_ = Data().Abundances().StartValues();
  return read;
}

int WorkingCopy::Fill(const Config &configure, const CompoundSystem *compound, const LevelInputs &inputs) {
  if (SetAbundanceSamples(configure.outStream, configure, *compound) == -1) return -1;
  if (ReportUnreadDataFiles(configure, true)) return kFillRefused;
  // Every refusal of the set-up prints its ERROR line before it returns: that line ends the run (WorkingCopy.h).
  return this->Fill(configure, compound, inputs, configure.project->built.data) == 0 ? 0 : kFillRefused;
}

int WorkingCopy::MakePoints(const Config &configure, const CompoundSystem *compound, const LevelInputs &inputs) {
  if (SetAbundanceSamples(configure.outStream, configure, *compound) == -1) return -1;
  ReportUnreadDataFiles(configure, false);
  return this->MakePoints(configure, compound, inputs, configure.project->built.extrapolations) == 0 ? 0 : kFillRefused;
}



bool ProjectLevelSection(const Config &configure, azure2::io::LevelSection &out) {
  if (!configure.project->levelSection) return false;
  out = *configure.project->levelSection;
  return true;
}

namespace {
// The pairs a reaction exits into: its leading term's and those of its components (the other terms of a SUM or a
// RATIO).  A capture pair that is only a component needs its external capture as much as a leading one: until
// 2026-11-11 only the leading terms were listed, and "elastic + 1000 x capture" ran without the external capture of
// its second term unless another data set had the capture as its leading term (tests/component_mapping (d)).
template <class F> void ForEachTermPair(const azure2::experiment::ReactionSpec &reaction, F &&f) {
  f(reaction.entranceKey, reaction.exitKey);
  for (const azure2::experiment::ComponentSpec &component : reaction.components) f(component.entranceKey, component.exitKey);
}
}  // namespace

std::set<int> ProjectExitPairKeys(const Config &configure) {
  std::set<int> keys;
  auto add = [&keys](int, int exitKey) { keys.insert(exitKey); };
  for (const auto &item : configure.project->built.data.datasets.items)
    if (item.spec.active) ForEachTermPair(item.spec.reaction, add);
  for (const auto &item : configure.project->built.extrapolations.extrapolations.items)
    if (item.spec.active) ForEachTermPair(item.spec.reaction, add);
  if (configure.paramMask & Config::CALCULATE_REACTION_RATE) keys.insert(configure.rateParams.exitPair);
  return keys;
}

std::vector<std::pair<int, int>> ProjectReactionPairs(const Config &configure, bool data) {
  std::vector<std::pair<int, int>> pairs;
  auto add = [&pairs](int entranceKey, int exitKey) { pairs.push_back({entranceKey, exitKey}); };
  if (data) {
    for (const auto &item : configure.project->built.data.datasets.items)
      if (item.spec.active) ForEachTermPair(item.spec.reaction, add);
  } else {
    for (const auto &item : configure.project->built.extrapolations.extrapolations.items)
      if (item.spec.active) ForEachTermPair(item.spec.reaction, add);
  }
  return pairs;
}

namespace {
// Every entrance pair key of the active items of `items` that names a pair of `compound`, requested as an entrance.
template <class Items> void RequestEntrancesOf(CompoundSystem &compound, const Items &items) {
  auto add = [&compound](int entranceKey, int) {
    if (compound.IsPairKey(entranceKey)) compound.RequestEntrance(compound.GetPairNumFromKey(entranceKey));
  };
  for (const auto &item : items)
    if (item.spec.active) ForEachTermPair(item.spec.reaction, add);
}
}  // namespace

void RequestDataEntrances(CompoundSystem &compound, const Config &configure) {
  RequestEntrancesOf(compound, configure.project->built.data.datasets.items);
}

void RequestExtrapolationEntrances(CompoundSystem &compound, const Config &configure) {
  RequestEntrancesOf(compound, configure.project->built.extrapolations.extrapolations.items);
}

namespace {
// Every entrance pair key of the active items of `items` that names a pair of `compound`, asked for in `requests`.
template <class Items>
void AddEntrancesOf(const CompoundSystem &compound, const Items &items, CompoundSetUpRequests &requests) {
  auto add = [&compound, &requests](int entranceKey, int) {
    if (compound.IsPairKey(entranceKey)) requests.AddEntrance(compound.GetPairNumFromKey(entranceKey));
  };
  for (const auto &item : items)
    if (item.spec.active) ForEachTermPair(item.spec.reaction, add);
}
}  // namespace

void RequestDataEntrances(const CompoundSystem &compound, const Config &configure, CompoundSetUpRequests &requests) {
  AddEntrancesOf(compound, configure.project->built.data.datasets.items, requests);
}

void RequestExtrapolationEntrances(const CompoundSystem &compound, const Config &configure,
                                   CompoundSetUpRequests &requests) {
  AddEntrancesOf(compound, configure.project->built.extrapolations.extrapolations.items, requests);
}

bool ProjectHasThmData(const Config &configure) {
  for (const auto &item : configure.project->built.data.datasets.items)
    if (item.spec.active && item.spec.reaction.thm) return true;
  return false;
}
