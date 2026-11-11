// The set-up input of a project, from either producer of specs (docs/upgrade/project_format3.md step 3b): a format-3
// project was read once by Config::ReadConfigFile (records -> Validate -> <levels> tables -> Build, Config::format3);
// an old (format 1/2) project is read here from its sections through convert/'s old-line adaptor
// (SpecsFromOldLines.h).  The old branch is TEMPORARY: step 3d deletes it with the adaptor.  Either way the engine's
// set-up (WorkingCopy::Fill, MakePoints, CNuc::SetAbundances) consumes the same specs.
#include "ProjectSetup.h"

#include <istream>
#include <memory>

#include "CNuc.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "azure2/convert/SpecsFromOldLines.h"
#include "azure2/experiment/ProjectReader.h"

namespace {

// The old <targetInt> section; a project without one gives a failed list (the old engine stopped when a segment
// needed it).
azure2::experiment::SpecList<azure2::experiment::TargetEffectSpec> OldTargetEffects(const Config &configure) {
  std::unique_ptr<std::istream> section = configure.SectionStream("targetInt");
  if (!section) {
    azure2::experiment::SpecList<azure2::experiment::TargetEffectSpec> none;
    none.failed = true;
    return none;
  }
  return azure2::convert::TargetEffectSpecsFromOldLines(*section);
}

azure2::experiment::DataSpecs OldDataSpecs(const Config &configure, bool &present) {
  azure2::experiment::DataSpecs specs;
  std::unique_ptr<std::istream> section = configure.SectionStream("segmentsData");
  present = section != nullptr;
  if (!section) return specs;
  specs.datasets = azure2::convert::DatasetSpecsFromOldLines(*section, configure.relativistic);
  if (!specs.datasets.failed) {
    std::unique_ptr<std::istream> settings = configure.SectionStream("parameterSettings");
    specs.priorCentreError = azure2::convert::PriorCentresFromOldRows(settings.get(), specs.datasets);
  }
  specs.targetEffects = OldTargetEffects(configure);
  return specs;
}

azure2::experiment::ExtrapolationSpecs OldExtrapolationSpecs(const Config &configure, bool &present) {
  azure2::experiment::ExtrapolationSpecs specs;
  std::unique_ptr<std::istream> section = configure.SectionStream("segmentsTest");
  present = section != nullptr;
  if (!section) return specs;
  specs.extrapolations = azure2::convert::ExtrapolationSpecsFromOldLines(*section, configure.relativistic);
  specs.targetEffects = OldTargetEffects(configure);
  return specs;
}

// The data files of a format-3 project that could not be read (Config::ReadFormat3Project builds it for a run without
// data): a data run refuses them as E_DATA_FILE and refuses the project (Config::ReportFormat3Refused, the line of a
// project refused when it is read), a run without data warns.  True when there were any.
bool ReportUnreadDataFiles(const Config &configure, bool dataRun) {
  bool any = false;
  for (azure2::experiment::Diagnostic d : configure.format3->diagnostics) {
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
  if (configure.format3) {
    if (ReportUnreadDataFiles(configure, true)) return kFillRefused;
    return this->Fill(configure, theCNuc, configure.format3->built.data);
  }
  bool present = false;
  const azure2::experiment::DataSpecs specs = OldDataSpecs(configure, present);
  return present ? this->Fill(configure, theCNuc, specs) : -1;
}

int WorkingCopy::MakePoints(const Config &configure, CNuc *theCNuc) {
  if (configure.format3) {
    ReportUnreadDataFiles(configure, false);
    return this->MakePoints(configure, theCNuc, configure.format3->built.extrapolations);
  }
  bool present = false;
  const azure2::experiment::ExtrapolationSpecs specs = OldExtrapolationSpecs(configure, present);
  return present ? this->MakePoints(configure, theCNuc, specs) : -1;
}

int CNuc::ReadAbundances(const Config &configure) {
  if (configure.format3) return this->SetAbundances(configure, configure.format3->built.abundances);
  std::unique_ptr<std::istream> section = configure.SectionStream("abundances");
  return this->SetAbundances(configure, azure2::convert::AbundanceSpecsFromOldLines(section.get()));
}

bool ProjectLevelSection(const Config &configure, azure2::io::LevelSection &out) {
  if (configure.format3) {
    if (!configure.format3->levelSection) return false;
    out = *configure.format3->levelSection;
    return true;
  }
  std::unique_ptr<std::istream> section = configure.SectionStream("levels");
  if (!section) return false;
  out = azure2::io::ReadLevelSection(*section);
  return true;
}

std::set<int> ProjectExitPairKeys(const Config &configure) {
  std::set<int> keys;
  bool present = false;
  const azure2::experiment::DataSpecs data =
      configure.format3 ? configure.format3->built.data : OldDataSpecs(configure, present);
  const azure2::experiment::ExtrapolationSpecs extrapolations =
      configure.format3 ? configure.format3->built.extrapolations : OldExtrapolationSpecs(configure, present);
  for (const auto &item : data.datasets.items)
    if (item.spec.active) keys.insert(item.spec.reaction.exitKey);
  for (const auto &item : extrapolations.extrapolations.items)
    if (item.spec.active) keys.insert(item.spec.reaction.exitKey);
  if (configure.paramMask & Config::CALCULATE_REACTION_RATE) keys.insert(configure.rateParams.exitPair);
  return keys;
}

bool ProjectReactionPairs(const Config &configure, bool data, std::vector<std::pair<int, int>> &pairs) {
  bool present = configure.format3 != nullptr;
  if (data) {
    const azure2::experiment::DataSpecs specs = configure.format3 ? configure.format3->built.data : OldDataSpecs(configure, present);
    for (const auto &item : specs.datasets.items)
      if (item.spec.active) pairs.push_back({item.spec.reaction.entranceKey, item.spec.reaction.exitKey});
  } else {
    const azure2::experiment::ExtrapolationSpecs specs =
        configure.format3 ? configure.format3->built.extrapolations : OldExtrapolationSpecs(configure, present);
    for (const auto &item : specs.extrapolations.items)
      if (item.spec.active) pairs.push_back({item.spec.reaction.entranceKey, item.spec.reaction.exitKey});
  }
  return present;
}

bool ProjectHasThmData(const Config &configure) {
  bool present = false;
  const azure2::experiment::DataSpecs specs = configure.format3 ? configure.format3->built.data : OldDataSpecs(configure, present);
  for (const auto &item : specs.datasets.items)
    if (item.spec.active && item.spec.reaction.thm) return true;
  return false;
}
