// TEMPORARY (step 3d of project_format3.md deletes this file): the set-up of an old (format 1/2) project.  Each
// function reads its sections through convert/'s old-line adaptor (SpecsFromOldLines.h) and hands the specs to the
// set-up that a format-3 project reaches through Schema::Build.
#include <istream>
#include <memory>

#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "azure2/convert/SpecsFromOldLines.h"

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

}  // namespace

int EData::Fill(const Config &configure, CNuc *theCNuc) {
  std::unique_ptr<std::istream> section = configure.SectionStream("segmentsData");
  if (!section) return -1;
  azure2::experiment::DataSpecs specs;
  specs.datasets = azure2::convert::DatasetSpecsFromOldLines(*section, configure.relativistic);
  if (!specs.datasets.failed) {
    std::unique_ptr<std::istream> settings = configure.SectionStream("parameterSettings");
    specs.priorCentreError = azure2::convert::PriorCentresFromOldRows(settings.get(), specs.datasets);
  }
  specs.targetEffects = OldTargetEffects(configure);
  return this->Fill(configure, theCNuc, specs);
}

int EData::MakePoints(const Config &configure, CNuc *theCNuc) {
  std::unique_ptr<std::istream> section = configure.SectionStream("segmentsTest");
  if (!section) return -1;
  azure2::experiment::ExtrapolationSpecs specs;
  specs.extrapolations = azure2::convert::ExtrapolationSpecsFromOldLines(*section, configure.relativistic);
  specs.targetEffects = OldTargetEffects(configure);
  return this->MakePoints(configure, theCNuc, specs);
}

int CNuc::ReadAbundances(const Config &configure) {
  std::unique_ptr<std::istream> section = configure.SectionStream("abundances");
  return this->SetAbundances(configure, azure2::convert::AbundanceSpecsFromOldLines(section.get()));
}
