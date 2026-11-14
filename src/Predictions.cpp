// The expression trees of the segments (Predictions.h; docs/upgrade/data_model.md step S10).
#include "Predictions.h"

#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include "CompoundSystem.h"
#include "CompoundSystemSetUp.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "PPair.h"
#include "azure2/experiment/ProjectReader.h"

namespace azure2 {

namespace {

int NucleusOf(const CompoundSystem &compound, int entranceKey) {
  if (!compound.IsPairKey(entranceKey)) return 0;
  return compound.GetPair(compound.GetPairNumFromKey(entranceKey))->GetNucleus();
}

/// The names of a format-3 project's records of one kind ("dataset", "extrapolation") in file order: the record
/// with key k (Schema::Build numbers them 1, 2, ... in that order, inactive ones counted) is names[k - 1].  Empty for
/// an old project, whose data sets have no names.
std::vector<std::string> RecordNames(const Config &configure, const std::string &kind) {
  std::vector<std::string> names;
  for (const experiment::ValidatedRecord &r : configure.project->validated.records)
    if (r.kind == kind) names.push_back(r.name);
  return names;
}

/// The trees of one WorkingCopy as lines "<what> <label>: <tree>", the label the data set's name where the project gives
/// one (format 3), else its key.
void AppendLines(WorkingCopy &data, const char *what, const std::vector<std::string> &names, std::vector<std::string> &lines) {
  for (int s = 1; s <= data.Data().NumSegments(); s++) {
    const Segment &segment = *data.Data().GetSegment(s);
    const experiment::Prediction *p = segment.GetPrediction();
    const int key = segment.GetSegmentKey();
    const bool named = key >= 1 && key <= (int)names.size() && !names[key - 1].empty();
    lines.push_back(std::string(what) + " " + (named ? names[key - 1] : std::to_string(key)) + ": " +
                    (p ? experiment::PredictionText(*p) : std::string("(no tree)")));
  }
}

}  // namespace

experiment::PredictionSetup PredictionSetupOf(const Segment &segment, const CompoundSystem &compound) {
  experiment::PredictionSetup s;
  s.reaction = segment.GetReaction();
  s.absolute = segment.IsCoefficientAbsolute();
  s.maxAngDistOrder = segment.GetMaxAngDistOrder();
  s.upos = segment.IsUPOS();
  if (s.upos) {
    s.uposL = segment.GetSecondaryDecayL();
    s.uposIc = segment.GetIc();
    s.uposDelta = segment.GetDelta();
  }
  experiment::PredictionTerm base;
  base.entranceKey = segment.GetEntranceKey();
  base.exitKey = segment.GetExitKey();
  base.nucleus = NucleusOf(compound, base.entranceKey);
  base.total = segment.IsTotalTerm();
  s.terms.push_back(base);
  // The leaves: the line's components first (with their fixed angles), then the further capture pairs of
  // a total capture, then a yield's total term (WorkingCopy::Fill, FinishTransformedSegment).
  const std::vector<experiment::ComponentSpec> &line = segment.GetComponents();
  for (int leaf = 1; leaf < segment.NumTerms(); leaf++) {
    const std::size_t k = leaf - 1;
    experiment::PredictionTerm t;
    t.entranceKey = segment.TermEntranceKey(leaf);
    t.exitKey = segment.TermExitKey(leaf);
    t.scaling = segment.LeafScaling(leaf);
    if (k < line.size() && line[k].hasAngle) {
      t.hasAngle = true;
      t.angle = line[k].angle;
    }
    t.nucleus = NucleusOf(compound, t.entranceKey);
    t.total = segment.IsTermTotal(leaf);
    t.yieldTotal = segment.IsYieldTotalLeaf(leaf);
    s.terms.push_back(t);
  }
  if (segment.IsNatural()) {
    s.sample = segment.GetNaturalSample();
    s.sampleName = segment.GetNaturalToken();
  }
  s.thickness = segment.GetNominalThickness();
  return s;
}

std::vector<std::string> ProjectPredictions(const std::string &projectFile, std::ostream &messages) {
  std::vector<std::string> lines;
  Config configure(messages);
  configure.configfile = projectFile;
  if (configure.ReadConfigFile() < 0) {
    lines.push_back("project: cannot be read");
    return lines;
  }
  // The trees need the compound's structure only, never the finished compound (CompoundSystemSetUp.h).
  bool levels = true;
  const bool read = InspectCompoundStructure(configure, [&](const CompoundSystem &structure, const LevelState &, const LevelInputs &inputs) {
    if (structure.NumPairs() == 0) {
      levels = false;
      return;
    }
    WorkingCopy data, extrapolations;
    if (data.Fill(configure, &structure, inputs) != 0)
      lines.push_back("data: set-up failed");
    else
      AppendLines(data, "data", RecordNames(configure, "dataset"), lines);
    if (extrapolations.MakePoints(configure, &structure, inputs) != 0)
      lines.push_back("extrapolation: set-up failed");
    else
      AppendLines(extrapolations, "extrapolation", RecordNames(configure, "extrapolation"), lines);
  });
  if (!read || !levels) lines.push_back("project: <levels> cannot be set up");
  return lines;
}

}  // namespace azure2

bool WorkingCopy::BuildPredictions(const CompoundSystem *compound, const Config &configure) {
  for (int s = 1; s <= this->Draft().NumSegments(); s++) {
    Segment *segment = this->Draft().GetSegment(s);
    try {
      segment->SetPrediction(std::make_shared<const azure2::experiment::Prediction>(
          azure2::experiment::BuildPrediction(azure2::PredictionSetupOf(*segment, *compound))));
    } catch (const std::exception &e) {
      configure.outStream << "ERROR: segment " << segment->GetSegmentKey() << ": no model expression: " << e.what()
                          << std::endl;
      return false;
    }
  }
  return true;
}
