// The expression trees of the segments (Predictions.h; docs/upgrade/data_model.md step S10).
#include "Predictions.h"

#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "ESegment.h"
#include "PPair.h"

namespace azure2 {

namespace {

int NucleusOf(CNuc &compound, int entranceKey) {
  if (!compound.IsPairKey(entranceKey)) return 0;
  return compound.GetPair(compound.GetPairNumFromKey(entranceKey))->GetNucleus();
}

/// The trees of one EData as lines "<what> <key>: <tree>".
void AppendLines(EData &data, const char *what, std::vector<std::string> &lines) {
  for (int s = 1; s <= data.NumSegments(); s++) {
    const ESegment &segment = *data.GetSegment(s);
    const experiment::Prediction *p = segment.GetPrediction();
    lines.push_back(std::string(what) + " " + std::to_string(segment.GetSegmentKey()) + ": " +
                    (p ? experiment::PredictionText(*p) : std::string("(no tree)")));
  }
}

}  // namespace

experiment::PredictionSetup PredictionSetupOf(const ESegment &segment, CNuc &compound) {
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
  // The component segments: the line's components first (with their fixed angles), then the further capture pairs of
  // a total capture, then a yield's total term (EData::Fill, FinishTransformedSegment).
  const std::vector<experiment::ComponentSpec> &line = segment.GetComponents();
  const std::vector<ESegment *> &components = segment.GetComponentSegments();
  for (std::size_t k = 0; k < components.size(); k++) {
    const ESegment *c = components[k];
    if (!c) throw std::logic_error("PredictionSetupOf: segment " + std::to_string(segment.GetSegmentKey()) +
                                   " has a null component");
    experiment::PredictionTerm t;
    t.entranceKey = c->GetEntranceKey();
    t.exitKey = c->GetExitKey();
    t.scaling = c->GetComponentScaling();
    if (k < line.size() && line[k].hasAngle) {
      t.hasAngle = true;
      t.angle = line[k].angle;
    }
    t.nucleus = NucleusOf(compound, t.entranceKey);
    t.total = c->IsTotalTerm();
    t.yieldTotal = c->IsYieldTotalTerm();
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
  CNuc compound;
  if (compound.Fill(configure) == -1 || compound.NumPairs() == 0) {
    lines.push_back("project: <levels> cannot be set up");
    return lines;
  }
  EData data, extrapolations;
  if (data.Fill(configure, &compound) == -1)
    lines.push_back("data: set-up failed");
  else
    AppendLines(data, "data", lines);
  if (extrapolations.MakePoints(configure, &compound) == -1)
    lines.push_back("extrapolation: set-up failed");
  else
    AppendLines(extrapolations, "extrapolation", lines);
  return lines;
}

}  // namespace azure2

bool EData::BuildPredictions(CNuc *compound, const Config &configure) {
  for (int s = 1; s <= this->NumSegments(); s++) {
    ESegment *segment = this->GetSegment(s);
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
