// The set-up of the compound system and of the data, once (CompoundSystemSetUp.h).
#include "CompoundSystemSetUp.h"

#include "CompoundSystem.h"
#include "Config.h"
#include "GSLException.h"

#include <ostream>

void CompoundSetUpRequests::AddEntrance(int pair) {
  for (int r : entrancePairs)
    if (r == pair) return;
  entrancePairs.push_back(pair);
}

void CompoundSetUpRequests::AddCaptureAy(int entrancePair, int exitPair) {
  const std::pair<int, int> request(entrancePair, exitPair);
  for (const std::pair<int, int> &r : captureAyPairs)
    if (r == request) return;
  captureAyPairs.push_back(request);
}

CompoundSetUp SetUpCompoundSystem(const Config &configure, const DataSetUpStep &dataSetUp,
                                  const std::map<int, double> &radii) {
  CompoundSetUp out;
  // Made in place and handed over at the end: the data's set-up may keep the structure's address, which is the
  // finished object's.
  std::unique_ptr<CompoundSystem> building(new CompoundSystem());
  if (building->Fill(configure, out.state, out.inputs, radii) == -1) {
    out.refusal = CompoundSetUp::Refusal::kUnreadable;
    return out;
  } else if (building->NumPairs() == 0 || building->NumJGroups() == 0) {
    out.refusal = CompoundSetUp::Refusal::kNoNuclearData;
    return out;
  }
  const std::optional<CompoundSetUpRequests> requests = dataSetUp(*building, out.inputs);
  if (!requests) {
    out.refusal = CompoundSetUp::Refusal::kData;
    return out;
  }
  try {
    building->Initialize(configure, out.state, out.inputs, *requests);
  } catch (const GSLException &e) {
    out.refusal = CompoundSetUp::Refusal::kFinish;
    out.finishError = e.what();
    return out;
  }
  out.compound.reset(building.release());
  return out;
}

void PrintSetUpRefusal(const CompoundSetUp &setUp, std::ostream &out) {
  switch (setUp.refusal) {
    case CompoundSetUp::Refusal::kNone: return;
    case CompoundSetUp::Refusal::kUnreadable: out << "Could not fill compound nucleus from file." << std::endl; return;
    case CompoundSetUp::Refusal::kNoNuclearData:
      out << "No nuclear data exists. Calculation not possible." << std::endl;
      return;
    case CompoundSetUp::Refusal::kData: return;
    case CompoundSetUp::Refusal::kFinish:
      out << setUp.finishError << std::endl << std::endl << "Calculation was aborted." << std::endl;
      return;
  }
}

bool InspectCompoundStructure(
    const Config &configure,
    const std::function<void(const CompoundSystem &structure, const LevelState &state, const LevelInputs &inputs)> &inspect) {
  CompoundSystem structure;
  LevelState state;
  LevelInputs inputs;
  if (structure.Fill(configure, state, inputs, std::map<int, double>()) == -1) return false;
  inspect(structure, state, inputs);
  return true;
}
