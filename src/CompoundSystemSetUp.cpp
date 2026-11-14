// The set-up of the compound system and of the data, once (CompoundSystemSetUp.h).
#include "CompoundSystemSetUp.h"

#include "CompoundSystem.h"
#include "Config.h"
#include "GSLException.h"

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
  for (int pair : requests->entrancePairs) building->RequestEntrance(pair);
  for (const std::pair<int, int> &r : requests->captureAyPairs) building->RequestCaptureAy(r.first, r.second);
  try {
    building->Initialize(configure, out.state, out.inputs);
  } catch (const GSLException &e) {
    out.refusal = CompoundSetUp::Refusal::kFinish;
    out.finishError = e.what();
    return out;
  }
  out.compound.reset(building.release());
  return out;
}

bool InspectCompoundStructure(const Config &configure,
                              const std::function<void(const CompoundSystem &structure, const LevelInputs &inputs)> &inspect) {
  CompoundSystem structure;
  LevelState state;
  LevelInputs inputs;
  if (structure.Fill(configure, state, inputs) == -1) return false;
  inspect(structure, inputs);
  return true;
}
