// The tests' set-up of a project: the compound system and the data made by the engine's one set-up
// (CompoundSystemSetUp.h), without the program's messages.
#ifndef SETUPFORTESTS_H
#define SETUPFORTESTS_H

#include <map>
#include <optional>

#include "CaptureAyRequests.h"
#include "CompoundSystem.h"
#include "CompoundSystemSetUp.h"
#include "Config.h"
#include "ProjectSetup.h"
#include "WorkingCopy.h"

namespace test_setup {

/// Whether the data's capture-A_y terms are asked for, as the program asks for them (RequestCaptureAyTerms).
enum class CaptureAy { kRequested, kNotRequested };

/// The compound with `data` filled from the project's data sets in its middle and the data sets' entrances asked for.
inline CompoundSetUp WithData(const Config &cfg, WorkingCopy &data, CaptureAy captureAy = CaptureAy::kNotRequested) {
  return SetUpCompoundSystem(cfg, [&](const CompoundSystem &structure,
                                      const LevelInputs &inputs) -> std::optional<CompoundSetUpRequests> {
    if (data.Fill(cfg, &structure, inputs) != 0) return std::nullopt;
    CompoundSetUpRequests requests;
    RequestDataEntrances(structure, cfg, requests);
    if (captureAy == CaptureAy::kRequested) RequestCaptureAyTerms(structure, data, requests);
    return requests;
  });
}

/// The compound alone: no data, nothing asked for (no entrance pathways).
inline CompoundSetUp WithoutData(const Config &cfg, const std::map<int, double> &radii = std::map<int, double>()) {
  return SetUpCompoundSystem(
      cfg, [](const CompoundSystem &, const LevelInputs &) { return std::optional<CompoundSetUpRequests>(CompoundSetUpRequests()); },
      radii);
}

/// The levels as the compound's structure is made, before its last part: the file's input values (not transformed)
/// and the state's zero-width flags.  false when the structure cannot be made.
inline bool LevelsBeforeTheLastPart(const Config &cfg, LevelState &state, LevelInputs &inputs) {
  return InspectCompoundStructure(cfg, [&](const CompoundSystem &, const LevelState &s, const LevelInputs &in) {
    state = s;
    inputs = in;
  });
}

}  // namespace test_setup

#endif
