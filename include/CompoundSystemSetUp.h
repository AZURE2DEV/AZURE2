#ifndef COMPOUNDSYSTEMSETUP_H
#define COMPOUNDSYSTEMSETUP_H

#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "LevelInputs.h"
#include "LevelState.h"

class CompoundSystem;
class Config;

/*!
 * The set-up of a project's compound system and of its data, in the one order the code needs
 * (src/CompoundSystemSetUp.cpp): the compound's structure from the <levels> rows (CompoundSystem::Fill), then the
 * caller's data set-up, which reads that structure and hands back what the compound's last part must build, then
 * the last part (CompoundSystem::Initialize: B_c, the input transformation, the pathway table, the external-capture
 * set-up).  The program's run, the sampler and the session call it; nobody else makes a CompoundSystem.
 */

/// What the data's set-up asks of the compound's last part: the pathways of each entrance pair, and the capture
/// analyzing-power terms of each (entrance, exit) pair.  Pair numbers are 1-based positions in the compound's pairs.
struct CompoundSetUpRequests {
  std::vector<int> entrancePairs;                      ///< in the order first asked, each once
  std::vector<std::pair<int, int> > captureAyPairs;    ///< (entrance, exit), in the order first asked, each once
  /// Ask for the pathways of entrance pair `pair`; a pair already asked for is not added again.
  void AddEntrance(int pair);
  /// Ask for the capture analyzing-power terms of (entrance pair, exit pair); a repeat is not added again.
  void AddCaptureAy(int entrancePair, int exitPair);
};

/// The caller's data set-up, run between the compound's structure and its last part.  It reads the structure (pairs,
/// J^pi groups, channels, levels' set-up data, model data) and the levels' input values as the file gives them, sets
/// up the caller's data from them, and returns what the last part must build; or nothing when it refused (it has said
/// why).  What the last part writes (B_c, the final-state energies, the pathway table, the external-capture blocks,
/// the capture masks) cannot be read here: those accessors throw std::logic_error.
using DataSetUpStep =
    std::function<std::optional<CompoundSetUpRequests>(const CompoundSystem &structure, const LevelInputs &inputs)>;

/// The outcome of SetUpCompoundSystem.
struct CompoundSetUp {
  /// Which step refused; the caller says it in its own words (the session is silent, the program is not).
  enum class Refusal {
    kNone,           ///< set up: `compound` is the finished system
    kUnreadable,     ///< the <levels> rows could not be made into a structure (CompoundSystem::Fill said why)
    kNoNuclearData,  ///< no pair or no J^pi group
    kData,           ///< the caller's data set-up refused
    kFinish          ///< the last part refused (B_c, the transformation): `finishError` is its message
  };
  Refusal refusal = Refusal::kNone;
  std::string finishError;
  /// The finished compound system, constant; null unless refusal is kNone.  One owner (the run, the sampler's run,
  /// the session); every evaluation, thread and walker borrows it as `const CompoundSystem *` or `&`.
  std::unique_ptr<const CompoundSystem> compound;
  /// The owner's levels: shaped for the compound's levels, the input values (transformed into the fit's basis by
  /// the last part when the project asks for it) and the zero-width flags.
  LevelState state;
  LevelInputs inputs;
};

/// Set up the compound system of the project of `configure` and, in the middle, the caller's data (`dataSetUp`).
/// `radii` overrides channel radii (pair number, 1-based -> fm), as the session's SetRadius does.  Exceptions other
/// than the last part's GSLException (a refusal, kFinish) pass through.
CompoundSetUp SetUpCompoundSystem(const Config &configure, const DataSetUpStep &dataSetUp,
                                  const std::map<int, double> &radii = std::map<int, double>());

/// The compound's structure only, never finished: `inspect` reads it with the levels as the structure's set-up makes
/// them (the input values as the file gives them), for a caller that needs the structure and no evaluation (the
/// segments' trees, Predictions.cpp; tests of the input values).  false when the <levels> rows could not be made into
/// a structure (`inspect` is then not called).
bool InspectCompoundStructure(
    const Config &configure,
    const std::function<void(const CompoundSystem &structure, const LevelState &state, const LevelInputs &inputs)> &inspect);

#endif
