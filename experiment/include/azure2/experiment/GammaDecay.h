#ifndef AZURE2_EXPERIMENT_GAMMADECAY_H
#define AZURE2_EXPERIMENT_GAMMADECAY_H

#include <string>
#include <vector>

#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/LevelTables.h"

namespace azure2 {
namespace experiment {

struct ValidationResult;  // Schema.h (which includes Build.h, which includes this header)

/*!
 * The radiative decay of a level written as its total radiative width, branchings and mixing ratios (a `gamma_decay`
 * record of <gammaDecays>, docs/upgrade/gamma_parametrisation.md): the validated record resolved against <levels>.
 * Plain data in engine units; no engine headers.  The engine (src/GammaDecayMap.cpp) turns it into partial widths
 * and parameters; Schema::Build calls BuildGammaDecays for the refusals of --validate and pyazr.
 */

/// One number of the record: `<value> fixed|fit [prior <centre> <width>[%]]`, in the engine unit (MeV for the width,
/// 1 for a branching or a mixing ratio).
struct GammaNumber {
  double value = 0.0;
  bool fit = false;
  bool hasPrior = false;
  double priorCentre = 0.0, priorWidth = 0.0;  ///< the width in the engine unit (a percentage already applied)
};

/// One final state of the level: the photon pair and its one or two multipolarities.
struct GammaFinalState {
  int pairKey = 0;
  std::string label;            ///< the pair as the record writes it, or PairLabel's text for a final state not listed
  std::vector<int> channels;    ///< the channels of the level's J^pi group (1-based), by increasing L
  std::vector<int> L;           ///< their multipolarities
  enum class Role { Listed, Remainder, NotListed } role = Role::NotListed;
  GammaNumber branching;        ///< b_f (Listed; NotListed: 0, fixed; Remainder: unused)
  bool hasMixing = false;       ///< a `mixing` line (only with two multipolarities)
  GammaNumber mixing;           ///< delta_f, Krane-Steffen (0, fixed without a line)
  int sign = 1;                 ///< sigma_f, the interference sign of the L amplitude
};

struct GammaDecaySpec {
  std::string name;
  int line = 0;
  int levelId = 0;
  int group = 0, level = 0;          ///< the engine's J^pi group and level in it (1-based)
  double levelEnergy = 0.0;          ///< MeV, E_x as written
  GammaNumber totalWidth;            ///< Gamma_gamma, MeV
  std::vector<GammaFinalState> finals;  ///< every photon final state of the level's group, in channel order
  int remainder = -1;                ///< index of the final state that takes 1 - sum of the others
  bool hasRemainderPrior = false;
  double remainderPriorCentre = 0.0, remainderPriorWidth = 0.0;
};

struct GammaDecayBuild {
  std::vector<GammaDecaySpec> specs;
  std::vector<Diagnostic> diagnostics;
  bool Ok() const;
};

/// Tolerance of the sum of fixed branchings (absolute).
constexpr double kBranchingSumTolerance = 1e-6;

/*!
 * The gamma_decay records of `validated` against <levels> (`tables`: its levels, pairs and row slots).  Refusals:
 * E_UNKNOWN_LEVEL, E_LEVEL_JPI (as a datum), E_UNKNOWN_PAIR / E_AMBIGUOUS_PAIR / E_NOT_A_CHANNEL (a final state that is no
 * photon channel of the level), E_GAMMA_MIXTURE (a photon width of the level written in <levels>), E_BRANCHING_SUM
 * (fixed branchings that do not sum to 1, or a free one without a remainder), E_BRANCHING_REMAINDER (a negative
 * remainder), E_MIXING_ONE_MULTIPOLE, E_GAMMA_RECORD (a level with two records, a final state named twice, a
 * remainder prior that is not two numbers, a level without photon channels).
 */
GammaDecayBuild BuildGammaDecays(const ValidationResult &validated, const LevelTables &tables);

}  // namespace experiment
}  // namespace azure2

#endif
