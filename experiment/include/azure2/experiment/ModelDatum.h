#ifndef AZURE2_EXPERIMENT_MODELDATUM_H
#define AZURE2_EXPERIMENT_MODELDATUM_H

#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/Schema.h"

namespace azure2 {
namespace experiment {

/*!
 * Model data (docs/upgrade/project_format3.md section 3.10, docs/upgrade/data_model.md section 4, decision M10): a
 * literature value of a level property or of a scattering length, compared with a function of the level parameters
 * rather than with the model at a beam energy.  A validated `datum` record (Schema::Validate) becomes a ModelDatum
 * once its level and pairs are found in <levels>: BuildModelData, which takes <levels> as the plain tables below, so
 * that it needs no engine (the engine fills them from its compound nucleus; the tests write them by hand).
 *
 * A ModelDatum is plain data, immutable after it is built.  Values are in the engine unit of the quantity's dimension
 * (MeV for energies and widths, 1, fm, fm^-1/2: Units.h), with the unit and factor as written kept for the output.
 */

/// The quantities a datum can be, each evaluated by the engine (src/ModelData.cpp).  ModelQuantities() maps the
/// table's names (Quantities.cpp, record kind datum) to them; repo_consistency check 27 fails for a datum quantity of
/// the table without a row there and a test named after it.
enum class ModelQuantity {
  ResonanceStrength,  ///< omega gamma = omega Gamma_a Gamma_b / Gamma (data_kinds.md eq. 8)
  ResonanceEnergy,    ///< the observed level energy E_lambda, as E_x, E_cm or E_lab of a pair
  TotalWidth,         ///< Gamma = sum_c |Gamma_c|
  PartialWidth,       ///< Gamma of one pair: the sum over its channels
  BranchingRatio,     ///< Gamma_b / Gamma
  Anc,                ///< |C| of a closed channel (pair, l, s)
  ScatteringLength    ///< Re and Im of the s-wave scattering length a_0 (include/ScatteringLength.h)
};

/// (table name, quantity) for every datum quantity the engine evaluates.
const std::vector<std::pair<std::string, ModelQuantity>> &ModelQuantities();

/// What BuildModelData needs of a pair of <levels> (the engine's PPair).
struct PairInfo {
  int key = 0;   ///< the pair key of <levels>
  int kind = 0;  ///< 0 particle pair, 10 photon (capture to a level of the compound nucleus), 20 beta
  int lightZ = 0, lightA = 0;  ///< the light particle (a photon: 0, 0)
  int heavyZ = 0, heavyA = 0;  ///< the heavy particle (photon pair: the compound nucleus)
  double lightExcitation = 0.0;  ///< MeV
  double heavyExcitation = 0.0;  ///< MeV (photon pair: E_x of the final level)
  double threshold = 0.0;  ///< MeV, the E_x of the compound nucleus at which the pair opens (S + e)
  int twoJLight = 0, twoJHeavy = 0;  ///< 2j of the two particles
  bool identical = false;            ///< the two particles are identical
  std::string label;                 ///< for messages ("p+15N")
};

/// One channel of a level's J^pi group.
struct LevelChannelInfo {
  int pairKey = 0;
  int l = 0;
  int twoS = 0;
};

/// What BuildModelData needs of a level of <levels> (the engine's ALevel and its J^pi group).
struct LevelInfo {
  std::vector<int> ids;  ///< the level IDs (column 9) of the <levels> lines read into this level
  int twoJ = 0;
  int parity = 0;
  double energy = 0.0;  ///< MeV, E_x as written
  /// MeV, the total width from the widths written in <levels> (the start of a fit), for the warnings; < 0: unknown
  double width = -1.0;
  int group = 0, level = 0;  ///< the engine's indices of the J^pi group and of the level in it (1-based)
  std::vector<LevelChannelInfo> channels;  ///< the channels of its J^pi group
};

struct ModelDatum {
  std::string name;
  std::string source;
  int line = 0;  ///< the record's header line
  ModelQuantity quantity = ModelQuantity::ResonanceStrength;
  std::string quantityName;
  double value = 0.0;  ///< in the engine unit of the quantity's dimension
  double error = 0.0;  ///< > 0, same unit
  bool hasImag = false;  ///< scattering length: Im a_0 given
  double imag = 0.0, imagError = 0.0;
  std::string unit;  ///< the unit as written (for the output)
  double unitFactor = 1.0;  ///< the written unit in engine units
  // Level quantities: the level.
  int levelId = 0;
  int group = 0, level = 0;  ///< the engine's indices (LevelInfo)
  std::string jpi;           ///< as written (3-)
  // Pairs, as <levels> keys (0 when the quantity has none).
  int entrancePair = 0;  ///< resonance_strength: the entrance pair a + A
  int exitPair = 0;      ///< resonance_strength, branching_ratio: the exit pair b + B
  int pair = 0;          ///< resonance_energy (E_cm, E_lab), partial_width, anc, scattering_length
  std::string frame;     ///< resonance_energy: E_x, E_cm or E_lab
  int l = -1;            ///< anc
  int twoS = -1;         ///< anc, scattering_length (2s)
  bool coherent = false; ///< scattering_length: the spin-weighted sum over the channel spins
};

struct ModelDataBuild {
  std::vector<ModelDatum> data;          ///< in record order
  std::vector<Diagnostic> diagnostics;   ///< refusals (errors) and warnings, in record order
  bool Ok() const;                       ///< no diagnostic of severity error
};

/*!
 * The model data of validated records (the `datum` records of `validated`; records that did not validate are left
 * out, their diagnostics are the validator's) against <levels>: each record's level by its ID and J^pi (F13 (a)), its
 * reaction or pair by the pairs' particles and states.  Refusals: E_UNKNOWN_LEVEL (no level has the ID, or two of the
 * same J^pi), E_LEVEL_JPI (the level's J^pi is not the datum's), E_UNKNOWN_PAIR (no pair has the particles and
 * states), E_AMBIGUOUS_PAIR (two pairs have them), E_NOT_A_CHANNEL (the pair, or (pair, l, s), is no channel of the
 * level; a scattering length's s is not an s-wave channel spin of the pair), E_OUT_OF_RANGE (a strength of a level
 * below the entrance threshold), E_FIELD_NOT_APPLICABLE (an ANC of an open channel; a coherent ANC or the coherent
 * scattering length of identical particles), E_REACTION_NOT_ALLOWED (a [total] exit).  Warnings for a strength
 * (data_kinds.md section 4.7, where eq. (8) and the area of the resonance differ): W_BROAD_LEVEL (Gamma / E_r > 1 %),
 * W_NEARBY_LEVEL (another level of the same J^pi within 10 Gamma), from the widths written in <levels>.
 */
ModelDataBuild BuildModelData(const ValidationResult &validated, const std::vector<LevelInfo> &levels,
                              const std::vector<PairInfo> &pairs);

}  // namespace experiment
}  // namespace azure2

#endif
