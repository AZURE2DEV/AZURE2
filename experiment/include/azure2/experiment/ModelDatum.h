#ifndef AZURE2_EXPERIMENT_MODELDATUM_H
#define AZURE2_EXPERIMENT_MODELDATUM_H

#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/LevelTables.h"
#include "azure2/experiment/Schema.h"

namespace azure2 {
namespace experiment {

/*!
 * Model data (docs/upgrade/project_format3.md section 3.10, docs/upgrade/data_model.md section 4, decision M10): a
 * literature value of a level property or of a scattering length, compared with a function of the level parameters
 * rather than with the model at a beam energy.  A validated `datum` record (Schema::Validate) becomes a ModelDatum
 * (SetupSpecs.h) once its level and pairs are found in <levels>: BuildModelData, which takes <levels> as the plain
 * tables of LevelTables.h, so that it needs no engine (they are made from the <levels> rows, LevelTablesOfRows; the
 * tests write them by hand).  Schema::Build calls it for the <modelData> records of a format-3 project.
 */

/// (table name, quantity) for every datum quantity the engine evaluates.
const std::vector<std::pair<std::string, ModelQuantity>> &ModelQuantities();

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
