#ifndef AZURE2_EXPERIMENT_MODELQUANTITY_H
#define AZURE2_EXPERIMENT_MODELQUANTITY_H

namespace azure2 {
namespace experiment {

/*!
 * The quantities a model datum can be (docs/upgrade/data_model.md section 11): the key between a row of the quantity
 * table (Quantities.cpp, `q.datum = Datum(ModelQuantity::..., ...)`: its name, fields, units and build rule) and the
 * engine's value (`LevelValue`, src/ModelData.cpp, a switch on this enum that the compiler checks: -Werror=switch).
 * A new value here is named by exactly one row [repo_consistency check 27] and has its case there.
 */
enum class ModelQuantity {
  ResonanceStrength,  ///< omega gamma = omega Gamma_a Gamma_b / Gamma (data_kinds.md eq. 8)
  ResonanceEnergy,    ///< the observed level energy E_lambda, as E_x, E_cm or E_lab of a pair
  TotalWidth,         ///< Gamma = sum_c |Gamma_c|
  PartialWidth,       ///< Gamma of one pair: the sum over its channels
  BranchingRatio,     ///< Gamma_b / Gamma
  Anc,                ///< |C| of a closed channel (pair, l, s)
  ScatteringLength    ///< Re and Im of the s-wave scattering length a_0 (include/ScatteringLength.h)
};

}  // namespace experiment
}  // namespace azure2

#endif
