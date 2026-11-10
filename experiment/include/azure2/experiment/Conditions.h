#ifndef AZURE2_EXPERIMENT_CONDITIONS_H
#define AZURE2_EXPERIMENT_CONDITIONS_H

#include <string>
#include <utility>
#include <vector>

namespace azure2 {
namespace experiment {

/*!
 * When a field applies, is required, or a unit or norm mode is allowed (docs/upgrade/project_format3.md section 7.1,
 * "conditional"): a closed set of conditions on one record, evaluated by the validator with the record's parsed
 * fields and the records it references.  Exported in the schema as {"if": <kind>, ...} objects so that a GUI can show
 * and hide fields with the same rules.
 */
struct Condition {
  enum Kind {
    Always,
    FieldIs,        ///< `field`'s word (an enum value or a structured field's leading keyword) is one of `values`
    FieldPresent,   ///< `field` is given
    UnitDimension,  ///< the record's `unit` has one of the dimensions `values`
    RatioReaction,  ///< the reaction expression is a ratio
    SingleReaction, ///< the reaction expression is one reaction (no sum, no ratio)
    AngleAxis,      ///< the record has an angle axis: a theta column (data set) or an angle grid (extrapolation)
    GammaExit,      ///< every reaction term has a photon exit
    ThmExperiment,  ///< `thm` names an experiment of <thm> (not `global`)
    Fold,           ///< the record is folded: a target that integrates or broadens, a beam, or a detector resolution
  };
  Kind kind = Always;
  std::string field;
  std::vector<std::string> values;
  bool negate = false;

  static Condition Is(const std::string &field, std::vector<std::string> values) {
    return {FieldIs, field, std::move(values), false};
  }
  static Condition Present(const std::string &field) { return {FieldPresent, field, {}, false}; }
  static Condition Absent(const std::string &field) { return {FieldPresent, field, {}, true}; }
  static Condition Dimension(std::vector<std::string> dims) { return {UnitDimension, "unit", std::move(dims), false}; }
  static Condition Of(Kind k, bool negate = false) { return {k, "", {}, negate}; }
};

const char *ConditionKindName(Condition::Kind k);
/// The condition in words, for messages: "kind is layer or sample", "the unit is per_charge", "no fold".
std::string Describe(const Condition &c);

}  // namespace experiment
}  // namespace azure2

#endif
