#ifndef AZURE2_EXPERIMENT_QUANTITIES_H
#define AZURE2_EXPERIMENT_QUANTITIES_H

#include <string>
#include <vector>

#include "azure2/experiment/Conditions.h"

namespace azure2 {
namespace experiment {

/*!
 * The quantity table of format 3 (docs/upgrade/project_format3.md sections 7.1, 7.2): one row per measured or model
 * quantity -- its units, the column roles it accepts, the reactions it is defined for, the fields it requires and
 * allows, its norm modes, its tree, its EXFOR codes and its refusals.  The validator (Schema::Validate) reads the rows
 * and nothing else; a new quantity is one row here, a tree builder (step 3) and a test.
 */

/// A unit, or every unit of a dimension, allowed for the value column when `when` holds.
struct UnitRule {
  std::string unit;       ///< a unit token, or empty
  std::string dimension;  ///< a dimension (Units.h), or empty
  Condition when;
};

/// What a quantity's reaction expression must be (section 3.5; checked term by term).  `flags` is a sum of the
/// bits below; each failing bit is reported with the message of the Refusal whose `rule` names it.
struct ReactionRule {
  enum Flag : unsigned {
    kSum = 1u << 0,                ///< a sum of reactions (and scaled terms) is allowed       (rule "composite")
    kRatio = 1u << 1,              ///< a ratio of two reactions is allowed                    (rule "composite")
    kElastic = 1u << 2,            ///< every term elastic                                     (rule "inelastic")
    kChargedEntrance = 1u << 3,    ///< both entrance particles charged                        (rule "neutral")
    kNeutralEntrance = 1u << 4,    ///< an uncharged entrance particle                         (rule "charged_entrance")
    kIdentical = 1u << 5,          ///< projectile and target the same nucleus                 (rule "non_identical")
    kNotIdentical = 1u << 6,       ///<                                                        (rule "identical")
    kNoChargedElastic = 1u << 7,   ///< no charged-particle elastic term                       (rule "charged_elastic")
    kNoTotal = 1u << 8,            ///< no [total] exit                                        (rule "total")
    kSpinHalfProjectile = 1u << 9, ///< projectile n, p, t or h; other letters refused        (rule "spin")
    kNeutralPair = 1u << 10,       ///< the datum's pair is uncharged                          (rule "charged_pair")
  };
  unsigned flags = 0;
  std::string text;  ///< the rule in words (schema, messages)
  bool Has(Flag f) const { return (flags & f) != 0; }
};

/// A field the quantity requires or allows, when `when` holds (fields not listed are refused for the quantity).
struct QuantityField {
  std::string field;
  bool required = false;
  Condition when;
};

/// Norm modes allowed when `when` holds (checked in order; the first that applies decides).
struct NormRule {
  Condition when;
  std::vector<std::string> modes;  ///< of fixed, fit, profiled
  bool prior = true;               ///< a prior is allowed
  std::string message;
};

struct Refusal {
  std::string code;       ///< the diagnostic code it produces
  std::string rule;       ///< a rule of ReactionRule ("identical", ...) or the condition in words
  std::string message;
  std::string checkedBy;  ///< "validate" (here), "build" (step 3: needs <levels> or the data file), "setup" (the
                          ///< engine's set-up: needs the stopping on the sub-point grid) or "importer"
};

/*!
 * One row of the table.  Quantities.cpp fills each row member by member, in one function per quantity; a member a
 * row does not set keeps the default below (no angle axis, no order or window column, no reaction constraint, no
 * field, no target kind, no norm rule beyond GeneralNormRules, no EXFOR code, no refusal).
 */
struct Quantity {
  std::string name;
  std::vector<std::string> aliases;
  std::string title;
  std::string doc;
  std::vector<std::string> recordKinds;  ///< dataset, extrapolation | datum
  std::string engineUnit;
  std::vector<UnitRule> units;
  std::vector<std::string> energyRoles;  ///< column roles of the energy axis
  std::vector<std::string> angleRoles;
  enum AngleAxis { NoAngle, OptionalAngle, RequiredAngle } angle = NoAngle;
  bool orderRole = false;    ///< the Legendre order column
  bool windowRoles = false;  ///< window_lo window_hi (a beam-profile kernel)
  ReactionRule reaction;
  std::vector<QuantityField> fields;
  std::vector<std::string> targetKinds;  ///< target kinds allowed when a target is referenced
  std::vector<NormRule> norm;
  std::string tree;
  std::vector<std::string> exfor;
  std::vector<Refusal> refusals;
  /// Non-empty for a quantity that is no measurable observable: a data set of it gets W_NOT_MEASURABLE with this text.
  std::string datasetWarning;
};

const std::vector<Quantity> &Quantities();
/// The row of a name or alias, or nullptr.
const Quantity *FindQuantity(const std::string &nameOrAlias);
/// The names of the rows that a record kind may use.
std::vector<std::string> QuantityNames(const std::string &recordKind);

/// The norm rules that hold for every quantity (unit arbitrary), checked after the row's own.
const std::vector<NormRule> &GeneralNormRules();

}  // namespace experiment
}  // namespace azure2

#endif
