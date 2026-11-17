// Helpers of the typed values (Values.h).
#include "azure2/experiment/Values.h"

namespace azure2 {
namespace experiment {

const char *ValueTypeName(ValueType t) {
  switch (t) {
    case ValueType::Number: return "number";
    case ValueType::Integer: return "integer";
    case ValueType::Unit: return "unit";
    case ValueType::NumberWithUnit: return "number_with_unit";
    case ValueType::Range: return "range";
    case ValueType::Nuisance: return "nuisance";
    case ValueType::Enum: return "enum";
    case ValueType::Reference: return "reference";
    case ValueType::Name: return "name";
    case ValueType::Path: return "path";
    case ValueType::Text: return "text";
    case ValueType::Formula: return "formula";
    case ValueType::NumberList: return "number_list";
    case ValueType::Spin: return "spin";
    case ValueType::SpinParity: return "spin_parity";
    case ValueType::SpinOperator: return "spin_operator";
    case ValueType::Nucleus: return "nucleus";
    case ValueType::Pair: return "pair";
    case ValueType::Composition: return "composition";
    case ValueType::Reaction: return "reaction_expr";
    case ValueType::Columns: return "columns";
    case ValueType::Grid: return "grid";
    case ValueType::Form: return "structured";
  }
  return "?";
}

const char *NuisanceModeName(NuisanceMode m) {
  switch (m) {
    case NuisanceMode::Fixed: return "fixed";
    case NuisanceMode::Fit: return "fit";
    case NuisanceMode::Profiled: return "profiled";
  }
  return "?";
}

bool SingleReaction::Elastic() const {
  const bool sameState = (state.kind == ResidualState::Ground && targetState.kind == ResidualState::Ground) ||
                         (state.kind == ResidualState::Excited && targetState.kind == ResidualState::Excited &&
                          state.energy.value == targetState.energy.value);
  return projectile.Same(ejectile) && target.Same(residual) && sameState && !projectile.excited && !ejectile.excited;
}

bool SingleReaction::ChargedEntrance() const { return projectile.Z > 0 && target.Z > 0; }

bool SingleReaction::IdenticalEntrance() const { return projectile.Z == target.Z && projectile.A == target.A; }

bool SingleReaction::GammaExit() const { return ejectile.light == "g"; }

std::string SingleReaction::Canonical() const {
  auto particle = [](const Particle &p) {
    std::string s = p.light.empty() ? p.text : p.light;
    if (p.excited) s += "[" + p.excitation.text + "]";
    return s;
  };
  std::string s = target.text + (targetState.kind == ResidualState::Excited ? "[" + targetState.energy.text + "]" : "") + "(" +
                  particle(projectile) + "," + particle(ejectile) + ")" + residual.text;
  if (state.kind == ResidualState::Total) s += "[total]";
  if (state.kind == ResidualState::Excited) s += "[" + state.energy.text + "]";
  return s;
}

std::string ReactionExpr::Tree() const {
  auto term = [](const ReactionTerm &t) {
    std::string s = t.reaction.Canonical();
    if (t.hasAngle) s += " @ " + t.angle.text + " " + t.angleUnit->token + " " + t.angleFrame;
    return t.hasScale ? "Scale(" + t.scale.text + ", " + s + ")" : s;
  };
  if (op == Single && terms.size() == 1) return term(terms[0]);
  std::string s = op == Ratio ? "Ratio(" : "Sum(";
  for (std::size_t i = 0; i < terms.size(); i++) s += (i ? ", " : "") + term(terms[i]);
  return s + ")";
}

const ColumnRole *Columns::Find(const std::string &role) const {
  for (const ColumnRole &r : roles)
    if (r.role == role) return &r;
  return nullptr;
}

const ColumnRole *Columns::EnergyRole() const {
  for (const ColumnRole &r : roles)
    if (r.role.compare(0, 2, "E_") == 0) return &r;
  return nullptr;
}

const ColumnRole *Columns::AngleRole() const {
  for (const ColumnRole &r : roles)
    if (r.role.compare(0, 6, "theta_") == 0) return &r;
  return nullptr;
}

const Value *Value::Get(const std::string &name) const {
  for (const ValuePart &p : parts)
    if (p.name == name) return &p.value;
  return nullptr;
}

}  // namespace experiment
}  // namespace azure2
