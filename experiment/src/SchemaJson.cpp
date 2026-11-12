// Schema::ToJson: the table as one JSON document (docs/upgrade/project_format3.md section 6.1), generated from the
// rows of Schema.cpp, Quantities.cpp, Units.cpp and Diagnostics.cpp and nothing else.  The committed copy
// tests/schema/schema_v3.json must equal it (check 21); docs/source/reference/format3.rst is generated from it
// (tools/gen_schema_docs.py, check 18).
#include <cmath>
#include <cstdio>
#include <string>

#include "azure2/experiment/Json.h"
#include "azure2/experiment/Parsers.h"
#include "azure2/experiment/Schema.h"
#include "azure2/physics/Elements.h"

namespace azure2 {
namespace experiment {

namespace {

// A ratio of two table factors to 15 digits: the factors are exact decimal powers, their quotient in binary is not.
std::string JsonNumber15(double v) {
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.15g", v);
  return buf;
}

void WriteCondition(JsonWriter &j, const Condition &c) {
  j.BeginObject();
  j.Key("if").String(ConditionKindName(c.kind));
  if (!c.field.empty() && c.kind != Condition::UnitDimension) j.Key("field").String(c.field);
  if (!c.values.empty()) j.Key("in").StringArray(c.values);
  if (c.negate) j.Key("not").Bool(true);
  j.Key("text").String(Describe(c));
  j.EndObject();
}

const char *TypeDoc(ValueType t) {
  switch (t) {
    case ValueType::Number: return "a real number, kept as written";
    case ValueType::Integer: return "an integer";
    case ValueType::Unit: return "one unit token of the unit table";
    case ValueType::NumberWithUnit: return "a number and its unit";
    case ValueType::Range: return "a low and a high number and their unit";
    case ValueType::Nuisance: return "a value with its mode; 'fit' without 'prior' is a free parameter; a width with % is relative "
                                     "to the centre; grid_max is the largest value the sub-point grids must span";
    case ValueType::Enum: return "one word of a closed list";
    case ValueType::Reference: return "the name of another record of the project";
    case ValueType::Name: return "a name [A-Za-z][A-Za-z0-9_-]* that is not a record";
    case ValueType::Path: return "a file path, relative to the project; quoted when it holds blanks";
    case ValueType::Text: return "free text, never read by the physics";
    case ValueType::Formula: return "a quoted expression in x and its parameters";
    case ValueType::NumberList: return "one or more numbers";
    case ValueType::Spin: return "a spin: 0, 1, 3/2, ...";
    case ValueType::SpinParity: return "a spin and parity: 0+, 3-, 1/2+, ...";
    case ValueType::Nucleus: return "a nuclide <A><Symbol>";
    case ValueType::Pair: return "a particle pair, light particle first";
    case ValueType::Composition: return "element symbols with their atoms per formula unit";
    case ValueType::Reaction: return "reactions <A><X>[<E_x>](<a>,<b>)<A><Y>[<E_x>|total] (the target's state for an excited entrance pair), a sum or a ratio of two, scaled terms "
                                     "<c>*<reaction>, terms at a fixed angle '@ <angle> <unit> lab|cm'";
    case ValueType::Columns: return "the role of each column of the data file, energy and angle roles with their unit";
    case ValueType::Grid: return "an energy or angle grid; range: E_k = lo + k step (F10)";
    case ValueType::Form: return "a structured value: a leading keyword chooses the alternative, then its parts";
  }
  return "";
}

void WriteSpecDetail(JsonWriter &j, const ValueSpec &s);

void WritePart(JsonWriter &j, const PartSpec &p) {
  j.BeginObject();
  j.Key("name").String(p.name);
  j.Key("keyword").Bool(p.keyword);
  j.Key("type").String(ValueTypeName(p.spec.type));
  j.Key("grammar").String(Grammar(p.spec));
  j.Key("optional").Bool(p.optional);
  if (!p.defaultText.empty()) j.Key("default").String(p.defaultText);
  WriteSpecDetail(j, p.spec);
  if (!p.doc.empty()) j.Key("doc").String(p.doc);
  j.EndObject();
}

void WriteSpecDetail(JsonWriter &j, const ValueSpec &s) {
  if (s.type == ValueType::Enum) {
    j.Key("values").StringArray(s.values);
    if (!s.aliases.empty()) {
      j.Key("aliases").BeginObject();
      for (const auto &a : s.aliases) j.Key(a.first).String(a.second);
      j.EndObject();
    }
  }
  if (s.type == ValueType::Reference) {
    j.Key("refers_to").StringArray(s.values);
    if (s.epoch) j.Key("epoch").Bool(true);
  }
  if (!s.refKinds.empty()) j.Key("or_refers_to").StringArray(s.refKinds);
  if (!s.dimension.empty()) {
    j.Key("dimension").String(s.dimension);
    if (s.dimension != "@quantity") j.Key("units").StringArray(UnitsOfDimension(s.dimension));
  }
  if (s.type == ValueType::Integer && s.min != LONG_MIN) j.Key("min").Integer(s.min);
  if (s.positive) j.Key("positive").Bool(true);
  if (s.profiled) j.Key("profiled").Bool(true);
  if (s.type == ValueType::Form) {
    j.Key("alternatives").BeginArray();
    for (const Alternative &a : s.form->alternatives) {
      j.BeginObject();
      j.Key("keyword").String(a.keyword);
      j.Key("parts").BeginArray();
      for (const PartSpec &p : a.parts) WritePart(j, p);
      j.EndArray();
      if (!a.doc.empty()) j.Key("doc").String(a.doc);
      j.EndObject();
    }
    j.EndArray();
  }
}

void WriteField(JsonWriter &j, const FieldDef &f) {
  j.Key(f.key).BeginObject();
  j.Key("title").String(f.title);
  j.Key("type").String(ValueTypeName(f.spec.type));
  j.Key("grammar").String(Grammar(f.spec));
  if (f.byQuantity) {
    j.Key("by_quantity").Bool(true);
  } else {
    if (f.applies.kind != Condition::Always || f.applies.negate) WriteCondition(j.Key("applies_when"), f.applies);
    if (f.required.kind == Condition::Always)
      j.Key("required").Bool(!f.required.negate);
    else
      WriteCondition(j.Key("required_when"), f.required);
  }
  if (f.repeatable) j.Key("repeatable").Bool(true);
  if (!f.defaultText.empty()) j.Key("default").String(f.defaultText);
  WriteSpecDetail(j, f.spec);
  j.Key("doc").String(f.doc);
  j.EndObject();
}

void WriteQuantity(JsonWriter &j, const Quantity &q) {
  j.Key(q.name).BeginObject();
  j.Key("title").String(q.title);
  j.Key("doc").String(q.doc);
  if (!q.aliases.empty()) j.Key("aliases").StringArray(q.aliases);
  j.Key("record_kinds").StringArray(q.recordKinds);
  j.Key("engine_unit").String(q.engineUnit);
  const Unit *engine = FindUnit(q.engineUnit);
  j.Key("units").BeginObject();
  for (const UnitRule &r : q.units) {
    const std::vector<std::string> tokens = r.dimension.empty() ? std::vector<std::string>{r.unit} : UnitsOfDimension(r.dimension);
    for (const std::string &t : tokens) {
      const Unit *u = FindUnit(t);
      j.Key(t);
      if (engine && std::string(u->dimension) == engine->dimension && std::string(u->dimension) != "arbitrary")
        j.Number(std::stod(JsonNumber15(u->factor / engine->factor)));  // 1e-3 / 1e-6 is 1000, not 1000.0000000000001
      else
        j.Null();  // per-charge and count units need q and the efficiency; arbitrary units have no factor
    }
  }
  j.EndObject();
  bool conditional = false;
  for (const UnitRule &r : q.units) conditional = conditional || r.when.kind != Condition::Always;
  if (conditional) {
    j.Key("unit_conditions").BeginArray();
    for (const UnitRule &r : q.units) {
      if (r.when.kind == Condition::Always) continue;
      j.BeginObject();
      j.Key(r.dimension.empty() ? "unit" : "dimension").String(r.dimension.empty() ? r.unit : r.dimension);
      WriteCondition(j.Key("when"), r.when);
      j.EndObject();
    }
    j.EndArray();
  }
  j.Key("columns").BeginObject();
  j.Key("energy").StringArray(q.energyRoles);
  j.Key("angle").StringArray(q.angleRoles);
  j.Key("angle_axis").String(q.angle == Quantity::NoAngle ? "none" : q.angle == Quantity::OptionalAngle ? "optional" : "required");
  j.Key("order").Bool(q.orderRole);
  j.Key("windows").Bool(q.windowRoles);
  j.EndObject();
  j.Key("reaction").BeginObject();
  j.Key("rule").String(q.reaction.text);
  static const std::pair<ReactionRule::Flag, const char *> kFlags[] = {
      {ReactionRule::kSum, "sum"}, {ReactionRule::kRatio, "ratio"}, {ReactionRule::kElastic, "elastic"},
      {ReactionRule::kChargedEntrance, "charged_entrance"}, {ReactionRule::kNeutralEntrance, "neutral_entrance"},
      {ReactionRule::kIdentical, "identical"}, {ReactionRule::kNotIdentical, "not_identical"},
      {ReactionRule::kNoChargedElastic, "no_charged_elastic"}, {ReactionRule::kNoTotal, "no_total"},
      {ReactionRule::kSpinHalfProjectile, "spin_half_projectile"}, {ReactionRule::kNeutralPair, "neutral_pair"}};
  std::vector<std::string> flags;
  for (const auto &f : kFlags)
    if (q.reaction.Has(f.first)) flags.push_back(f.second);
  j.Key("flags").StringArray(flags);
  j.EndObject();
  j.Key("fields").BeginObject();
  for (const QuantityField &f : q.fields) {
    j.Key(f.field).BeginObject();
    j.Key("required").Bool(f.required);
    if (f.when.kind != Condition::Always || f.when.negate) WriteCondition(j.Key("when"), f.when);
    j.EndObject();
  }
  j.EndObject();
  if (!q.targetKinds.empty()) j.Key("target_kinds").StringArray(q.targetKinds);
  if (!q.norm.empty()) {
    j.Key("norm_rules").BeginArray();
    for (const NormRule &r : q.norm) {
      j.BeginObject();
      WriteCondition(j.Key("when"), r.when);
      j.Key("modes").StringArray(r.modes);
      j.Key("prior").Bool(r.prior);
      j.Key("message").String(r.message);
      j.EndObject();
    }
    j.EndArray();
  }
  j.Key("tree").String(q.tree);
  j.Key("exfor").StringArray(q.exfor);
  j.Key("refusals").BeginArray();
  for (const Refusal &r : q.refusals) {
    j.BeginObject();
    j.Key("code").String(r.code);
    j.Key("rule").String(r.rule);
    j.Key("message").String(r.message);
    j.Key("checked_by").String(r.checkedBy);
    j.EndObject();
  }
  j.EndArray();
  if (!q.datasetWarning.empty()) j.Key("dataset_warning").String(q.datasetWarning);
  j.EndObject();
}

}  // namespace

std::string Schema::ToJson() {
  JsonWriter j(1);
  j.BeginObject();
  j.Key("schema_version").String(kVersion);
  j.Key("sections").BeginObject();
  for (const auto &s : Sections()) j.Key(s.first).StringArray(s.second);
  j.EndObject();

  j.Key("value_types").BeginObject();
  for (int t = 0; t <= static_cast<int>(ValueType::Form); t++) {
    ValueSpec s(static_cast<ValueType>(t));
    const bool dimensioned = s.type == ValueType::Unit || s.type == ValueType::NumberWithUnit || s.type == ValueType::Range ||
                             s.type == ValueType::Nuisance || s.type == ValueType::Grid;
    if (dimensioned) s.dimension = "*";  // any dimension: the field says which
    j.Key(ValueTypeName(s.type)).BeginObject();
    if (s.type != ValueType::Form) j.Key("grammar").String(Grammar(s));
    j.Key("doc").String(TypeDoc(s.type));
    j.EndObject();
  }
  j.EndObject();

  j.Key("records").BeginObject();
  for (const RecordKind &k : Kinds()) {
    j.Key(k.kind).BeginObject();
    j.Key("title").String(k.title);
    j.Key("doc").String(k.doc);
    j.Key("sections").StringArray(k.sections);
    j.Key("fields").BeginObject();
    for (const FieldDef &f : k.fields) WriteField(j, f);
    j.EndObject();
    if (!k.oneOf.empty()) {
      j.Key("one_of").BeginArray();
      for (const auto &g : k.oneOf) j.StringArray(g);
      j.EndArray();
    }
    if (!k.refused.empty()) {
      j.Key("refused").BeginObject();
      for (const auto &r : k.refused) j.Key(r.first).String(r.second);
      j.EndObject();
    }
    if (!k.deleted.empty()) {
      j.Key("deleted").BeginObject();
      for (const auto &r : k.deleted) j.Key(r.first).String(r.second);
      j.EndObject();
    }
    j.EndObject();
  }
  j.EndObject();

  j.Key("quantities").BeginObject();
  for (const Quantity &q : Quantities()) WriteQuantity(j, q);
  j.EndObject();
  j.Key("general_norm_rules").BeginArray();
  for (const NormRule &r : GeneralNormRules()) {
    j.BeginObject();
    WriteCondition(j.Key("when"), r.when);
    j.Key("modes").StringArray(r.modes);
    j.Key("prior").Bool(r.prior);
    j.Key("message").String(r.message);
    j.EndObject();
  }
  j.EndArray();

  j.Key("units").BeginObject();
  for (const Unit &u : Units()) {
    j.Key(u.token).BeginObject();
    j.Key("dimension").String(u.dimension);
    j.Key("factor").Number(u.factor);
    j.Key("engine_unit").String(FindDimension(u.dimension)->base);
    if (!u.aliases.empty()) j.Key("aliases").StringArray(u.aliases);
    j.EndObject();
  }
  j.EndObject();
  j.Key("dimensions").BeginObject();
  for (const Dimension &d : Dimensions()) {
    j.Key(d.name).BeginObject();
    j.Key("engine_unit").String(d.base);
    j.Key("doc").String(d.doc);
    j.EndObject();
  }
  j.EndObject();

  j.Key("column_roles").BeginObject();
  for (const auto &r : AxisRoles()) j.Key(r.first).BeginObject().Key("unit").String(r.second).EndObject();
  for (const std::string &r : PlainRoles()) j.Key(r).BeginObject().Key("unit").Null().EndObject();
  j.EndObject();

  j.Key("lists").BeginObject();
  j.Key("frames").StringArray({"lab", "cm"});
  j.Key("nuisance_modes").StringArray({"fixed", "fit", "profiled"});
  j.Key("light_particles").StringArray({"n", "p", "d", "t", "h", "a", "g", "e-", "b-", "b+"});
  j.Key("grid_kinds").StringArray({"range", "list", "file", "adaptive", "none"});
  std::vector<std::string> elements;
  for (int z = 0; z <= physics::kMaxElementZ; z++) elements.push_back(physics::ElementSymbol(z));
  j.Key("elements").StringArray(elements);
  j.EndObject();

  j.Key("diagnostic_codes").BeginObject();
  for (const DiagnosticCode &c : DiagnosticCodes()) {
    j.Key(c.code).BeginObject();
    j.Key("severity").String(SeverityName(c.severity));
    j.Key("meaning").String(c.meaning);
    j.EndObject();
  }
  j.EndObject();
  j.EndObject();
  return j.Str();
}

}  // namespace experiment
}  // namespace azure2
