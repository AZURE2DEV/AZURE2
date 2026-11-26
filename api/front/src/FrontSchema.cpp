// The table of project format 3 as a front end's forms read it (Engine.h): the record kinds and their fields in the
// schema's order, the sections, the quantities and their units, the diagnostic codes.  Read from experiment/'s table
// (Schema::Kinds, Quantities, DiagnosticCodes), the one AZURE2 --schema writes; nothing is spelled here.
#include "azure2/front/Engine.h"
#include "azure2/support/ClassicNumeric.h"

#include <algorithm>
#include <climits>
#include <set>

#include "azure2/experiment/Conditions.h"
#include "azure2/experiment/Diagnostics.h"
#include "azure2/experiment/Quantities.h"
#include "azure2/experiment/Schema.h"
#include "azure2/experiment/Units.h"
#include "azure2/experiment/Values.h"
#include "azure2/physics/Widths.h"

namespace azure2 {
namespace front {

namespace ex = azure2::experiment;

static ValueInfo InfoOfSpec(const ex::ValueSpec &spec) {
  ValueInfo v;
  v.type = ex::ValueTypeName(spec.type);
  v.grammar = ex::Grammar(spec);
  v.dimension = spec.dimension;
  if (spec.dimension == "@quantity")
    v.unitOfQuantity = true;
  else if (!spec.dimension.empty())
    v.units = ex::UnitsOfDimension(spec.dimension);
  if (spec.type == ex::ValueType::Enum || spec.type == ex::ValueType::Reference) v.values = spec.values;
  v.positive = spec.positive;
  v.hasMin = spec.type == ex::ValueType::Integer && spec.min != LONG_MIN;
  v.min = v.hasMin ? spec.min : 0;
  v.epoch = spec.epoch;
  v.profiled = spec.profiled;
  v.none = spec.none;
  v.refKinds = spec.refKinds;
  if (spec.form)
    for (const ex::Alternative &a : spec.form->alternatives) {
      AlternativeInfo alt;
      alt.keyword = a.keyword;
      alt.doc = a.doc;
      for (const ex::PartSpec &p : a.parts) {
        PartInfo part;
        part.name = p.name;
        part.keyword = p.keyword;
        part.optional = p.optional;
        part.defaultText = p.defaultText;
        part.doc = p.doc;
        part.value = InfoOfSpec(p.spec);
        alt.parts.push_back(part);
      }
      v.alternatives.push_back(alt);
    }
  return v;
}

namespace {

/// The field a field belongs with on a form: the one its requirement names, or the one whose presence it needs, or
/// the structured field whose alternative it depends on (FieldInfo::group).
std::string GroupOf(const ex::FieldDef &f, const ex::RecordKind &k) {
  for (const ex::Condition *c : {&f.required, &f.applies}) {
    if (c->negate || c->field.empty() || c->field == f.key) continue;
    if (c->kind == ex::Condition::FieldPresent) return c->field;
    if (c->kind == ex::Condition::FieldIs) {
      const ex::FieldDef *named = k.Find(c->field);
      if (named && named->spec.type == ex::ValueType::Form) return c->field;
    }
  }
  return std::string();
}

/// The columns of a table of records (KindInfo::summaryFields): every field a user sets often, but free text and the
/// column roles; of a kind with many such fields, only the required ones, the references to other records, the
/// windows and the grids (a field required only under a condition is a column of a small kind only).  Required first.
void SetSummary(KindInfo &info, const ex::RecordKind &k) {
  std::vector<const ex::FieldDef *> candidates;
  for (const ex::FieldDef &f : k.fields) {
    if (f.rare || f.key == info.includeField) continue;
    if (f.spec.type == ex::ValueType::Text || f.spec.type == ex::ValueType::Columns) continue;
    if (f.byQuantity && f.spec.type != ex::ValueType::Reference) continue;
    candidates.push_back(&f);
  }
  const bool few = candidates.size() <= 8;
  std::vector<std::string> required, others;
  for (const ex::FieldDef *f : candidates) {
    const bool always = f->required.kind == ex::Condition::Always && !f->required.negate && !f->byQuantity &&
                        f->applies.kind == ex::Condition::Always && !f->applies.negate;
    const bool kept = few || always || f->spec.type == ex::ValueType::Reference ||
                      f->spec.type == ex::ValueType::Range || f->spec.type == ex::ValueType::Grid;
    if (!kept) continue;
    (always ? required : others).push_back(f->key);
  }
  info.summaryFields = required;
  info.summaryFields.insert(info.summaryFields.end(), others.begin(), others.end());
  for (FieldInfo &x : info.fields)
    x.summary = std::find(info.summaryFields.begin(), info.summaryFields.end(), x.key) != info.summaryFields.end();
}

}  // namespace

std::vector<KindInfo> RecordKinds() {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::vector<KindInfo> out;
  for (const ex::RecordKind &k : ex::Schema::Kinds()) {
    KindInfo info;
    info.kind = k.kind;
    info.title = k.title;
    info.doc = k.doc;
    info.sections = k.sections;
    for (const ex::FieldDef &f : k.fields) {
      FieldInfo x;
      x.key = f.key;
      x.title = f.title;
      x.type = ex::ValueTypeName(f.spec.type);
      x.grammar = ex::Grammar(f.spec);
      x.doc = f.doc;
      x.byQuantity = f.byQuantity;
      if (!f.byQuantity) {
        if (f.applies.kind != ex::Condition::Always || f.applies.negate) x.appliesWhen = ex::Describe(f.applies);
        if (f.required.kind == ex::Condition::Always)
          x.required = !f.required.negate;
        else
          x.requiredWhen = ex::Describe(f.required);
      }
      x.rare = f.rare;
      x.repeatable = f.repeatable;
      x.defaultText = f.defaultText;
      if (f.spec.type == ex::ValueType::Enum || f.spec.type == ex::ValueType::Reference) x.values = f.spec.values;
      if (f.spec.dimension == "@quantity")
        x.unitOfQuantity = true;
      else if (!f.spec.dimension.empty())
        x.units = ex::UnitsOfDimension(f.spec.dimension);
      x.value = InfoOfSpec(f.spec);
      x.group = GroupOf(f, k);
      // the on/off switch of a record: an optional yes/no word whose default is its first word
      if (info.includeField.empty() && f.spec.type == ex::ValueType::Enum && f.spec.values.size() == 2 &&
          f.required.kind == ex::Condition::Always && f.required.negate && f.applies.kind == ex::Condition::Always &&
          f.defaultText == f.spec.values[0]) {
        info.includeField = f.key;
        info.includeOn = f.spec.values[0];
        info.includeOff = f.spec.values[1];
      }
      info.fields.push_back(x);
    }
    info.oneOf = k.oneOf;
    SetSummary(info, k);
    out.push_back(info);
  }
  return out;
}

std::vector<std::pair<std::string, std::vector<std::string>>> RecordSections() { return ex::Schema::Sections(); }

std::vector<std::string> QuantityNames(const std::string &kind) { return ex::QuantityNames(kind); }

std::vector<QuantityInfo> Quantities(const std::string &kind) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::vector<QuantityInfo> out;
  for (const std::string &name : ex::QuantityNames(kind)) {
    const ex::Quantity *q = ex::FindQuantity(name);
    if (!q) continue;
    out.push_back({q->name, q->title, q->doc, QuantityUnits(name)});
  }
  return out;
}

std::vector<std::string> QuantityUnits(const std::string &quantity) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::vector<std::string> out;
  const ex::Quantity *q = ex::FindQuantity(quantity);
  if (!q) return out;
  std::set<std::string> seen;
  for (const ex::UnitRule &rule : q->units) {
    std::vector<std::string> units;
    if (!rule.unit.empty()) units.push_back(rule.unit);
    if (!rule.dimension.empty())
      for (const std::string &u : ex::UnitsOfDimension(rule.dimension)) units.push_back(u);
    for (const std::string &u : units)
      if (seen.insert(u).second) out.push_back(u);
  }
  return out;
}

std::vector<std::pair<std::string, std::string>> DiagnosticCodes() {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  std::vector<std::pair<std::string, std::string>> out;
  for (const ex::DiagnosticCode &c : ex::DiagnosticCodes()) out.push_back({c.code, c.meaning});
  return out;
}

double SingleParticleWidth(int z1, int z2, double mu_u, int l, double a_fm, double E_MeV) {
  const support::ClassicNumericScope numeric;  // the numbers of files in the "C" locale (ClassicNumeric.h)
  azure2::physics::ChannelKinematics kin;
  kin.z1 = z1;
  kin.z2 = z2;
  kin.mu_u = mu_u;
  return azure2::physics::single_particle_width(kin, l, a_fm, E_MeV);
}

}  // namespace front
}  // namespace azure2
