// The table of project format 3 as a front end's forms read it (Engine.h): the record kinds and their fields in the
// schema's order, the sections, the quantities and their units, the diagnostic codes.  Read from experiment/'s table
// (Schema::Kinds, Quantities, DiagnosticCodes), the one AZURE2 --schema writes; nothing is spelled here.
#include "azure2/front/Engine.h"

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

std::vector<KindInfo> RecordKinds() {
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
      info.fields.push_back(x);
    }
    out.push_back(info);
  }
  return out;
}

std::vector<std::pair<std::string, std::vector<std::string>>> RecordSections() { return ex::Schema::Sections(); }

std::vector<std::string> QuantityNames(const std::string &kind) { return ex::QuantityNames(kind); }

std::vector<std::string> QuantityUnits(const std::string &quantity) {
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
  std::vector<std::pair<std::string, std::string>> out;
  for (const ex::DiagnosticCode &c : ex::DiagnosticCodes()) out.push_back({c.code, c.meaning});
  return out;
}

double SingleParticleWidth(int z1, int z2, double mu_u, int l, double a_fm, double E_MeV) {
  azure2::physics::ChannelKinematics kin;
  kin.z1 = z1;
  kin.z2 = z2;
  kin.mu_u = mu_u;
  return azure2::physics::single_particle_width(kin, l, a_fm, E_MeV);
}

}  // namespace front
}  // namespace azure2
