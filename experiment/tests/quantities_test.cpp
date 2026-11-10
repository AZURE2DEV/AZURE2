// Every quantity row of the table (docs/upgrade/project_format3.md section 7.2) exercised through the validator:
//  1. each row has at least one valid record here (a new row without one fails this test), and each validates;
//  2. units: every unit token of the unit table is tried on every record; the ones the row allows (its unit rules,
//     their conditions evaluated here independently of the validator) give no diagnostic, every other exactly
//     E_UNIT_NOT_ALLOWED; per-charge units need charge_state, count units also the efficiency;
//  3. required fields: each field the row requires, removed, gives exactly E_FIELD_REQUIRED naming it;
//  4. fields the row does not list are refused with E_FIELD_NOT_APPLICABLE.
//   experiment_quantities_test
#include <algorithm>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "azure2/experiment/Schema.h"
#include "azure2/io/Document.h"
#include "azure2/io/RecordReader.h"

using namespace azure2::experiment;
using azure2::io::Document;

namespace {

int failures = 0;
int checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

std::vector<Diagnostic> Diagnose(const std::string &text) {
  std::vector<Diagnostic> out;
  const Document doc = Document::Parse(text);
  std::vector<SectionRecords> sections;
  for (const Document::Section &s : doc.Sections()) {
    if (s.depth != 0 || Schema::KindsOfSection(s.name).empty()) continue;
    const azure2::io::RecordReadResult r = azure2::io::ReadRecords(doc, s, Schema::KindsOfSection(s.name));
    for (const azure2::io::RecordDiagnostic &d : r.diagnostics) out.push_back(FromRecordDiagnostic(d, s.name));
    sections.push_back({s.name, r.set.records});
  }
  const ValidationResult v = Schema::Validate(sections);
  // Errors only: a quantity's data-set warning (W_NOT_MEASURABLE) is tested in diagnostics_test.
  for (const Diagnostic &d : v.diagnostics)
    if (d.severity == Severity::Error) out.push_back(d);
  return out;
}

const std::string kApparatus =
    "<targets>\n"
    "target t1\n    kind layer\n    active 13C\n    areal_density 1e18 atoms/cm2 fixed\n"
    "    stopping formula \"a0\" a0=1e-20 MeV*cm2/atom\n    straggling none\n"
    "target ts\n    kind sample\n    active 20Ne\n    areal_density 0.01 atoms/b fixed\n"
    "beam b1\n    spread 1 keV lab\n"
    "detector d1\n    attenuation 1 0.99\n"
    "</targets>\n";

// One record of a quantity, with what the unit conditions of its row depend on.
struct Variant {
  Variant(std::string k, std::vector<std::pair<std::string, std::string>> f, bool r = false, std::string e = "", bool a = false)
      : kind(std::move(k)), fields(std::move(f)), ratio(r), expansion(std::move(e)), angle(a) {}
  std::string kind;                                         // dataset or datum
  std::vector<std::pair<std::string, std::string>> fields;  // key, value
  bool ratio;
  std::string expansion;
  bool angle;
};

std::string Text(const Variant &v, const std::string &skip = "", const std::string &unit = "",
                 const std::vector<std::pair<std::string, std::string>> &extra = {}) {
  std::string body = v.kind + " r1\n";
  for (const auto &f : v.fields) {
    if (f.first == skip) continue;
    std::string value = f.second;
    if (!unit.empty() && f.first == "unit") value = unit;
    if (!unit.empty() && (f.first == "value" || f.first == "error")) value = value.substr(0, value.find(' ')) + " " + unit;
    body += "    " + f.first + " " + value + "\n";
  }
  for (const auto &f : extra) body += "    " + f.first + " " + f.second + "\n";
  const std::string section = v.kind == "datum" ? "modelData" : "data";
  return kApparatus + "<" + section + ">\n" + body + "</" + section + ">\n";
}

std::vector<std::pair<std::string, std::string>> DS(const std::string &reaction, const std::string &quantity,
                                                    const std::string &unit, const std::string &columns,
                                                    std::vector<std::pair<std::string, std::string>> more = {},
                                                    const std::string &norm = "1 fit") {
  std::vector<std::pair<std::string, std::string>> f = {{"reaction", reaction}, {"quantity", quantity}, {"unit", unit},
                                                        {"file", "f.dat"},      {"columns", columns},   {"norm", norm}};
  f.insert(f.end(), more.begin(), more.end());
  return f;
}

bool Holds(const Condition &c, const Variant &v, const Unit &u) {
  bool r = true;
  switch (c.kind) {
    case Condition::Always: r = true; break;
    case Condition::RatioReaction: r = v.ratio; break;
    case Condition::AngleAxis: r = v.angle; break;
    case Condition::FieldIs:
      r = c.field == "expansion" && std::find(c.values.begin(), c.values.end(), v.expansion) != c.values.end();
      break;
    case Condition::UnitDimension:
      r = std::find(c.values.begin(), c.values.end(), u.dimension) != c.values.end();
      break;
    default:
      check(std::string("condition kind ") + ConditionKindName(c.kind) + " is known to this test", false);
  }
  return r != c.negate;
}

}  // namespace

int main() {
  const std::string M = "E_lab[MeV] value error", MA = "E_lab[MeV] theta_lab[deg] value error";
  std::map<std::string, std::vector<Variant>> variants = {
      {"sigma",
       {{"dataset", DS("15N(p,a)12C", "sigma", "b", M)},
        {"dataset", DS("15N(p,a)12C[4.439] / 15N(p,a)12C", "sigma", "1", M), true}}},
      {"sigma_nuclear", {{"dataset", DS("12C(p,p)12C", "sigma_nuclear", "b", M)}}},
      {"sigma_total", {{"dataset", DS("20Ne(n,n)20Ne", "sigma_total", "b", M)}}},
      {"dsigma_domega",
       {{"dataset", DS("15N(p,a)12C", "dsigma_domega", "b/sr", MA), false, "", true},
        {"dataset", DS("12C(p,p)12C / 12C(p,p)12C @ 90 deg cm", "dsigma_domega", "1", MA), true, "", true}}},
      {"sfactor", {{"dataset", DS("12C(p,g)13N", "sfactor", "keV*b", "E_cm[keV] value error")}}},
      {"ratio_to_rutherford", {{"dataset", DS("12C(a,a)12C", "ratio_to_rutherford", "1", MA), false, "", true}}},
      {"ratio_to_mott", {{"dataset", DS("4He(a,a)4He", "ratio_to_mott", "1", MA), false, "", true}}},
      {"legendre", {}},
      {"analyzing_power", {{"dataset", DS("12C(p,p)12C", "analyzing_power", "1", MA), false, "", true}}},
      {"phase_shift", {{"dataset", DS("12C(p,p)12C", "phase_shift", "deg", M, {{"J", "1/2"}, {"l", "0"}, {"s", "1/2"}})}}},
      {"yield",
       {{"dataset", DS("13C(p,g)14N", "yield", "1/particle", M, {{"target", "t1"}, {"quadrature", "points 20"}})},
        {"dataset", DS("13C(p,g)14N", "yield", "1/particle/sr", MA, {{"target", "t1"}, {"quadrature", "points 20"}}), false, "",
         true}}},
      {"transmission", {{"dataset", DS("20Ne(n,n)20Ne", "transmission", "1", M, {{"target", "ts"}})}}},
      {"neutron_yield", {{"dataset", DS("20Ne(n,g)21Ne", "neutron_yield", "1", M, {{"target", "ts"}})}}},
      {"beta_delayed_spectrum", {{"dataset", DS("16N(b-,a)12C", "beta_delayed_spectrum", "arbitrary", "E_decay[keV] value error")}}},
      {"resonance_strength",
       {{"datum", {{"quantity", "resonance_strength"}, {"level", "5"}, {"J", "3-"}, {"reaction", "15N(p,a)12C"},
                   {"value", "2.63 keV"}, {"error", "0.15 keV"}}}}},
      {"resonance_energy",
       {{"datum", {{"quantity", "resonance_energy"}, {"level", "5"}, {"J", "3-"}, {"frame", "E_cm"}, {"pair", "p+15N"},
                   {"value", "1.0 MeV"}, {"error", "1 MeV"}}}}},
      {"total_width", {{"datum", {{"quantity", "total_width"}, {"level", "5"}, {"J", "3-"}, {"value", "14 keV"}, {"error", "1 keV"}}}}},
      {"partial_width",
       {{"datum", {{"quantity", "partial_width"}, {"level", "5"}, {"J", "3-"}, {"pair", "a+12C"}, {"value", "12 keV"},
                   {"error", "1 keV"}}}}},
      {"branching_ratio",
       {{"datum", {{"quantity", "branching_ratio"}, {"level", "5"}, {"J", "3-"}, {"reaction", "15N(p,a)12C"}, {"value", "88 %"},
                   {"error", "1 %"}}}}},
      {"anc",
       {{"datum", {{"quantity", "anc"}, {"level", "1"}, {"J", "1/2-"}, {"pair", "p+12C"}, {"l", "1"}, {"s", "1/2"},
                   {"value", "1.6 fm^-1/2"}, {"error", "0.1 fm^-1/2"}}}}},
      {"scattering_length",
       {{"datum", {{"quantity", "scattering_length"}, {"pair", "n+20Ne"}, {"s", "1/2"}, {"value", "4.6 fm"}, {"error", "0.1 fm"}}}}},
  };
  for (const char *e : {"a_over_a0", "rsl", "absolute", "2l2", "4pi"}) {
    const bool normalised = std::string(e) == "a_over_a0" || std::string(e) == "rsl";
    const std::string unit = normalised ? "1" : std::string(e) == "4pi" ? "b" : "b/sr";
    variants["legendre"].push_back({"dataset", DS("12C(a,a)12C", "legendre", unit, "E_cm[MeV] order value error",
                                                  {{"expansion", e}}, normalised ? "1 fixed" : "1 fit"),
                                    false, e});
  }

  // Plausible values of every quantity-dependent field, to try them where the row does not list them.
  const std::map<std::string, std::string> plausible = {
      {"target", "t1"}, {"beam", "b1"}, {"detector", "d1"}, {"J", "1/2"}, {"l", "0"}, {"s", "1/2"},
      {"expansion", "absolute"}, {"multipole", "E1"}, {"secondary_gamma", "L 1 Ic 1 delta 0"}, {"charge_state", "1"},
      {"efficiency", "0.1 fixed"}, {"thm", "global"}, {"level", "1"}, {"reaction", "15N(p,a)12C"}, {"frame", "E_x"},
      {"pair", "p+15N"}, {"imag", "1 fm"}, {"imag_error", "1 fm"}, {"orders", "0 2"}};

  for (const Quantity &q : Quantities()) {
    const auto it = variants.find(q.name);
    check(q.name + ": has a record in this test", it != variants.end() && !it->second.empty());
    if (it == variants.end()) continue;
    for (const Variant &v : it->second) {
      const std::string label = q.name + (v.expansion.empty() ? "" : " " + v.expansion) + (v.ratio ? " (ratio)" : "") +
                                (v.angle ? " (angle)" : "");
      const std::vector<Diagnostic> base = Diagnose(Text(v));
      check(label + ": validates", base.empty());
      for (const Diagnostic &d : base) std::printf("        %s\n", d.ToText().c_str());

      // 2. Units.
      for (const Unit &u : Units()) {
        bool allowed = false;
        for (const UnitRule &r : q.units) {
          const bool matches = r.dimension.empty() ? r.unit == u.token : r.dimension == u.dimension;
          allowed = allowed || (matches && Holds(r.when, v, u));
        }
        std::vector<std::pair<std::string, std::string>> extra;
        const std::string dim = u.dimension;
        if (allowed && (dim == "per_charge" || dim == "counts_per_charge")) {
          const std::vector<Diagnostic> without = Diagnose(Text(v, "", u.token));
          check(label + " " + u.token + ": charge_state required",
                !without.empty() && without[0].code == std::string("E_FIELD_REQUIRED"));
          extra.push_back({"charge_state", "1"});
          if (dim == "counts_per_charge") extra.push_back({"efficiency", "0.02 fit prior 0.02 3%"});
        }
        const std::vector<Diagnostic> d = Diagnose(Text(v, "", u.token, extra));
        if (allowed) {
          check(label + ": unit " + u.token + " allowed", d.empty());
        } else {
          check(label + ": unit " + u.token + " refused with E_UNIT_NOT_ALLOWED",
                !d.empty() && d[0].code == std::string("E_UNIT_NOT_ALLOWED"));
          for (const Diagnostic &x : d)
            if (x.code != std::string("E_UNIT_NOT_ALLOWED") && x.code != std::string("E_MODE_NOT_ALLOWED"))
              check(label + ": unit " + u.token + " gives only unit problems (" + x.code + ")", false);
        }
      }

      // 3. Required fields.
      for (const QuantityField &f : q.fields) {
        if (!f.required || f.when.kind != Condition::Always) continue;
        bool present = false;
        for (const auto &x : v.fields) present = present || x.first == f.field;
        if (!present) continue;
        const std::vector<Diagnostic> d = Diagnose(Text(v, f.field));
        const bool one = d.size() == 1 && d[0].code == std::string("E_FIELD_REQUIRED") && d[0].field == f.field;
        check(label + ": without " + f.field + " exactly E_FIELD_REQUIRED", one);
        if (!one)
          for (const Diagnostic &x : d) std::printf("        %s\n", x.ToText().c_str());
      }

      // 4. Fields the row does not list.
      const RecordKind *kind = Schema::FindKind(v.kind);
      for (const FieldDef &f : kind->fields) {
        if (!f.byQuantity) continue;
        bool listed = false;
        for (const QuantityField &qf : q.fields) listed = listed || qf.field == f.key;
        if (listed) continue;
        const auto p = plausible.find(f.key);
        check(f.key + ": a plausible value in this test", p != plausible.end());
        if (p == plausible.end()) continue;
        const std::vector<Diagnostic> d = Diagnose(Text(v, "", "", {{f.key, p->second}}));
        bool refused = false;
        for (const Diagnostic &x : d) refused = refused || (x.code == std::string("E_FIELD_NOT_APPLICABLE") && x.field == f.key);
        check(label + ": " + f.key + " refused", refused);
      }
    }
  }
  std::printf("%d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
