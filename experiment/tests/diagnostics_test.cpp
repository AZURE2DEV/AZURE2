// One case per diagnostic code (docs/upgrade/project_format3.md section 6.2, check 24): each project text below
// differs from a valid one in one place and must give exactly one diagnostic, with that code.  The list covers
// every code of DiagnosticCodes() (checked here, and by repo_consistency check 24 by name).  Also: the syntax-only
// fixes (nearest key, '?' placeholder, unit spelling, nearest name), the JSON form of a diagnostic, and that a
// text with many problems reports all of them.
//   experiment_diagnostics_test
#include <cstdio>
#include <set>
#include <string>
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
  out.insert(out.end(), v.diagnostics.begin(), v.diagnostics.end());
  return out;
}

// A valid sigma data set, a target layer, an infinite target, a beam, and a model datum.
const std::string kDataset =
    "dataset d1\n"
    "    reaction       15N(p,a)12C\n"
    "    quantity       sigma\n"
    "    unit           b\n"
    "    file           f.dat\n"
    "    columns        E_lab[MeV] value error\n"
    "    norm           1 fixed\n";
const std::string kTargets =
    "<targets>\n"
    "target t1\n"
    "    kind           layer\n"
    "    active         15N\n"
    "    areal_density  1e18 atoms/cm2 fixed\n"
    "    stopping       formula \"a0\" a0=1e-20 MeV*cm2/atom\n"
    "    straggling     none\n"
    "target tinf\n"
    "    kind           infinite\n"
    "    active         15N\n"
    "    stopping       formula \"a0\" a0=1e-20 MeV*cm2/atom\n"
    "    straggling     none\n"
    "beam b1\n"
    "    spread         1 keV lab\n"
    "</targets>\n";
const std::string kDatum =
    "datum x1\n"
    "    quantity       resonance_strength\n"
    "    level          5\n"
    "    J              3-\n"
    "    reaction       15N(p,a)12C\n"
    "    value          2.63 keV\n"
    "    error          0.15 keV\n";

// `base` with the line of field `key` replaced by `line` ("" removes it).
std::string Set(const std::string &base, const std::string &key, const std::string &line) {
  std::string out;
  std::size_t a = 0;
  bool done = false;
  while (a < base.size()) {
    const std::size_t b = base.find('\n', a);
    const std::string l = base.substr(a, b - a);
    if (!done && l.compare(0, 4 + key.size() + 1, "    " + key + " ") == 0) {
      if (!line.empty()) out += line + "\n";
      done = true;
    } else {
      out += l + "\n";
    }
    a = b + 1;
  }
  return out;
}
std::string Data(const std::string &records) { return "<data>\n" + records + "</data>\n"; }
std::string Folded(const std::string &dataset) { return dataset + "    beam           b1\n    quadrature     points 20\n"; }
std::string Model(const std::string &records) { return "<modelData>\n" + records + "</modelData>\n"; }

struct Case {
  const char *code;
  std::string text;
};

}  // namespace

int main() {
  check("the valid base texts give no diagnostic",
        Diagnose(Data(kDataset)).empty() && Diagnose(kTargets + Data(Folded(kDataset))).empty() &&
            Diagnose(Model(kDatum)).empty());

  // Constructs the old engine reads and the table took from step 3a on (the converter's survey of the test tree):
  // each validates with no diagnostic.
  check("thm global with norm profiled (ESegment::IsProfiledNorm profiles every fitted THM norm)",
        Diagnose(Data(Set(Set(kDataset, "norm", "    norm profiled"), "unit", "    unit arbitrary") + "    thm global\n")).empty());
  check("legendre data folded with a beam",
        Diagnose(kTargets + Data(Folded(Set(Set(Set(Set(kDataset, "quantity", "    quantity legendre"), "unit", "    unit 1"),
                                                "columns", "    columns E_lab[MeV] order value error"),
                                            "reaction", "    reaction 15N(p,p)15N") +
                                        "    expansion a_over_a0\n")))
            .empty());
  check("dsigma_domega with the window columns of a beam profile",
        Diagnose("<targets>\nbeam bp\n    profile skewed_gaussian xi 0.3 MeV omega 0.03 MeV alpha 1 weight 1\n"
                 "    truncation 0\n    detailed_balance no\n</targets>\n" +
                 Data(Set(Set(Set(kDataset, "quantity", "    quantity dsigma_domega"), "unit", "    unit b/sr"), "columns",
                          "    columns E_lab[MeV] theta_lab[deg] value error window_lo window_hi") +
                      "    beam bp\n    quadrature points 20\n"))
            .empty());
  check("beta_delayed_spectrum with a fixed norm (the feeding amplitudes carry the scale)",
        Diagnose(Data(Set(Set(Set(Set(kDataset, "quantity", "    quantity beta_delayed_spectrum"), "unit", "    unit arbitrary"),
                              "columns", "    columns E_decay[MeV] value error"),
                          "reaction", "    reaction 16N(b-,a)12C")))
            .empty());
  check("scattering_length of a charged pair (Coulomb-modified)",
        Diagnose(Model("datum a1\n    quantity scattering_length\n    pair p+7Be\n    s 2\n    value -7 fm\n    error 1 fm\n")).empty());

  const std::vector<Case> cases = {
      // The reader
      {"E_UNTERMINATED_QUOTE", Data(Set(kDataset, "file", "    file \"f.dat"))},
      {"E_QUOTE_IN_TOKEN", Data(Set(kDataset, "file", "    file f\"x"))},
      {"E_UNKNOWN_KIND", Data("datset d0\n" + kDataset)},
      {"E_UNKNOWN_KIND", Data(kDataset + "datset d2\n")},  // read as a field of d1; the validator says what it is
      {"E_HEADER_NAME", Data("dataset d1 extra\n" + kDataset.substr(kDataset.find('\n') + 1))},
      {"E_KEY_SYNTAX", Data(kDataset + "    9x 1\n")},
      // Records and names
      {"E_BAD_NAME", Data("dataset 1d\n" + kDataset.substr(kDataset.find('\n') + 1))},
      {"E_DUPLICATE_NAME", Data(kDataset + "normalization d1\n    value 1 fixed\n")},
      // Fields
      {"E_UNKNOWN_FIELD", Data(kDataset + "    taget t1\n")},
      {"E_FIELD_TWICE", Data(kDataset + "    unit b\n")},
      {"E_FIELD_REQUIRED", Data(Set(kDataset, "file", ""))},
      {"E_FIELD_NOT_APPLICABLE", Data(kDataset + "    charge_state 1\n")},
      {"E_VALUE_MISSING", Data(kDataset + "    covariance\n")},
      {"E_PLACEHOLDER", Data(kDataset + "    covariance ?\n")},
      // Values
      {"E_TOO_FEW_TOKENS", Data(kDataset + "    energy_window 1\n")},
      {"E_TOO_MANY_TOKENS", Data(Set(kDataset, "unit", "    unit b extra"))},
      {"E_NOT_A_NUMBER", Data(kDataset + "    energy_window a 2 MeV\n")},
      {"E_NOT_AN_INTEGER", Model(Set(kDatum, "level", "    level 5.5"))},
      {"E_OUT_OF_RANGE", Model(Set(kDatum, "error", "    error -0.15 keV"))},
      {"E_EMPTY_RANGE", Data(kDataset + "    energy_window 2 1 MeV\n")},
      {"E_NOT_IN_LIST", Data(Set(kDataset, "quantity", "    quantity sigmaa"))},
      {"E_UNKNOWN_UNIT", Data(Set(kDataset, "unit", "    unit MB"))},
      {"E_UNIT_MISSING", Data(kDataset + "    energy_window 1 2\n")},
      {"E_UNIT_NOT_ALLOWED", Data(kDataset + "    energy_window 1 2 deg\n")},
      {"E_UNIT_NOT_ALLOWED", Data(Set(kDataset, "unit", "    unit b/sr"))},
      {"E_NUISANCE_MODE", Data(Set(kDataset, "norm", "    norm 1"))},
      {"E_MODE_NOT_ALLOWED", Data(Set(kDataset, "norm", "    norm profiled"))},
      {"E_MODE_NOT_ALLOWED", Data(Set(Set(kDataset, "unit", "    unit arbitrary"), "norm", "    norm 1 fit prior 1 5%"))},
      {"E_PRIOR_SYNTAX", Data(Set(kDataset, "norm", "    norm 1 fit prior 1"))},
      {"E_PART_MISSING", "<targets>\nbeam b1\n    profile skewed_gaussian xi 0.3 MeV omega 0.03 MeV alpha 1\n"
                         "    truncation 0\n    detailed_balance no\n</targets>\n"},
      {"E_PART_TWICE", "<targets>\nbeam b1\n    profile skewed_gaussian xi 0.3 MeV xi 0.3 MeV omega 0.03 MeV alpha 1 weight 1\n"
                       "    truncation 0\n    detailed_balance no\n</targets>\n"},
      {"E_FORMULA_SYNTAX", "<targets>\nbeam b1\n    spread formula a0*x a0=1 MeV lab\n</targets>\n"},
      {"E_SPIN_SYNTAX", Model(Set(kDatum, "J", "    J 3"))},
      {"E_NUCLEUS_SYNTAX", Set(kTargets, "active", "    active C13")},
      {"E_UNKNOWN_ELEMENT", Set(kTargets, "active", "    active 13Xx")},
      {"E_PAIR_SYNTAX", Model("datum a1\n    quantity scattering_length\n    pair n20Ne\n    s coherent\n"
                              "    value 4.6 fm\n    error 0.1 fm\n")},
      {"E_WINDOW_SYNTAX", kTargets + Data(Folded(kDataset) + "    fold_windows wide MeV\n")},
      // Reactions
      {"E_REACTION_SYNTAX", Data(Set(kDataset, "reaction", "    reaction 15N(p,a"))},
      {"E_REACTION_ORDINAL", Data(Set(kDataset, "reaction", "    reaction 15N(p,a1)12C"))},
      {"E_REACTION_BALANCE", Data(Set(kDataset, "reaction", "    reaction 15N(p,a)13C"))},
      {"E_REACTION_NOT_ALLOWED", Data(Set(kDataset, "reaction", "    reaction 12C(p,p)12C"))},
      // Column roles
      {"E_UNKNOWN_ROLE", Data(Set(kDataset, "columns", "    columns E_lab[MeV] valu error"))},
      {"E_ROLE_UNIT", Data(Set(kDataset, "columns", "    columns E_lab value error"))},
      {"E_ROLE_TWICE", Data(Set(kDataset, "columns", "    columns E_lab[MeV] value value error"))},
      {"E_ROLE_REQUIRED", Data(Set(kDataset, "columns", "    columns E_lab[MeV] value"))},
      {"E_ROLE_NOT_ALLOWED", Data(Set(kDataset, "columns", "    columns E_lab[MeV] theta_lab[deg] value error"))},
      // E_x is no column role (step 3a): no quantity takes an excitation-energy column; 'frame E_x' of a datum is a
      // field value
      {"E_UNKNOWN_ROLE", Data(Set(kDataset, "columns", "    columns E_x[MeV] value error"))},
      // Between records
      {"E_UNKNOWN_REFERENCE", kTargets + Data(kDataset + "    beam b2\n    quadrature points 20\n")},
      {"E_REFERENCE_KIND", Data(Set(kDataset, "norm", "    norm d1"))},
      {"E_UNKNOWN_EPOCH", kTargets + Data(Folded(kDataset) + "    target t1@run9\n")},
      {"E_TARGET_KIND", kTargets + Data(Folded(kDataset) + "    target tinf\n")},
      {"E_DEGENERATE_NUISANCE",
       kTargets + Data(Set(Set(Set(Folded(kDataset), "unit", "    unit counts/uC\n    charge_state 1\n    efficiency 0.02 fit"),
                               "quantity", "    quantity yield"),
                           "norm", "    norm 1 fit") +
                       "    target t1\n")},
  };

  std::set<std::string> tested;
  for (const Case &c : cases) {
    const std::vector<Diagnostic> d = Diagnose(c.text);
    const bool one = d.size() == 1 && d[0].code == c.code;
    check(std::string(c.code) + ": exactly this one diagnostic", one);
    if (!one)
      for (const Diagnostic &x : d) std::printf("        got %s\n", x.ToText().c_str());
    if (one) {
      check(std::string(c.code) + ": a line and a message naming what is wrong", d[0].line > 0 && !d[0].message.empty());
      tested.insert(c.code);
    }
  }
  {
    // A record in the wrong section cannot be written in a file (the reader takes the section's kinds), but a front
    // end that builds records (pyazr's add_*) can put one there.
    azure2::io::Record r;
    r.kind = "normalization";
    r.name = "N";
    r.line = 3;
    r.col = 1;
    azure2::io::Field f;
    f.key = "value";
    f.line = 4;
    f.col = 5;
    f.tokens = {{"1", 11, false}, {"fixed", 13, false}};
    r.fields.push_back(f);
    const ValidationResult v = Schema::Validate({{"targets", {r}}});
    const bool one = v.diagnostics.size() == 1 && v.diagnostics[0].code == std::string("E_WRONG_SECTION");
    check("E_WRONG_SECTION: exactly this one diagnostic", one);
    if (one) tested.insert("E_WRONG_SECTION");
  }
  for (const DiagnosticCode &c : DiagnosticCodes()) check(std::string(c.code) + " has a case here", tested.count(c.code) == 1);

  // Fixes change syntax only.
  {
    const std::vector<Diagnostic> d = Diagnose(Data(kDataset + "    taget t1\n"));
    check("unknown field: fix replace_key target", d.size() == 1 && d[0].hasFix && d[0].fix.action == "replace_key" &&
                                                       d[0].fix.text == "target" && d[0].fix.line == 9);
    check("unknown field: JSON", d.size() == 1 && d[0].ToJson().find("\"code\": \"E_UNKNOWN_FIELD\"") != std::string::npos &&
                                     d[0].ToJson().find("\"fix\": {\"action\": \"replace_key\", \"line\": 9") != std::string::npos);
  }
  {
    const std::vector<Diagnostic> d = Diagnose(Data(Set(kDataset, "file", "")));
    check("missing field: fix inserts '?' after the record's last line",
          d.size() == 1 && d[0].hasFix && d[0].fix.action == "insert_field" && d[0].fix.afterLine == 7 &&
              d[0].fix.text == "    file           ?" && d[0].field == "file" && d[0].record == "d1");
  }
  {
    const std::vector<Diagnostic> d = Diagnose(Data(Set(kDataset, "unit", "    unit MB")));
    check("unit spelling: fix replace_token mb", d.size() == 1 && d[0].hasFix && d[0].fix.text == "mb" && d[0].fix.col == 10);
    const std::vector<Diagnostic> e = Diagnose(Data(Set(kDataset, "unit", "    unit MEV")));
    check("no fix where the case decides (MeV, meV)", e.size() == 1 && !e[0].hasFix);
  }
  {
    const std::vector<Diagnostic> d = Diagnose(kTargets + Data(kDataset + "    beam b2\n    quadrature points 20\n"));
    check("unknown name: fix the nearest beam", d.size() == 1 && d[0].hasFix && d[0].fix.text == "b1" &&
                                                    d[0].allowed == std::vector<std::string>{"b1"});
  }
  {
    // Every problem, not the first.
    const std::string text = Data("dataset d1\n    reaction 15N(p,a)13C\n    quantity sigma\n    unit MeV\n"
                                  "    columns E_lab value\n    norm 1\n    taget x\n    covariance ?\n");
    const std::vector<Diagnostic> d = Diagnose(text);
    std::set<std::string> codes;
    for (const Diagnostic &x : d) codes.insert(x.code);
    check("all problems reported (balance, unit, role unit, mode, unknown field, placeholder, missing file)",
          codes.count("E_REACTION_BALANCE") && codes.count("E_UNIT_NOT_ALLOWED") && codes.count("E_ROLE_UNIT") &&
              codes.count("E_NUISANCE_MODE") && codes.count("E_UNKNOWN_FIELD") && codes.count("E_PLACEHOLDER") &&
              codes.count("E_FIELD_REQUIRED") && d.size() == 7);
  }
  {
    // A data set without its quantity: one diagnostic, not one per field that depends on the quantity.
    const std::vector<Diagnostic> d = Diagnose(kTargets + Data(Set(Folded(kDataset), "quantity", "") + "    target t1\n"));
    check("missing quantity: only E_FIELD_REQUIRED", d.size() == 1 && d[0].field == "quantity");
  }
  {
    // A neutron sample is not a fold: no quadrature.
    const std::string ts = "<targets>\ntarget ts\n    kind sample\n    active 20Ne\n    areal_density 0.01 atoms/b fixed\n</targets>\n";
    check("sample: quadrature refused", Diagnose(ts + Data(Set(Set(Set(kDataset, "reaction", "    reaction 20Ne(n,n)20Ne"),
                                                                       "quantity", "    quantity transmission"), "unit", "    unit 1") +
                                                              "    target ts\n    quadrature points 20\n")).size() == 1);
  }
  {
    // Story (b) with identical particles: the note's message (section 9.2).
    const std::vector<Diagnostic> d =
        Diagnose(Data("dataset aC_RTH\n    reaction 4He(a,a)4He\n    quantity ratio_to_rutherford\n    unit 1\n"
                      "    file f.dat\n    columns E_lab[MeV] theta_lab[deg] value error\n    norm 1 fit prior 1 3%\n"));
    check("ratio_to_rutherford with identical particles: use ratio_to_mott",
          d.size() == 1 && d[0].code == std::string("E_REACTION_NOT_ALLOWED") &&
              d[0].message.find("ratio_to_mott") != std::string::npos && d[0].field == "reaction");
  }
  {
    // Story (a) without charge_state: the error of section 4.2.
    const std::string ya = kTargets + Data(Set(Set(Folded(kDataset), "unit", "    unit counts/uC\n    efficiency 0.02 fit prior 0.02 3%"),
                                               "quantity", "    quantity yield") +
                                           "    target t1\n");
    const std::vector<Diagnostic> d = Diagnose(ya);
    check("per-charge unit without charge_state", d.size() == 1 && d[0].code == std::string("E_FIELD_REQUIRED") &&
                                                       d[0].field == "charge_state" &&
                                                       d[0].message.find("charge state") != std::string::npos);
  }

  std::printf("%d checks, %d failures (%zu cases)\n", checks, failures, cases.size());
  return failures == 0 ? 0 : 1;
}
