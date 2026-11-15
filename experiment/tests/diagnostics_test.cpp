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

#include "azure2/experiment/GammaDecay.h"
#include "azure2/experiment/LayerChecks.h"
#include "azure2/experiment/ModelDatum.h"
#include "azure2/experiment/ParameterFileNames.h"
#include "azure2/experiment/Schema.h"
#include "azure2/experiment/StoppingChecks.h"
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
  check("thm global with norm profiled (Segment::IsProfiledNorm profiles every fitted THM norm)",
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

  // step 3a, director's D1-D3: sigma_nuclear in an extrapolation (no warning); a pair key on one reaction; a target state
  check("sigma_nuclear in an extrapolation validates",
        Diagnose("<extrapolations>\nextrapolation x1\n    reaction 12C(p,p)12C\n    quantity sigma_nuclear\n    unit b\n"
                 "    energies range 0.5 1 MeV lab step 0.1\n</extrapolations>\n")
            .empty());
  check("entrance_pair and exit_pair on one reaction validate",
        Diagnose(Data(kDataset + "    entrance_pair 1\n    exit_pair 2\n")).empty());
  check("entrance_pair and exit_pair on a sum validate (the build checks them against the terms' labels)",
        Diagnose(Data(Set(kDataset, "reaction", "    reaction 15N(p,a)12C + 15N(p,a)12C[4.439]") +
                      "    entrance_pair 1\n    exit_pair 2\n")).empty());
  check("a target state: 16O[2.143](n,n)16O[2.143] is elastic (sigma_total validates)",
        Diagnose(Data(Set(Set(kDataset, "reaction", "    reaction 16O[2.143](n,n)16O[2.143]"), "quantity", "    quantity sigma_total")))
            .empty());
  {
    const std::vector<Diagnostic> d =
        Diagnose(Data(Set(Set(kDataset, "reaction", "    reaction 16O[2.143](n,n)16O"), "quantity", "    quantity sigma_total")));
    check("a target state different from the residual's: not elastic",
          d.size() == 1 && d[0].code == std::string("E_REACTION_NOT_ALLOWED"));
  }

  {
    // a deleted field: refused, never read; the message names the field and says what to do
    const std::vector<Diagnostic> d = Diagnose(kTargets + Data(Folded(kDataset) + "    fold_windows auto\n"));
    check("fold_windows (deleted): one E_FIELD_DELETED that names the field and says to remove it",
          d.size() == 1 && d[0].code == std::string("E_FIELD_DELETED") && d[0].message.find("'fold_windows'") != std::string::npos &&
              d[0].message.find("remove it") != std::string::npos);
    if (!d.empty()) std::printf("        %s\n", d[0].ToText().c_str());
  }

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
      {"E_FIELD_DELETED", kTargets + Data(Folded(kDataset) + "    fold_windows 0.2-0.4 MeV blend 2 keV\n")},
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
      {"E_PRIOR_CENTRE", Data(Set(kDataset, "norm", "    norm 1 fit prior 0 15%"))},
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
      // a data set of a quantity that is not measurable
      {"W_NOT_MEASURABLE", Data(Set(Set(kDataset, "reaction", "    reaction 15N(p,p)15N"), "quantity", "    quantity sigma_nuclear"))},
      // a target state that is not a number
      {"E_REACTION_SYNTAX", Data(Set(kDataset, "reaction", "    reaction 16O[total](n,n)16O"))},
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
  {
    // A parameter file's names against the run's parameters (ParameterFileNames.h).
    ParameterFileNames n;
    n.ofRun = {"energy_1", "width_1_1", "d1.norm", "x1.shift.offset"};
    n.optionalOfRun = {"x1.shift.offset"};
    n.ofOtherRuns = {"d2.norm"};
    n.inFile = {"energy_1", "width_1_1", "d1.norm", "d2.norm"};
    check("every name of the run (an optional one left out, another run's name given): no diagnostic",
          MatchParameterFileNames("p.sav", n).empty());
    n.inFile = {"energy_1", "width_1_1", "segment_1_norm"};
    const std::vector<Diagnostic> partial = MatchParameterFileNames("p.sav", n);
    const bool warn = partial.size() == 1 && partial[0].code == "W_PARAMETER_FILE_NAMES" &&
                      partial[0].severity == Severity::Warning && partial[0].file == "p.sav" &&
                      partial[0].message == "2 of its 3 names applied; names this project has no parameter for (skipped): "
                                            "segment_1_norm; parameters it does not name (the project's values kept): d1.norm";
    check("W_PARAMETER_FILE_NAMES: a partial match, the names of both sides", warn);
    if (warn) tested.insert("W_PARAMETER_FILE_NAMES");
    n.inFile = {"segment_1_norm", "d2.norm"};
    const std::vector<Diagnostic> none = MatchParameterFileNames("p.sav", n);
    n.inFile = {};
    const std::vector<Diagnostic> empty = MatchParameterFileNames("p.sav", n);
    const bool refused = none.size() == 1 && none[0].code == "E_PARAMETER_FILE_NO_MATCH" &&
                         none[0].severity == Severity::Error &&
                         none[0].message.find("(segment_1_norm, d2.norm)") != std::string::npos && empty.size() == 1 &&
                         empty[0].code == "E_PARAMETER_FILE_NO_MATCH";
    check("E_PARAMETER_FILE_NO_MATCH: no name of the run in the file (or no name at all)", refused);
    if (refused) tested.insert("E_PARAMETER_FILE_NO_MATCH");
  }
  {
    // <parameterSettings> rows against the project's parameter names (ParameterFileNames.h).
    const std::vector<std::string> known = {"energy_1", "d1.norm", "x1.shift.offset"};
    check("<parameterSettings> rows that name parameters: no diagnostic",
          UnknownParameterSettings("p.azr", {{"d1.norm", 40}, {"x1.shift.offset", 41}}, known).empty());
    const std::vector<Diagnostic> u = UnknownParameterSettings("p.azr", {{"d1.norm", 40}, {"screening_1", 41}}, known);
    const bool warn = u.size() == 1 && u[0].code == "W_PARAMETER_SETTING_UNKNOWN" && u[0].severity == Severity::Warning &&
                      u[0].line == 41 && u[0].found == "screening_1" &&
                      u[0].ToText() == "p.azr:41: warning W_PARAMETER_SETTING_UNKNOWN: <parameterSettings> names "
                                       "'screening_1', which is no parameter of this project; the row is ignored";
    check("W_PARAMETER_SETTING_UNKNOWN: the row's line and name", warn);
    if (warn) tested.insert("W_PARAMETER_SETTING_UNKNOWN");
  }
  {
    // The model data against <levels> (BuildModelData): T2's pairs and the 3- levels as plain tables
    // (tests/upgrade/T2/README.md section 1); each case changes one thing of a valid datum or of the tables.
    auto pair = [](int key, int lightZ, int lightA, int heavyZ, int heavyA, double heavyE, double threshold, int twoJ1,
                   int twoJ2) {
      PairInfo p;
      p.key = key;
      p.lightZ = lightZ;
      p.lightA = lightA;
      p.heavyZ = heavyZ;
      p.heavyA = heavyA;
      p.heavyExcitation = heavyE;
      p.threshold = threshold;
      p.twoJLight = twoJ1;
      p.twoJHeavy = twoJ2;
      return p;
    };
    std::vector<PairInfo> pairs = {pair(1, 1, 1, 7, 15, 0.0, 12.1274, 1, 1), pair(2, 2, 4, 6, 12, 0.0, 7.16192, 0, 0),
                                   pair(3, 2, 4, 6, 12, 4.43891, 7.16192 + 4.43891, 0, 4)};
    const std::vector<LevelChannelInfo> channels = {{1, 2, 2}, {1, 4, 2}, {2, 3, 0}};  // no a1 channel
    std::vector<LevelInfo> levels(3);
    levels[0] = {{4}, 6, -1, 13.0734, 1e-3, 2, 1, channels};
    levels[1] = {{5}, 6, -1, 13.2652, 1e-3, 2, 2, channels};
    levels[2] = {{6}, 6, -1, 11.0, 1e-3, 2, 3, channels};  // below the p + 15N threshold
    auto build = [&](const std::string &records, const std::vector<LevelInfo> &lv, const std::vector<PairInfo> &pr) {
      const std::vector<Diagnostic> v = Diagnose(Model(records));
      if (!v.empty()) return v;
      const Document doc = Document::Parse(Model(records));
      const azure2::io::RecordReadResult r = azure2::io::ReadRecords(doc, doc.Sections()[0], Schema::KindsOfSection("modelData"));
      return BuildModelData(Schema::Validate({{"modelData", r.set.records}}), lv, pr).diagnostics;
    };
    check("model data: the valid datum builds with no diagnostic", build(kDatum, levels, pairs).empty());
    std::vector<LevelInfo> broad = levels, nearby = levels;
    broad[1].width = 0.01462;  // Gamma / E_r = 14.62 keV / 1.1378 MeV = 1.3 %
    nearby[0].energy = 13.26;  // 5.2 keV from level 5, Gamma = 1 keV
    std::vector<PairInfo> twice = pairs;
    twice.push_back(pairs[1]);
    twice.back().key = 4;
    struct BuildCase {
      const char *code;
      std::string records;
      const std::vector<LevelInfo> *levels;
      const std::vector<PairInfo> *pairs;
    };
    const std::vector<BuildCase> buildCases = {
        {"E_UNKNOWN_LEVEL", Set(kDatum, "level", "    level 9"), &levels, &pairs},
        {"E_LEVEL_JPI", Set(kDatum, "J", "    J 2+"), &levels, &pairs},
        {"E_UNKNOWN_PAIR", Set(kDatum, "reaction", "    reaction 15N(p,d)14N"), &levels, &pairs},
        {"E_AMBIGUOUS_PAIR", kDatum, &levels, &twice},
        {"E_NOT_A_CHANNEL", Set(kDatum, "reaction", "    reaction 15N(p,a)12C[4.439]"), &levels, &pairs},
        {"E_OUT_OF_RANGE", Set(kDatum, "level", "    level 6"), &levels, &pairs},
        {"E_FIELD_NOT_APPLICABLE",
         "datum c1\n    quantity anc\n    level 5\n    J 3-\n    pair p+15N\n    l 2\n    s 1\n"
         "    value 1 fm^-1/2\n    error 0.1 fm^-1/2\n",
         &levels, &pairs},
        {"W_BROAD_LEVEL", kDatum, &broad, &pairs},
        {"W_NEARBY_LEVEL", kDatum, &nearby, &pairs},
    };
    for (const BuildCase &c : buildCases) {
      const std::vector<Diagnostic> d = build(c.records, *c.levels, *c.pairs);
      const bool one = d.size() == 1 && d[0].code == c.code;
      check(std::string(c.code) + " (model data): exactly this one diagnostic", one);
      if (!one)
        for (const Diagnostic &x : d) std::printf("        got %s\n", x.ToText().c_str());
      if (one) {
        check(std::string(c.code) + ": a line and a message naming what is wrong", d[0].line > 0 && !d[0].message.empty());
        tested.insert(c.code);
      }
    }
  }
  {
    // A level's radiative decay against <levels> (BuildGammaDecays): a 1/2+ level of 15O with a proton channel and
    // photon channels to the ground state (E1), to 6.7931 MeV (M1 and E2) and to 5.183 MeV (M1); each case changes one
    // thing of a valid record or of the tables.
    auto pair = [](int key, int kind, int lightZ, int lightA, int heavyZ, int heavyA, double heavyE) {
      PairInfo p;
      p.key = key;
      p.kind = kind;
      p.lightZ = lightZ;
      p.lightA = lightA;
      p.heavyZ = heavyZ;
      p.heavyA = heavyA;
      p.heavyExcitation = heavyE;
      return p;
    };
    LevelTables tables;
    tables.pairs = {pair(1, 0, 1, 1, 7, 14, 0.0), pair(2, 10, 0, 0, 8, 15, 0.0), pair(3, 10, 0, 0, 8, 15, 6.7931),
                    pair(4, 10, 0, 0, 8, 15, 5.183)};
    tables.levels.resize(1);
    tables.levels[0] = {{5}, 1, 1, 7.5567, 1e-3, 1, 1, {{1, 0, 1}, {2, 1, 1}, {3, 1, 3}, {3, 2, 3}, {4, 1, 1}}};
    for (int c = 1; c <= 5; c++) tables.rowSlots.push_back({1, 1, c, c == 1 ? 975.0 : 0.0});
    LevelTables mixture = tables;
    mixture.rowSlots[4].value = 0.0065;
    const std::string kDecay =
        "gamma_decay d5\n    level 5\n    J 1/2+\n    total_width 0.0385 eV fit\n    branching g+15O[6.7931] 0.3 fit\n"
        "    branching g+15O 0.1 fixed\n    remainder g+15O[5.183] prior 0.6 0.05\n    mixing g+15O[6.7931] 0.1 fit\n";
    auto build = [&](const std::string &records, const LevelTables &t) {
      const std::string text = "<gammaDecays>\n" + records + "</gammaDecays>\n";
      const std::vector<Diagnostic> v = Diagnose(text);
      if (!v.empty()) return v;
      const Document doc = Document::Parse(text);
      const azure2::io::RecordReadResult r = azure2::io::ReadRecords(doc, doc.Sections()[0], Schema::KindsOfSection("gammaDecays"));
      return BuildGammaDecays(Schema::Validate({{"gammaDecays", r.set.records}}), t).diagnostics;
    };
    check("gamma_decay: the valid record builds with no diagnostic", build(kDecay, tables).empty());
    struct DecayCase {
      const char *code;
      std::string records;
      const LevelTables *tables;
    };
    const std::vector<DecayCase> decayCases = {
        {"E_GAMMA_MIXTURE", kDecay, &mixture},
        {"E_BRANCHING_SUM", Set(kDecay, "remainder", ""), &tables},
        {"E_BRANCHING_REMAINDER", Set(kDecay, "branching", "    branching g+15O[6.7931] 0.95 fit"), &tables},
        {"E_MIXING_ONE_MULTIPOLE", Set(kDecay, "mixing", "    mixing g+15O[5.183] 0.1 fit"), &tables},
        {"E_GAMMA_RECORD", Set(kDecay, "remainder", "    remainder g+15O[6.7931]"), &tables},
    };
    for (const DecayCase &c : decayCases) {
      const std::vector<Diagnostic> d = build(c.records, *c.tables);
      const bool one = d.size() == 1 && d[0].code == c.code;
      check(std::string(c.code) + " (gamma_decay): exactly this one diagnostic", one);
      if (!one)
        for (const Diagnostic &x : d) std::printf("        got %s\n", x.ToText().c_str());
      if (one) tested.insert(c.code);
    }
  }
  // the set-up's refusals of a yield on a layer (LayerChecks.h): the engine evaluates the stopping on its grid and
  // calls these; each alone, naming the target's field and saying what to write instead
  {
    const LayerPoint p{"C13A", "data/y.dat", 0.35};
    const Diagnostic s = LayerStoppingUndefined(p, 0.21, -1e-21);
    const bool okS = s.code == "E_LAYER_STOPPING" && s.record == "C13A" && s.field == "stopping" && s.found == "-1e-21" &&
                     s.message.find("0.21 MeV") != std::string::npos && s.message.find("positive") != std::string::npos;
    check("E_LAYER_STOPPING: the target's stopping, the energy and the value", okS);
    if (okS) tested.insert(s.code);
    const Diagnostic e = LayerExitBelowModel(p, 3e19, 0.00107);
    const bool okE = e.code == "E_LAYER_EXIT_BELOW_MODEL" && e.record == "C13A" && e.field == "areal_density" &&
                     e.found == "3e+19" && e.message.find("kind infinite") != std::string::npos;
    check("E_LAYER_EXIT_BELOW_MODEL: the target's areal density, and the infinite target to use instead", okE);
    if (okE) tested.insert(e.code);
  }
  // the set-up's refusals of a stopping from the library (StoppingChecks.h): each alone, naming the target's stopping
  {
    const Diagnostic m = StoppingDataMissing("TP", "stopping: no SRIM data");
    const bool okM = m.code == "E_STOPPING_DATA" && m.record == "TP" && m.field == "stopping" &&
                     m.message.find("--srim-data-dir=DIR") != std::string::npos &&
                     m.message.find("AZURE2_SRIM_DATA_DIR") != std::string::npos;
    check("E_STOPPING_DATA: the target's stopping, and where the files are read from", okM);
    if (okM) tested.insert(m.code);
    // the table's row (the schema and the format-3 page copy it) names both ways the refusal names
    std::string meaning;
    for (const DiagnosticCode &c : DiagnosticCodes())
      if (std::string(c.code) == "E_STOPPING_DATA") meaning = c.meaning;
    check("E_STOPPING_DATA: the table's row names --srim-data-dir=DIR and AZURE2_SRIM_DATA_DIR, as the refusal",
          meaning.find("--srim-data-dir=DIR") != std::string::npos &&
              meaning.find("AZURE2_SRIM_DATA_DIR") != std::string::npos);
    const Diagnostic r = StoppingOutOfRange("TP", "1H", 0.85, "outside the table's 0.93 - 1.5 MeV");
    const bool okR = r.code == "E_STOPPING_RANGE" && r.record == "TP" && r.found == "0.85" &&
                     r.message.find("0.85 MeV") != std::string::npos && r.message.find("1H") != std::string::npos;
    check("E_STOPPING_RANGE: the energy, the projectile and the library's reason", okR);
    if (okR) tested.insert(r.code);
  }
  // the build checks need <levels> pairs and data files: their cases are in build_test.cpp (Schema::Build);
  // E_UNKNOWN_PAIR and E_AMBIGUOUS_PAIR have a case in both (a datum's pair above, a data set's reaction there)
  // (E_THM_OLD_LIST: build_test.cpp, the reading path; E_LAB_JACOBIAN_AT_LIMIT: frames_test.cpp, the frames of the set-up)
  const std::set<std::string> buildCodes = {"E_REFERENCE_TABLE", "E_PAIR_KEY",    "E_DATA_FILE",   "W_DATA_FILE_NOT_READ", "E_DATA_ROW", "W_DATA_ERROR",
                                            "E_COVARIANCE_FILE", "E_FINAL_STATE_NOT_LEVEL", "E_FINAL_STATE_SHARED", "E_FINAL_STATE_UNBOUND",
                                            "E_NOT_BUILT",    "E_THM_RELATIVISTIC", "E_THM_COVARIANCE", "E_FOLD_CONFLICT",
                                            "E_THM_EXPERIMENT_UNKNOWN", "E_THM_OLD_LIST", "E_LAB_JACOBIAN_AT_LIMIT",
                                            "E_POTENTIAL_TWICE", "E_POTENTIAL_PAIR_KIND", "E_POTENTIAL_RELATIVISTIC",
                                            "E_THM_HYBRID_POTENTIAL", "E_OLD_POTENTIAL_BLOCK", "E_CONFIG_KEY_REMOVED", "E_CONFIG_OPTION_VALUE",
                                            "E_UNKNOWN_SECTION", "E_UNREAD_LINE",
                                            "E_CONFIG_OPTIONS_CONTRADICT", "E_SO_TWO_SPINS", "E_SO_ZERO_SPIN",
                                            "E_SO_RELATIVISTIC", "E_SO_THM", "E_SO_NEEDS_J", "E_SO_NO_TERM",
                                            "E_STOPPING_PROJECTILE", "E_STOPPING_ACTIVE", "E_STOPPING_SCALE_CONFLICT"};
  for (const DiagnosticCode &c : DiagnosticCodes())
    if (!buildCodes.count(c.code)) check(std::string(c.code) + " has a case here", tested.count(c.code) == 1);

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
