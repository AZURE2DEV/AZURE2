// The format-3 table on the note's own examples (docs/upgrade/project_format3.md):
//  1. every fixture of tests/data/ -- T2 in format 3 (section 3.11), the user stories (a)-(d) of section 9, the
//     examples of sections 3.5-3.10 and 8.1 -- reads and validates with zero diagnostics;
//  2. the typed values of T2 and story (a) are what the text says (reaction trees, nuisances, roles, units, the
//     per-charge factor), so that step 3's Build reads values, not text;
//  3. Schema::ToJson() equals the committed tests/schema/schema_v3.json byte for byte (check 21) and has the
//     sections of section 6.1.
//   experiment_schema_test <tests/data dir> <tests/schema/schema_v3.json>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/experiment/Schema.h"
#include "azure2/io/Document.h"
#include "azure2/io/RecordReader.h"
#include "azure2/physics/Constants.h"

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

std::string Slurp(const std::string &path) {
  std::ifstream in(path, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

// Read every record section of a project text and validate it, as step 3's reader will.
struct Checked {
  std::vector<Diagnostic> diagnostics;
  ValidationResult result;
};

Checked CheckText(const std::string &text) {
  Checked c;
  const Document doc = Document::Parse(text);
  std::vector<SectionRecords> sections;
  for (const Document::Section &s : doc.Sections()) {
    if (s.depth != 0 || Schema::KindsOfSection(s.name).empty()) continue;
    const azure2::io::RecordReadResult r = azure2::io::ReadRecords(doc, s, Schema::KindsOfSection(s.name));
    for (const azure2::io::RecordDiagnostic &d : r.diagnostics) c.diagnostics.push_back(FromRecordDiagnostic(d, s.name));
    sections.push_back({s.name, r.set.records});
  }
  c.result = Schema::Validate(sections);
  c.diagnostics.insert(c.diagnostics.end(), c.result.diagnostics.begin(), c.result.diagnostics.end());
  return c;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::printf("usage: experiment_schema_test <tests/data dir> <schema_v3.json>\n");
    return 2;
  }
  const std::string dir = argv[1];
  const std::vector<std::string> fixtures = {"T2_format3.azr",         "story_a_yield.azr",      "story_b_rutherford.azr",
                                             "story_b_five_files.azr", "story_c_model_data.azr", "story_d_exfor.azr",
                                             "section_examples.azr"};
  for (const std::string &f : fixtures) {
    const std::string text = Slurp(dir + "/" + f);
    check(f + ": read", !text.empty());
    const Checked c = CheckText(text);
    check(f + ": zero diagnostics", c.diagnostics.empty());
    for (const Diagnostic &d : c.diagnostics) std::printf("        %s\n", d.ToText().c_str());
    check(f + ": records", !c.result.records.empty());
  }

  {
    const Checked t2 = CheckText(Slurp(dir + "/T2_format3.azr"));
    const ValidationResult &r = t2.result;
    check("T2: 13 data sets, 13 apparatus records, 2 shared nuisances", r.records.size() == 13 + 13 + 2);
    const ValidatedRecord *s09 = r.Find("s09_ratio");
    check("s09: ratio tree", s09 && s09->Get("reaction")->reaction.op == ReactionExpr::Ratio &&
                                 s09->Get("reaction")->reaction.Tree() == "Ratio(15N(p,a)12C[4.43891], 15N(p,a)12C)");
    check("s09: quantity row and unit", s09 && s09->quantity && s09->quantity->name == "sigma" &&
                                            std::string(s09->Get("unit")->unit->token) == "1");
    const ValidatedRecord *s08 = r.Find("s08_sum");
    check("s08: sum of two terms, excited state 4.43891 as written",
          s08 && s08->Get("reaction")->reaction.op == ReactionExpr::Sum && s08->Get("reaction")->reaction.terms.size() == 2 &&
              s08->Get("reaction")->reaction.terms[1].reaction.state.kind == ResidualState::Excited &&
              s08->Get("reaction")->reaction.terms[1].reaction.state.energy.text == "4.43891");
    const ValidatedRecord *s01 = r.Find("s01_int_thin");
    const Nuisance *n = s01 ? &s01->Get("norm")->nuisance : nullptr;
    check("s01: norm 1 fit prior 1 5%", n && n->hasValue && n->value.text == "1" && n->mode == NuisanceMode::Fit &&
                                          n->hasPrior && n->priorWidth.text == "5" && n->priorPercent);
    const Value *shift = s01 ? s01->Get("shift") : nullptr;
    check("s01: shift offset in MeV", shift && shift->word == "offset" && shift->Get("value") &&
                                          std::string(shift->Get("value")->nuisance.unit->token) == "MeV");
    const Value *cols = s01 ? s01->Get("columns") : nullptr;
    check("s01: columns E_lab[MeV] _ value error", cols && cols->columns.roles.size() == 4 &&
                                                       cols->columns.EnergyRole()->role == "E_lab" &&
                                                       std::string(cols->columns.EnergyRole()->unit->token) == "MeV" &&
                                                       !cols->columns.AngleRole());
    const ValidatedRecord *s02 = r.Find("s02_diff090");
    check("s02: norm and shift are references", s02 && s02->Get("norm")->type == ValueType::Reference &&
                                                    s02->Get("norm")->word == "D" && s02->Get("shift")->word.empty() &&
                                                    s02->Get("shift")->Get("scale")->word == "G");
    // the old line gives no multiplier and no points per width: the converter writes only the points
    check("s02: quadrature parts", s02 && s02->Get("quadrature") && s02->Get("quadrature")->Get("points")->integer == 20 &&
                                       !s02->Get("quadrature")->Get("per_width"));
    const ValidatedRecord *s11 = r.Find("s11_thm_hoes_a");
    check("s11: norm profiled, thm LC09", s11 && s11->Get("norm")->nuisance.mode == NuisanceMode::Profiled &&
                                             !s11->Get("norm")->nuisance.hasValue && s11->Get("thm")->word == "LC09");
    const ValidatedRecord *t5 = r.Find("tgt_s05_implanted");
    const Nuisance *a = t5 ? &t5->Get("areal_density")->nuisance : nullptr;
    check("tgt_s05: areal density text kept, grid_max", a && a->value.text == "2e+17" && a->value.value == 2e17 &&
                                                            a->hasGridMax && a->gridMax.text == "3e+17" &&
                                                            std::string(a->unit->token) == "atoms/cm2");
    const Value *strag = t5 ? t5->Get("straggling") : nullptr;
    check("tgt_s05: straggling coefficient 0.6 keV^1/2 cm",
          strag && strag->word == "coefficient" && strag->Get("value")->number.text == "0.6" &&
              std::string(strag->Get("value")->unit->token) == "keV^1/2" && strag->Get("frame")->word == "cm");
    const ValidatedRecord *t1 = r.Find("tgt_s01_int_thin");
    const Value *stop = t1 ? t1->Get("stopping") : nullptr;
    check("tgt_s01: stopping formula with two parameters",
          stop && stop->word == "formula" && stop->Get("expression")->formula.expression == "a0*x^a1" &&
              stop->Get("expression")->formula.params.size() == 2 &&
              stop->Get("expression")->formula.params[1].value.text == "-0.45" &&
              std::string(stop->Get("unit")->unit->token) == "MeV*cm2/atom");
    const ValidatedRecord *nat = r.Find("natN");
    check("natN: isotope 15N", nat && nat->Get("isotope")->Get("nucleus")->nucleus.Z == 7 &&
                                   nat->Get("isotope")->Get("nucleus")->nucleus.A == 15 &&
                                   nat->Get("isotope_sum")->word == "none");
    const ValidatedRecord *b7 = r.Find("beam_s07_profile");
    check("beam_s07: profile skewed_gaussian", b7 && b7->Get("profile")->word == "skewed_gaussian" &&
                                                   b7->Get("profile")->Get("omega")->number.text == "0.035");
    const ValidatedRecord *d2 = r.Find("det_targetInt2");
    check("det_targetInt2: five Q_k", d2 && d2->Get("attenuation")->numbers.size() == 5);
  }
  {
    const Checked a = CheckText(Slurp(dir + "/story_a_yield.azr"));
    const ValidatedRecord *y = a.result.Find("Y_curve");
    check("story a: per-charge unit, charge state 1, efficiency with prior",
          y && std::string(y->Get("unit")->unit->dimension) == "counts_per_charge" && y->Get("charge_state")->integer == 1 &&
              y->Get("efficiency")->nuisance.hasPrior);
    if (y) {
      const double perUC = ParticlesPerUnitCharge(*y->Get("unit")->unit, 1);
      check("story a: 1e-6 C / (q e) protons per uC", std::fabs(perUC / (1e-6 / azure2::physics::kElementaryCharge_C) - 1) < 1e-15);
    }
    const ValidatedRecord *c = a.result.Find("C13A");
    check("story a: stopping evaluation scale fit", c && c->Get("stopping")->word == "evaluation" &&
                                                        c->Get("stopping")->Get("scale")->nuisance.mode == NuisanceMode::Fit);
  }
  {
    const Checked d = CheckText(Slurp(dir + "/story_d_exfor.azr"));
    const ValidatedRecord *r = d.result.Find("aC_RTH");
    check("story d: the alias rth is read as ratio_to_rutherford", r && r->quantity && r->quantity->name == "ratio_to_rutherford" &&
                                                                       r->Get("quantity")->word == "ratio_to_rutherford");
  }
  {
    const Checked e = CheckText(Slurp(dir + "/section_examples.azr"));
    const ValidatedRecord *r = e.result.Find("ratio_fixed_angle");
    check("section 3.5: fixed-angle denominator",
          r && r->Get("reaction")->reaction.Tree() == "Ratio(12C(p,p)12C, 12C(p,p)12C @ 90 deg cm)");
    r = e.result.Find("sum_scaled");
    check("section 3.5: scaled term", r && r->Get("reaction")->reaction.Tree() == "Sum(20Ne(n,n)20Ne, Scale(0.5, 20Ne(n,a)17O))");
    r = e.result.Find("epoch_run2");
    check("section 3.7: target@epoch", r && r->Get("target")->word == "C13B" && r->Get("target")->epoch == "run2");
    r = e.result.Find("S_points");
    check("grid range points log", r && r->Get("energies")->grid.hasPoints && r->Get("energies")->grid.points == 200 &&
                                       r->Get("energies")->grid.log && r->Get("energies")->grid.frame == "cm");
    r = e.result.Find("S_adaptive");
    check("grid adaptive", r && r->Get("energies")->grid.maxPoints == 2000 && r->Get("energies")->grid.tolerance.text == "1e-3");
    r = e.result.Find("p12C");
    check("potential: pair, R_C, Woods-Saxon parts with units", r && r->kind == "potential" && r->section == "potentials" &&
                                                                    r->Get("pair")->pair.heavy.A == 12 &&
                                                                    r->Get("coulomb_radius")->number.value == 3.5 &&
                                                                    r->Get("central")->word == "woods_saxon" &&
                                                                    r->Get("central")->Get("diffuseness")->number.text == "0.6" &&
                                                                    std::string(r->Get("central")->Get("depth")->unit->token) == "MeV");
    r = e.result.Find("a12C");
    check("potential: an excited pair, a Gaussian in keV, matching_radius and step left out",
          r && r->Get("central")->word == "gaussian" && r->Get("central")->Get("width")->number.value == 2.5 &&
              std::string(r->Get("central")->Get("depth")->unit->token) == "keV" && !r->Has("matching_radius") && !r->Has("step"));
    r = e.result.Find("anc_gs");
    check("datum spin-parity 1/2-", r && r->Get("J")->spin.twice == 1 && r->Get("J")->spin.parity == -1);
    r = e.result.Find("width_a1");
    check("datum pair a+12C[4.439]", r && r->Get("pair")->pair.light.light == "a" && r->Get("pair")->pair.heavy.A == 12 &&
                                         r->Get("pair")->pair.state.energy.text == "4.439");
  }

  // Within a section, record kinds and the field keys of its kinds are disjoint (section 3.1), so the reader can
  // tell a header from a field by its first token.
  for (const auto &section : Schema::Sections())
    for (const std::string &kind : section.second)
      for (const FieldDef &f : Schema::FindKind(kind)->fields)
        for (const std::string &other : section.second)
          check("<" + section.first + ">: field " + kind + "." + f.key + " is not a kind of the section", f.key != other);

  // The schema document and its committed copy.
  const std::string json = Schema::ToJson();
  for (const char *key : {"\"schema_version\": \"3.0\"", "\"sections\"", "\"value_types\"", "\"records\"", "\"quantities\"",
                          "\"units\"", "\"diagnostic_codes\"", "\"ratio_to_rutherford\"", "\"charge_state\""})
    check(std::string("schema has ") + key, json.find(key) != std::string::npos);
  const std::string golden = Slurp(argv[2]);
  check("Schema::ToJson() + newline == tests/schema/schema_v3.json (regenerate with AZURE2 --schema > "
        "tests/schema/schema_v3.json and bump schema_version if the table changed: check 21)",
        golden == json + "\n");

  std::printf("%d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
