// The value-type parsers (Parsers.h) and the tables they read, one type at a time:
//  numbers keep their text; integers; units, aliases and the case-only fix; nuisances (modes, priors with and without
//  %, grid_max, references, profiled); spins with and without parity; nuclides and light particles, ordinals refused;
//  pairs; compositions; formulas; reaction expressions (trees, scales, fixed angles, aliases, [total], balance);
//  column roles; grids; windows; the unit table (dimensions, factors, per-charge particles) and the JSON writer.
//   experiment_parsers_test
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "azure2/experiment/Json.h"
#include "azure2/experiment/Parsers.h"
#include "azure2/io/RecordReader.h"
#include "azure2/physics/Constants.h"

using namespace azure2::experiment;

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

struct Parsed {
  bool ok = false;
  Value value;
  std::vector<Diagnostic> diagnostics;
};

// Parse `text` (one field's value) with `spec`; ok only if it parsed and nothing was left.
Parsed Parse(const ValueSpec &spec, const std::string &text) {
  Parsed p;
  std::string comment;
  std::vector<azure2::io::RecordDiagnostic> problems;
  const std::vector<azure2::io::Token> tokens = azure2::io::Tokenize(text, 1, comment, problems);
  ParseContext ctx;
  ctx.out = &p.diagnostics;
  ctx.base.recordKind = "test";
  ctx.base.record = "r";
  TokenCursor cur;
  cur.tokens = &tokens;
  cur.line = 1;
  cur.endCol = static_cast<int>(text.size()) + 1;
  p.ok = ParseValue(spec, cur, ctx, p.value) && cur.AtEnd() && problems.empty();
  return p;
}

std::string Code(const Parsed &p) { return p.diagnostics.empty() ? "" : p.diagnostics[0].code; }

ValueSpec S(ValueType t, const std::string &dimension = "") { return ValueSpec(t, dimension); }

}  // namespace

int main() {
  // Numbers and integers.
  {
    Parsed p = Parse(S(ValueType::Number), "2e+17");
    check("number: text kept", p.ok && p.value.number.text == "2e+17" && p.value.number.value == 2e17);
    check("number: 'nan', 'inf', '1x' refused",
          Code(Parse(S(ValueType::Number), "nan")) == "E_NOT_A_NUMBER" && Code(Parse(S(ValueType::Number), "inf")) == "E_NOT_A_NUMBER" &&
              Code(Parse(S(ValueType::Number), "1x")) == "E_NOT_A_NUMBER");
    ValueSpec i = S(ValueType::Integer);
    i.min = 1;
    check("integer", Parse(i, "3").ok && Parse(i, "3").value.integer == 3);
    check("integer: below min", Code(Parse(i, "0")) == "E_OUT_OF_RANGE");
    check("integer: 1.0 is not one", Code(Parse(i, "1.0")) == "E_NOT_AN_INTEGER");
  }
  // Units.
  {
    check("unit: alias µb -> ub", Parse(S(ValueType::Unit, "cross_section"), "\xC2\xB5" "b").ok &&
                                      std::string(Parse(S(ValueType::Unit, "cross_section"), "\xC2\xB5" "b").value.unit->token) == "ub");
    check("unit: wrong dimension", Code(Parse(S(ValueType::Unit, "energy"), "deg")) == "E_UNIT_NOT_ALLOWED");
    check("number with unit", Parse(S(ValueType::NumberWithUnit, "energy"), "0.3 keV").ok);
    check("number with unit: missing unit", Code(Parse(S(ValueType::NumberWithUnit, "energy"), "0.3")) == "E_UNIT_MISSING");
    check("range", Parse(S(ValueType::Range, "angle"), "90 90 deg").ok);
  }
  // Nuisances.
  {
    ValueSpec n = S(ValueType::Nuisance, "areal_density");
    Parsed p = Parse(n, "2e+17 atoms/cm2 fit prior 2e+17 10% grid_max 3e+17");
    check("nuisance: full form", p.ok && p.value.nuisance.mode == NuisanceMode::Fit && p.value.nuisance.hasPrior &&
                                     p.value.nuisance.priorPercent && p.value.nuisance.priorWidth.text == "10" &&
                                     p.value.nuisance.gridMax.text == "3e+17");
    check("nuisance: unit required", Code(Parse(n, "2e17 fixed")) == "E_UNIT_MISSING");
    check("nuisance: mode required", Code(Parse(n, "2e17 atoms/cm2")) == "E_NUISANCE_MODE");
    check("nuisance: absolute prior width", Parse(S(ValueType::Nuisance, "energy"), "0 MeV fit prior 0 0.003").ok &&
                                                !Parse(S(ValueType::Nuisance, "energy"), "0 MeV fit prior 0 0.003").value.nuisance.priorPercent);
    check("nuisance: dimensionless takes no unit", Code(Parse(S(ValueType::Nuisance), "1 b fit")) == "E_NUISANCE_MODE");
    check("nuisance: profiled only where allowed", Code(Parse(S(ValueType::Nuisance), "profiled")) == "E_MODE_NOT_ALLOWED");
    ValueSpec norm = S(ValueType::Nuisance);
    norm.profiled = true;
    norm.refKinds = {"norm"};
    check("nuisance: profiled alone", Parse(norm, "profiled").ok && !Parse(norm, "profiled").value.nuisance.hasValue);
    check("nuisance: a name is a reference", Parse(norm, "D").ok && Parse(norm, "D").value.type == ValueType::Reference);
    check("nuisance: zero prior width refused", Code(Parse(norm, "1 fit prior 1 0")) == "E_PRIOR_SYNTAX");
    check("nuisance: 'fit' without prior is free", Parse(norm, "1 fit").value.nuisance.FreeWithoutPrior());
  }
  // Spins, nuclei, pairs, compositions, formulas.
  {
    check("spin 3/2", Parse(S(ValueType::Spin), "3/2").ok && Parse(S(ValueType::Spin), "3/2").value.spin.twice == 3);
    check("spin 1.5 refused (one spelling)", Code(Parse(S(ValueType::Spin), "1.5")) == "E_SPIN_SYNTAX");
    check("spin 2/2 refused", Code(Parse(S(ValueType::Spin), "2/2")) == "E_SPIN_SYNTAX");
    check("spin-parity 0+", Parse(S(ValueType::SpinParity), "0+").ok && Parse(S(ValueType::SpinParity), "0+").value.spin.parity == 1);
    check("nucleus 13C", Parse(S(ValueType::Nucleus), "13C").value.nucleus.Z == 6);
    check("nucleus: A < Z", Code(Parse(S(ValueType::Nucleus), "3C")) == "E_NUCLEUS_SYNTAX");
    check("nucleus: case matters", Code(Parse(S(ValueType::Nucleus), "4HE")) == "E_UNKNOWN_ELEMENT");
    check("pair a+12C[4.439]", Parse(S(ValueType::Pair), "a+12C[4.439]").ok);
    check("pair p+15N", Parse(S(ValueType::Pair), "p+15N").value.pair.light.Z == 1);
    check("composition C 1 O 2", Parse(S(ValueType::Composition), "C 1 O 2").value.composition.size() == 2);
    check("composition: zero atoms", Code(Parse(S(ValueType::Composition), "C 0")) == "E_OUT_OF_RANGE");
    Parsed f = Parse(S(ValueType::Formula), "\"a0*x^a1\" a0=1.6e-20 a1=-0.45");
    check("formula with parameters", f.ok && f.value.formula.params.size() == 2 && f.value.formula.params[0].value.text == "1.6e-20");
    check("formula: bad parameter", Code(Parse(S(ValueType::Formula), "\"x\" a0=")) == "E_FORMULA_SYNTAX");
  }
  // Reactions.
  {
    const ValueSpec r = S(ValueType::Reaction);
    check("reaction: alias a0 is the ground state", Parse(r, "15N(p,a0)12C").ok &&
                                                        Parse(r, "15N(p,a0)12C").value.reaction.Tree() == "15N(p,a)12C");
    check("reaction: [0] is the ground state", Parse(r, "15N(p,a)12C[0]").value.reaction.terms[0].reaction.state.kind ==
                                                   ResidualState::Ground);
    check("reaction: [total] with a photon", Parse(r, "12C(p,g)13N[total]").ok);
    check("reaction: [total] with a particle", Code(Parse(r, "15N(p,a)12C[total]")) == "E_REACTION_SYNTAX");
    check("reaction: 4He inside the parentheses is 'a'", Code(Parse(r, "12C(4He,4He)12C")) == "E_REACTION_SYNTAX");
    check("reaction: heavy-ion projectile", Parse(r, "12C(12C,p)23Na").ok);
    check("reaction: scale with blanks", Parse(r, "20Ne(n,n)20Ne + 0.5 * 20Ne(n,a)17O").value.reaction.Tree() ==
                                             "Sum(20Ne(n,n)20Ne, Scale(0.5, 20Ne(n,a)17O))");
    check("reaction: a ratio of sums is refused", Code(Parse(r, "15N(p,a)12C + 15N(p,a)12C / 15N(p,a)12C")) == "E_REACTION_SYNTAX");
    check("reaction: three-term ratio refused", Code(Parse(r, "15N(p,a)12C / 15N(p,a)12C / 15N(p,a)12C")) == "E_REACTION_SYNTAX");
    check("reaction: '+' without blanks", Code(Parse(r, "15N(p,a)12C+15N(p,a)12C")) == "E_REACTION_SYNTAX");
    check("reaction: fixed angle needs a unit and frame", Code(Parse(r, "12C(p,p)12C @ 90")) == "E_TOO_FEW_TOKENS");
    check("reaction: unknown element", Code(Parse(r, "15Q(p,a)12C")) == "E_UNKNOWN_ELEMENT");
    check("reaction: particle excitation", Parse(r, "12C(a,a[3.0])12C").ok &&
                                               Parse(r, "12C(a,a[3.0])12C").value.reaction.terms[0].reaction.ejectile.excited);
    const SingleReaction e = Parse(r, "4He(a,a)4He").value.reaction.terms[0].reaction;
    check("reaction: elastic, charged, identical", e.Elastic() && e.ChargedEntrance() && e.IdenticalEntrance());
  }
  // Column roles and grids.
  {
    const ValueSpec c = S(ValueType::Columns);
    check("columns: from_file", Parse(c, "from_file").ok && Parse(c, "from_file").value.columns.fromFile);
    check("columns: rad", Parse(c, "E_cm[keV] theta_cm[rad] value error").ok);
    check("columns: unit on value refused", Code(Parse(c, "E_lab[MeV] value[b] error")) == "E_ROLE_UNIT");
    check("columns: wrong dimension", Code(Parse(c, "E_lab[deg] value error")) == "E_ROLE_UNIT");
    check("columns: '_' may repeat", Parse(c, "E_lab[MeV] _ _ value error").ok);
    ValueSpec g = S(ValueType::Grid, "energy");
    check("grid range step", Parse(g, "range 0.02 0.8 MeV cm step 0.002").ok);
    check("grid range needs step or points", Code(Parse(g, "range 0.02 0.8 MeV cm")) == "E_PART_MISSING");
    check("grid range empty", Code(Parse(g, "range 0.8 0.02 MeV cm step 0.01")) == "E_EMPTY_RANGE");
    check("grid frame", Code(Parse(g, "list 1 2 MeV lb")) == "E_NOT_IN_LIST");
    check("grid none only for angles", Code(Parse(g, "none")) == "E_NOT_IN_LIST");
    check("grid adaptive needs all three", Code(Parse(g, "adaptive 0 1 MeV cm tolerance 1e-3")) == "E_PART_MISSING");
    ValueSpec w = S(ValueType::Windows, "energy");
    check("windows: auto:1e-3:0.1", Parse(w, "auto:1e-3:0.1").ok && Parse(w, "auto:1e-3:0.1").value.windows.margin.text == "0.1");
    check("windows: exponents in ranges", Parse(w, "1e-3-2e-3 MeV").ok &&
                                              Parse(w, "1e-3-2e-3 MeV").value.windows.windows[0].hi.value == 2e-3);
    check("windows: empty", Code(Parse(w, "0.4-0.2 MeV")) == "E_EMPTY_RANGE");
  }
  // The unit table.
  {
    for (const Unit &u : Units()) {
      check(std::string(u.token) + ": dimension exists", FindDimension(u.dimension) != nullptr);
      check(std::string(u.token) + ": factor > 0", u.factor > 0);
    }
    for (const Dimension &d : Dimensions()) {
      const Unit *base = FindUnit(d.base);
      check(std::string(d.name) + ": base unit has factor 1", base && base->factor == 1.0 && std::string(base->dimension) == d.name);
    }
    check("keV = 1e-3 MeV, atoms/b = 1e24 atoms/cm2", FindUnit("keV")->factor == 1e-3 && FindUnit("atoms/b")->factor == 1e24);
    check("protons per mC at q = 2", std::fabs(ParticlesPerUnitCharge(*FindUnit("1/mC"), 2) /
                                               (1e-3 / (2 * azure2::physics::kElementaryCharge_C)) - 1) < 1e-15);
  }
  // JSON.
  {
    check("json number: shortest round trip", JsonNumber(1e-3) == "0.001" && JsonNumber(1e24) == "1e+24" && JsonNumber(0.1) == "0.1");
    check("json escape", JsonEscape("a\"b\\c\n") == "\"a\\\"b\\\\c\\n\"");
    JsonWriter j(0);
    j.BeginObject().Key("a").StringArray({"x", "y"}).Key("b").BeginArray().Integer(1).Null().EndArray().EndObject();
    check("json writer", j.Str() == "{\"a\": [\"x\", \"y\"], \"b\": [1, null]}");
  }
  std::printf("%d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
