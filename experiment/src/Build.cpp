// Schema::Build (Schema.h): validated format-3 records and the pairs of <levels> -> the engine's set-up specs
// (SetupSpecs.h).  docs/upgrade/project_format3.md section 8 is the specification of what each record means in the
// engine; this file resolves reactions (section 3.5), builds data sets (3.4), shared nuisances (3.6),
// extrapolations (3.8) and THM membership (3.9).  The effects of targets, beams and detectors: BuildEffects.cpp.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

#include "azure2/experiment/BuildContext.h"
#include "azure2/experiment/CovarianceFile.h"
#include "azure2/experiment/Frames.h"
#include "azure2/experiment/Grids.h"
#include "azure2/experiment/Measurement.h"
#include "azure2/experiment/ModelDatum.h"
#include "azure2/physics/Elements.h"
#include "azure2/physics/SpinObservable.h"

namespace azure2 {
namespace experiment {

using build::Context;
using build::Member;

std::vector<std::string> ParameterNamesOf(const BuiltProject &b) {
  std::vector<std::string> out;
  auto add = [&](const std::string &name) {
    if (!name.empty() && std::find(out.begin(), out.end(), name) == out.end()) out.push_back(name);
  };
  for (const auto &item : b.data.datasets.items) {
    const ParameterNames &n = item.spec.names;
    add(n.norm);
    add(n.energyShift);
    add(n.energyShiftSqrt);
    for (const std::string &f : n.shiftFunction) add(f);
    add(n.thickness);
  }
  for (const auto &item : b.extrapolations.extrapolations.items) {
    add(item.spec.names.norm);
    add(item.spec.names.energyShift);
    add(item.spec.names.energyShiftSqrt);
  }
  for (const IsospinNucleusEntry &e : b.isospin.nuclei)
    if (!e.reference) add(e.record + ".shift");
  for (const IsospinMixEntry &e : b.isospin.mixes) add(e.record + ".amplitude");
  for (const AbundanceRowSpec &row : b.abundances.abundances.rows)
    if (!row.isSum && !row.isotope.empty()) add(AbundanceParameterName(row.sample, row.isotope, row.z, row.a));
  for (const SpecList<TargetEffectSpec> *effects : {&b.data.targetEffects, &b.extrapolations.targetEffects})
    for (const auto &item : effects->items) {
      const TargetEffectSpec &e = item.spec;
      if (e.targetName.empty()) continue;
      if (e.screening) add(e.targetName + ".screening");
      if (e.densityFit) add(e.targetName + ".areal_density");
      if (e.dopplerFit) add(e.targetName + ".doppler");
      if (e.stopping.hasScale) add(e.stopping.scaleName);
    }
  for (const auto &source : b.parameterSources) add(source.first);
  return out;
}

bool BuiltProject::Ok() const {
  for (const Diagnostic &d : diagnostics)
    if (d.severity == Severity::Error) return false;
  return true;
}

namespace build {

void Context::Error(const char *code, const ValidatedRecord &record, const std::string &field, const std::string &message,
                    const std::vector<std::string> &allowed) {
  Diagnostic d;
  d.code = code;
  d.section = record.section;
  d.recordKind = record.kind;
  d.record = record.name;
  d.field = field;
  d.line = record.line;
  for (const ValidatedField &f : record.fields)
    if (f.key == field) {
      d.line = f.line;
      d.column = f.col;
      break;
    }
  d.allowed = allowed;
  d.message = record.kind + " " + record.name + (field.empty() ? "" : ", field " + field) + ": " + message;
  out_.diagnostics.push_back(d);
}

void Context::Warn(const char *code, const ValidatedRecord &record, const std::string &field, const std::string &message,
                   const std::string &found) {
  Error(code, record, field, message);
  out_.diagnostics.back().severity = Severity::Warning;
  out_.diagnostics.back().hasFound = !found.empty();
  out_.diagnostics.back().found = found;
}

const PairInfo *Context::Pair(int key) const { return PairByKey(key, pairs_); }

const ValidatedRecord *Context::Record(const std::string &name, const std::string &kind) const {
  for (const ValidatedRecord &r : v_.records)
    if (r.name == name && r.kind == kind) return &r;
  return nullptr;
}

int Context::Position(const ValidatedRecord *r) const { return r ? (int)(r - v_.records.data()) : -1; }

}  // namespace build

namespace {

// ---------------------------------------------------------------- reactions (section 3.5)

std::string Keys(const std::vector<int> &keys) {
  std::string s;
  for (int k : keys) s += (s.empty() ? "" : ", ") + std::to_string(k);
  return s;
}

}  // namespace

namespace build {

bool ResolvePair(Context &c, const ValidatedRecord &r, const std::string &labelField, const std::string &what,
                 const Particle &light, const Particle &heavy, const ResidualState &state, const std::string &text,
                 const std::string &keyField, int &key) {
  const std::vector<int> keys = PairsNamed(light, heavy, state, c.pairs_);
  const Value *given = keyField.empty() ? nullptr : r.Get(keyField);
  if (keys.empty()) {
    c.Error(codes::kUnknownPair, r, labelField, UnknownPairMessage(what), PairListing(c.pairs_));
    return false;
  }
  if (keys.size() == 1) {
    if (given) {
      c.Error(codes::kPairKey, r, keyField,
              keyField + " " + std::to_string(given->integer) + ": the label " + text + " names one pair (key " +
                  std::to_string(keys[0]) + ") already; a key is written only where the label names two pairs");
      return false;
    }
    key = keys[0];
    return true;
  }
  if (!given) {
    c.Error(codes::kAmbiguousPair, r, labelField,
            AmbiguousPairMessage(what, keys,
                                 keyField.empty() ? "a term of a sum or ratio cannot carry a pair key"
                                                  : "write " + keyField + " <key> to say which"),
            {Keys(keys)});
    return false;
  }
  if (std::find(keys.begin(), keys.end(), (int)given->integer) == keys.end()) {
    c.Error(codes::kPairKey, r, keyField,
            keyField + " " + std::to_string(given->integer) + " is not a pair the label " + text + " names (" + Keys(keys) + ")",
            {Keys(keys)});
    return false;
  }
  key = (int)given->integer;
  return true;
}

}  // namespace build

namespace {

// How a record's pair-key field applies to a term (section 3.5): a single reaction's key is written only where its
// label names several pairs; a sum's or a ratio's key applies to every term whose label names several pairs on that
// side, and a term whose label names one pair takes that pair whatever the key says.
enum class TermKeys { kSingleReaction, kTermOfComposite };

// One side of a reaction: the pair key, by the label and (where the label names several pairs) the pair-key field
// `keyField`.  `ambiguous` is set when the label names several pairs.  Returns false after an error.
bool ResolveSide(Context &c, const ValidatedRecord &r, const std::string &side, const Particle &light, const Particle &heavy,
                 const ResidualState &state, const std::string &text, const std::string &keyField, TermKeys use,
                 bool &ambiguous, int &key) {
  ambiguous = PairsNamed(light, heavy, state, c.pairs_).size() > 1;
  const std::string field = use == TermKeys::kTermOfComposite && !ambiguous ? "" : keyField;
  return build::ResolvePair(c, r, "reaction", "the " + side + " of " + text, light, heavy, state, text, field, key);
}

// Which sides of a reaction's terms named several pairs (the record's keys were needed there).
struct SidesNamingSeveral {
  bool entrance = false, exit = false;
};

// The entrance and exit pair keys of one reaction term (exit -1 for [total]).
bool ResolveTerm(Context &c, const ValidatedRecord &r, const SingleReaction &t, TermKeys use, SidesNamingSeveral &several,
                 int &entrance, int &exit) {
  bool ambiguous = false;
  const bool in = ResolveSide(c, r, "entrance", t.projectile, t.target, t.targetState, t.text, "entrance_pair", use,
                              ambiguous, entrance);
  several.entrance = several.entrance || ambiguous;
  if (!in) return false;
  if (t.state.kind == ResidualState::Total) {
    exit = -1;
    return true;
  }
  const bool out = ResolveSide(c, r, "exit", t.ejectile, t.residual, t.state, t.text, "exit_pair", use, ambiguous, exit);
  several.exit = several.exit || ambiguous;
  return out;
}

// A sum's or a ratio's pair key that no term needs: every term's label names one pair on that side, so the key would
// be a second source of truth (as on a single reaction, ResolvePair).
bool CompositeKeysNeeded(Context &c, const ValidatedRecord &r, const SidesNamingSeveral &several) {
  bool ok = true;
  for (const auto &side : {std::make_pair("entrance_pair", several.entrance), std::make_pair("exit_pair", several.exit)}) {
    const Value *given = r.Get(side.first);
    if (!given || side.second) continue;
    c.Error(codes::kPairKey, r, side.first,
            std::string(side.first) + " " + std::to_string(given->integer) + ": every term of the " +
                (r.Get("reaction")->reaction.op == ReactionExpr::Sum ? "sum" : "ratio") +
                " names one pair on that side already; a key is written only where a term's label names two pairs");
    ok = false;
  }
  return ok;
}

// ---------------------------------------------------------------- fields

std::string Word(const ValidatedRecord &r, const std::string &key, const std::string &absent = "") {
  const Value *v = r.Get(key);
  return v ? v->word : absent;
}

// The reaction, its terms and their pair keys; the observable is set by the caller.
bool ReadReaction(Context &c, const ValidatedRecord &r, const std::string &angleFrame, ReactionSpec &spec, Member &m) {
  const Value *v = r.Get("reaction");
  if (!v) return false;
  const ReactionExpr &e = v->reaction;
  const TermKeys use = e.op == ReactionExpr::Single ? TermKeys::kSingleReaction : TermKeys::kTermOfComposite;
  bool ok = true;
  SidesNamingSeveral several;
  for (std::size_t i = 0; i < e.terms.size(); i++) {
    const ReactionTerm &t = e.terms[i];
    int in = 0, out = 0;
    if (!ResolveTerm(c, r, t.reaction, use, several, in, out)) {
      ok = false;
      continue;
    }
    m.termEntranceKeys.push_back(in);
    if (t.hasAngle && !angleFrame.empty() && t.angleFrame != angleFrame) {
      c.Error(codes::kNotBuilt, r, "reaction",
              "the term at " + t.angle.text + " deg " + t.angleFrame + " is in another frame than the data set's angles (" +
                  angleFrame + "); this engine converts a fixed angle in the data set's frame only");
      ok = false;
    }
    if (i == 0) {
      spec.entranceKey = in;
      spec.exitKey = out;
      if (t.hasScale) {
        c.Error(codes::kNotBuilt, r, "reaction", "a factor on the first term: this engine scales the other terms only");
        ok = false;
      }
      continue;
    }
    ComponentSpec comp;
    comp.entranceKey = in;
    comp.exitKey = out;
    comp.hasAngle = t.hasAngle;
    comp.angle = t.hasAngle ? Context::InBase(t.angle, t.angleUnit) : 0.0;
    comp.scaling = t.hasScale ? t.scale.value : 1.0;
    spec.components.push_back(comp);
  }
  if (ok && use == TermKeys::kTermOfComposite) ok = CompositeKeysNeeded(c, r, several);
  if (e.op == ReactionExpr::Sum) spec.combination = ReactionSpec::Combination::Sum;
  if (e.op == ReactionExpr::Ratio) spec.combination = ReactionSpec::Combination::Ratio;
  m.entranceKey = spec.entranceKey;
  return ok;
}

// The phase-shift partial wave and its check against the entrance pair's spins (a channel spin s couples j1 and j2,
// and J couples l and s).
bool ReadPhaseShift(Context &c, const ValidatedRecord &r, ReactionSpec &spec) {
  const Value *J = r.Get("J"), *l = r.Get("l"), *s = r.Get("s");
  if (!J || !l || !s) return false;
  spec.phaseJ = J->spin.twice / 2.0;
  spec.phaseL = (int)l->integer;
  spec.phaseS = s->spin.twice / 2.0;
  const PairInfo *p = c.Pair(spec.entranceKey);
  if (!p) return true;
  const double sv = spec.phaseS, jv = spec.phaseJ, lv = spec.phaseL;
  const double j1 = p->twoJLight / 2.0, j2 = p->twoJHeavy / 2.0;
  const bool sOk = sv >= std::fabs(j1 - j2) - 1e-6 && sv <= j1 + j2 + 1e-6;
  const bool jOk = jv >= std::fabs(lv - sv) - 1e-6 && jv <= lv + sv + 1e-6;
  if (!sOk || !jOk) {
    c.Error(codes::kFieldNotApplicable, r, sOk ? "J" : "s",
            "(J, l, s) = (" + J->spin.text + ", " + l->number.text + ", " + s->spin.text + ") is not a partial wave of pair " +
                PairLabel(*p) + " (spins " + std::to_string(j1) + " and " + std::to_string(j2) + ")");
    return false;
  }
  return true;
}

// The S factor S = E sigma exp(2 pi eta) of the angle-integrated cross section (data_kinds.md section 4.1): a factor
// on the model at the point's c.m. energy, defined for a charged entrance pair; the modified S factor multiplies it by
// exp(g E), g the record's 'shape' (the validator required it and checked its unit).
bool SFactorOf(Context &c, const ValidatedRecord &r, ReactionSpec &spec) {
  const PairInfo *p = c.Pair(spec.entranceKey);
  if (p && (p->kind != 0 || p->lightZ * p->heavyZ == 0)) {
    c.Error(codes::kReactionNotAllowed, r, "reaction",
            "pair " + PairLabel(*p) + " has an uncharged particle: the S factor is defined for charged particles (the Gamow "
            "factor exp(2 pi eta)); give quantity sigma");
    return false;
  }
  if (r.Has("thm")) {
    c.Error(codes::kNotBuilt, r, "thm",
            "a THM data set as an S factor: THM cross sections are in arbitrary units and reach below the entrance threshold, "
            "where exp(2 pi eta) is not defined; give quantity sigma");
    return false;
  }
  spec.scale.sFactor = true;
  if (const Value *g = r.Get("shape")) spec.scale.sFactorShape_perMeV = ToBase(g->number.value, g->unit);
  return true;
}

// A reaction given with entrance_pair / exit_pair: its label names several pairs and the keys choose; elastic means
// the same pair in and out, so the quantity's elastic rules are checked here on the keys of every term (the validator
// left them).
bool PairsDecideElastic(Context &c, const ValidatedRecord &r, const ReactionSpec &spec) {
  if (!r.quantity || (!r.Has("entrance_pair") && !r.Has("exit_pair"))) return true;
  const ReactionRule &rule = r.quantity->reaction;
  std::vector<std::pair<int, int>> terms = {{spec.entranceKey, spec.exitKey}};
  for (const ComponentSpec &t : spec.components) terms.push_back({t.entranceKey, t.exitKey});
  for (const std::pair<int, int> &t : terms) {
    const bool elastic = t.first == t.second;
    const PairInfo *p = c.Pair(t.first);
    const std::string pairs = "entrance pair " + std::to_string(t.first) + ", exit pair " + std::to_string(t.second);
    if (rule.Has(ReactionRule::kElastic) && !elastic) {
      c.Error(codes::kReactionNotAllowed, r, "exit_pair",
              "quantity '" + r.quantity->name + "' is defined for elastic scattering (one pair in and out), not " + pairs);
      return false;
    }
    if (rule.Has(ReactionRule::kNoChargedElastic) && elastic && p && p->kind == 0 && p->lightZ * p->heavyZ != 0) {
      c.Error(codes::kReactionNotAllowed, r, "exit_pair",
              "quantity '" + r.quantity->name + "' is not defined for charged-particle elastic scattering (" + pairs + ")");
      return false;
    }
  }
  return true;
}

// The observable of a data set or an extrapolation from its quantity and its angle frame ("" angle-integrated).
bool ReadObservable(Context &c, const ValidatedRecord &r, const std::string &angleFrame, ReactionSpec &spec) {
  // capture into a broad unbound final state (the table allows the field for a gamma exit): keV for the engine
  if (const Value *w = r.Get("line_window")) spec.lineWindow_keV = ToBase(w->number.value, w->unit) * 1.0e3;
  const std::string q = r.quantity ? r.quantity->name : "";
  const bool data = r.kind == "dataset";
  using O = Observable;
  if (q == "sigma" || q == "sigma_nuclear" || q == "sigma_total" || q == "beta_delayed_spectrum" || q == "transmission" ||
      q == "neutron_yield" || q == "sfactor" || q == "modified_sfactor" || q == "ratio_to_reference" || (q == "yield" && angleFrame.empty())) {
    const std::string multipole = Word(r, "multipole");
    spec.observable = multipole == "E1" ? O::CrossSectionE1 : multipole == "E2" ? O::CrossSectionE2 : O::CrossSection;
    if (spec.exitKey == -1 && spec.combination == ReactionSpec::Combination::None) spec.observable = O::TotalCapture;
    spec.total = q == "sigma_total";
    if (!PairsDecideElastic(c, r, spec)) return false;
    if (q == "sfactor" || q == "modified_sfactor") return SFactorOf(c, r, spec);
    if (q == "beta_delayed_spectrum") {
      const PairInfo *p = c.Pair(spec.entranceKey);
      if (p && p->kind != 20) {
        c.Error(codes::kReactionNotAllowed, r, "reaction", "a beta-delayed spectrum needs a beta entrance pair (type 20)");
        return false;
      }
    }
    return true;
  }
  if (q == "dsigma_domega" || q == "yield") {
    spec.observable = angleFrame == "cm" ? O::DifferentialCM : O::DifferentialLab;
    if (q == "yield" && angleFrame == "cm") {
      c.Error(codes::kNotBuilt, r, "", "a yield per steradian at a c.m. angle: this engine folds lab-angle yields only");
      return false;
    }
    return true;
  }
  if (q == "ratio_to_rutherford" || q == "ratio_to_mott") {
    // dsigma/dOmega over the Coulomb cross section, both at the c.m. angle (data_kinds.md section 4.3); lab angles are
    // converted at load as for dsigma_domega, the value is not (it is the same in both frames).
    const bool mott = q == "ratio_to_mott";
    spec.observable = angleFrame == "cm" ? O::DifferentialCM : O::DifferentialLab;
    spec.coulombReference = mott ? CoulombReference::Mott : CoulombReference::Rutherford;
    // The validator read the reaction's labels; here the pairs of <levels> they resolved to decide.
    const PairInfo *p = c.Pair(spec.entranceKey);
    if (spec.combination != ReactionSpec::Combination::None || spec.exitKey != spec.entranceKey) {
      c.Error(codes::kReactionNotAllowed, r, "reaction",
              "the ratio to the " + std::string(mott ? "Mott" : "Rutherford") +
                  " cross section is defined for one elastic reaction (entrance and exit the same pair of <levels>): "
                  "give dsigma_domega");
      return false;
    }
    if (p && (p->kind != 0 || p->lightZ * p->heavyZ == 0)) {
      c.Error(codes::kReactionNotAllowed, r, "reaction",
              "pair " + PairLabel(*p) + " has no Coulomb scattering (a neutral or non-particle pair): give dsigma_domega");
      return false;
    }
    if (p && p->identical != mott) {
      c.Error(codes::kReactionNotAllowed, r, "quantity",
              mott ? "pair " + PairLabel(*p) + " has no identical particles: use quantity 'ratio_to_rutherford'"
                   : "pair " + PairLabel(*p) + " has identical particles; the Coulomb reference is Mott-symmetrised: use "
                                               "quantity 'ratio_to_mott'");
      return false;
    }
    return true;
  }
  if (q == "analyzing_power") {
    spec.observable = O::AnalyzingPower;
    // The engine's analysing power is that of a polarised spin-1/2 projectile on a target; a beta decay or a photon
    // entrance has none (the label rule "spin" lets a beta through: its letter is no d, a or g).
    const PairInfo *entrance = c.Pair(spec.entranceKey);
    if (entrance && entrance->kind != 0) {
      c.Error(codes::kReactionNotAllowed, r, "reaction",
              "an analyzing power needs a particle entrance (a polarized beam on a target); entrance pair " +
                  std::to_string(spec.entranceKey) + " is " + (entrance->kind == 20 ? "a beta decay" : "a photon pair"));
      return false;
    }
    if (angleFrame != "cm") {
      c.Error(codes::kNotBuilt, r, "", "an analyzing power at a lab angle: this engine computes it at c.m. angles");
      return false;
    }
    return true;
  }
  if (q == "spin_observable") {
    // A_y's route with the operator's trace in place of A_y's (physics/SpinObservable.h): the amplitude in the product
    // basis of the four particles' spin projections, so a reaction between two particle pairs.
    spec.observable = O::AnalyzingPower;
    const Value *op = r.Get("operator");
    if (!op) return false;  // the validator refused it
    physics::SpinObservable observable;
    std::string why;
    if (!observable.Parse(op->word, why)) return false;  // refused by the value's parser
    spec.spinObservable = op->word;
    if (spec.combination != ReactionSpec::Combination::None || spec.exitKey < 0) {
      c.Error(codes::kReactionNotAllowed, r, "reaction", "a spin observable is that of one reaction: one entrance and one exit pair");
      return false;
    }
    const PairInfo *entrance = c.Pair(spec.entranceKey), *exit = c.Pair(spec.exitKey);
    if (!entrance || !exit) return false;  // the pairs were refused already
    if (entrance->kind != 0 || exit->kind != 0) {
      c.Error(codes::kReactionNotAllowed, r, "reaction",
              "a spin observable needs particle pairs on both sides; " + std::string(entrance->kind != 0 ? "entrance" : "exit") +
                  " pair " + std::to_string(entrance->kind != 0 ? spec.entranceKey : spec.exitKey) + " is " +
                  ((entrance->kind != 0 ? entrance->kind : exit->kind) == 20 ? "a beta decay" : "a photon pair") +
                  " (for capture only A_y, quantity analyzing_power, is computed)");
      return false;
    }
    if (!observable.Check(0.5 * entrance->twoJLight, 0.5 * entrance->twoJHeavy, 0.5 * exit->twoJLight, 0.5 * exit->twoJHeavy,
                          why)) {
      c.Error(codes::kReactionNotAllowed, r, "operator", "operator " + op->word + ": " + why);
      return false;
    }
    if (angleFrame != "cm") {
      c.Error(codes::kNotBuilt, r, "", "a spin observable at a lab angle: this engine computes it at c.m. angles");
      return false;
    }
    return true;
  }
  if (q == "phase_shift") {
    spec.observable = O::PhaseShift;
    return PairsDecideElastic(c, r, spec) && ReadPhaseShift(c, r, spec);
  }
  if (q == "legendre") {
    spec.observable = data ? O::AngularCoefficients : O::AngularDistribution;
    spec.legendre = LegendreExpansionOf(Word(r, "expansion", "a_over_a0"));
    return true;
  }
  c.Error(codes::kNotBuilt, r, "quantity", "quantity " + q + ": this engine cannot set it up yet (a later step)");
  return false;
}

// A data file's path as the engine opens it: relative to the project directory ("" = the working directory).  Messages
// name the file as the project writes it, never this path.
std::string DataPath(const Context &c, const std::string &file) { return c.ProjectPath(file); }

// The reference table of a ratio_to_reference record (ReferenceTable.h): read once, its energies taken from the data
// set's energy column (an extrapolation: its grid) -- `toMeV` MeV per unit, c.m. of the entrance pair when `cm` --
// to the pair's lab frame in MeV, its values to b.  E_REFERENCE_TABLE for a file that cannot be read or is not a table.
bool ReferenceInto(Context &c, const ValidatedRecord &r, double toMeV, bool cm, ReactionSpec &spec) {
  const Value *ref = r.Get("reference");
  if (!ref) return true;
  const Value *file = ref->Get("file");
  const Value *unit = ref->Get("unit");
  if (!file) return false;  // the validator refused it
  const std::string path = DataPath(c, file->word);
  std::ifstream in(path);
  std::vector<double> energy, value;
  const std::string why = in ? ReadReferenceRows(in, energy, value) : "cannot be read";
  if (!why.empty()) {
    c.Error(codes::kReferenceTable, r, "reference", "reference table " + file->word + ": " + why);
    return false;
  }
  const PairInfo *p = c.Pair(spec.entranceKey);
  if (cm && !p) return false;  // the pair was refused already
  Partition entrance;
  if (p) {
    entrance.m1 = p->lightMass;
    entrance.m2 = p->heavyMass;
  }
  ReferenceTable &t = spec.reference;
  t.file = file->word;
  const double toBarn = unit && unit->unit ? unit->unit->factor : 1.0;
  for (std::size_t i = 0; i < energy.size(); i++) {
    const double e = energy[i] * toMeV;
    t.energyLab_MeV.push_back(cm ? BeamLabFromCM(entrance, e, c.opt_.relativistic) : e);
    t.sigma_b.push_back(value[i] * toBarn);
  }
  return true;
}

// E_REFERENCE_TABLE for each lab energy (MeV) of `energies` outside the reference table (none: true).
bool InsideReference(Context &c, const ValidatedRecord &r, const ReferenceTable &t, const std::vector<double> &energies,
                     const std::string &what) {
  bool ok = true;
  for (double e : energies) {
    if (!std::isnan(ReferenceAt(t, e))) continue;
    std::ostringstream m;
    m << what << " at E_lab = " << e << " MeV is outside the reference table " << t.file << " (E_lab " << t.energyLab_MeV.front()
      << " - " << t.energyLab_MeV.back() << " MeV): the table is not extrapolated";
    c.Error(codes::kReferenceTable, r, "reference", m.str());
    ok = false;
  }
  return ok;
}

// The data file exists and its first row holds every column the roles name (a run that may not read the data).
bool CheckFirstDataRow(Context &c, const ValidatedRecord &r, const std::string &file, std::ifstream &in, std::size_t columns) {
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream row(line);
    std::string token;
    std::size_t n = 0;
    if (!(row >> token) || token[0] == '#') continue;
    n = 1;
    while (row >> token) n++;
    if (n < columns) {
      c.Error(codes::kDataFile, r, "file",
              file + ": its first row has " + std::to_string(n) + " columns, the roles name " + std::to_string(columns));
      return false;
    }
    return true;
  }
  c.Error(codes::kDataFile, r, "file", file + " holds no data row");
  return false;
}

// An error of `r` at `field` with the text found.
void ErrorFound(Context &c, const char *code, const ValidatedRecord &r, const std::string &field, const std::string &message,
                const std::string &found) {
  c.Error(code, r, field, message);
  c.out_.diagnostics.back().hasFound = true;
  c.out_.diagnostics.back().found = found;
}

// The covariance file of a data set read by its one reader (CovarianceFile.h), and the data set's file among its
// "file" lines, as the set-up requires; null after an error.
std::shared_ptr<CovarianceFile> CheckCovarianceFile(Context &c, const ValidatedRecord &r, const DatasetSpec &d,
                                                    const std::string &dataPath) {
  CovarianceFileAccess access;
  access.resolve = [&c](const std::string &name) { return DataPath(c, name); };
  std::string err;
  std::shared_ptr<CovarianceFile> f = ReadCovarianceFile(d.covarianceFile, access, err);
  if (!f) {
    ErrorFound(c, codes::kCovarianceFile, r, "covariance", err, d.covarianceFile);
    return nullptr;
  }
  if (f->FileIndex(CanonicalPath(dataPath)) < 0) {
    ErrorFound(c, codes::kCovarianceFile, r, "covariance",
               "the data file " + d.dataFile + " is not among the \"file\" lines of the covariance file " + d.covarianceFile +
                   " (" + std::to_string(f->files.size()) + " file(s), " + std::to_string(f->N) + " rows)",
               d.covarianceFile);
    return nullptr;
  }
  return f;
}

// The data file of an active data set as the run reads it.  A run that may not use the data (DataFiles::NotUsed: the
// engine's own build, which runs before the menu) checks the first row's columns only.  A run that uses them
// (--validate, --describe, pyazr.validate) reads every row as the set-up does and lists every row it would refuse:
// one that does not parse, has fewer columns than the roles, an angular-coefficient row without a Legendre order, and
// the covariance file: readable, covering the data file with its row count, a positive variance on every row inside
// the window.  A non-positive uncertainty inside the window (no covariance file) is a warning: the run keeps the row.
bool CheckDataFile(Context &c, const ValidatedRecord &r, const DatasetSpec &d, std::size_t columns) {
  const std::string path = DataPath(c, d.dataFile);
  std::ifstream in(path);
  if (!in) {
    if (c.opt_.dataFiles == DataFiles::NotUsed) {  // an extrapolation run never opens it
      c.Warn(codes::kDataFileNotRead, r, "file", "cannot read " + d.dataFile + "; not read: this run does not use data",
             d.dataFile);
      return true;
    }
    c.Error(codes::kDataFile, r, "file", "cannot read " + d.dataFile);
    return false;
  }
  // A ratio to a reference table reads every row in every run: the table must cover each kept energy.
  if (c.opt_.dataFiles == DataFiles::NotUsed && d.reaction.reference.Empty()) return CheckFirstDataRow(c, r, d.dataFile, in, columns);
  std::vector<double> keptEnergies;
  bool ok = true;
  std::shared_ptr<CovarianceFile> cov;
  if (!d.covarianceFile.empty()) {
    cov = CheckCovarianceFile(c, r, d, path);
    ok = cov != nullptr;
  }
  const std::string canonicalData = cov ? CanonicalPath(path) : "";
  auto checkKept = [&](const DataRow &row) -> std::string {
    keptEnergies.push_back(row.energy);
    if (cov) {  // the set-up takes sqrt(C_ii) as the row's uncertainty (Segment::Fill)
      const int fileRow = cov->RowOf(canonicalData, row.record);
      const double v = fileRow < 0 ? 0.0 : cov->Variance(fileRow, row.value);
      if (v > 0.0) return "";
      std::ostringstream m;
      m << "row " << fileRow + 1 << " of the covariance file " << d.covarianceFile << " has a variance " << v << " <= 0";
      return m.str();
    }
    if (row.error > 0.0) return "";
    return row.error == 0.0 ? "the uncertainty is 0: the chi^2 skips the point (it is still counted in N)"
                            : "the uncertainty is negative: the chi^2 squares it as if it were positive";
  };
  const bool orders = d.reaction.observable == Observable::AngularCoefficients;
  const DataRowCheck check = CheckDataRows(in, d.dataFile, WindowOf(d), d.columns, orders, columns, checkKept);
  for (const DataRowProblem &p : check.problems) {
    const bool kept = p.kind == DataRowProblem::Kind::RefusedKept;
    if (kept && !cov) {  // the run computes the chi^2 without the point (chi_squared.rst): a warning
      c.Warn(codes::kDataError, r, "file", p.message, p.found);
      continue;
    }
    const char *code = kept ? codes::kCovarianceFile
                       : p.kind == DataRowProblem::Kind::TooFewColumns ? codes::kDataFile : codes::kDataRow;
    ErrorFound(c, code, r, kept ? "covariance" : "file", p.message, p.found);
    ok = false;
  }
  if (check.problemsNotListed > 0) {
    c.Error(codes::kDataRow, r, "file", std::to_string(check.problemsNotListed) + " further problem(s) in " + d.dataFile + " not listed");
    ok = false;
  }
  if (check.provenance.rowsRead == 0 && check.problems.empty()) {
    c.Error(codes::kDataFile, r, "file", d.dataFile + " holds no data row");
    ok = false;
  }
  if (!d.reaction.reference.Empty()) ok = InsideReference(c, r, d.reaction.reference, keptEnergies, "the row") && ok;
  return ok;
}

// ---------------------------------------------------------------- shared nuisances (section 3.6)

void NormInto(const Nuisance &n, DatasetSpec &d) {
  if (n.mode == NuisanceMode::Profiled) {
    // A THM experiment's profiled scale (Segment::IsProfiledNorm): a fitted norm, eliminated in closed form.
    d.norm = 1.0;
    d.varyNorm = true;
    return;
  }
  d.norm = n.value.value;
  d.varyNorm = n.mode == NuisanceMode::Fit;
  if (n.hasPrior) {
    // The engine's prior width is in per cent of the nominal norm (AZURECalc: n0 sigma / 100).
    d.normError = n.priorPercent ? n.priorWidth.value : 100.0 * n.priorWidth.value / n.priorCentre.value;
    if (n.priorCentre.value != n.value.value) {
      d.normCentre.written = true;
      d.normCentre.value = n.priorCentre.value;
      d.normCentre.source = "norm prior " + n.priorCentre.text;
    }
  }
}

// One shift term (offset or sqrt) into the data set; `what` names it in messages.
void OffsetInto(const Nuisance &n, DatasetSpec &d) {
  d.energyShift = Context::InBase(n.value, n.unit);
  d.varyEnergyShift = n.mode == NuisanceMode::Fit;
  if (n.hasPrior) {
    d.energyShiftError = Context::InBase(n.priorWidth, n.unit);
    const double centre = Context::InBase(n.priorCentre, n.unit);
    if (centre != 0.0) {
      d.shiftCentre.written = true;
      d.shiftCentre.value = centre;
      d.shiftCentre.source = "shift offset prior " + n.priorCentre.text;
    }
  }
}

void SqrtInto(const Nuisance &n, DatasetSpec &d) {
  d.hasSqrtShift = true;
  d.energyShiftSqrt = Context::InBase(n.value, n.unit);
  d.varyEnergyShiftSqrt = n.mode == NuisanceMode::Fit;
  if (n.hasPrior) d.energyShiftSqrtError = Context::InBase(n.priorWidth, n.unit);
  d.shiftSqrtCentre = build::PriorCentreOf(n, Context::InBase(n.priorCentre, n.unit), "shift sqrt");
}

void FunctionInto(const ValidatedRecord &scale, DatasetSpec &d) {
  const Value *f = scale.Get("function");
  if (!f) return;
  const Value *expr = f->Get("expression");
  d.shiftFunction = expr ? expr->formula.expression : "";
  for (const Value *p : scale.GetAll("param")) {
    const Value *value = p->Get("value");
    if (!value) continue;
    const Nuisance &n = value->nuisance;
    d.shiftFunctionValues.push_back(n.value.value);
    d.shiftFunctionErrors.push_back(n.hasPrior ? n.priorWidth.value : 0.0);
    d.shiftFunctionVary.push_back(n.mode == NuisanceMode::Fit);
    d.shiftFunctionCentres.push_back(build::PriorCentreOf(n, n.priorCentre.value, "shift function parameter"));
  }
}

// ---------------------------------------------------------------- records

bool BuildDataset(Context &c, const ValidatedRecord &r, DatasetSpec &d, Member &m,
                  std::map<std::string, int> &firstActiveOfScale) {
  d.active = m.active;
  bool ok = true;
  const Value *cols = r.Get("columns");
  const ColumnRole *angleRole = cols ? cols->columns.AngleRole() : nullptr;
  const std::string angleFrame = angleRole ? (angleRole->role == "theta_cm" ? "cm" : "lab") : "";
  ok = ReadReaction(c, r, angleFrame, d.reaction, m) && ok;
  if (!ok) return false;
  ok = ReadObservable(c, r, angleFrame, d.reaction) && ok;
  ok = build::ValueUnit(c, r, d.reaction, &d.columns) && ok;
  d.absolute = !LegendreExpansionNormalised(d.reaction.legendre);
  if (const Value *w = r.Get("angle_window")) {
    d.minA = Context::InBase(w->lo, w->unit);
    d.maxA = Context::InBase(w->hi, w->unit);
  } else {  // every angle (an angle-integrated observable has none; the engine reports the window, normalizations.out)
    d.minA = 0.0;
    d.maxA = 180.0;
  }
  if (const Value *f = r.Get("file")) d.dataFile = f->word;
  if (cols) {
    const std::string second = angleRole ? angleRole->role : (r.quantity && r.quantity->name == "legendre") ? "order" : "_";
    ok = build::DataColumns(c, r, cols->columns, second, d) && ok;
    if (ok) ok = ReferenceInto(c, r, d.columns.energy, d.columns.energyCM, d.reaction) && ok;
  }
  build::EnergyWindow(r, d);
  // norm: its own nuisance, or a shared normalization record (its values in every member; the engine takes the group's
  // from its first loaded member, WorkingCopy::ResolveNormGroups)
  if (const Value *n = r.Get("norm")) {
    if (n->type == ValueType::Reference) {
      const ValidatedRecord *g = c.Record(n->word, "normalization");
      const Value *gv = g ? g->Get("value") : nullptr;
      if (gv) NormInto(gv->nuisance, d);
      d.normGroup = n->word;
    } else {
      NormInto(n->nuisance, d);
    }
  }
  // shift terms: offset, sqrt, a shared energy_scale (summed by the engine in that order, section 3.6)
  int offsets = 0, sqrts = 0, scales = 0;
  for (const Value *s : r.GetAll("shift")) {
    if (s->word == "offset") {
      offsets++;
      if (const Value *v = s->Get("value")) OffsetInto(v->nuisance, d);
    } else if (s->word == "sqrt") {
      sqrts++;
      if (const Value *v = s->Get("value")) SqrtInto(v->nuisance, d);
    } else if (const Value *ref = s->Get("scale")) {
      scales++;
      const ValidatedRecord *g = c.Record(ref->word, "energy_scale");
      if (!g) continue;
      const bool sharesOffset = g->Has("offset");
      if (sharesOffset) {
        // shiftgroup: the scale's offset is the group's constant shift; its function belongs to the group's first
        // active member (WorkingCopy::ResolveShiftFunctions refuses a function on a follower)
        d.shiftGroup = g->name;
        if (const Value *o = g->Get("offset")) OffsetInto(o->nuisance, d);
        if (const Value *q = g->Get("sqrt")) {
          sqrts++;
          SqrtInto(q->nuisance, d);
        }
      } else {
        d.shiftFunctionGroup = g->name;
      }
      if (m.active && !firstActiveOfScale.count(g->name)) firstActiveOfScale[g->name] = d.key;
      if (m.active && firstActiveOfScale[g->name] == d.key) FunctionInto(*g, d);
    }
  }
  if (offsets > 1 || sqrts > 1 || scales > 1 || (offsets > 0 && !d.shiftGroup.empty())) {
    c.Error(codes::kFieldTwice, r, "shift",
            "two terms of one kind (two offsets, two sqrt terms, two energy_scale records, or an offset and a scale "
            "that holds one); the shift has one term of each kind");
    ok = false;
  }
  if (const Value *cov = r.Get("covariance")) d.covarianceFile = cov->word;
  if (const Value *g = r.Get("secondary_gamma")) {
    d.upos = true;
    if (const Value *L = g->Get("L")) d.secondaryDecayL = (int)L->integer;
    if (const Value *Ic = g->Get("Ic")) d.secondaryFinalSpin = Ic->spin.twice / 2.0;
    if (const Value *delta = g->Get("delta")) d.secondaryMixingRatio = delta->number.value;
  }
  if (const Value *thm = r.Get("thm")) {
    d.reaction.thm = true;
    if (c.opt_.relativistic) {
      c.Error(codes::kThmRelativistic, r, "thm", "THM data sets are not supported with relativistic=1 (the THM vertex, line shape and "
                                                 "kinematics are non-relativistic); remove relativistic=1 from "
                                                 "<config>, or the thm field of the data set");
      ok = false;
    }
    if (!d.covarianceFile.empty()) {
      c.Error(codes::kThmCovariance, r, "covariance", "a THM data set's chi^2 is its experiment's (profiled norm), not a sum a "
                                                      "covariance could amend; remove the covariance field");
      ok = false;
    }
    if (thm->word != "global") {
      auto &xs = c.out_.thm.experiments;
      auto it = std::find_if(xs.begin(), xs.end(), [&](const ThmMembershipSpec::Experiment &x) { return x.name == thm->word; });
      if (it == xs.end()) {
        xs.push_back({thm->word, {}});
        it = xs.end() - 1;
      }
      if (m.active) it->datasetKeys.push_back(d.key);
    }
    if (const Value *w = r.Get("thm_weight")) c.out_.thm.datasetWeights.push_back({d.key, w->word});
  }
  // the data file last: the window, the column conversion and the covariance file are known
  if (ok && cols && m.active && c.opt_.readDataFiles) ok = CheckDataFile(c, r, d, cols->columns.roles.size());
  return ok;
}

// The grid of an extrapolation (Grids.h).  `lo` and `hi` are its bounds as the segment reports them.
bool ReadGrid(Context &c, const ValidatedRecord &r, const std::string &field, const Grid &g, const PairInfo *entrance,
              std::vector<double> &values, double &lo, double &hi, std::string *file = nullptr) {
  GridSetting setting;
  setting.energies = field == "energies";
  setting.entrance = entrance;
  setting.relativistic = c.opt_.relativistic;
  setting.projectDir = c.opt_.projectDir;
  GridValues v = ReadGridValues(g, setting);
  if (!v.ok) {
    c.Error(v.code, r, field, v.error);
    return false;
  }
  values = std::move(v.values);
  lo = v.lo;
  hi = v.hi;
  if (file) *file = v.file;
  return true;
}

bool BuildExtrapolation(Context &c, const ValidatedRecord &r, ExtrapolationSpec &x, Member &m) {
  x.active = m.active;
  const Value *angles = r.Get("angles");
  const bool hasAngles = angles && angles->grid.kind != "none";
  const std::string angleFrame = hasAngles ? angles->grid.frame : "";
  bool ok = ReadReaction(c, r, angleFrame, x.reaction, m);
  if (!ok) return false;
  ok = ReadObservable(c, r, angleFrame, x.reaction) && ok;
  ok = build::ValueUnit(c, r, x.reaction, nullptr) && ok;
  if (const Value *e = r.Get("energies")) {
    ok = ReadGrid(c, r, "energies", e->grid, c.Pair(x.reaction.entranceKey), x.energies, x.minE, x.maxE, &x.gridFile) && ok;
    if (ok && r.Has("reference")) {
      ok = ReferenceInto(c, r, e->grid.unit ? e->grid.unit->factor : 1.0, e->grid.frame == "cm", x.reaction);
      const std::vector<double> ends = x.energies.empty() ? std::vector<double>{x.minE, x.maxE} : x.energies;
      if (ok) ok = InsideReference(c, r, x.reaction.reference, ends, "the grid point");
    }
    if (e->grid.kind == "adaptive") {
      x.adaptive.on = true;
      x.adaptive.tolerance = e->grid.tolerance.value;
      x.adaptive.maxPoints = e->grid.maxPoints;
    }
  }
  if (hasAngles) ok = ReadGrid(c, r, "angles", angles->grid, nullptr, x.angles, x.minA, x.maxA) && ok;
  else x.angles = {0.0};
  if (const Value *o = r.Get("orders")) {
    const Value *first = o->Get("first"), *last = o->Get("last");
    if (first && first->integer != 0) {
      c.Error(codes::kNotBuilt, r, "orders", "this engine computes the orders 0 ... K (orders 0 <K>)");
      ok = false;
    }
    if (last) x.maxAngDistOrder = (int)last->integer;
  }
  if (r.Has("thm")) {
    x.reaction.thm = true;
    if (c.opt_.relativistic) {
      c.Error(codes::kThmRelativistic, r, "thm", "THM extrapolations are not supported with relativistic=1 (the THM vertex is "
                                                 "non-relativistic); remove relativistic=1 from <config>, or the "
                                                 "thm field of the extrapolation");
      ok = false;
    }
    if (const Value *w = r.Get("thm_weight")) c.out_.thm.extrapolationWeights.push_back({x.key, w->word});
  }
  return ok;
}

// ---------------------------------------------------------------- parameter names (SetupSpecs.h, ParameterNames)

// The names of the parameters the engine makes for each data set and extrapolation: the data set's own, a shared
// record's for the first active member of its group (the engine's leader rule, WorkingCopy::Resolve*Groups), and for a
// later member, whose parameter the engine keeps fixed, <record>.<field>[<data set>].
void NameParameters(Context &c) {
  std::map<std::string, int> normLeader, shiftLeader, thicknessLeader;
  auto shared = [](std::map<std::string, int> &leaders, const std::string &label, int key, const std::string &field,
                   const std::string &dataset, const std::string &own) {
    if (label.empty()) return dataset + "." + own;
    const int leader = leaders.emplace(label, key).first->second;
    return label + "." + field + (leader == key ? "" : "[" + dataset + "]");
  };
  std::map<std::string, ParameterSource> &sources = c.out_.parameterSources;
  // The occurrence of a data set's shift field with the alternative `word` (offset, sqrt), or -1.
  auto shiftOccurrence = [](const ValidatedRecord &r, const std::string &word) {
    const std::vector<const Value *> shifts = r.GetAll("shift");
    for (std::size_t k = 0; k < shifts.size(); k++)
      if (shifts[k]->word == word) return (int)k;
    return -1;
  };
  for (std::size_t i = 0; i < c.data_.size(); i++) {
    const Member &m = c.data_[i];
    DatasetSpec &d = c.out_.data.datasets.items[i].spec;
    const std::string &name = m.record->name;
    if (!m.active) continue;
    d.names.norm = shared(normLeader, d.normGroup, d.key, "value", name, "norm");
    d.names.energyShift = shared(shiftLeader, d.shiftGroup, d.key, "offset", name, "shift.offset");
    d.names.energyShiftSqrt = name + ".shift.sqrt";
    const std::string &scale = d.shiftFunctionGroup.empty() ? d.shiftGroup : d.shiftFunctionGroup;
    for (std::size_t j = 0; j < d.shiftFunctionValues.size(); j++) d.names.shiftFunction.push_back(scale + ".a" + std::to_string(j));
    // where each value is written (a follower's parameter of a group is fixed and written nowhere)
    if (d.normGroup.empty())
      sources[d.names.norm] = {name, "norm", 0, "", ""};
    else if (normLeader[d.normGroup] == d.key)
      sources[d.names.norm] = {d.normGroup, "value", 0, "", ""};
    const int offset = shiftOccurrence(*m.record, "offset"), sqrt = shiftOccurrence(*m.record, "sqrt");
    if (!d.shiftGroup.empty()) {
      if (shiftLeader[d.shiftGroup] == d.key) sources[d.names.energyShift] = {d.shiftGroup, "offset", 0, "", ""};
    } else if (offset >= 0) {
      sources[d.names.energyShift] = {name, "shift", offset, "value", ""};
    }
    if (sqrt >= 0) sources[d.names.energyShiftSqrt] = {name, "shift", sqrt, "value", ""};
    for (std::size_t j = 0; j < d.names.shiftFunction.size(); j++)
      sources[d.names.shiftFunction[j]] = {scale, "param", (int)j, "value", ""};
    if (d.thicknessFit && !d.thicknessGroup.empty()) {
      const int leader = thicknessLeader.emplace(d.thicknessGroup, d.key).first->second;
      if (leader != d.key) d.names.thickness += "[" + name + "]";
    }
  }
  for (std::size_t i = 0; i < c.extrapolations_.size(); i++) {
    ExtrapolationSpec &x = c.out_.extrapolations.extrapolations.items[i].spec;
    const std::string &name = c.extrapolations_[i].record->name;
    x.names.norm = name + ".norm";
    x.names.energyShift = name + ".shift.offset";
    x.names.energyShiftSqrt = name + ".shift.sqrt";
  }
}

}  // namespace

BuiltProject Schema::Build(const ValidationResult &validated, const LevelTables &levels, const BuildOptions &options) {
  BuiltProject out;
  Context c(validated, levels.pairs, options, out);
  std::map<std::string, int> firstActiveOfScale;
  for (const ValidatedRecord &r : validated.records) {
    if (r.kind != "dataset" && r.kind != "extrapolation") continue;
    Member m;
    m.record = &r;
    m.isData = r.kind == "dataset";
    m.active = options.allRecordsActive || Word(r, "active", "yes") != "no";
    if (m.isData) {
      SpecList<DatasetSpec>::Item item;
      item.spec.key = m.key = (int)c.data_.size() + 1;
      BuildDataset(c, r, item.spec, m, firstActiveOfScale);
      out.data.datasets.items.push_back(item);
      c.data_.push_back(m);
    } else {
      SpecList<ExtrapolationSpec>::Item item;
      item.spec.key = m.key = (int)c.extrapolations_.size() + 1;
      BuildExtrapolation(c, r, item.spec, m);
      out.extrapolations.extrapolations.items.push_back(item);
      c.extrapolations_.push_back(m);
    }
  }
  for (auto &x : out.thm.experiments) std::sort(x.datasetKeys.begin(), x.datasetKeys.end());
  build::BuildEffects(c);
  build::BuildPotentials(c);
  build::BuildChannelRadii(c);
  build::BuildIsospin(c, levels);
  if (options.setUpChecks) build::CheckFinalStates(c, levels);
  build::CheckUnboundFinalStates(c, levels);
  NameParameters(c);
  ModelDataBuild modelData = BuildModelData(validated, levels.levels, levels.pairs);
  out.modelData = std::move(modelData.data);
  out.diagnostics.insert(out.diagnostics.end(), modelData.diagnostics.begin(), modelData.diagnostics.end());
  GammaDecayBuild gammaDecays = BuildGammaDecays(validated, levels);
  out.gammaDecays = std::move(gammaDecays.specs);
  out.diagnostics.insert(out.diagnostics.end(), gammaDecays.diagnostics.begin(), gammaDecays.diagnostics.end());
  std::stable_sort(out.diagnostics.begin(), out.diagnostics.end(),
                   [](const Diagnostic &a, const Diagnostic &b) { return a.line < b.line; });
  if (!out.Ok()) {
    out.data.datasets.failed = out.extrapolations.extrapolations.failed = true;
    out.abundances.failed = true;
  }
  return out;
}

}  // namespace experiment
}  // namespace azure2
