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
#include "azure2/experiment/ModelDatum.h"
#include "azure2/physics/Elements.h"

namespace azure2 {
namespace experiment {

using build::Context;
using build::Member;

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

// One side of a reaction: the pair key, by the label and (where the label names several pairs) the pair-key field.
// `keyField` is "" for a term of a sum or ratio, which cannot carry a key.  Returns false after an error.
bool ResolveSide(Context &c, const ValidatedRecord &r, const std::string &side, const Particle &light, const Particle &heavy,
                 const ResidualState &state, const std::string &text, const std::string &keyField, int &key) {
  const std::vector<int> keys = PairsNamed(light, heavy, state, c.pairs_);
  const Value *given = keyField.empty() ? nullptr : r.Get(keyField);
  if (keys.empty()) {
    c.Error(codes::kUnknownPair, r, "reaction", UnknownPairMessage("the " + side + " of " + text), PairListing(c.pairs_));
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
    c.Error(codes::kAmbiguousPair, r, "reaction",
            AmbiguousPairMessage("the " + side + " of " + text, keys,
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

// The entrance and exit pair keys of one reaction term (exit -1 for [total]).
bool ResolveTerm(Context &c, const ValidatedRecord &r, const SingleReaction &t, bool single, int &entrance, int &exit) {
  if (!ResolveSide(c, r, "entrance", t.projectile, t.target, t.targetState, t.text, single ? "entrance_pair" : "", entrance))
    return false;
  if (t.state.kind == ResidualState::Total) {
    exit = -1;
    return true;
  }
  return ResolveSide(c, r, "exit", t.ejectile, t.residual, t.state, t.text, single ? "exit_pair" : "", exit);
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
  const bool single = e.op == ReactionExpr::Single;
  bool ok = true;
  for (std::size_t i = 0; i < e.terms.size(); i++) {
    const ReactionTerm &t = e.terms[i];
    int in = 0, out = 0;
    if (!ResolveTerm(c, r, t.reaction, single, in, out)) {
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
  if (e.op == ReactionExpr::Sum) spec.combination = ReactionSpec::Combination::Sum;
  if (e.op == ReactionExpr::Ratio) spec.combination = ReactionSpec::Combination::Ratio;
  m.entranceKey = spec.entranceKey;
  return ok;
}

// The unit of the values: the engine compares the numbers as written, in its own unit.
bool EngineUnit(Context &c, const ValidatedRecord &r) {
  const Value *u = r.Get("unit");
  if (u && u->unit && u->unit->factor != 1.0) {
    c.Error(codes::kNotBuilt, r, "unit",
            std::string("values in ") + u->unit->token + ": this engine reads the values in the unit of the quantity (" +
                (r.quantity ? r.quantity->engineUnit : std::string("?")) + "); a conversion is a later step");
    return false;
  }
  return true;
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

// The observable of a data set or an extrapolation from its quantity and its angle frame ("" angle-integrated).
bool ReadObservable(Context &c, const ValidatedRecord &r, const std::string &angleFrame, ReactionSpec &spec) {
  const std::string q = r.quantity ? r.quantity->name : "";
  const bool data = r.kind == "dataset";
  using O = Observable;
  if (q == "sigma" || q == "sigma_nuclear" || q == "sigma_total" || q == "beta_delayed_spectrum" || q == "transmission" ||
      q == "neutron_yield" || (q == "yield" && angleFrame.empty())) {
    const std::string multipole = Word(r, "multipole");
    spec.observable = multipole == "E1" ? O::CrossSectionE1 : multipole == "E2" ? O::CrossSectionE2 : O::CrossSection;
    if (spec.exitKey == -1 && spec.combination == ReactionSpec::Combination::None) spec.observable = O::TotalCapture;
    spec.total = q == "sigma_total";
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
    if (angleFrame != "cm") {
      c.Error(codes::kNotBuilt, r, "", "an analyzing power at a lab angle: this engine computes it at c.m. angles");
      return false;
    }
    return true;
  }
  if (q == "phase_shift") {
    spec.observable = O::PhaseShift;
    return ReadPhaseShift(c, r, spec);
  }
  if (q == "legendre") {
    spec.observable = data ? O::AngularCoefficients : O::AngularDistribution;
    const std::string expansion = Word(r, "expansion", "a_over_a0");
    if (expansion != "a_over_a0" && (!data || expansion != "absolute")) {
      c.Error(codes::kNotBuilt, r, "expansion", "expansion " + expansion + ": this engine has a_over_a0" +
                                                    std::string(data ? " and absolute" : "") + " only");
      return false;
    }
    return true;
  }
  c.Error(codes::kNotBuilt, r, "quantity", "quantity " + q + ": this engine cannot set it up yet (a later step)");
  return false;
}

// The column layout the engine reads: E, angle (or order, or unused), value, error[, window_lo, window_hi], unused.
bool CheckColumns(Context &c, const ValidatedRecord &r, const Columns &cols, const std::string &secondRole) {
  std::vector<std::string> want = {r.quantity && r.quantity->name == "beta_delayed_spectrum" ? "E_decay" : "E_lab", secondRole,
                                   "value", "error"};
  bool ok = !cols.fromFile && cols.roles.size() >= 4;
  for (std::size_t i = 0; ok && i < cols.roles.size(); i++) {
    const ColumnRole &role = cols.roles[i];
    const std::string expected = i < 4 ? want[i] : i == 4 ? "window_lo" : i == 5 ? "window_hi" : "_";
    if (role.role != expected && !(i >= 4 && role.role == "_") && !(i == 1 && expected == "_" && role.role == "_")) ok = false;
    if (role.unit && role.unit->factor != 1.0) ok = false;
  }
  if (!ok) {
    std::string layout;
    for (const std::string &w : want) layout += (layout.empty() ? "" : " ") + w;
    c.Error(codes::kNotBuilt, r, "columns",
            "this engine reads the columns " + layout + "[ window_lo window_hi] in MeV and deg, in this order (other "
                                                        "layouts, units and 'from_file' are a later step)");
  }
  return ok;
}

// The data file exists and its first row holds every column the roles name.
bool CheckDataFile(Context &c, const ValidatedRecord &r, const std::string &file, std::size_t columns) {
  const bool absolute = !file.empty() && (file[0] == '/' || file[0] == '\\');
  const std::string path = absolute || c.opt_.projectDir.empty() ? file : c.opt_.projectDir + "/" + file;
  std::ifstream in(path);
  if (!in) {
    if (c.opt_.dataFiles == DataFiles::NotUsed) {  // an extrapolation run never opens it
      c.Warn(codes::kDataFileNotRead, r, "file", "cannot read " + path + "; not read: this run does not use data", path);
      return true;
    }
    c.Error(codes::kDataFile, r, "file", "cannot read " + path);
    return false;
  }
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
              path + ": its first row has " + std::to_string(n) + " columns, the roles name " + std::to_string(columns));
      return false;
    }
    return true;
  }
  c.Error(codes::kDataFile, r, "file", path + " holds no data row");
  return false;
}

// ---------------------------------------------------------------- shared nuisances (section 3.6)

void NormInto(const Nuisance &n, DatasetSpec &d) {
  if (n.mode == NuisanceMode::Profiled) {
    // A THM experiment's profiled scale (ESegment::IsProfiledNorm): a fitted norm, eliminated in closed form.
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
      d.hasNormCentre = true;
      d.normCentre = n.priorCentre.value;
      d.normCentreSource = "norm prior " + n.priorCentre.text;
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
      d.hasShiftCentre = true;
      d.shiftCentre = centre;
      d.shiftCentreSource = "shift offset prior " + n.priorCentre.text;
    }
  }
}

void SqrtInto(const Nuisance &n, DatasetSpec &d) {
  d.hasSqrtShift = true;
  d.energyShiftSqrt = n.value.value;
  d.varyEnergyShiftSqrt = n.mode == NuisanceMode::Fit;
  if (n.hasPrior) d.energyShiftSqrtError = n.priorWidth.value;
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
  ok = EngineUnit(c, r) && ok;
  d.absolute = Word(r, "expansion") == "absolute";
  if (const Value *w = r.Get("energy_window")) {
    d.minE = Context::InBase(w->lo, w->unit);
    d.maxE = Context::InBase(w->hi, w->unit);
  } else {
    d.minE = 0.0;
    d.maxE = std::numeric_limits<double>::max();
  }
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
    ok = CheckColumns(c, r, cols->columns, second) && ok;
    if (ok && m.active && c.opt_.readDataFiles) ok = CheckDataFile(c, r, d.dataFile, cols->columns.roles.size()) && ok;
  }
  // norm: its own nuisance, or a shared normalization record (its values in every member; the engine takes the group's
  // from its first loaded member, EData::ResolveNormGroups)
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
        // active member (EData::ResolveShiftFunctions refuses a function on a follower)
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
                                                 "kinematics are non-relativistic)");
      ok = false;
    }
    if (!d.covarianceFile.empty()) {
      c.Error(codes::kThmCovariance, r, "covariance", "a THM data set's chi^2 is its experiment's (profiled norm), not a sum a "
                                                      "covariance could amend");
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
  return ok;
}

// The grid of an extrapolation (decision F10: range lo hi step s is E_k = lo + k s).  `bounds` are its lo and hi.
bool ReadGrid(Context &c, const ValidatedRecord &r, const std::string &field, const Grid &g, std::vector<double> &values,
              double &lo, double &hi) {
  if (g.frame == "cm" && field == "energies") {
    c.Error(codes::kNotBuilt, r, field, "an energy grid in the c.m. frame: this engine takes lab energies (a later step)");
    return false;
  }
  if (g.kind == "list") {
    for (const Number &n : g.values) values.push_back(Context::InBase(n, g.unit));
  } else if (g.kind == "range" && g.hasStep) {
    lo = Context::InBase(g.lo, g.unit);
    hi = Context::InBase(g.hi, g.unit);
    const double step = Context::InBase(g.step, g.unit);
    const long n = step > 0 ? (long)std::floor((hi - lo) / step + 1e-9) : 0;
    for (long k = 0; k <= n; k++) values.push_back(lo + k * step);
    return true;
  } else if (g.kind == "range" && g.hasPoints) {
    lo = Context::InBase(g.lo, g.unit);
    hi = Context::InBase(g.hi, g.unit);
    for (long k = 0; k < g.points; k++) {
      const double t = g.points > 1 ? (double)k / (double)(g.points - 1) : 0.0;
      values.push_back(g.log ? lo * std::pow(hi / lo, t) : lo + t * (hi - lo));
    }
    return true;
  } else {
    c.Error(codes::kNotBuilt, r, field, "a grid of kind " + g.kind + ": this engine takes range and list grids (a later step)");
    return false;
  }
  if (!values.empty()) {
    lo = values.front();
    hi = values.back();
  }
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
  ok = EngineUnit(c, r) && ok;
  if (const Value *e = r.Get("energies")) ok = ReadGrid(c, r, "energies", e->grid, x.energies, x.minE, x.maxE) && ok;
  if (hasAngles) ok = ReadGrid(c, r, "angles", angles->grid, x.angles, x.minA, x.maxA) && ok;
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
      c.Error(codes::kThmRelativistic, r, "thm", "THM extrapolations are not supported with relativistic=1");
      ok = false;
    }
    if (const Value *w = r.Get("thm_weight")) c.out_.thm.extrapolationWeights.push_back({x.key, w->word});
  }
  return ok;
}

// ---------------------------------------------------------------- parameter names (SetupSpecs.h, ParameterNames)

// The names of the parameters the engine makes for each data set and extrapolation: the data set's own, a shared
// record's for the first active member of its group (the engine's leader rule, EData::Resolve*Groups), and for a
// later member, whose parameter the engine keeps fixed, <record>.<field>[<data set>].
void NameParameters(Context &c) {
  std::map<std::string, int> normLeader, shiftLeader, thicknessLeader;
  auto shared = [](std::map<std::string, int> &leaders, const std::string &label, int key, const std::string &field,
                   const std::string &dataset, const std::string &own) {
    if (label.empty()) return dataset + "." + own;
    const int leader = leaders.emplace(label, key).first->second;
    return label + "." + field + (leader == key ? "" : "[" + dataset + "]");
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
    m.active = Word(r, "active", "yes") != "no";
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
  NameParameters(c);
  ModelDataBuild modelData = BuildModelData(validated, levels.levels, levels.pairs);
  out.modelData = std::move(modelData.data);
  out.diagnostics.insert(out.diagnostics.end(), modelData.diagnostics.begin(), modelData.diagnostics.end());
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
