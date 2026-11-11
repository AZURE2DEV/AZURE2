// Schema::Build, part 2 (BuildContext.h): the effects of targets, beams and detectors (project_format3.md section 3.7)
// as the engine's TargetEffectSpecs, the sample thicknesses of transmission and neutron-yield data, and the abundances
// of the targets' isotopes.
//
// The engine sets up one effect per distinct combination of target, beam, detector, quadrature and fold windows; every
// data set (and extrapolation) naming that combination gets a copy of it (EData::SetUpTargetEffects), and the effect's
// parameters (screening U_e, Doppler temperature) are shared by the copies.  A target is one set of parameters
// (section 3.7): its areal density is shared through the engine's density group when several effects use it; a
// target whose screening or temperature is free must be used through one combination.  The effects are ordered by
// the first of their records in the file, which is the order of the old <targetInt> lines the converter wrote them from.
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include "azure2/experiment/BuildContext.h"

namespace azure2 {
namespace experiment {
namespace build {

namespace {

std::string Text(const Value *v) {
  if (!v) return "";
  if (v->type == ValueType::Form) {
    std::string s = v->word;
    for (const ValuePart &p : v->parts) s += " " + p.name + "=" + Text(&p.value);
    return s;
  }
  if (v->type == ValueType::Windows) {
    std::string s = v->windows.automatic ? "auto" : "";
    for (const auto &w : v->windows.windows) s += " " + w.lo.text + "-" + w.hi.text;
    return s + ":" + v->windows.tolerance.text + ":" + v->windows.margin.text + (v->windows.unit ? v->windows.unit->token : "");
  }
  if (v->type == ValueType::Integer || v->type == ValueType::Number) return v->number.text;
  if (v->type == ValueType::NumberWithUnit) return v->number.text + (v->unit ? v->unit->token : "");
  return v->word + (v->epoch.empty() ? "" : "@" + v->epoch);
}

struct Combination {
  const ValidatedRecord *target = nullptr, *beam = nullptr, *detector = nullptr;
  const Value *quadrature = nullptr, *windows = nullptr;
  std::string signature;
  int position = 0;  ///< the first of its records in the file
  std::vector<int> keys;
};

bool Folds(const ValidatedRecord *t) {
  return t && (t->Has("doppler") || t->Has("screening") ||
               (t->Get("kind") && (t->Get("kind")->word == "layer" || t->Get("kind")->word == "infinite")));
}

bool FreeNuisance(const ValidatedRecord *t, const std::string &field) {
  const Value *v = t ? t->Get(field) : nullptr;
  if (v && v->type == ValueType::Form && !v->parts.empty()) v = &v->parts[0].value;
  return v && (v->nuisance.mode == NuisanceMode::Fit || v->nuisance.hasPrior || v->nuisance.hasGridMax);
}

// The text of a nuisance's prior width as the old reader kept it ("10%" or "0.5"), in the engine's unit of its
// dimension: the text as written for a percentage or a value in that unit, else the converted number.
std::string WidthText(const Nuisance &n) {
  if (!n.hasPrior) return "";
  if (n.priorPercent) return n.priorWidth.text + "%";
  if (!n.unit || n.unit->factor == 1.0) return n.priorWidth.text;
  std::ostringstream s;
  s.precision(17);
  s << ToBase(n.priorWidth.value, n.unit);
  return s.str();
}

// A formula's value (x in MeV) in the engine's unit: the expression as written for the engine's unit, else
// (expression)*factor.
std::string FormulaInBase(const std::string &expression, const Unit *unit) {
  if (!unit || unit->factor == 1.0) return expression;
  std::ostringstream s;
  s.precision(17);
  s << "(" << expression << ")*" << unit->factor;
  return s.str();
}

void TargetInto(Context &c, const ValidatedRecord &t, const ValidatedRecord &user, TargetEffectSpec &e) {
  const std::string kind = t.Get("kind") ? t.Get("kind")->word : "none";
  if (kind == "layer" || kind == "infinite") {
    e.targetIntegration = true;
    if (const Value *d = t.Get("areal_density")) {
      const Nuisance &n = d->nuisance;
      e.density = ToBase(n.value.value, n.unit);
      if (n.mode == NuisanceMode::Fit || n.hasPrior || n.hasGridMax) {
        e.densityFit = true;
        e.densityVary = n.mode == NuisanceMode::Fit;
        e.densitySigmaText = n.hasPrior ? WidthText(n) : "0";
        if (n.hasGridMax) e.densityMax = ToBase(n.gridMax.value, n.unit);
        if (n.hasPrior && n.priorCentre.value != n.value.value)
          c.Error(codes::kNotBuilt, t, "areal_density", "a prior centred away from the value: this engine centres it on the value");
      }
    }
    if (const Value *s = t.Get("stopping")) {
      if (s->word != "formula") {
        c.Error(codes::kNotBuilt, t, "stopping", "stopping " + s->word + ": this engine takes a formula (Phase 5)");
      } else if (const Value *f = s->Get("expression")) {
        const Value *u = s->Get("unit");
        e.stoppingPowerEquation = FormulaInBase(f->formula.expression, u ? u->unit : nullptr);
        for (std::size_t i = 0; i < f->formula.params.size(); i++) {
          if (f->formula.params[i].name != "a" + std::to_string(i))
            c.Error(codes::kNotBuilt, t, "stopping", "this engine names the formula's parameters a0, a1, ... in order");
          e.stoppingPowerParameters.push_back(f->formula.params[i].value.value);
        }
      }
    }
    if (const Value *s = t.Get("straggling")) {
      e.hasStragglingFlag = true;
      e.straggling = s->word != "none";
      e.hasStragglingCoefficient = true;
      const Value *v = s->Get("value");
      // the engine's c_s is in keV^1/2 (TargetEffect: sigma_cm = c_s sqrt(dE_cm [keV]))
      e.stragglingCoefficient = s->word == "coefficient" && v ? ToUnit(v->number.value, v->unit, "keV^1/2") : 0.0;
    }
  }
  if (const Value *d = t.Get("doppler")) {
    const Value *temp = d->Get("temperature");
    if (temp) {
      const Nuisance &n = temp->nuisance;
      // the effective temperature of the free-gas model (K) is positive; T <= 0 has no Doppler width (TargetEffect)
      if (!(n.value.value > 0.0))
        c.Error(codes::kOutOfRange, t, "doppler", "free_gas: the effective temperature must be > 0 K, not " + n.value.text);
      e.doppler = true;
      e.dopplerTemperature = n.value.value;
      if (n.mode == NuisanceMode::Fit || n.hasPrior || n.hasGridMax) {
        e.dopplerFit = true;
        e.dopplerVaryText = n.mode == NuisanceMode::Fit ? "1" : "0";
        e.dopplerSigmaText = n.hasPrior ? n.priorWidth.text : "0";
        e.dopplerMaxText = n.hasGridMax ? n.gridMax.text : "";
      }
    }
  }
  if (const Value *s = t.Get("screening")) {
    const Nuisance &n = s->nuisance;
    // the engine's U_e and its prior width are in eV (TargetEffect)
    e.screening = true;
    e.screeningUe = ToUnit(n.value.value, n.unit, "eV");
    e.screeningVary = n.mode == NuisanceMode::Fit;
    e.screeningSigma = n.hasPrior ? (n.priorPercent ? n.priorWidth.value / 100.0 * e.screeningUe : ToUnit(n.priorWidth.value, n.unit, "eV")) : 0.0;
  }
  (void)user;
}

void BeamInto(Context &c, const ValidatedRecord &b, TargetEffectSpec &e) {
  if (const Value *s = b.Get("spread")) {
    const Value *frame = s->Get("frame");
    const bool lab = frame && frame->word == "lab";
    if (s->word == "formula") {
      const Value *f = s->Get("expression");
      e.convolutionCoefficients = true;
      const Value *u = s->Get("unit");
      e.convolutionEquation = f ? FormulaInBase(f->formula.expression, u ? u->unit : nullptr) : "";
      if (f)
        for (const auto &p : f->formula.params) e.convolutionCoefficientValues.push_back(p.value.value);
      e.sigmaLab = lab;
    } else if (const Value *sigma = s->Get("sigma")) {
      if (!lab) c.Error(codes::kNotBuilt, b, "spread", "a constant sigma in the c.m. frame: this engine takes it in the lab");
      e.convolution = true;
      e.sigma = Context::InBase(sigma->number, sigma->unit);
    }
  }
  for (const Value *p : b.GetAll("profile")) {
    e.beamProfile = true;
    TargetEffectSpec::ProfileComponent pc;
    if (const Value *v = p->Get("xi")) pc.xi = Context::InBase(v->number, v->unit);
    if (const Value *v = p->Get("omega")) pc.omega = Context::InBase(v->number, v->unit);
    if (const Value *v = p->Get("alpha")) pc.alpha = v->number.value;
    if (const Value *v = p->Get("weight")) pc.weight = v->number.value;
    e.beamProfileComponents.push_back(pc);
  }
  if (const Value *v = b.Get("truncation")) e.beamTruncation = v->number.value;
  if (const Value *v = b.Get("detailed_balance")) e.beamDetailedBalance = v->word == "yes";
}

void DetectorInto(Context &c, const ValidatedRecord &d, TargetEffectSpec &e) {
  if (const Value *q = d.Get("attenuation")) {
    e.qCoefficients = true;
    for (const Number &n : q->numbers) e.qCoefficientValues.push_back(n.value);
  }
  if (const Value *r = d.Get("energy_resolution")) {
    const Value *sigma = r->Get("sigma");
    const Value *frame = r->Get("frame");
    if (!e.beamProfile || (frame && frame->word != "lab"))
      c.Error(codes::kNotBuilt, d, "energy_resolution",
              "this engine applies a detector resolution (lab) only with a beam profile (the TPC sigma of the profile)");
    if (sigma) e.beamTpcSigma = Context::InBase(sigma->number, sigma->unit);
  }
  if (const Value *u = d.Get("udr")) {
    e.udr = true;
    if (const Value *v = u->Get("file")) e.udrFile = v->word;
    if (const Value *v = u->Get("flight_path")) e.udrFlightPath = Context::InBase(v->number, v->unit);
    if (const Value *v = u->Get("burst")) e.udrBurstFwhm = Context::InBase(v->number, v->unit);
    if (const Value *v = u->Get("channel")) e.udrChannelWidth = Context::InBase(v->number, v->unit);
    if (const Value *v = u->Get("centred")) e.udrCentred = v->word == "yes";
  }
}

void QuadratureInto(const Value *q, const Value *w, TargetEffectSpec &e) {
  if (q) {
    if (const Value *v = q->Get("points")) e.numIntegrationPoints = (int)v->integer;
    // given or not (left out, the engine's default: 20, and a pure Gaussian convolution's smaller core)
    if (const Value *v = q->Get("width_multiplier")) {
      e.hasResonanceWidthMultiplier = true;
      e.resonanceWidthMultiplier = v->number.value;
    }
    if (const Value *v = q->Get("per_width")) {
      e.hasPointsPerWidth = true;
      e.pointsPerWidth = v->number.value;
    }
  }
  if (w) {
    e.hasRanges = true;
    const Value *win = w->Get("windows");
    if (win) {
      const FoldWindows &f = win->windows;
      e.autoWindows = f.automatic;
      e.hasAutoWindowTolerance = f.hasTolerance;
      e.autoWindowTolerance = f.tolerance.value;
      e.hasAutoWindowMargin = f.hasMargin;
      e.autoWindowMargin = f.margin.value;
      for (const auto &x : f.windows) e.ranges.push_back({Context::InBase(x.lo, f.unit), Context::InBase(x.hi, f.unit)});
    }
    if (const Value *v = w->Get("blend")) e.transitionWidth = Context::InBase(v->number, v->unit);
    if (const Value *v = w->Get("auto_tolerance")) e.autoTolerance = v->number.value;
  }
}

/*!
 * The folds a ratio to the Rutherford or Mott cross section cannot take (data_kinds.md section 4.3): the engine folds
 * the Coulomb cross section with the numerator's kernel on the numerator's sub-points (src/CoulombReference.cpp), which
 * covers the energy kernels (beam spread, target layer, beam profile, UDR) and the edge blend.  Not covered:
 *  - an automatic tolerance, which picks the folded or the bare value of each observable from its own sub-point
 *    values, so numerator and denominator could take different choices;
 *  - free-gas Doppler broadening, whose fitted temperature has an analytic derivative through the numerator only;
 *  - electron screening, a factor on the nuclear cross section, not part of a ratio to Coulomb scattering;
 *  - detector attenuation coefficients, which attenuate the Legendre terms of the nuclear amplitude and leave the
 *    Coulomb amplitude, which has no finite Legendre series (sigma_R ~ theta^-4), as it is.
 */
void CoulombRatioFolds(Context &c, const ValidatedRecord &r, const Combination &x) {
  const std::string q = r.quantity ? r.quantity->name : "";
  if (q != "ratio_to_rutherford" && q != "ratio_to_mott") return;
  const std::string what = "a ratio to the Coulomb cross section ";
  if (x.target && x.target->Has("doppler"))
    c.Error(codes::kNotBuilt, r, "target",
            what + "with a Doppler-broadened target: this engine folds the Coulomb cross section with the energy kernels only");
  if (x.target && x.target->Has("screening"))
    c.Error(codes::kNotBuilt, r, "target",
            what + "with electron screening: screening multiplies the nuclear cross section; the ratio is defined without it");
  if (x.detector && x.detector->Has("attenuation"))
    c.Error(codes::kNotBuilt, r, "detector",
            what + "with attenuation coefficients: they attenuate the Legendre terms of the nuclear amplitude only, and the "
                   "Coulomb cross section has no finite Legendre series");
  if (x.windows && x.windows->Get("auto_tolerance"))
    c.Error(codes::kNotBuilt, r, "fold_windows",
            what + "with an automatic tolerance: it chooses the folded or the bare value of each observable from that "
                   "observable's sub-points, which numerator and denominator could choose differently");
}

// The effects of the data sets (isData) or of the extrapolations.
SpecList<TargetEffectSpec> Effects(Context &c, const std::vector<Member> &members) {
  std::vector<Combination> combos;
  for (const Member &m : members) {
    if (!m.active) continue;
    const ValidatedRecord &r = *m.record;
    Combination x;
    const Value *t = r.Get("target"), *b = r.Get("beam"), *d = r.Get("detector");
    x.target = t ? c.Record(t->word, "target") : nullptr;
    if (t && !t->epoch.empty()) c.Error(codes::kNotBuilt, r, "target", "a target epoch: this engine has one areal density per target");
    if (!Folds(x.target)) x.target = nullptr;
    x.beam = b ? c.Record(b->word, "beam") : nullptr;
    x.detector = d ? c.Record(d->word, "detector") : nullptr;
    if (!x.target && !x.beam && !x.detector) continue;
    x.quadrature = r.Get("quadrature");
    x.windows = r.Get("fold_windows");
    CoulombRatioFolds(c, r, x);
    x.signature = (x.target ? x.target->name : "") + "|" + (x.beam ? x.beam->name : "") + "|" +
                  (x.detector ? x.detector->name : "") + "|" + Text(x.quadrature) + "|" + Text(x.windows);
    x.position = 1 << 30;
    for (const ValidatedRecord *p : {x.target, x.beam, x.detector})
      if (p) x.position = std::min(x.position, c.Position(p));
    auto it = std::find_if(combos.begin(), combos.end(), [&](const Combination &y) { return y.signature == x.signature; });
    if (it == combos.end()) {
      combos.push_back(x);
      it = combos.end() - 1;
    }
    it->keys.push_back(m.key);
  }
  std::stable_sort(combos.begin(), combos.end(), [](const Combination &a, const Combination &b) { return a.position < b.position; });
  // a target is one set of parameters: through several combinations, only its areal density can be shared
  std::map<const ValidatedRecord *, int> usesOfTarget;
  for (const Combination &x : combos)
    if (x.target) usesOfTarget[x.target]++;
  SpecList<TargetEffectSpec> out;
  std::set<const ValidatedRecord *> reported;
  for (const Combination &x : combos) {
    SpecList<TargetEffectSpec>::Item item;
    TargetEffectSpec &e = item.spec;
    e.active = true;
    e.datasetKeys = x.keys;
    e.source = "target " + (x.target ? x.target->name : std::string("-")) + ", beam " + (x.beam ? x.beam->name : std::string("-")) +
               ", detector " + (x.detector ? x.detector->name : std::string("-"));
    if (x.target) {
      e.targetName = x.target->name;
      TargetInto(c, *x.target, *x.target, e);
      if (usesOfTarget[x.target] > 1) {
        if ((x.target->Has("screening") || FreeNuisance(x.target, "doppler")) && !reported.count(x.target) &&
            reported.insert(x.target).second)
          c.Error(codes::kFoldConflict, *x.target, "",
                  "its data sets differ in beam, detector, quadrature or fold windows, and its screening or temperature is a "
                  "parameter: this engine shares such a parameter only among data sets folded alike");
        if (e.densityFit) e.densityGroup = x.target->name;
      }
    }
    if (x.beam) BeamInto(c, *x.beam, e);
    if (x.detector) DetectorInto(c, *x.detector, e);
    QuadratureInto(x.quadrature, x.windows, e);
    out.items.push_back(item);
  }
  return out;
}

// Neutron-sample thicknesses (transmission, neutron_yield): the areal density of a target of kind sample, in atoms/b.
void SampleThickness(Context &c, std::vector<Member> &members, bool isData) {
  std::map<std::string, int> users;
  for (const Member &m : members)
    if (m.active && m.record->Get("target")) users[m.record->Get("target")->word]++;
  for (std::size_t i = 0; i < members.size(); i++) {
    const Member &m = members[i];
    const std::string q = m.record->quantity ? m.record->quantity->name : "";
    if (q != "transmission" && q != "neutron_yield") continue;
    const Value *t = m.record->Get("target");
    const ValidatedRecord *target = t ? c.Record(t->word, "target") : nullptr;
    const Value *d = target ? target->Get("areal_density") : nullptr;
    if (!d) continue;
    const Nuisance &n = d->nuisance;
    // a thickness of 0 would make the engine read the data set as a plain cross section (ESegment: n > 0 marks it)
    if (!(n.value.value > 0.0)) {
      c.Error(codes::kOutOfRange, *target, "areal_density",
              "the sample thickness of a " + q + " must be > 0 atoms/b, not " + n.value.text);
      continue;
    }
    // the engine's sample thickness and its prior width are in atoms/b (ESegment)
    ReactionSpec &r = isData ? c.out_.data.datasets.items[i].spec.reaction : c.out_.extrapolations.extrapolations.items[i].spec.reaction;
    (q == "transmission" ? r.transmission : r.yieldThickness) = ToUnit(n.value.value, n.unit, "atoms/b");
    if (!isData) continue;
    DatasetSpec &ds = c.out_.data.datasets.items[i].spec;
    ds.thicknessFit = n.mode == NuisanceMode::Fit || n.hasPrior;
    ds.thicknessVary = n.mode == NuisanceMode::Fit;
    ds.thicknessSigma = !n.hasPrior ? 0.0
                        : n.priorPercent ? n.priorWidth.value / 100.0 * (q == "transmission" ? r.transmission : r.yieldThickness)
                                         : ToUnit(n.priorWidth.value, n.unit, "atoms/b");
    if (ds.thicknessFit && users[target->name] > 1) ds.thicknessGroup = target->name;
    if (ds.thicknessFit) ds.names.thickness = target->name + ".areal_density";  // a follower's: Build's NameParameters
  }
}

// The isotopes of a target: one abundance row per compound nucleus its data sets reach from that isotope.
void Abundances(Context &c) {
  std::vector<const ValidatedRecord *> targets;
  for (const ValidatedRecord &r : c.v_.records)
    if (r.kind == "target" && r.Has("isotope")) targets.push_back(&r);
  if (targets.empty()) return;
  AbundancesSpec &a = c.out_.abundances.abundances;
  a.present = true;
  for (const ValidatedRecord *t : targets) {
    std::vector<const Member *> users;
    for (const std::vector<Member> *list : {&c.data_, &c.extrapolations_})
      for (const Member &m : *list)
        if (m.record->Get("target") && m.record->Get("target")->word == t->name) users.push_back(&m);
    for (const Member *m : users) {
      ReactionSpec &r = m->isData ? c.out_.data.datasets.items[m->key - 1].spec.reaction
                                  : c.out_.extrapolations.extrapolations.items[m->key - 1].spec.reaction;
      r.sample = t->name;
    }
    const std::vector<const Value *> isotopes = t->GetAll("isotope");
    for (std::size_t k = 0; k < isotopes.size(); k++) {
      const Value *iso = isotopes[k];
      const Value *nucleus = iso->Get("nucleus"), *abundance = iso->Get("abundance");
      if (!nucleus || !abundance) continue;
      // where the abundance is written, under the parameter's name (the engine's, CNuc::AbundanceParamName)
      c.out_.parameterSources[AbundanceParameterName(t->name, nucleus->nucleus.text, 0, 0)] = {t->name, "isotope", (int)k, "abundance"};
      std::set<std::pair<int, int>> compounds;
      std::set<std::pair<int, int>> lights;  // the light particles of the target's data sets
      for (const Member *m : users)
        for (int key : m->termEntranceKeys)
          if (const PairInfo *p = c.Pair(key)) {
            lights.insert({p->lightZ, p->lightA});
            if (p->heavyZ == nucleus->nucleus.Z && p->heavyA == nucleus->nucleus.A)
              compounds.insert({p->lightZ + p->heavyZ, p->lightA + p->heavyA});
          }
      // an isotope no data set reaches is still part of the composition (the sum-to-one prior): its compound nucleus is
      // the one the pairs of <levels> form with the data sets' light particle
      if (compounds.empty())
        for (const PairInfo &p : c.pairs_)
          if (p.kind == 0 && lights.count({p.lightZ, p.lightA}) && p.heavyZ == nucleus->nucleus.Z && p.heavyA == nucleus->nucleus.A)
            compounds.insert({p.lightZ + p.heavyZ, p.lightA + p.heavyA});
      for (const auto &zA : compounds) {
        AbundanceRowSpec row;
        row.sample = t->name;
        row.isotope = nucleus->nucleus.text;
        row.z = zA.first;
        row.a = zA.second;
        const Nuisance &n = abundance->nuisance;
        row.value = n.value.value;
        row.vary = n.mode == NuisanceMode::Fit;
        row.sigma = n.hasPrior ? n.priorWidth.value : 0.0;
        a.rows.push_back(row);
      }
    }
    if (const Value *s = t->Get("isotope_sum")) {
      AbundanceRowSpec row;
      row.sample = t->name;
      row.isSum = true;
      const Value *w = s->Get("width");
      row.sumSigma = s->word == "prior" && w ? w->number.value : 0.0;
      a.rows.push_back(row);
    }
  }
}

}  // namespace

void BuildEffects(Context &c) {
  c.out_.data.targetEffects = Effects(c, c.data_);
  c.out_.extrapolations.targetEffects = Effects(c, c.extrapolations_);
  SampleThickness(c, c.data_, true);
  SampleThickness(c, c.extrapolations_, false);
  Abundances(c);
}

}  // namespace build

std::string AbundanceParameterName(const std::string &sample, const std::string &isotope, int z, int a) {
  if (!isotope.empty()) return sample + ".isotope." + isotope;
  return "abundance_" + (sample.empty() ? std::string() : sample + "_") + std::to_string(z) + "_" + std::to_string(a);
}

}  // namespace experiment
}  // namespace azure2
