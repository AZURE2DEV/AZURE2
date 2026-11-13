// Schema::Build, part 2 (BuildContext.h): the effects of targets, beams and detectors (project_format3.md section 3.7)
// as the engine's TargetEffectSpecs, the sample thicknesses of transmission and neutron-yield data, and the abundances
// of the targets' isotopes.
//
// The engine sets up one effect per distinct combination of target, beam, detector and quadrature; every
// data set (and extrapolation) naming that combination gets a copy of it (WorkingCopy::SetUpTargetEffects), and the effect's
// parameters (screening U_e, Doppler temperature) are shared by the copies.  A target is one set of parameters
// (section 3.7): its areal density is shared through the engine's density group when several effects use it; a
// target whose screening or temperature is free must be used through one combination.  The effects are ordered by
// the first of their records in the file, (for a converted project the order of the old lines).
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include "azure2/experiment/BuildContext.h"
#include "azure2/experiment/Kernels.h"
#include "azure2/physics/Elements.h"

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
  if (v->type == ValueType::Integer || v->type == ValueType::Number) return v->number.text;
  if (v->type == ValueType::NumberWithUnit) return v->number.text + (v->unit ? v->unit->token : "");
  return v->word + (v->epoch.empty() ? "" : "@" + v->epoch);
}

struct Combination {
  const ValidatedRecord *target = nullptr, *beam = nullptr, *detector = nullptr;
  const Value *quadrature = nullptr;
  bool layerYield = false;  ///< a yield on a layer: per incident particle, not the target-averaged cross section
  bool decay = false;  ///< its data sets are beta-delayed spectra (no beam: the energy is the particle's in the decay)
  const PairInfo *projectile = nullptr;  ///< the entrance pair whose light particle crosses the target (a scaled stopping)
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

// The composition as a formula ("Ta 2 O 5" -> "Ta2O5"; a count of 1 left out), the label of a stopping scale.
std::string MaterialFormula(const Value &composition) {
  std::string out;
  for (const auto &part : composition.composition) out += part.first + (part.second.value == 1.0 ? "" : part.second.text);
  return out;
}

// The stopping record's source, composition, phase and compound correction, and its scale nuisance (SetupSpecs.h
// StoppingSpec).  The projectile, the key and the name of the scale are the combination's (Effects).
void StoppingInto(Context &c, const ValidatedRecord &t, const Value &stopping, StoppingSpec &s) {
  using Source = StoppingSpec::Source;
  s.source = stopping.word == "evaluation" ? Source::Evaluation : stopping.word == "table" ? Source::Table : Source::Formula;
  if (const Value *comp = t.Get("composition")) {
    for (const auto &part : comp->composition) s.composition.push_back({physics::ElementZ(part.first), part.second.value});
    s.material = MaterialFormula(*comp);
  }
  if (const Value *a = t.Get("active")) s.activeZ = a->nucleus.Z;
  if (s.source != Source::Formula) {
    bool present = false;
    for (const auto &zn : s.composition) present = present || (zn.first == s.activeZ && zn.second > 0.0);
    if (!present)
      c.Error(codes::kStoppingActive, t, "composition",
              "the active nucleus's element is not in the composition: the stopping per active nucleus divides the "
              "material's stopping by its atom fraction (energy_grid.md section 6)");
  }
  if (const Value *p = t.Get("phase")) s.phase = p->word;
  if (const Value *cc = t.Get("compound_correction")) s.compoundCorrection = cc->number.value;
  if (s.source == Source::Table) {
    if (const Value *f = stopping.Get("file")) s.tableFile = f->word;  // opened relative to the project (TargetStopping)
    if (const Value *u = stopping.Get("energy_unit")) s.tableEnergyUnit = u->word;
    if (const Value *u = stopping.Get("stopping_unit")) s.tableStoppingUnit = u->word;
  }
  const Value *scale = stopping.Get("scale");
  // left out: a formula has no scale; an evaluation or a table has one, fixed at 1
  s.hasScale = scale || s.source != Source::Formula;
  if (!scale) return;
  const Nuisance &n = scale->nuisance;
  s.scale = n.value.value;
  if (!(s.scale > 0.0)) c.Error(codes::kOutOfRange, t, "stopping", "the stopping scale must be > 0, not " + n.value.text);
  s.scaleVary = n.mode == NuisanceMode::Fit;
  s.scaleHasPrior = n.hasPrior;
  if (n.hasPrior) {
    if (n.priorCentre.value != n.value.value)
      c.Error(codes::kNotBuilt, t, "stopping", "a stopping-scale prior centred away from the value: this engine centres it on the value");
    s.scalePriorWidth = n.priorPercent ? n.priorWidth.value / 100.0 * n.priorCentre.value : n.priorWidth.value;
  }
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
      StoppingInto(c, t, *s, e.stopping);
      if (s->word != "formula") {
        // evaluation, table: the engine's stopping/ library (TargetStopping.h)
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

// A table kernel's parameter written as a nuisance (KernelParameter::fittable): its entry of the effect's apparatus
// parameters (SetupSpecs.h ApparatusParameterSpec), shared by every effect whose data sets name this record, and where
// its fitted value is written back.  Nothing here depends on the kernel.
void ApparatusParameterInto(Context &c, const ValidatedRecord &record, const std::string &field, const KernelParameter &p,
                            int slot, const Nuisance &n, TargetEffectSpec &e) {
  ApparatusParameterSpec a;
  a.row = "kernel";
  a.slot = slot;
  a.key = record.name + "|" + field + "|" + p.name;
  a.owner = {{"record", record.name}, {"field", field}, {"parameter", p.name}};
  a.unit = p.unit;
  a.value = ToBase(n.value.value, n.unit);
  a.vary = n.mode == NuisanceMode::Fit;
  if (!(a.value > 0.0)) c.Error(codes::kOutOfRange, record, field, "the " + p.name + " must be > 0, not " + n.value.text);
  if (n.hasPrior) {
    if (n.priorCentre.value != n.value.value)
      c.Error(codes::kNotBuilt, record, field, "a prior centred away from the value: this engine centres it on the value");
    a.priorWidth = n.priorPercent ? n.priorWidth.value / 100.0 * a.value : ToBase(n.priorWidth.value, n.unit);
  }
  if (n.hasGridMax) {
    a.gridMax = ToBase(n.gridMax.value, n.unit);
    if (!(a.gridMax >= a.value))
      c.Error(codes::kOutOfRange, record, field, "grid_max " + n.gridMax.text + " is below the value " + n.value.text);
  }
  c.out_.parameterSources[a.owner[0].second + "." + field + "." + p.name] = {record.name, field, 0, p.name};
  e.apparatusParameters.push_back(a);
}

// A detector's energy resolution: with a beam profile the profile's TPC resolution; on a beta-delayed spectrum (no
// beam) a Gaussian fold of the model in the observed particle energy, the engine's convolution kernel on the lab energy
// (for a decay at rest the energy the data give); elsewhere not built (a beam's spread states a Gaussian there).
void DetectorInto(Context &c, const ValidatedRecord &d, bool decay, TargetEffectSpec &e) {
  if (const Value *q = d.Get("attenuation")) {
    e.qCoefficients = true;
    for (const Number &n : q->numbers) e.qCoefficientValues.push_back(n.value);
  }
  const Value *r = d.Get("energy_resolution");
  // A table kernel (Kernels.h): its name and its parameters, in the order of its row, in their engine units.
  if (const Kernel *k = r ? KernelOfForm("detector", "energy_resolution", r->word) : nullptr) {
    e.kernel = k->name;
    for (const KernelParameter &p : k->parameters) {
      const Value *v = r->Get(p.name);
      // fitted, or with a prior or a grid_max (as an areal density): a parameter; a fixed value is a constant of the kernel
      const Nuisance *n = (v && v->type == ValueType::Nuisance) ? &v->nuisance : nullptr;
      if (n && (n->mode == NuisanceMode::Fit || n->hasPrior || n->hasGridMax)) {
        ApparatusParameterInto(c, d, "energy_resolution", p, (int)e.kernelParameters.size(), *n, e);
        e.kernelParameters.push_back(e.apparatusParameters.back().value);
      } else if (n) {
        e.kernelParameters.push_back(ToBase(n->value.value, n->unit));
        if (!(e.kernelParameters.back() > 0.0)) c.Error(codes::kOutOfRange, d, "energy_resolution", "the " + p.name + " must be > 0, not " + n->value.text);
      } else {
        e.kernelParameters.push_back(v ? Context::InBase(v->number, v->unit) : 0.0);
      }
    }
  } else if (r) {
    const Value *sigma = r->Get("sigma");
    const Value *frame = r->Get("frame");
    if (frame && frame->word != "lab")
      c.Error(codes::kNotBuilt, d, "energy_resolution", "a detector resolution in the c.m. frame: this engine takes it in the lab");
    else if (decay && !e.beamProfile) {
      e.convolution = true;
      if (sigma) e.sigma = Context::InBase(sigma->number, sigma->unit);
      return;
    } else if (!e.beamProfile)
      c.Error(codes::kNotBuilt, d, "energy_resolution",
              "this engine applies a detector resolution (lab) only with a beam profile (the TPC sigma of the profile) or "
              "on a beta-delayed spectrum; a Gaussian on the beam energy is the beam's spread");
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

void QuadratureInto(const Value *q, TargetEffectSpec &e) {
  if (!q) return;
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

/*!
 * The folds a ratio to the Rutherford or Mott cross section cannot take (data_kinds.md section 4.3): the engine folds
 * the Coulomb cross section with the numerator's kernel on the numerator's sub-points (src/CoulombReference.cpp), which
 * covers the energy kernels (beam spread, target layer, beam profile, UDR).  Not covered:
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
}

// A target whose stopping depends on the projectile through the engine: an evaluation, a table, or any stopping with a
// scale (one scale per material and projectile).
bool ScaledStopping(const ValidatedRecord *t) {
  const Value *s = t ? t->Get("stopping") : nullptr;
  return s && (s->word != "formula" || s->Get("scale"));
}

// The entrance pair whose light particle crosses the target, for a data set: one particle for every term, an ion the
// stopping library computes (Z = 1..92; stopping/'s SRIM tables).  nullptr after a refusal.
const PairInfo *StoppingProjectile(Context &c, const Member &m) {
  const PairInfo *first = nullptr;
  for (int key : m.termEntranceKeys) {
    const PairInfo *p = c.Pair(key);
    if (!p) return nullptr;  // an unknown pair key is refused where the reaction is resolved
    if (p->kind != 0 || p->lightZ < 1 || p->lightZ > 92) {
      c.Error(codes::kStoppingProjectile, *m.record, "target",
              "the stopping of its target needs a charged projectile of Z = 1..92 (SRIM's tables), and the entrance pair " +
                  std::to_string(key) + " has none (Z = " + std::to_string(p->lightZ) + ")");
      return nullptr;
    }
    if (first && (first->lightZ != p->lightZ || first->lightA != p->lightA)) {
      c.Error(codes::kStoppingProjectile, *m.record, "target",
              "its terms have different projectiles; the stopping of the target is one projectile's: give one data set per "
              "projectile");
      return nullptr;
    }
    if (!first) first = p;
  }
  return first;
}

// The projectile of an effect and the key and name of its stopping scale (SetupSpecs.h StoppingSpec): an evaluation is
// shared by every target of the same composition, phase and compound correction; a table by the targets naming the same
// file; a formula is its target's own.
void ScaleOf(const ValidatedRecord &t, const PairInfo &p, StoppingSpec &s) {
  s.projectileZ = p.lightZ;
  s.projectileA = p.lightA;
  s.projectileMass = p.lightMass;
  const std::string projectile = std::to_string(p.lightA) + physics::ElementSymbol(p.lightZ);
  std::ostringstream key;
  key.precision(17);
  std::string material;
  switch (s.source) {
    case StoppingSpec::Source::Evaluation:
      key << "evaluation|" << s.material << "|" << s.phase << "|" << s.compoundCorrection;
      material = s.material + (s.phase == "gas" ? "_gas" : "");
      if (s.compoundCorrection != 1.0) {
        std::ostringstream cc;
        cc << s.compoundCorrection;
        material += "_C" + cc.str();
      }
      break;
    case StoppingSpec::Source::Table:
      key << "table|" << s.tableFile << "|" << s.material << "|" << s.tableEnergyUnit << "|" << s.tableStoppingUnit;
      material = t.name;
      break;
    case StoppingSpec::Source::Formula:
      key << "formula|" << t.name;
      material = t.name;
      break;
  }
  key << "|" << p.lightZ << "," << p.lightA;
  s.scaleKey = key.str();
  s.scaleMaterial = material;
  s.scaleProjectile = projectile;
  s.scaleName = "stopping_scale." + material + "." + projectile;
}

// The targets that share a stopping scale must write it alike (one parameter, one start value, one prior); where its
// fitted value is written: the first of them.
void SharedScales(Context &c, const std::vector<Combination> &combos, SpecList<TargetEffectSpec> &out) {
  std::map<std::string, std::size_t> first;
  for (std::size_t i = 0; i < out.items.size(); i++) {
    const StoppingSpec &s = out.items[i].spec.stopping;
    if (!s.hasScale || s.scaleKey.empty()) continue;
    auto it = first.emplace(s.scaleKey, i).first;
    if (it->second == i) {
      c.out_.parameterSources[s.scaleName] = {combos[i].target->name, "stopping", 0, "scale"};
      continue;
    }
    const StoppingSpec &a = out.items[it->second].spec.stopping;
    if (a.scale != s.scale || a.scaleVary != s.scaleVary || a.scaleHasPrior != s.scaleHasPrior ||
        a.scalePriorWidth != s.scalePriorWidth)
      c.Error(codes::kStoppingScaleConflict, *combos[i].target, "stopping",
              "its stopping scale is " + s.scaleName + ", shared with target " + combos[it->second].target->name +
                  " (one scale per material and projectile, decision D6), and the two write it differently: write the "
                  "same scale field on both");
  }
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
    x.decay = r.quantity && r.quantity->name == "beta_delayed_spectrum";
    CoulombRatioFolds(c, r, x);
    // A yield on a layer folds the same sub-points to another quantity (reactions per incident particle): an effect of
    // its own, whose areal density is shared with the target's other effects as below.
    x.layerYield = x.target && x.target->Get("kind") && x.target->Get("kind")->word == "layer" && r.quantity &&
                   r.quantity->name == "yield";
    // A stopping that depends on the projectile (an evaluation, a table, a scale): one projectile per effect.
    if (ScaledStopping(x.target) && !(x.projectile = StoppingProjectile(c, m))) continue;
    x.signature = (x.target ? x.target->name : "") + "|" + (x.beam ? x.beam->name : "") + "|" +
                  (x.detector ? x.detector->name : "") + "|" + Text(x.quadrature) +
                  (x.layerYield ? "|yield" : "") + (x.decay ? "|decay" : "") +
                  (x.projectile ? "|" + std::to_string(x.projectile->lightZ) + "," + std::to_string(x.projectile->lightA) : "");
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
               ", detector " + (x.detector ? x.detector->name : std::string("-")) +
               (x.layerYield ? "; yield per incident particle on the layer: 1e-24 int_{E_out}^{E_0} <sigma>/epsilon dE, "
                               "int_{E_out}^{E_0} dE/epsilon = N" : "");
    if (x.target) {
      e.targetName = x.target->name;
      TargetInto(c, *x.target, *x.target, e);
      if (x.projectile) ScaleOf(*x.target, *x.projectile, e.stopping);
      e.layerYield = x.layerYield;
      if (usesOfTarget[x.target] > 1) {
        if ((x.target->Has("screening") || FreeNuisance(x.target, "doppler")) && !reported.count(x.target) &&
            reported.insert(x.target).second)
          c.Error(codes::kFoldConflict, *x.target, "",
                  "its data sets differ in beam, detector or quadrature, and its screening or temperature is a "
                  "parameter: this engine shares such a parameter only among data sets folded alike");
        if (e.densityFit) e.densityGroup = x.target->name;
      }
    }
    if (x.beam) BeamInto(c, *x.beam, e);
    if (x.detector) DetectorInto(c, *x.detector, x.decay, e);
    QuadratureInto(x.quadrature, e);
    out.items.push_back(item);
  }
  SharedScales(c, combos, out);
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
    // a thickness of 0 would make the engine read the data set as a plain cross section (Segment: n > 0 marks it)
    if (!(n.value.value > 0.0)) {
      c.Error(codes::kOutOfRange, *target, "areal_density",
              "the sample thickness of a " + q + " must be > 0 atoms/b, not " + n.value.text);
      continue;
    }
    // the engine's sample thickness and its prior width are in atoms/b (Segment)
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

// The total cross section of a target with isotopes (quantity sigma_total with a sum): the abundance-weighted sum of
// each nucleus's total by the optical theorem, one elastic term per nucleus.  A sum means nothing else: refused
// without isotopes on the target, with two terms on one pair, or with a factor on a term (the abundance is the
// weight).  A term is elastic by its label (a term of a sum carries no pair keys); a single reaction with keys by
// its keys (Build.cpp PairsDecideElastic).  An inactive record is checked when it is
// switched on (BuildOptions::allRecordsActive): the converter writes an old inactive line without its sample.
void NaturalTotals(Context &c, const std::vector<Member> &members) {
  for (const Member &m : members) {
    const ValidatedRecord &r = *m.record;
    if (!m.active || !r.quantity || r.quantity->name != "sigma_total") continue;
    const ReactionSpec &spec = m.isData ? c.out_.data.datasets.items[m.key - 1].spec.reaction
                                        : c.out_.extrapolations.extrapolations.items[m.key - 1].spec.reaction;
    std::vector<ComponentSpec> terms = {ComponentSpec{spec.entranceKey, spec.exitKey}};
    terms.insert(terms.end(), spec.components.begin(), spec.components.end());
    if (spec.combination != ReactionSpec::Combination::Sum) continue;
    const Value *name = r.Get("target");
    const ValidatedRecord *target = name ? c.Record(name->word, "target") : nullptr;
    if (!target || !target->Has("isotope")) {
      c.Error(codes::kReactionNotAllowed, r, "reaction",
              "a sum of total cross sections is the total of a target with isotopes (their abundances weight the terms): "
              "give the data set a target with 'isotope' fields, or one reaction");
      continue;
    }
    std::set<int> pairs;
    for (const ComponentSpec &t : terms) {
      if (!pairs.insert(t.entranceKey).second) {
        c.Error(codes::kReactionNotAllowed, r, "reaction",
                "two terms on pair " + std::to_string(t.entranceKey) + ": a target's total has one term per nucleus");
        break;
      }
      if (t.scaling != 1.0) {
        c.Error(codes::kReactionNotAllowed, r, "reaction",
                "a factor on a term of a target's total: the abundance of the term's nucleus is its weight");
        break;
      }
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
  NaturalTotals(c, c.data_);
  NaturalTotals(c, c.extrapolations_);
}

}  // namespace build

std::string AbundanceParameterName(const std::string &sample, const std::string &isotope, int z, int a) {
  if (!isotope.empty()) return sample + ".isotope." + isotope;
  return "abundance_" + (sample.empty() ? std::string() : sample + "_") + std::to_string(z) + "_" + std::to_string(a);
}

}  // namespace experiment
}  // namespace azure2
