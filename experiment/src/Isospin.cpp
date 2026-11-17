// The records of <isospin> (Isospin.h; BuildContext.h, BuildIsospin): each names its nucleus, pairs and levels as the
// project writes them; here they are resolved against the pairs and level IDs of <levels> and checked against one
// another.  What needs the channels and the parameters is the set-up's (the engine's IsospinSet).
#include "azure2/experiment/Isospin.h"

#include <cmath>
#include <set>

#include "azure2/experiment/BuildContext.h"

namespace azure2 {
namespace experiment {
namespace build {

namespace {

std::string PairText(const ParticlePair &p) {
  std::string s = p.light.text + (p.light.excited ? "[" + p.light.excitation.text + "]" : "") + "+" + p.heavy.text;
  if (p.state.kind == ResidualState::Excited) s += "[" + p.state.energy.text + "]";
  return s;
}

/// The key of the pair field `field` of `r`, or 0 after an error.
int PairKeyOf(Context &c, const ValidatedRecord &r, const std::string &field) {
  const Value *pair = r.Get(field);
  if (!pair) return 0;
  const std::string text = PairText(pair->pair);
  int key = 0;
  if (!ResolvePair(c, r, field, "the pair " + text, pair->pair.light, pair->pair.heavy, pair->pair.state, text, "pair_key", key))
    return 0;
  return key;
}

/// A nuisance as the specs carry it; false after an error (a profiled value, a prior without a width).
bool NumberOf(Context &c, const ValidatedRecord &r, const std::string &field, IsospinNumber &out) {
  const Value *v = r.Get(field);
  if (!v) return false;
  const Nuisance &n = v->nuisance;
  if (n.mode == NuisanceMode::Profiled || !n.hasValue) {
    c.Error(codes::kIsospinRecord, r, field, "the " + field + " of an " + r.kind + " record is fixed or fit, not profiled");
    return false;
  }
  auto base = [&](double x) { return n.unit ? ToBase(x, n.unit) : x; };
  out.value = base(n.value.value);
  out.fit = n.mode == NuisanceMode::Fit;
  if (n.hasPrior) {
    out.hasPrior = true;
    out.priorCentre = base(n.priorCentre.value);
    out.priorWidth = n.priorPercent ? std::fabs(out.priorCentre) * n.priorWidth.value / 100.0 : base(n.priorWidth.value);
    if (!(out.priorWidth > 0.0)) {
      c.Error(codes::kIsospinRecord, r, field,
              "the prior of the " + field + " has no width (a percentage of a centre 0, or a width that is not > 0)");
      return false;
    }
  }
  return true;
}

bool HasLevel(const LevelTables &tables, int id) {
  for (const LevelInfo &level : tables.levels)
    for (const int i : level.ids)
      if (i == id) return true;
  return false;
}

}  // namespace

void BuildIsospin(Context &c, const LevelTables &tables) {
  IsospinSpec &spec = c.out_.isospin;
  const ValidatedRecord *firstRecord = nullptr, *reference = nullptr;
  for (const ValidatedRecord &r : c.v_.records) {
    if (r.kind.compare(0, 8, "isospin_") != 0) continue;
    if (!firstRecord) firstRecord = &r;
    if (r.kind == "isospin_reference" || r.kind == "isospin_shift") {
      const Value *nucleus = r.Get("nucleus");
      if (!nucleus || !r.Get("threshold")) continue;  // reported by the validator
      IsospinNucleusEntry e;
      e.Z = nucleus->nucleus.Z;
      e.A = nucleus->nucleus.A;
      e.reference = r.kind == "isospin_reference";
      e.record = r.name;
      e.line = r.line;
      e.pairKey = PairKeyOf(c, r, "threshold");
      if (!e.pairKey) continue;
      const PairInfo *info = c.Pair(e.pairKey);
      if (info && (info->kind != 0 || info->lightZ + info->heavyZ != e.Z || info->lightA + info->heavyA != e.A)) {
        c.Error(codes::kIsospinPair, r, "threshold",
                "pair " + std::to_string(e.pairKey) + " (" + PairLabel(*info) + ") is not a particle pair of " +
                    NuclideLabel(e.Z, e.A) + ": the energies of the nucleus's isospin levels are measured from the "
                    "threshold of one of its particle pairs");
        continue;
      }
      if (!e.reference && !NumberOf(c, r, "shift", e.shift)) continue;
      const IsospinNucleusEntry *first = nullptr;
      for (const IsospinNucleusEntry &o : spec.nuclei)
        if (o.Z == e.Z && o.A == e.A) first = &o;
      if (first) {
        c.Error(codes::kIsospinRecord, r, "nucleus",
                NuclideLabel(e.Z, e.A) + " has a record already (" + first->record + ", line " +
                    std::to_string(first->line) + "); one isospin_reference or isospin_shift per nucleus");
        continue;
      }
      if (e.reference) {
        if (reference) {
          c.Error(codes::kIsospinRecord, r, "nucleus",
                  "a second isospin_reference (the first: " + reference->name + ", line " +
                      std::to_string(reference->line) + "); one nucleus is the reference, the others have isospin_shift records");
          continue;
        }
        reference = &r;
      }
      spec.nuclei.push_back(e);
    } else if (r.kind == "isospin_level") {
      const Value *level = r.Get("level");
      const Value *T = r.Get("T");
      if (!level || !T) continue;
      IsospinLevelEntry e;
      e.levelId = (int)level->integer;
      e.twoT = T->spin.twice;
      e.record = r.name;
      e.line = r.line;
      if (T->spin.parity != 0 || e.twoT < 0) {
        c.Error(codes::kIsospinRecord, r, "T", "the isospin T is a non-negative integer or half-integer, without a parity");
        continue;
      }
      const Value *sE = r.Get("soft_energy");
      const Value *sG = r.Get("soft_amplitude");
      if ((sE != nullptr) != (sG != nullptr)) {
        c.Error(codes::kIsospinRecord, r, sE ? "soft_energy" : "soft_amplitude",
                "a soft tie has both widths: soft_energy (of the level energies) and soft_amplitude (of the reduced-width "
                "amplitudes, MeV^1/2)");
        continue;
      }
      if (sE) {
        e.soft = true;
        e.sigmaE = Context::InBase(sE->number, sE->unit);
        e.sigmaGamma = sG->number.value;
      }
      if (spec.Declares(e.levelId)) {
        c.Error(codes::kIsospinRecord, r, "level", "level ID " + std::to_string(e.levelId) + " has an isospin_level record already");
        continue;
      }
      if (!HasLevel(tables, e.levelId)) {
        c.Error(codes::kIsospinRecord, r, "level", "no active level of <levels> has the ID " + std::to_string(e.levelId));
        continue;
      }
      spec.levels.push_back(e);
    } else if (r.kind == "isospin_mix") {
      const Value *level = r.Get("level");
      const Value *l = r.Get("l");
      const Value *s = r.Get("s");
      if (!level || !l || !s || !r.Get("pair")) continue;
      IsospinMixEntry e;
      e.levelId = (int)level->integer;
      e.l = (int)l->integer;
      e.twoS = s->spin.twice;
      e.record = r.name;
      e.line = r.line;
      e.pairKey = PairKeyOf(c, r, "pair");
      if (!e.pairKey) continue;
      if (!NumberOf(c, r, "amplitude", e.amplitude)) continue;
      spec.mixes.push_back(e);
    } else if (r.kind == "isospin_factor") {
      const Value *level = r.Get("level");
      const Value *factor = r.Get("factor");
      if (!level || !factor || !r.Get("pair")) continue;
      IsospinFactorEntry e;
      e.levelId = (int)level->integer;
      e.value = factor->number.value;
      e.record = r.name;
      e.line = r.line;
      e.pairKey = PairKeyOf(c, r, "pair");
      if (!e.pairKey) continue;
      spec.factors.push_back(e);
    } else if (r.kind == "isospin_pair") {
      const Value *t1 = r.Get("t_light");
      const Value *t2 = r.Get("t_heavy");
      if (!t1 || !t2 || !r.Get("pair")) continue;
      IsospinPairEntry e;
      e.twoT1 = t1->spin.twice;
      e.twoT2 = t2->spin.twice;
      e.record = r.name;
      e.line = r.line;
      if (t1->spin.parity != 0 || t2->spin.parity != 0 || e.twoT1 < 0 || e.twoT2 < 0) {
        c.Error(codes::kIsospinRecord, r, t1->spin.parity != 0 || e.twoT1 < 0 ? "t_light" : "t_heavy",
                "an isospin is a non-negative integer or half-integer, without a parity");
        continue;
      }
      e.pairKey = PairKeyOf(c, r, "pair");
      if (!e.pairKey) continue;
      const PairInfo *info = c.Pair(e.pairKey);
      if (info && info->kind != 0) {
        c.Error(codes::kIsospinPair, r, "pair",
                "pair " + std::to_string(e.pairKey) + " (" + PairLabel(*info) + ") is a " + (info->kind == 10 ? "photon" : "beta") +
                    " pair; the isospin relations are between the channels of particle pairs");
        continue;
      }
      bool twice = false;
      for (const IsospinPairEntry &o : spec.pairs) twice = twice || o.pairKey == e.pairKey;
      if (twice) {
        c.Error(codes::kIsospinRecord, r, "pair", "pair " + std::to_string(e.pairKey) + " has an isospin_pair record already");
        continue;
      }
      spec.pairs.push_back(e);
    }
  }
  if (!firstRecord) return;
  // The records against one another (each refusal at the record that cannot stand).
  if (!reference && !spec.levels.empty())
    c.Error(codes::kIsospinRecord, *firstRecord, "",
            "<isospin> has no isospin_reference record: one nucleus is the reference, whose level energies and "
            "amplitudes are the free parameters");
  for (const ValidatedRecord &r : c.v_.records) {
    if (r.kind != "isospin_mix" && r.kind != "isospin_factor") continue;
    const Value *level = r.Get("level");
    if (!level) continue;
    const int id = (int)level->integer;
    const IsospinLevelEntry *declared = nullptr;
    for (const IsospinLevelEntry &e : spec.levels)
      if (e.levelId == id) declared = &e;
    bool built = false;
    for (const IsospinMixEntry &e : spec.mixes) built = built || (e.record == r.name && r.kind == "isospin_mix");
    for (const IsospinFactorEntry &e : spec.factors) built = built || (e.record == r.name && r.kind == "isospin_factor");
    if (!built) continue;  // refused above
    if (!declared)
      c.Error(codes::kIsospinRecord, r, "level", "level ID " + std::to_string(id) + " has no isospin_level record");
    else if (declared->soft && r.kind == "isospin_mix")
      c.Error(codes::kIsospinRecord, r, "level",
              "level ID " + std::to_string(id) + " has a soft tie (isospin_level " + declared->record +
                  "), which frees its forbidden channels already; an isospin_mix belongs to a hard tie");
  }
}

}  // namespace build
}  // namespace experiment
}  // namespace azure2
