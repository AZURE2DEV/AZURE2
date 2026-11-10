// Model data: a validated datum record resolved against <levels> (ModelDatum.h; project_format3.md section 3.10).
#include "azure2/experiment/ModelDatum.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace azure2 {
namespace experiment {

const std::vector<std::pair<std::string, ModelQuantity>> &ModelQuantities() {
  static const std::vector<std::pair<std::string, ModelQuantity>> rows = {
      {"resonance_strength", ModelQuantity::ResonanceStrength},
      {"resonance_energy", ModelQuantity::ResonanceEnergy},
      {"total_width", ModelQuantity::TotalWidth},
      {"partial_width", ModelQuantity::PartialWidth},
      {"branching_ratio", ModelQuantity::BranchingRatio},
      {"anc", ModelQuantity::Anc},
      {"scattering_length", ModelQuantity::ScatteringLength},
  };
  return rows;
}

bool ModelDataBuild::Ok() const {
  for (const Diagnostic &d : diagnostics)
    if (d.severity == Severity::Error) return false;
  return true;
}

namespace {

std::string Fmt(double x) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.6g", x);
  return buf;
}

std::string JpiText(int twoJ, int parity) {
  const std::string j = twoJ % 2 ? std::to_string(twoJ) + "/2" : std::to_string(twoJ / 2);
  return j + (parity < 0 ? "-" : "+");
}

/// An excitation energy of <levels> (MeV) against one written in a label: equal to the digits written (4.439 is
/// 4.43891); not written means the ground state (0).
bool SameExcitation(double levelsMeV, bool written, const experiment::Number &n) {
  if (!written) return std::fabs(levelsMeV) < 1e-9;
  const std::string &t = n.text;
  const std::size_t dot = t.find('.');
  int decimals = 0;
  if (t.find_first_of("eE") != std::string::npos) return std::fabs(levelsMeV - n.value) <= 1e-9 * std::fabs(n.value);
  if (dot != std::string::npos) decimals = (int)(t.size() - dot - 1);
  return std::fabs(levelsMeV - n.value) <= 0.5 * std::pow(10.0, -decimals) * (1.0 + 1e-12);
}

bool SameParticle(int Z, int A, const Particle &p) { return p.Z == Z && p.A == A; }

/// The keys of the pairs with these particles and states: light particle `light`, heavy `heavy` in state `state`.
std::vector<int> PairsOf(const Particle &light, const Particle &heavy, const ResidualState &state,
                         const std::vector<PairInfo> &pairs) {
  std::vector<int> keys;
  const bool photon = light.light == "g";
  for (const PairInfo &p : pairs) {
    if (photon != (p.kind == 10) || p.kind == 20) continue;
    if (!photon && (!SameParticle(p.lightZ, p.lightA, light) ||
                    !SameExcitation(p.lightExcitation, light.excited, light.excitation)))
      continue;
    if (!SameParticle(p.heavyZ, p.heavyA, heavy)) continue;
    if (!SameExcitation(p.heavyExcitation, state.kind == ResidualState::Excited, state.energy)) continue;
    keys.push_back(p.key);
  }
  return keys;
}

const PairInfo *PairByKey(int key, const std::vector<PairInfo> &pairs) {
  for (const PairInfo &p : pairs)
    if (p.key == key) return &p;
  return nullptr;
}

class Builder {
 public:
  Builder(const std::vector<LevelInfo> &levels, const std::vector<PairInfo> &pairs, ModelDataBuild &out)
      : levels_(levels), pairs_(pairs), out_(out) {}

  void Add(const ValidatedRecord &r) {
    record_ = &r;
    ok_ = true;
    ModelDatum d;
    d.name = r.name;
    d.line = r.line;
    const Value *q = r.Get("quantity");
    const Value *value = r.Get("value");
    const Value *error = r.Get("error");
    if (!q || !value || !error || !value->unit || !error->unit) return;  // the validator said why
    d.quantityName = q->word;
    bool known = false;
    for (const auto &row : ModelQuantities())
      if (row.first == q->word) {
        d.quantity = row.second;
        known = true;
      }
    if (!known) return;  // a datum quantity of the table without an evaluation: repo_consistency check 27
    d.unit = value->unit->token;
    d.unitFactor = value->unit->factor;
    d.value = value->number.value * value->unit->factor;
    d.error = error->number.value * error->unit->factor;
    if (const Value *s = r.Get("source")) d.source = s->word;
    if (const Value *im = r.Get("imag")) {
      const Value *ie = r.Get("imag_error");
      if (!im->unit || !ie || !ie->unit) return;
      d.hasImag = true;
      d.imag = im->number.value * im->unit->factor;
      d.imagError = ie->number.value * ie->unit->factor;
    }
    if (d.quantity == ModelQuantity::ScatteringLength) {
      ScatteringLengthPair(d);
    } else {
      const LevelInfo *level = Level(d);
      if (!level) return;
      LevelPairs(d, *level);
    }
    if (ok_) out_.data.push_back(d);
  }

 private:
  void Report(Severity severity, const char *code, const std::string &field, const std::string &message) {
    Diagnostic x;
    x.severity = severity;
    x.code = code;
    x.section = "modelData";
    x.recordKind = "datum";
    x.record = record_->name;
    x.field = field;
    x.line = record_->line;
    for (const ValidatedField &f : record_->fields)
      if (f.key == field) {
        x.line = f.line;
        x.column = f.col;
      }
    x.message = "datum '" + record_->name + "'" + (field.empty() ? "" : ", field '" + field + "'") + ": " + message;
    out_.diagnostics.push_back(x);
    if (severity == Severity::Error) ok_ = false;
  }

  /// The level of the datum: its ID and J^pi (F13 (a)).
  const LevelInfo *Level(ModelDatum &d) {
    const Value *id = record_->Get("level");
    const Value *jpi = record_->Get("J");
    if (!id || !jpi) {
      ok_ = false;
      return nullptr;
    }
    d.levelId = (int)id->integer;
    d.jpi = jpi->spin.text;
    std::vector<const LevelInfo *> withId, matching;
    for (const LevelInfo &l : levels_)
      for (int i : l.ids)
        if (i == d.levelId) {
          withId.push_back(&l);
          if (l.twoJ == jpi->spin.twice && l.parity == jpi->spin.parity) matching.push_back(&l);
        }
    if (withId.empty()) {
      Report(Severity::Error, codes::kUnknownLevel, "level",
             "no active level of <levels> has the ID " + std::to_string(d.levelId));
      return nullptr;
    }
    if (matching.empty()) {
      const LevelInfo &l = *withId.front();
      Report(Severity::Error, codes::kLevelJpi, "J",
             "level " + std::to_string(d.levelId) + " is " + JpiText(l.twoJ, l.parity) + " (E_x = " + Fmt(l.energy) +
                 " MeV); the datum says " + jpi->spin.text);
      return nullptr;
    }
    if (matching.size() > 1) {
      Report(Severity::Error, codes::kUnknownLevel, "level",
             "two levels of J^pi " + jpi->spin.text + " have the ID " + std::to_string(d.levelId) + " (E_x = " +
                 Fmt(matching[0]->energy) + " and " + Fmt(matching[1]->energy) + " MeV): give them distinct IDs");
      return nullptr;
    }
    d.group = matching[0]->group;
    d.level = matching[0]->level;
    return matching[0];
  }

  /// One pair key from a list of candidates, or 0 with a refusal.
  int OnePair(const std::vector<int> &keys, const std::string &field, const std::string &what) {
    if (keys.empty()) {
      Report(Severity::Error, codes::kUnknownPair, field, "no pair of <levels> is " + what);
      return 0;
    }
    if (keys.size() > 1) {
      Report(Severity::Error, codes::kAmbiguousPair, field,
             "pairs " + std::to_string(keys[0]) + " and " + std::to_string(keys[1]) + " of <levels> are both " + what +
                 "; the label does not tell them apart");
      return 0;
    }
    return keys[0];
  }

  bool IsChannel(const LevelInfo &level, int pairKey, int l = -1, int twoS = -1) const {
    for (const LevelChannelInfo &c : level.channels)
      if (c.pairKey == pairKey && (l < 0 || c.l == l) && (twoS < 0 || c.twoS == twoS)) return true;
    return false;
  }

  void NeedChannel(const LevelInfo &level, int pairKey, const std::string &field, int l = -1, int twoS = -1) {
    if (!ok_ || IsChannel(level, pairKey, l, twoS)) return;
    const PairInfo *p = PairByKey(pairKey, pairs_);
    std::string what = "pair " + std::to_string(pairKey) + (p ? " (" + p->label + ")" : "");
    if (l >= 0) what += " with l = " + std::to_string(l) + " and s = " + (twoS % 2 ? std::to_string(twoS) + "/2" : std::to_string(twoS / 2));
    Report(Severity::Error, codes::kNotAChannel, field,
           what + " is no channel of level " + std::to_string(record_->Get("level")->integer) + " (" +
               JpiText(level.twoJ, level.parity) + ")");
  }

  int ReactionPair(const SingleReaction &r, bool entrance) {
    const std::string what = entrance ? "the entrance " + r.projectile.text + " + " + r.target.text
                                      : "the exit " + r.ejectile.text + " + " + r.residual.text;
    if (entrance) return OnePair(PairsOf(r.projectile, r.target, r.targetState, pairs_), "reaction", what);
    return OnePair(PairsOf(r.ejectile, r.residual, r.state, pairs_), "reaction", what);
  }

  void LevelPairs(ModelDatum &d, const LevelInfo &level) {
    switch (d.quantity) {
      case ModelQuantity::ResonanceStrength:
      case ModelQuantity::BranchingRatio: {
        const Value *reaction = record_->Get("reaction");
        if (!reaction || reaction->reaction.terms.empty()) {
          ok_ = false;
          return;
        }
        const SingleReaction &r = reaction->reaction.terms[0].reaction;
        if (reaction->reaction.op != ReactionExpr::Single || r.state.kind == ResidualState::Total) {
          Report(Severity::Error, codes::kReactionNotAllowed, "reaction",
                 "a " + d.quantityName + " names one entrance and one exit pair (no sum, ratio or [total])");
          return;
        }
        d.entrancePair = ReactionPair(r, true);
        d.exitPair = ReactionPair(r, false);
        if (!ok_) return;
        NeedChannel(level, d.entrancePair, "reaction");
        NeedChannel(level, d.exitPair, "reaction");
        if (!ok_ || d.quantity != ModelQuantity::ResonanceStrength) return;
        const PairInfo *in = PairByKey(d.entrancePair, pairs_);
        const double Er = level.energy - in->threshold;  // MeV, c.m. of the entrance pair
        if (Er <= 0.0) {
          Report(Severity::Error, codes::kOutOfRange, "level",
                 "level " + std::to_string(d.levelId) + " lies " + Fmt(-Er) + " MeV below the " + in->label +
                     " threshold (E_r = " + Fmt(Er) + " MeV): a subthreshold level has no strength; use an anc");
          return;
        }
        Warnings(level, Er);
        return;
      }
      case ModelQuantity::ResonanceEnergy: {
        const Value *frame = record_->Get("frame");
        if (!frame) {
          ok_ = false;
          return;
        }
        d.frame = frame->word;
        if (d.frame == "E_x") return;
        PairField(d, level);
        return;
      }
      case ModelQuantity::TotalWidth:
        return;
      case ModelQuantity::PartialWidth:
        PairField(d, level);
        return;
      case ModelQuantity::Anc: {
        const Value *l = record_->Get("l");
        const Value *s = record_->Get("s");
        if (!l || !s) {
          ok_ = false;
          return;
        }
        if (s->word == "coherent") {
          Report(Severity::Error, codes::kFieldNotApplicable, "s", "an ANC belongs to one channel spin, not 'coherent'");
          return;
        }
        PairField(d, level);
        if (!ok_) return;
        d.l = (int)l->integer;
        d.twoS = s->Get("spin") ? s->Get("spin")->spin.twice : -1;
        NeedChannel(level, d.pair, "s", d.l, d.twoS);
        if (!ok_) return;
        const PairInfo *p = PairByKey(d.pair, pairs_);
        if (level.energy >= p->threshold)
          Report(Severity::Error, codes::kFieldNotApplicable, "pair",
                 "level " + std::to_string(d.levelId) + " lies " + Fmt(level.energy - p->threshold) + " MeV above the " +
                     p->label + " threshold: an ANC belongs to a closed channel (give a partial_width)");
        return;
      }
      case ModelQuantity::ScatteringLength:
        return;
    }
  }

  void PairField(ModelDatum &d, const LevelInfo &level) {
    const Value *pair = record_->Get("pair");
    if (!pair) {
      ok_ = false;
      return;
    }
    d.pair = OnePair(PairsOf(pair->pair.light, pair->pair.heavy, pair->pair.state, pairs_), "pair",
                     pair->pair.light.text + " + " + pair->pair.heavy.text);
    if (ok_) NeedChannel(level, d.pair, "pair");
  }

  void Warnings(const LevelInfo &level, double Er) {
    if (level.width < 0.0) return;
    if (level.width > 0.01 * Er)
      Report(Severity::Warning, codes::kBroadLevel, "level",
             "Gamma / E_r = " + Fmt(100.0 * level.width / Er) + " % > 1 % (Gamma = " + Fmt(1e3 * level.width) +
                 " keV, E_r = " + Fmt(1e3 * Er) + " keV, from <levels>): omega gamma of eq. (8) and the area of the "
                 "resonance differ (data_kinds.md section 4.7)");
    for (const LevelInfo &o : levels_)
      if (o.group == level.group && o.level != level.level && std::fabs(o.energy - level.energy) < 10.0 * level.width)
        Report(Severity::Warning, codes::kNearbyLevel, "level",
               "another level of J^pi " + JpiText(level.twoJ, level.parity) + " lies at E_x = " + Fmt(o.energy) +
                   " MeV, within 10 Gamma (" + Fmt(1e3 * level.width) + " keV) of this one at " + Fmt(level.energy) +
                   " MeV: the two interfere, and omega gamma of eq. (8) is not the area of the resonance");
  }

  void ScatteringLengthPair(ModelDatum &d) {
    const Value *pair = record_->Get("pair");
    const Value *s = record_->Get("s");
    if (!pair || !s) {
      ok_ = false;
      return;
    }
    d.pair = OnePair(PairsOf(pair->pair.light, pair->pair.heavy, pair->pair.state, pairs_), "pair",
                     pair->pair.light.text + " + " + pair->pair.heavy.text);
    if (!ok_) return;
    const PairInfo *p = PairByKey(d.pair, pairs_);
    if (p->kind != 0) {
      Report(Severity::Error, codes::kNotAChannel, "pair", "a scattering length belongs to a particle pair");
      return;
    }
    if (s->word == "coherent") {
      d.coherent = true;
      if (p->identical)
        Report(Severity::Error, codes::kFieldNotApplicable, "s",
               "the coherent combination is not defined for a pair of identical particles");
      return;
    }
    d.twoS = s->Get("spin") ? s->Get("spin")->spin.twice : -1;
    // s-wave channel spins: |j1 - j2| ... j1 + j2; identical particles, s even (l + s even).
    const int lo = std::abs(p->twoJLight - p->twoJHeavy), hi = p->twoJLight + p->twoJHeavy;
    const bool allowed = d.twoS >= lo && d.twoS <= hi && (d.twoS - lo) % 2 == 0 && (!p->identical || d.twoS % 4 == 0);
    if (!allowed)
      Report(Severity::Error, codes::kNotAChannel, "s",
             "pair " + std::to_string(p->key) + " (" + p->label + ") has no s-wave channel with channel spin " +
                 s->Get("spin")->spin.text);
  }

  const std::vector<LevelInfo> &levels_;
  const std::vector<PairInfo> &pairs_;
  ModelDataBuild &out_;
  const ValidatedRecord *record_ = nullptr;
  bool ok_ = true;
};

}  // namespace

ModelDataBuild BuildModelData(const ValidationResult &validated, const std::vector<LevelInfo> &levels,
                              const std::vector<PairInfo> &pairs) {
  ModelDataBuild out;
  Builder b(levels, pairs, out);
  for (const ValidatedRecord &r : validated.records)
    if (r.kind == "datum") b.Add(r);
  return out;
}

}  // namespace experiment
}  // namespace azure2
