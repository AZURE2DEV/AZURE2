// A gamma_decay record resolved against <levels> (GammaDecay.h; docs/upgrade/gamma_parametrisation.md section 3.1).
#include "azure2/experiment/GammaDecay.h"

#include <cmath>
#include <cstdio>

#include "azure2/experiment/Schema.h"
#include "azure2/experiment/Units.h"

namespace azure2 {
namespace experiment {

bool GammaDecayBuild::Ok() const {
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

std::string PairText(const ParticlePair &p) {
  std::string s = p.light.text + "+" + p.heavy.text;
  if (p.state.kind == ResidualState::Excited) s += "[" + p.state.energy.text + "]";
  return s;
}

/// A nuisance value in the engine unit (dimensionless numbers have no unit).
GammaNumber NumberOf(const Nuisance &n) {
  GammaNumber g;
  auto base = [&](double x) { return n.unit ? ToBase(x, n.unit) : x; };
  g.value = base(n.value.value);
  g.fit = n.mode == NuisanceMode::Fit;
  if (n.hasPrior) {
    g.hasPrior = true;
    g.priorCentre = base(n.priorCentre.value);
    g.priorWidth = n.priorPercent ? std::fabs(g.priorCentre) * n.priorWidth.value / 100.0 : base(n.priorWidth.value);
  }
  return g;
}

class Builder {
 public:
  Builder(const LevelTables &tables, GammaDecayBuild &out) : tables_(tables), out_(out) {}

  void Add(const ValidatedRecord &r) {
    record_ = &r;
    ok_ = true;
    GammaDecaySpec s;
    s.name = r.name;
    s.line = r.line;
    const LevelInfo *level = Level(s);
    if (!level) return;
    for (const GammaDecaySpec &other : out_.specs)
      if (other.group == s.group && other.level == s.level)
        Report(codes::kGammaRecord, "level", "level " + std::to_string(s.levelId) + " already has the record '" + other.name +
                                                 "': one gamma_decay record per level");
    // The level's photon final states, in the order of the group's channels.
    for (std::size_t c = 0; c < level->channels.size(); c++) {
      const PairInfo *p = PairByKey(level->channels[c].pairKey, tables_.pairs);
      if (!p || p->kind != 10) continue;
      GammaFinalState *f = nullptr;
      for (GammaFinalState &x : s.finals)
        if (x.pairKey == p->key) f = &x;
      if (!f) {
        s.finals.emplace_back();
        f = &s.finals.back();
        f->pairKey = p->key;
        f->label = PairLabel(*p);
      }
      f->channels.push_back((int)c + 1);
      f->L.push_back(level->channels[c].l);
    }
    if (s.finals.empty()) {
      Report(codes::kGammaRecord, "level", "level " + std::to_string(s.levelId) + " has no photon channel");
      return;
    }
    for (GammaFinalState &f : s.finals)
      if (f.L.size() == 2 && f.L[0] > f.L[1]) {
        std::swap(f.L[0], f.L[1]);
        std::swap(f.channels[0], f.channels[1]);
      }
    Mixture(s);
    if (const Value *w = r.Get("total_width")) s.totalWidth = NumberOf(w->nuisance);
    if (s.totalWidth.value < 0.0) Report(codes::kOutOfRange, "total_width", "the total radiative width is negative");
    for (const Value *b : r.GetAll("branching")) {
      GammaFinalState *f = Final(s, b->Get("final"), "branching");
      if (!f) continue;
      f->role = GammaFinalState::Role::Listed;
      f->branching = NumberOf(b->Get("value")->nuisance);
      if (f->branching.value < 0.0 || f->branching.value > 1.0)
        Report(codes::kOutOfRange, "branching", "the branching to " + f->label + " is " + Fmt(f->branching.value) +
                                                    ", outside [0, 1]");
    }
    if (const Value *rem = r.Get("remainder")) {
      GammaFinalState *f = Final(s, rem->Get("final"), "remainder");
      if (f) {
        f->role = GammaFinalState::Role::Remainder;
        s.remainder = (int)(f - s.finals.data());
      }
      if (const Value *prior = rem->Get("prior")) {
        if (prior->numbers.size() != 2 || !(prior->numbers[1].value > 0.0))
          Report(codes::kGammaRecord, "remainder", "the remainder's prior is '<centre> <width>' with a width > 0");
        else {
          s.hasRemainderPrior = true;
          s.remainderPriorCentre = prior->numbers[0].value;
          s.remainderPriorWidth = prior->numbers[1].value;
        }
      }
    }
    for (const Value *m : r.GetAll("mixing")) {
      GammaFinalState *f = Final(s, m->Get("final"), "mixing", /*allowRepeat=*/true);
      if (!f) continue;
      if (f->L.size() != 2) {
        Report(codes::kMixingOneMultipole, "mixing", f->label + " has one multipolarity (L = " + std::to_string(f->L[0]) +
                                                         ") in level " + std::to_string(s.levelId) +
                                                         "'s channels: a mixing ratio needs L and L+1");
        continue;
      }
      if (f->hasMixing) Report(codes::kGammaRecord, "mixing", f->label + " has two mixing lines");
      f->hasMixing = true;
      f->mixing = NumberOf(m->Get("value")->nuisance);
    }
    for (const Value *g : r.GetAll("sign")) {
      GammaFinalState *f = Final(s, g->Get("final"), "sign", /*allowRepeat=*/true);
      if (f) f->sign = g->Get("sign")->word == "-" ? -1 : 1;
    }
    Branchings(s);
    if (ok_) out_.specs.push_back(s);
  }

 private:
  void Report(const char *code, const std::string &field, const std::string &message) {
    Diagnostic x;
    x.severity = Severity::Error;
    x.code = code;
    x.section = "gammaDecays";
    x.recordKind = "gamma_decay";
    x.record = record_->name;
    x.field = field;
    x.line = record_->line;
    for (const ValidatedField &f : record_->fields)
      if (f.key == field) {
        x.line = f.line;
        x.column = f.col;
      }
    x.message = "gamma_decay '" + record_->name + "'" + (field.empty() ? "" : ", field '" + field + "'") + ": " + message;
    out_.diagnostics.push_back(x);
    ok_ = false;
  }

  /// The level: its ID and J^pi, as a datum names it (F13 (a)).
  const LevelInfo *Level(GammaDecaySpec &s) {
    const Value *id = record_->Get("level");
    const Value *jpi = record_->Get("J");
    if (!id || !jpi) return nullptr;
    s.levelId = (int)id->integer;
    std::vector<const LevelInfo *> withId, matching;
    for (const LevelInfo &l : tables_.levels)
      for (int i : l.ids)
        if (i == s.levelId) {
          withId.push_back(&l);
          if (l.twoJ == jpi->spin.twice && l.parity == jpi->spin.parity) matching.push_back(&l);
        }
    if (withId.empty()) {
      Report(codes::kUnknownLevel, "level", "no active level of <levels> has the ID " + std::to_string(s.levelId));
      return nullptr;
    }
    if (matching.size() != 1) {
      Report(matching.empty() ? codes::kLevelJpi : codes::kUnknownLevel, "J",
             matching.empty() ? "level " + std::to_string(s.levelId) + " (E_x = " + Fmt(withId[0]->energy) +
                                    " MeV) is not " + jpi->spin.text
                              : "two levels of J^pi " + jpi->spin.text + " have the ID " + std::to_string(s.levelId));
      return nullptr;
    }
    s.group = matching[0]->group;
    s.level = matching[0]->level;
    s.levelEnergy = matching[0]->energy;
    return matching[0];
  }

  /// A level written with this record carries no photon width of its own in <levels> (one source per number).
  void Mixture(const GammaDecaySpec &s) {
    for (std::size_t row = 0; row < tables_.rowSlots.size(); row++) {
      const LevelRowSlot &slot = tables_.rowSlots[row];
      if (slot.group != s.group || slot.level != s.level || slot.value == 0.0) continue;
      for (const GammaFinalState &f : s.finals)
        for (std::size_t k = 0; k < f.channels.size(); k++)
          if (f.channels[k] == slot.channel)
            Report(codes::kGammaMixture, "level",
                   "row " + std::to_string(row + 1) + " of <levels> gives level " + std::to_string(s.levelId) + " the L = " +
                       std::to_string(f.L[k]) + " width " + Fmt(slot.value) + " eV to " + f.label +
                       ", and this record gives the level's radiative widths: write 0 in that row's width column "
                       "(column 12) or delete this record");
    }
  }

  GammaFinalState *Final(GammaDecaySpec &s, const Value *pair, const std::string &field, bool allowRepeat = false) {
    if (!pair) return nullptr;
    const std::string text = PairText(pair->pair);
    const std::vector<int> keys = PairsNamed(pair->pair.light, pair->pair.heavy, pair->pair.state, tables_.pairs);
    if (keys.empty()) {
      Report(codes::kUnknownPair, field, UnknownPairMessage("the final state " + text));
      return nullptr;
    }
    if (keys.size() > 1) {
      Report(codes::kAmbiguousPair, field, AmbiguousPairMessage("the final state " + text, keys, "the label does not tell them apart"));
      return nullptr;
    }
    for (GammaFinalState &f : s.finals)
      if (f.pairKey == keys[0]) {
        if (!allowRepeat && f.role != GammaFinalState::Role::NotListed) {
          Report(codes::kGammaRecord, field, text + " is named twice (as a branching or the remainder)");
          return nullptr;
        }
        f.label = text;
        return &f;
      }
    Report(codes::kNotAChannel, field, text + " is no photon channel of level " + std::to_string(s.levelId));
    return nullptr;
  }

  /// The constraint sum b_f = 1 (gamma_parametrisation.md section 2.6).
  void Branchings(GammaDecaySpec &s) {
    if (!ok_) return;
    bool anyFit = false;
    double sum = 0.0;
    for (std::size_t f = 0; f < s.finals.size(); f++)
      if ((int)f != s.remainder) {
        sum += s.finals[f].branching.value;
        anyFit = anyFit || s.finals[f].branching.fit;
      }
    if (s.finals.size() == 1) {
      // One final state: b = 1, whatever is written.
      if (s.finals[0].role == GammaFinalState::Role::Listed && std::fabs(s.finals[0].branching.value - 1.0) > kBranchingSumTolerance)
        Report(codes::kBranchingSum, "branching", "the level has one final state, " + s.finals[0].label + ": its branching is 1");
      s.finals[0].role = GammaFinalState::Role::Remainder;
      s.remainder = 0;
      return;
    }
    if (s.remainder < 0) {
      if (anyFit) {
        Report(codes::kBranchingSum, "branching",
               "a branching is fitted, so one final state must take the rest: add 'remainder <final state>'");
        return;
      }
      if (std::fabs(sum - 1.0) > kBranchingSumTolerance) {
        Report(codes::kBranchingSum, "branching", "the fixed branchings sum to " + Fmt(sum) + ", not 1");
        return;
      }
      // All fixed and summing to 1: the last one listed is the remainder.
      for (std::size_t f = 0; f < s.finals.size(); f++)
        if (s.finals[f].role == GammaFinalState::Role::Listed) s.remainder = (int)f;
      if (s.remainder >= 0) s.finals[s.remainder].role = GammaFinalState::Role::Remainder;
      return;
    }
    if (sum > 1.0 + kBranchingSumTolerance)
      Report(codes::kBranchingRemainder, "remainder",
             "the other branchings sum to " + Fmt(sum) + ", so the remainder to " + s.finals[s.remainder].label +
                 " would be " + Fmt(1.0 - sum) + " < 0");
  }

  const LevelTables &tables_;
  GammaDecayBuild &out_;
  const ValidatedRecord *record_ = nullptr;
  bool ok_ = true;
};

}  // namespace

GammaDecayBuild BuildGammaDecays(const ValidationResult &validated, const LevelTables &tables) {
  GammaDecayBuild out;
  Builder b(tables, out);
  for (const ValidatedRecord &r : validated.records)
    if (r.kind == "gamma_decay") b.Add(r);
  return out;
}

}  // namespace experiment
}  // namespace azure2
