// The `channel_radius` records of <channelRadii> (ChannelRadii.h; BuildContext.h, BuildChannelRadii): each names a
// particle pair of <levels> by its label, an l and optionally a channel spin s and a total J, and gives those
// channels a radius of their own.  docs/upgrade/dev_features.md section 2.
#include "azure2/experiment/ChannelRadii.h"

#include <map>

#include "azure2/experiment/BuildContext.h"

namespace azure2 {
namespace experiment {

namespace {
std::string HalfInteger(int twice) {
  return twice % 2 == 0 ? std::to_string(twice / 2) : std::to_string(twice) + "/2";
}
}  // namespace

std::string ChannelRadiusEntry::Text() const {
  std::string s = "pair " + std::to_string(pair_key) + ", l = " + std::to_string(l);
  if (has_s) s += ", s = " + HalfInteger(two_s);
  if (has_J) s += ", J = " + HalfInteger(two_J);
  return s;
}

const ChannelRadiusEntry *ChannelRadiusFor(const std::vector<ChannelRadiusEntry> &entries, int key, int l, int twoS,
                                           int twoJ, bool *ambiguous) {
  const ChannelRadiusEntry *best = nullptr;
  bool tie = false;
  for (const ChannelRadiusEntry &e : entries) {
    if (!e.Matches(key, l, twoS, twoJ)) continue;
    if (!best || e.Specificity() > best->Specificity()) {
      best = &e;
      tie = false;
    } else if (e.Specificity() == best->Specificity()) {
      // Two matching entries of one specificity are an s-only and a J-only entry (a repeated selection is refused
      // when the records are built); a later, more specific entry clears the tie.
      tie = true;
    }
  }
  if (ambiguous) *ambiguous = tie;
  return best;
}

namespace build {

namespace {
std::string PairText(const ParticlePair &p) {
  std::string s = p.light.text + (p.light.excited ? "[" + p.light.excitation.text + "]" : "") + "+" + p.heavy.text;
  if (p.state.kind == ResidualState::Excited) s += "[" + p.state.energy.text + "]";
  return s;
}
}  // namespace

void BuildChannelRadii(Context &c) {
  for (const ValidatedRecord &r : c.v_.records) {
    if (r.kind != "channel_radius") continue;
    const Value *pair = r.Get("pair");
    const Value *l = r.Get("l");
    const Value *radius = r.Get("radius");
    if (!pair || !l || !radius) continue;  // reported by the validator
    const std::string text = PairText(pair->pair);
    int key = 0;
    if (!ResolvePair(c, r, "pair", "the pair " + text, pair->pair.light, pair->pair.heavy, pair->pair.state, text, "pair_key", key))
      continue;
    const PairInfo *info = c.Pair(key);
    if (info && info->kind != 0) {
      c.Error(codes::kChannelRadiusPairKind, r, "pair",
              "pair " + std::to_string(key) + " (" + PairLabel(*info) + ") is a " + (info->kind == 10 ? "photon" : "beta") +
                  " pair; a channel radius belongs to a channel of a particle pair");
      continue;
    }
    ChannelRadiusEntry e;
    e.pair_key = key;
    e.l = (int)l->integer;
    e.radius_fm = Context::InBase(radius->number, radius->unit);
    if (const Value *s = r.Get("s")) {
      e.has_s = true;
      e.two_s = s->spin.twice;
    }
    if (const Value *J = r.Get("J")) {
      e.has_J = true;
      e.two_J = J->spin.twice;
    }
    e.record = r.name;
    e.line = r.line;
    const ChannelRadiusEntry *first = nullptr;
    for (const ChannelRadiusEntry &o : c.out_.channelRadii)
      if (o.SameSelection(e)) first = &o;
    if (first) {
      c.Error(codes::kChannelRadiusTwice, r, "l",
              e.Text() + " has a channel radius already (channel_radius " + first->record + ", line " +
                  std::to_string(first->line) + "); one record per selection of channels");
      continue;
    }
    c.out_.channelRadii.push_back(e);
  }

  // The THM vertex pieces are built per l of the entrance pair, at the pair's radius: a THM data set or extrapolation
  // whose entrance pair has a channel with a radius of its own is refused until the vertex carries the channel's.
  for (const std::vector<Member> *members : {&c.data_, &c.extrapolations_})
    for (const Member &m : *members) {
      if (!m.record->Has("thm")) continue;
      for (const ChannelRadiusEntry &e : c.out_.channelRadii) {
        if (e.pair_key != m.entranceKey) continue;
        const PairInfo *info = c.Pair(m.entranceKey);
        c.Error(codes::kChannelRadiusThm, *m.record, "thm",
                "the entrance pair " + std::to_string(m.entranceKey) + " (" + (info ? PairLabel(*info) : std::string("?")) +
                    ") has a channel with a radius of its own (channel_radius " + e.record + ", line " +
                    std::to_string(e.line) + "); the THM vertex is built at the pair's radius: remove the record, or "
                    "the THM data set");
        break;
      }
    }
}

}  // namespace build
}  // namespace experiment
}  // namespace azure2
