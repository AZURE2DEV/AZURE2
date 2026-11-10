// The `potential` records of hybrid-potential settings (PotentialRecords.h).
#include "azure2/convert/PotentialRecords.h"

#include <algorithm>
#include <cctype>

#include "OldProject.h"

namespace azure2 {
namespace convert {

namespace {

void AddField(io::Record &r, const std::string &key, const std::vector<std::string> &words) {
  io::Field f;
  f.key = key;
  for (const std::string &w : words) {
    io::Token t;
    t.text = w;
    f.tokens.push_back(t);
  }
  r.fields.push_back(f);
}

/// Another pair with the same particles and states: the label alone does not name `p` (section 3.5).
bool Ambiguous(const experiment::PairInfo &p, const std::vector<experiment::PairInfo> &pairs) {
  for (const experiment::PairInfo &o : pairs)
    if (o.key != p.key && o.kind == p.kind && o.lightZ == p.lightZ && o.lightA == p.lightA && o.heavyZ == p.heavyZ &&
        o.heavyA == p.heavyA && o.lightExcitation == p.lightExcitation && o.heavyExcitation == p.heavyExcitation)
      return true;
  return false;
}

/// "p+12C" -> "p12C", "a[3.0]+12C[4.439]" -> "a_3_0_12C_4_439": a format-3 name ([A-Za-z][A-Za-z0-9_-]*).
std::string NameOf(const std::string &label) {
  std::string s;
  for (char c : label) {
    if (c == '+' || c == ']') continue;
    s += std::isalnum((unsigned char)c) ? c : '_';
  }
  if (s.empty() || !std::isalpha((unsigned char)s[0])) s = "p" + s;
  return s;
}

}  // namespace

PotentialRecordsResult PotentialRecords(const experiment::HybridPotentialSettings &settings,
                                        const std::vector<experiment::PairInfo> &pairs, const std::set<std::string> &taken,
                                        const std::map<int, std::string> &namesOfPairs) {
  PotentialRecordsResult out;
  std::set<std::string> names = taken;
  std::set<int> keys;
  for (const experiment::PairInfo &p : pairs) keys.insert(p.key);
  for (const auto &own : settings.pairs) {
    if (!keys.count(own.first))
      out.notes.push_back("pair=" + std::to_string(own.first) + " names no pair of <levels>; its section is not written");
    else if (!own.second.enabled)
      out.notes.push_back("pair " + std::to_string(own.first) + ": useHybridPotential=0, no potential record");
  }
  if (!settings.defaults.enabled) out.notes.push_back("the default section has useHybridPotential=0: no record from it");
  bool defaultUsed = false;
  std::vector<const experiment::PairInfo *> ordered;
  for (const experiment::PairInfo &p : pairs) ordered.push_back(&p);
  std::sort(ordered.begin(), ordered.end(), [](const experiment::PairInfo *a, const experiment::PairInfo *b) { return a->key < b->key; });
  for (const experiment::PairInfo *p : ordered) {
    const bool own = settings.pairs.count(p->key) > 0;
    const experiment::HybridPotentialSetting &s = settings.For(p->key);
    if (!s.enabled) continue;
    const std::string label = experiment::PairLabel(*p);
    if (p->kind != 0) {
      out.notes.push_back(std::string(own ? "the section of" : "the default section applies to") + " pair " +
                          std::to_string(p->key) + " (" + label + "), a " + (p->kind == 10 ? "photon" : "beta") +
                          " pair: no record (the engine computes no channel function of a potential for it)");
      continue;
    }
    if (!own) defaultUsed = true;
    io::Record r;
    r.kind = "potential";
    const auto own_name = namesOfPairs.find(p->key);
    const std::string base = own_name != namesOfPairs.end() ? own_name->second : NameOf(label);
    r.name = base;
    for (int i = 2; names.count(r.name); i++) r.name = base + "_" + std::to_string(i);
    names.insert(r.name);
    AddField(r, "pair", {label});
    if (Ambiguous(*p, pairs)) AddField(r, "pair_key", {std::to_string(p->key)});
    AddField(r, "coulomb_radius", {ShortestText(s.coulomb_radius_fm), "fm"});
    if (s.shape == physics::FormShape::Gaussian)
      AddField(r, "central", {"gaussian", "depth", ShortestText(s.V0), "MeV", "width", ShortestText(s.r0), "fm"});
    else
      AddField(r, "central", {"woods_saxon", "depth", ShortestText(s.V0), "MeV", "radius", ShortestText(s.R), "fm",
                              "diffuseness", ShortestText(s.a), "fm"});
    AddField(r, "matching_radius", {ShortestText(s.matching_radius_fm), "fm"});
    AddField(r, "step", {ShortestText(s.step_fm), "fm"});
    out.records.records.push_back(r);
  }
  if (settings.defaults.enabled && !defaultUsed)
    out.notes.push_back("the default section applies to no particle pair: no record from it");
  return out;
}

}  // namespace convert
}  // namespace azure2
