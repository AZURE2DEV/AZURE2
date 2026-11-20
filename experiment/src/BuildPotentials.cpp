// The `potential` records of <potentials> (BuildContext.h, BuildPotentials): docs/upgrade/project_format3.md section
// 3.11 and docs/upgrade/spin_orbit_potential.md section 3.3.  Each record names a particle pair of <levels> by its
// label (the pair-key rule of the data sets, section 3.5) and becomes that pair's HybridPotentialSetting; numbers are
// converted to fm and MeV with the unit table (bit-exact for fm and MeV).  A pair without a record has no potential.
#include <algorithm>
#include <map>

#include "azure2/experiment/BuildContext.h"
#include "azure2/experiment/SpinOrbitChecks.h"

namespace azure2 {
namespace experiment {
namespace build {

namespace {

/// The pair as written: p+12C, a[3.0]+12C[4.439].
std::string PairText(const ParticlePair &p) {
  std::string s = p.light.text + (p.light.excited ? "[" + p.light.excitation.text + "]" : "") + "+" + p.heavy.text;
  if (p.state.kind == ResidualState::Excited) s += "[" + p.state.energy.text + "]";
  return s;
}

/// The value of a length or energy field or part in fm or MeV (the unit table's factor).
double InBase(const Value &v) { return Context::InBase(v.number, v.unit); }

/// The central term of the record into `s`.
void ReadCentral(const Value &central, HybridPotentialSetting &s) {
  const Value *depth = central.Get("depth");
  if (depth) s.V0 = InBase(*depth);
  if (central.word == "gaussian") {
    s.shape = physics::FormShape::Gaussian;
    if (const Value *w = central.Get("width")) s.r0 = InBase(*w);
    return;
  }
  s.shape = physics::FormShape::WoodsSaxon;
  if (const Value *r = central.Get("radius")) s.R = InBase(*r);
  if (const Value *a = central.Get("diffuseness")) s.a = InBase(*a);
}

/// The spin-orbit term of the record (spin_orbit woods_saxon depth radius diffuseness on light|heavy operator 2l.s).
physics::SpinOrbitTerm ReadSpinOrbit(const Value &so) {
  physics::SpinOrbitTerm t;
  if (const Value *d = so.Get("depth")) t.depth_MeV = InBase(*d);
  t.form.shape = physics::FormShape::WoodsSaxon;
  if (const Value *r = so.Get("radius")) t.form.radius_fm = InBase(*r);
  if (const Value *a = so.Get("diffuseness")) t.form.diffuseness_fm = InBase(*a);
  const Value *on = so.Get("on");
  t.acts_on = on && on->word == "heavy" ? physics::SpinOrbitParticle::Heavy : physics::SpinOrbitParticle::Light;
  return t;
}

/// "potential p12C: step not given, 0.0275 fm (the default)", the default as the schema row states it.
std::string DefaultNote(const ValidatedRecord &r, const std::string &field) {
  const RecordKind *kind = Schema::FindKind("potential");
  const FieldDef *f = kind ? kind->Find(field) : nullptr;
  return "potential " + r.name + ": " + field + " not given, " + (f ? f->defaultText : std::string("?")) + " (the default)";
}

}  // namespace

void BuildPotentials(Context &c) {
  std::map<int, const ValidatedRecord *> recordOfPair;
  std::map<std::pair<int, int>, const ValidatedRecord *> recordOfPairL;  // records with an l list
  std::map<int, const ValidatedRecord *> anyRecordOfPair;                // the first record of each pair, either kind
  for (const ValidatedRecord &r : c.v_.records) {
    if (r.kind != "potential") continue;
    const Value *pair = r.Get("pair");
    const Value *central = r.Get("central");
    const Value *coulomb = r.Get("coulomb_radius");
    if (!pair || !central || !coulomb) continue;  // reported by the validator
    const std::string text = PairText(pair->pair);
    int key = 0;
    if (!ResolvePair(c, r, "pair", "the pair " + text, pair->pair.light, pair->pair.heavy, pair->pair.state, text, "pair_key", key))
      continue;
    const PairInfo *info = c.Pair(key);
    if (info && info->kind != 0) {
      c.Error(codes::kPotentialPairKind, r, "pair",
              "pair " + std::to_string(key) + " (" + PairLabel(*info) + ") is a " + (info->kind == 10 ? "photon" : "beta") +
                  " pair; a hybrid potential bends the radial waves of two particles and belongs to a particle pair");
      continue;
    }
    if (const Value *depth = central->Get("depth"); depth && !(depth->number.value != 0.0)) {
      c.Error(codes::kOutOfRange, r, "central", "the depth must not be 0 (> 0 attractive, < 0 repulsive)");
      continue;
    }
    const Value *spinOrbit = r.Get("spin_orbit");
    // the l list: integers >= 0, none named by another record of the pair, and not with a spin-orbit term (whose keys
    // are the pair's)
    std::vector<int> ls;
    if (const Value *lv = r.Get("l")) {
      std::string problem;
      for (const Number &n : lv->numbers) {
        const int l = static_cast<int>(n.value);
        if (!(n.value >= 0.0) || n.value != l) problem = "'" + n.text + "' is not an integer >= 0";
        else if (std::find(ls.begin(), ls.end(), l) != ls.end()) problem = "l = " + n.text + " is named twice";
        else if (recordOfPairL.count({key, l}))
          problem = "l = " + n.text + " of pair " + std::to_string(key) + " (" + text + ") has a potential already (potential " +
                    recordOfPairL[{key, l}]->name + ", line " + std::to_string(recordOfPairL[{key, l}]->line) + ")";
        ls.push_back(l);
      }
      if (problem.empty() && spinOrbit) problem = "a record with an l list cannot have a spin-orbit term (give the term on the pair's record without l)";
      if (!problem.empty()) {
        c.Error(codes::kPotentialL, r, "l", problem);
        continue;
      }
    } else {
      const auto first = recordOfPair.find(key);
      if (first != recordOfPair.end()) {
        c.Error(codes::kPotentialTwice, r, "pair",
                "pair " + std::to_string(key) + " (" + text + ") has a potential already (potential " + first->second->name +
                    ", line " + std::to_string(first->second->line) + "); one record per pair, or records with different l");
        continue;
      }
      recordOfPair[key] = &r;
    }
    for (int l : ls) recordOfPairL[{key, l}] = &r;
    anyRecordOfPair.emplace(key, &r);
    if (spinOrbit && c.opt_.relativistic) {
      c.Error(codes::kSoRelativistic, r, "spin_orbit", "a spin-orbit term cannot be combined with relativistic=1 (<config>)");
      continue;
    }
    if (c.opt_.relativistic) {
      // The potential acts in a Schroedinger equation between a_c and the matching radius: no relativistic form.
      c.Error(codes::kPotentialRelativistic, r, "", "a hybrid potential cannot be combined with relativistic=1 (<config>)");
      continue;
    }
    HybridPotentialSetting s;
    s.enabled = true;
    ReadCentral(*central, s);
    if (spinOrbit) {
      s.spin_orbit = ReadSpinOrbit(*spinOrbit);
      const SpinOrbitRefusal refusal = info ? SpinOrbitPairRefusal(info->twoJLight, info->twoJHeavy, s.spin_orbit->acts_on,
                                                                   std::to_string(key) + " (" + PairLabel(*info) + ")")
                                            : SpinOrbitRefusal();
      if (refusal) {
        c.Error(refusal.code.c_str(), r, "spin_orbit", refusal.message);
        continue;
      }
    }
    if (const Value *b = r.Get("interior_boundary")) {
      if (info && info->twoJLight != 0 && info->twoJHeavy != 0) {
        c.Error(codes::kPotInteriorTwoSpins, r, "interior_boundary",
                "pair " + std::to_string(key) + " (" + PairLabel(*info) + ") has two particles with spin: the interior R-function of "
                    "its potential would couple the channel spins of one (l, J), which is not implemented "
                    "(docs/upgrade/interior_potential.md section 7)");
        continue;
      }
      s.interior_boundary = b->number.value;
    }
    s.coulomb_radius_fm = InBase(*coulomb);
    // r_m and the step: the table's defaults (HybridPotentialSetting's, the solver's) when absent, written to the log
    if (const Value *m = r.Get("matching_radius"))
      s.matching_radius_fm = InBase(*m);
    else
      c.out_.potentialNotes.push_back(DefaultNote(r, "matching_radius"));
    if (const Value *h = r.Get("step"))
      s.step_fm = InBase(*h);
    else
      c.out_.potentialNotes.push_back(DefaultNote(r, "step"));
    if (ls.empty()) {
      c.out_.potentials.pairs[key] = s;
      c.out_.potentialRecords[key] = {r.name, r.line};
    }
    for (int l : ls) {
      c.out_.potentials.byL[{key, l}] = s;
      c.out_.potentialRecordsByL[{key, l}] = {r.name, r.line};
    }
  }

  // The THM vertex uses point-Coulomb waves: a THM data set or extrapolation whose entrance pair has a potential is
  // refused (the engine's ThmVertexAllowsPair makes the same check on an old project).
  for (const std::vector<Member> *members : {&c.data_, &c.extrapolations_})
    for (const Member &m : *members) {
      if (!m.record->Has("thm")) continue;
      const auto p = anyRecordOfPair.find(m.entranceKey);
      if (p == anyRecordOfPair.end()) continue;
      const PairInfo *info = c.Pair(m.entranceKey);
      const bool spinOrbit = p->second->Has("spin_orbit");
      c.Error(spinOrbit ? codes::kSoThm : codes::kThmHybridPotential, *m.record, "thm",
              "the entrance pair " + std::to_string(m.entranceKey) + " (" + (info ? PairLabel(*info) : std::string("?")) +
                  ") has a hybrid potential" + (spinOrbit ? " with a spin-orbit term" : "") + " (potential " +
                  p->second->name + "); the THM vertex uses point-Coulomb waves and is not built with a potential: "
                  "remove the potential record of the pair, or the THM data set");
    }
}

}  // namespace build
}  // namespace experiment
}  // namespace azure2
