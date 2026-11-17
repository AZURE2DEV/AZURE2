// The checks of the final states of external capture and of capture into a final state above a threshold
// (BuildContext.h; the rules in FinalStates.h), made by Schema::Build.
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/experiment/BuildContext.h"
#include "azure2/experiment/FinalStates.h"

namespace azure2 {
namespace experiment {
namespace build {

// The refusals of the set-up's external capture (CompoundSystem::ParseExternalCapture takes them from the same rule): the final
// state of a photon pair with external capture is one level of <levels>; a pair that an active data set or
// extrapolation exits into is refused without one (the set-up only warns about a pair nothing exits into); two levels
// within 1 keV, or a level that is an earlier pair's final state, are refused in any case.  External capture is taken
// as on (the API sets it up always; the command line when a data set exits into such a pair).
void CheckFinalStates(Context &c, const LevelTables &levels) {
  if (levels.externalCapturePairs.empty()) return;
  // the record of the first active data set or extrapolation exiting into a pair (the diagnostic's place), by key
  std::map<int, const ValidatedRecord *> user;
  for (std::size_t i = 0; i < c.data_.size(); i++)
    if (c.data_[i].active) user.emplace(c.out_.data.datasets.items[i].spec.reaction.exitKey, c.data_[i].record);
  for (std::size_t i = 0; i < c.extrapolations_.size(); i++)
    if (c.extrapolations_[i].active)
      user.emplace(c.out_.extrapolations.extrapolations.items[i].spec.reaction.exitKey, c.extrapolations_[i].record);
  auto report = [&](const char *code, int key, int otherKey, const std::string &message) {
    auto u = user.find(key);
    if (u == user.end() && otherKey) u = user.find(otherKey);
    if (u != user.end()) {
      c.Error(code, *u->second, "reaction", message);
      return;
    }
    Diagnostic d;
    d.code = code;
    d.section = "levels";
    d.message = message;
    c.out_.diagnostics.push_back(d);
  };
  for (const FinalStateOfPair &f : ExternalCaptureFinalStates(levels, levels.externalCapturePairs)) {
    const PairInfo *p = PairByKey(f.pairKey, levels.pairs);
    std::ostringstream state;
    state << "J^pi = " << (f.twoJ % 2 ? std::to_string(f.twoJ) + "/2" : std::to_string(f.twoJ / 2)) << (f.parity > 0 ? "+" : "-")
          << " at E_x = " << f.energy << " MeV of " << NuclideLabel(p->heavyZ, p->heavyA);
    const std::string pair = "photon pair " + std::to_string(f.pairKey) + " (" + PairLabel(*p) + ")";
    if (f.levels.empty()) {
      if (!user.count(f.pairKey)) continue;  // the set-up warns and switches its external capture off
      report(codes::kFinalStateNotLevel, f.pairKey, 0,
             "the external capture into " + pair + " needs its final state, " + state.str() + ", as a level of <levels> (with " +
                 "its ANC and gamma widths); none is listed" + (f.groupExists ? " within 1 keV in that J^pi group" : ", and no level has that J^pi") +
                 ": add the level, or give the pair the energy of a level");
    } else if (f.levels.size() > 1) {
      report(codes::kFinalStateNotLevel, f.pairKey, 0,
             "the final state of " + pair + ", " + state.str() + ", is ambiguous: " + std::to_string(f.levels.size()) +
                 " levels of that J^pi lie within 1 keV of it; move the extra level(s) by more than 1 keV");
    } else if (f.sharedWithPairKey) {
      report(codes::kFinalStateShared, f.pairKey, f.sharedWithPairKey,
             "photon pairs " + std::to_string(f.sharedWithPairKey) + " and " + std::to_string(f.pairKey) +
                 " both have the final state " + state.str() + " (level " + std::to_string(levels.levels[f.levels[0]].level) +
                 " of that J^pi group); a photon pair has its own final state: remove the duplicated pair, or give it the "
                 "energy of another level");
    }
  }
}


// An external-capture final state at a particle threshold or too close above a charged one (FinalStates.h,
// UnboundFinalStates): refused in every reading (the run, --validate, pyazr), since its outgoing wave is not evaluated.
void CheckUnboundFinalStates(Context &c, const LevelTables &levels) {
  std::vector<CaptureRequest> requests;
  std::vector<const ValidatedRecord *> recordOf;
  auto ask = [&](const ReactionSpec &spec, const ValidatedRecord *record) {
    requests.push_back({spec.entranceKey, spec.exitKey, (int)recordOf.size()});
    for (const ComponentSpec &t : spec.components) requests.push_back({t.entranceKey, t.exitKey, (int)recordOf.size()});
    recordOf.push_back(record);
  };
  for (std::size_t i = 0; i < c.data_.size(); i++)
    if (c.data_[i].active) ask(c.out_.data.datasets.items[i].spec.reaction, c.data_[i].record);
  for (std::size_t i = 0; i < c.extrapolations_.size(); i++)
    if (c.extrapolations_[i].active) ask(c.out_.extrapolations.extrapolations.items[i].spec.reaction, c.extrapolations_[i].record);
  const std::vector<UnboundFinalState> found = UnboundFinalStates(levels, requests);
  if (found.empty()) return;
  auto energy = [](double e) {
    std::ostringstream o;
    o << e << " MeV";
    return o.str();
  };
  for (const UnboundFinalState &u : found) {
    const PairInfo &photon = *PairByKey(u.photonPairKey, levels.pairs);
    const PairInfo &particle = *PairByKey(u.particlePairKey, levels.pairs);
    const std::string state = "J^pi = " +
                              (photon.twoJHeavy % 2 ? std::to_string(photon.twoJHeavy) + "/2" : std::to_string(photon.twoJHeavy / 2)) +
                              (photon.heavyParity > 0 ? "+" : "-") + " of " + NuclideLabel(photon.heavyZ, photon.heavyA);
    const std::string pairs = "photon pair " + std::to_string(photon.key) + " (" + PairLabel(photon) + ")";
    const std::string threshold = "the threshold of pair " + std::to_string(particle.key) + " (" + PairLabel(particle) +
                                  "), E_x = " + energy(u.threshold);
    const LevelInfo &level = levels.levels[u.level];
    std::string ids;
    for (int id : level.ids) ids += (ids.empty() ? "" : ", ") + std::to_string(id);
    std::ostringstream how;
    if (u.finalEnergy > u.threshold)
      how << "lies only " << (u.finalEnergy - u.threshold) * 1e3 << " keV above " << threshold << " (eta = " << u.eta << " > "
          << kMaxOutgoingWaveEta << ")";
    else
      how << "lies at " << threshold;
    const std::string message =
        "the final state of " + pairs + ", " + state + " at E_x = " + energy(u.finalEnergy) + " (the level of <levels> with " +
        "level ID " + ids + "), " + how.str() + ", which has a channel in its J^pi group: an unbound final state is a narrow " +
        "resonance with an outgoing wave in that channel, which cannot be evaluated this close to the threshold. Lower the " +
        "final state below " + energy(u.threshold) + " (it is then bound), move it higher, or switch external capture into pair " +
        std::to_string(photon.key) + " off (column 31 = 0)";
    if (u.request >= 0) {  // on the first data set or extrapolation that captures into the pair
      c.Error(codes::kFinalStateUnbound, *recordOf[u.request], "reaction", message);
    } else {  // nothing asks it (a rate can): in <levels>
      Diagnostic d;
      d.code = codes::kFinalStateUnbound;
      d.section = "levels";
      d.message = message;
      c.out_.diagnostics.push_back(d);
    }
  }
}

}  // namespace build
}  // namespace experiment
}  // namespace azure2
