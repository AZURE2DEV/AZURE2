// The final states of external capture (FinalStates.h).
#include "azure2/experiment/FinalStates.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "azure2/physics/ChannelFunctions.h"

namespace azure2 {
namespace experiment {

std::vector<FinalStateOfPair> ExternalCaptureFinalStates(const LevelTables &tables, const std::vector<int> &pairKeys) {
  std::vector<FinalStateOfPair> out;
  std::map<int, int> ownerOfLevel;  // index into tables.levels -> the pair key that has it as its final state
  for (int key : pairKeys) {
    const PairInfo *pair = PairByKey(key, tables.pairs);
    if (!pair || pair->kind != 10) continue;
    FinalStateOfPair f;
    f.pairKey = key;
    f.twoJ = pair->twoJHeavy;
    f.parity = pair->heavyParity;
    f.energy = pair->heavyExcitation;
    for (std::size_t i = 0; i < tables.levels.size(); i++) {
      const LevelInfo &level = tables.levels[i];
      // the J^pi group of the compound nucleus (a photon pair's heavy particle is the compound nucleus)
      if (level.twoJ != f.twoJ || level.parity != f.parity || level.nucleusZ != pair->heavyZ || level.nucleusA != pair->heavyA)
        continue;
      f.groupExists = true;
      if (std::fabs(level.energy - f.energy) <= 1.0e-3) f.levels.push_back((int)i);
    }
    if (f.levels.size() == 1) {
      const auto owner = ownerOfLevel.emplace(f.levels[0], key);
      if (!owner.second) f.sharedWithPairKey = owner.first->second;
    }
    out.push_back(f);
  }
  return out;
}

const FinalStateOfPair *FinalStateOf(const std::vector<FinalStateOfPair> &states, int key) {
  for (const FinalStateOfPair &f : states)
    if (f.pairKey == key) return &f;
  return nullptr;
}

std::vector<UnboundFinalState> UnboundFinalStates(const LevelTables &tables, const std::vector<CaptureRequest> &requests) {
  std::vector<UnboundFinalState> out;
  const std::vector<FinalStateOfPair> finals = ExternalCaptureFinalStates(tables, tables.externalCapturePairs);
  // E_x of the final state of photon pair `photon`: E_gamma = 0 there (the photon channel's E_x - S - e2), and the
  // energy external capture takes for its level (none, or an ambiguous one: E_FINAL_STATE_NOT_LEVEL)
  auto finalEnergyOf = [&](const PairInfo &photon) {
    const FinalStateOfPair *f = FinalStateOf(finals, photon.key);
    if (!f || f->levels.size() != 1) return photon.threshold;
    return std::max(photon.threshold, tables.levels[f->levels[0]].energy);
  };
  // whether request r captures into photon pair `photon`: its exit, or every photon pair of its nucleus (total capture)
  auto asks = [&](const CaptureRequest &r, const PairInfo &photon) {
    const PairInfo *in = PairByKey(r.entranceKey, tables.pairs);
    if (!in || in->kind != 0 || photon.kind != 10) return false;
    if (r.exitKey >= 0) return photon.key == r.exitKey;
    return photon.heavyZ == in->lightZ + in->heavyZ && photon.heavyA == in->lightA + in->heavyA;
  };
  for (const FinalStateOfPair &f : finals) {
    // no final state, or an ambiguous one: E_FINAL_STATE_NOT_LEVEL, or its external capture switched off
    if (f.levels.size() != 1) continue;
    const PairInfo &photon = *PairByKey(f.pairKey, tables.pairs);
    const double finalEnergy = finalEnergyOf(photon);
    UnboundFinalState u;
    u.photonPairKey = f.pairKey;
    for (const LevelChannelInfo &ch : tables.levels[f.levels[0]].channels) {
      const PairInfo *p = PairByKey(ch.pairKey, tables.pairs);
      if (!p || p->kind != 0 || finalEnergy < p->threshold) continue;
      const double above = finalEnergy - p->threshold;
      double eta = 0.0;
      if (above > 0.0) {
        physics::ChannelKinematics kin;
        kin.z1 = p->lightZ;
        kin.z2 = p->heavyZ;
        kin.mu_u = p->lightMass * p->heavyMass / (p->lightMass + p->heavyMass);
        eta = physics::eta(kin, above);
        if (!(eta > kMaxOutgoingWaveEta)) continue;  // an unbound final state the engine computes
      }
      // at the threshold, or too close above it: the pair nearest its threshold
      if (u.particlePairKey == 0 || above < u.finalEnergy - u.threshold) {
        u.particlePairKey = p->key;
        u.threshold = p->threshold;
        u.finalEnergy = finalEnergy;
        u.eta = eta;
      }
    }
    if (u.particlePairKey == 0) continue;
    u.finalEnergy = finalEnergy;
    u.level = f.levels[0];
    for (const CaptureRequest &r : requests)
      if (asks(r, photon)) {
        u.request = r.request;
        break;
      }
    out.push_back(u);
  }
  return out;
}

}  // namespace experiment
}  // namespace azure2
