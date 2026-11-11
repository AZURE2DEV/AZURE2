// The final states of external capture (FinalStates.h).
#include "azure2/experiment/FinalStates.h"

#include <cmath>
#include <map>

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

}  // namespace experiment
}  // namespace azure2
