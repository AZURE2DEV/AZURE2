// The level of a <levels> line (LevelGrouping.h): the rule of CompoundSystem::Fill since 2026-10-24, moved here unchanged so
// that pyazr reads levels as the engine does.
#include "azure2/io/LevelGrouping.h"

#include <cmath>

namespace azure2 {
namespace io {

LevelGrouper::Slot LevelGrouper::Assign(int group, double energy, int id) {
  std::vector<Level> &levels = groups_[group];
  Slot slot;
  int sameEnergy = 0, nearby = 0;
  for (std::size_t lv = 1; lv <= levels.size() && !slot.level; lv++) {
    const double dE = std::fabs(levels[lv - 1].energy - energy);
    if (dE > 1e-3) continue;
    if (levels[lv - 1].ids.count(id)) slot.level = (int)lv;
    else if (dE == 0. && !sameEnergy) sameEnergy = (int)lv;
    else if (!nearby) nearby = (int)lv;
  }
  if (slot.level) {
    slot.kind = Kind::Joined;
  } else if (sameEnergy) {
    slot.level = sameEnergy;
    slot.kind = Kind::SameEnergy;
    slot.other = sameEnergy;
    slot.otherID = *levels[sameEnergy - 1].ids.begin();
  } else {
    levels.push_back(Level{energy, {}});
    slot.level = (int)levels.size();
    slot.kind = nearby ? Kind::NewNearby : Kind::NewLevel;
    if (nearby) {
      slot.other = nearby;
      slot.otherID = *levels[nearby - 1].ids.begin();
    }
  }
  levels[slot.level - 1].ids.insert(id);
  return slot;
}

}  // namespace io
}  // namespace azure2
