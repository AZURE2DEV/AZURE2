#include "LevelParameterLayout.h"

LevelParameterLayout::LevelParameterLayout(const std::vector<int> &numLevels, const std::vector<int> &numChannels) {
  int levelNumber = 1;
  for (std::size_t g = 0; g < numLevels.size(); g++) {
    const int j = (int)g + 1;
    first_.push_back((int)levels_.size());
    for (int la = 1; la <= numLevels[g]; la++, levelNumber++) {
      Level level;
      level.group = j;
      level.level = la;
      level.levelNumber = levelNumber;
      level.energy = (int)slots_.size();
      slots_.push_back({Slot::Kind::kEnergy, j, la, 0, levelNumber});
      for (int ch = 1; ch <= numChannels[g]; ch++) {
        level.amplitude.push_back((int)slots_.size());
        slots_.push_back({Slot::Kind::kAmplitude, j, la, ch, levelNumber});
      }
      levels_.push_back(std::move(level));
    }
  }
}

std::string LevelParameterLayout::Name(const Slot &slot) {
  if (slot.kind == Slot::Kind::kEnergy) return "energy_" + std::to_string(slot.levelNumber);
  return "width_" + std::to_string(slot.levelNumber) + "_" + std::to_string(slot.channel);
}
