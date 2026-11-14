#include "LevelInputs.h"

LevelInputs::LevelInputs(const std::vector<int> &numLevels, const std::vector<int> &numChannels) {
  groups_.resize(numLevels.size());
  for (std::size_t j = 0; j < numLevels.size(); j++) {
    Group &g = groups_[j];
    const std::size_t n = numLevels[j], nc = numChannels[j];
    g.numChannels = (int)nc;
    g.E.assign(n, 0.0);
    g.gamma.assign(n * nc, 0.0);
  }
}
