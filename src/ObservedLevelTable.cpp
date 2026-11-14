#include "ObservedLevelTable.h"

ObservedLevelTable::ObservedLevelTable(const std::vector<int> &numLevels, const std::vector<int> &numChannels) {
  levels_.resize(numLevels.size());
  for (std::size_t j = 0; j < numLevels.size(); j++) {
    ObservedLevel empty;
    const std::size_t nc = (std::size_t)numChannels.at(j);
    empty.gamma.assign(nc, 0.0);
    empty.external.assign(nc, complex(0.0, 0.0));
    empty.W.assign(nc, 0.0);
    empty.value.assign(nc, 0.0);
    empty.kind.assign(nc, 0);
    levels_[j].assign((std::size_t)numLevels[j], empty);
  }
}
