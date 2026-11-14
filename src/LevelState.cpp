#include "LevelState.h"
#include <atomic>
#include <limits>

namespace {

/// The version generator: names level states uniquely in the process (LevelState::Version).
std::uint64_t NewVersion() {
  static std::atomic<std::uint64_t> counter{0};
  return ++counter;
}

}  // namespace

void LevelState::SetGammaDecayFill(int r, GammaDecayFill fill) {
  if ((int)gammaDecayFills_.size() <= r) gammaDecayFills_.resize(r + 1);
  gammaDecayFills_[r] = std::move(fill);
}

LevelState::LevelState(const std::vector<int> &numLevels, const std::vector<int> &numChannels) {
  const double never = std::numeric_limits<double>::quiet_NaN();
  groups_.resize(numLevels.size());
  for (std::size_t j = 0; j < numLevels.size(); j++) {
    Group &g = groups_[j];
    const std::size_t n = numLevels[j], nc = numChannels[j];
    g.numChannels = (int)nc;
    g.E.assign(n, 0.0);
    g.inputE.assign(n, 0.0);
    g.sqrtNF.assign(n, 1.0);
    g.parkNorm.assign(n, 1.0);
    g.shiftEnergy.assign(n, never);
    g.d2Energy.assign(n, never);
    g.gamma.assign(n * nc, 0.0);
    g.ecGamma.assign(n * nc, 0.0);
    g.shift.assign(n * nc, 0.0);
    g.d1.assign(n * nc, 0.0);
    g.d2.assign(n * nc, 0.0);
    g.zero.assign(n * nc, 0);
    g.version.resize(n);
    for (std::uint64_t &v : g.version) v = NewVersion();
  }
}

void LevelState::Touch(int j, int la) {
  groups_[j - 1].version[la - 1] = NewVersion();
}

void LevelState::SetE(int j, int la, double energy) {
  Touch(j, la);
  groups_[j - 1].E[la - 1] = energy;
}

void LevelState::SetGamma(int j, int la, int ch, double gamma) {
  Touch(j, la);
  groups_[j - 1].gamma[At(j, la, ch)] = gamma;
}

void LevelState::SetShift(int j, int la, int ch, double value) {
  Touch(j, la);
  groups_[j - 1].shift[At(j, la, ch)] = value;
}

void LevelState::SetShiftSecondDerivative(int j, int la, int ch, double value) {
  Touch(j, la);
  groups_[j - 1].d2[At(j, la, ch)] = value;
}

void LevelState::SetParkNorm(int j, int la, double value) {
  Touch(j, la);
  groups_[j - 1].parkNorm[la - 1] = value;
}
