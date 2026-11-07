#include "azure2/physics/JpiBlock.h"

#include <stdexcept>
#include <string>
#include <utility>

namespace azure2 {
namespace physics {

JpiBlock::JpiBlock(int two_J, int parity, int nucleus, std::vector<int> channels, std::vector<double> boundary,
                   std::vector<double> energies_MeV, std::vector<double> gamma) :
  two_J_(two_J),
  parity_(parity),
  nucleus_(nucleus),
  channels_(std::move(channels)),
  B_(std::move(boundary)),
  E_(std::move(energies_MeV)),
  gamma_(std::move(gamma)) {
  const std::string what = "JpiBlock (2J = " + std::to_string(two_J) + "): ";
  if (two_J < 0) throw std::invalid_argument(what + "negative 2J");
  if (parity != 1 && parity != -1) throw std::invalid_argument(what + "parity must be +1 or -1, got " + std::to_string(parity));
  if (nucleus < 0) throw std::invalid_argument(what + "negative nucleus index");
  if (B_.size() != channels_.size())
    throw std::invalid_argument(what + std::to_string(channels_.size()) + " channels but " + std::to_string(B_.size()) + " boundary values");
  if (gamma_.size() != E_.size() * channels_.size())
    throw std::invalid_argument(what + std::to_string(E_.size()) + " levels x " + std::to_string(channels_.size()) +
                                " channels need as many amplitudes, got " + std::to_string(gamma_.size()));
  for (std::size_t k = 0; k < channels_.size(); k++) {
    if (channels_[k] < 0) throw std::invalid_argument(what + "negative channel index");
    for (std::size_t k2 = 0; k2 < k; k2++)
      if (channels_[k2] == channels_[k]) throw std::invalid_argument(what + "channel " + std::to_string(channels_[k]) + " listed twice");
  }
}

int JpiBlock::find_channel(int c) const {
  for (std::size_t k = 0; k < channels_.size(); k++)
    if (channels_[k] == c) return (int)k;
  return -1;
}

}  // namespace physics
}  // namespace azure2
