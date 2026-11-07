#include "azure2/physics/System.h"

#include <stdexcept>
#include <string>
#include <utility>

namespace azure2 {
namespace physics {

namespace {

// The partition kind a channel kind belongs to: particle channels to particle partitions, E/M multipoles to photon
// partitions, Fermi and Gamow-Teller transitions to beta partitions.
PartitionKind PartitionKindOf(ChannelKind kind) {
  switch (kind) {
    case ChannelKind::Particle: return PartitionKind::Particle;
    case ChannelKind::Electric:
    case ChannelKind::Magnetic: return PartitionKind::Photon;
    case ChannelKind::Fermi:
    case ChannelKind::GamowTeller: return PartitionKind::Beta;
  }
  throw std::invalid_argument("System: unknown channel kind");
}

}  // namespace

namespace {
Convention WithBasis(LevelBasis basis) {
  Convention c;
  c.basis = basis;
  return c;
}
}  // namespace

System::System(LevelBasis basis, std::vector<Nucleus> nuclei, std::vector<Partition> partitions,
               std::vector<Channel> channels, std::vector<JpiBlock> blocks) :
  System(WithBasis(basis), std::move(nuclei), std::move(partitions), std::move(channels), std::move(blocks)) {}

System::System(Convention convention, std::vector<Nucleus> nuclei, std::vector<Partition> partitions,
               std::vector<Channel> channels, std::vector<JpiBlock> blocks) :
  convention_(convention),
  nuclei_(std::move(nuclei)),
  partitions_(std::move(partitions)),
  channels_(std::move(channels)),
  blocks_(std::move(blocks)) {
  if ((convention_.electric_sign != 1 && convention_.electric_sign != -1) ||
      (convention_.magnetic_sign != 1 && convention_.magnetic_sign != -1))
    throw std::invalid_argument("System: the convention's photon-channel signs must be +1 or -1");
  const int nNuclei = (int)nuclei_.size();
  const int nPartitions = (int)partitions_.size();
  const int nChannels = (int)channels_.size();

  for (std::size_t a = 0; a < partitions_.size(); a++) {
    const Partition &p = partitions_[a];
    const std::string what = "System: partition " + std::to_string(a) + ": ";
    if (p.nucleus < 0 || p.nucleus >= nNuclei) throw std::invalid_argument(what + "nucleus index out of range");
    for (int i = 0; i < 2; i++) {
      if (p.two_j[i] < 0) throw std::invalid_argument(what + "negative spin");
      if (p.parity[i] != 1 && p.parity[i] != -1) throw std::invalid_argument(what + "parity must be +1 or -1");
    }
  }

  for (std::size_t c = 0; c < channels_.size(); c++) {
    const Channel &ch = channels_[c];
    const std::string what = "System: channel " + std::to_string(c) + ": ";
    if (ch.partition < 0 || ch.partition >= nPartitions) throw std::invalid_argument(what + "partition index out of range");
    if (PartitionKindOf(ch.kind) != partitions_[ch.partition].kind)
      throw std::invalid_argument(what + "channel kind does not match the kind of its partition");
    if (ch.l < 0 || ch.two_s < 0) throw std::invalid_argument(what + "negative l or s");
    for (std::size_t c2 = 0; c2 < c; c2++) {
      const Channel &o = channels_[c2];
      if (o.partition == ch.partition && o.kind == ch.kind && o.l == ch.l && o.two_s == ch.two_s)
        throw std::invalid_argument(what + "same (partition, kind, l, s) as channel " + std::to_string(c2));
    }
  }

  for (std::size_t b = 0; b < blocks_.size(); b++) {
    const JpiBlock &blk = blocks_[b];
    const std::string what = "System: block " + std::to_string(b) + ": ";
    if (blk.nucleus() >= nNuclei) throw std::invalid_argument(what + "nucleus index out of range");
    for (std::size_t k = 0; k < blk.num_channels(); k++) {
      const int c = blk.channel(k);
      if (c >= nChannels) throw std::invalid_argument(what + "channel index out of range");
      // A block is one compound nucleus's level matrix: its channels lead to partitions of that nucleus.
      if (partitions_[channels_[c].partition].nucleus != blk.nucleus())
        throw std::invalid_argument(what + "channel " + std::to_string(c) + " belongs to another compound nucleus");
    }
    for (std::size_t b2 = 0; b2 < b; b2++) {
      const JpiBlock &o = blocks_[b2];
      if (o.two_J() == blk.two_J() && o.parity() == blk.parity() && o.nucleus() == blk.nucleus())
        throw std::invalid_argument(what + "same J^pi and nucleus as block " + std::to_string(b2));
    }
  }
}

std::size_t System::num_levels() const {
  std::size_t n = 0;
  for (const JpiBlock &b : blocks_) n += b.num_levels();
  return n;
}

int System::find_channel(int partition, ChannelKind kind, int l, int two_s) const {
  for (std::size_t c = 0; c < channels_.size(); c++) {
    const Channel &ch = channels_[c];
    if (ch.partition == partition && ch.kind == kind && ch.l == l && ch.two_s == two_s) return (int)c;
  }
  return -1;
}

}  // namespace physics
}  // namespace azure2
