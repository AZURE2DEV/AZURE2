#include "CNucSystem.h"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "AZUREParams.h"
#include "CNuc.h"
#include "Config.h"

namespace phys = azure2::physics;

namespace {

// 2j as an integer.  The engine stores spins as doubles read from the file; a value that is not a multiple of 1/2
// would be stored rounded, so it is refused instead (no silent fallback).
int TwiceSpin(double j, const std::string &what) {
  const double twice = 2.0 * j;
  const double rounded = std::floor(twice + 0.5);
  if (twice != rounded || rounded < 0.0) {
    std::ostringstream msg;
    msg << "BuildSystemFromCNuc: " << what << " = " << j << " is not a non-negative multiple of 1/2";
    throw std::runtime_error(msg.str());
  }
  return (int)rounded;
}

phys::PartitionKind PartitionKindOf(int pType, int pairKey) {
  switch (pType) {
    case 0: return phys::PartitionKind::Particle;
    case 10: return phys::PartitionKind::Photon;
    case 20: return phys::PartitionKind::Beta;
  }
  throw std::runtime_error("BuildSystemFromCNuc: pair " + std::to_string(pairKey) + " has pair type " +
                           std::to_string(pType) + ", which the System has no partition kind for");
}

phys::ChannelKind ChannelKindOf(char radType) {
  switch (radType) {
    case 'P': return phys::ChannelKind::Particle;
    case 'E': return phys::ChannelKind::Electric;
    case 'M': return phys::ChannelKind::Magnetic;
    case 'F': return phys::ChannelKind::Fermi;
    case 'G': return phys::ChannelKind::GamowTeller;
  }
  throw std::runtime_error(std::string("BuildSystemFromCNuc: radiation type '") + radType +
                           "', which the System has no channel kind for");
}

phys::LevelBasis BasisOf(const Config &configure) {
  switch (AZUREParams::BasisForConfig(configure)) {
    case AZUREParams::kBasisStandard: return phys::LevelBasis::Formal;
    case AZUREParams::kBasisBrune: return phys::LevelBasis::Brune;
    case AZUREParams::kBasisPark: return phys::LevelBasis::Park;
  }
  throw std::runtime_error("BuildSystemFromCNuc: the configuration names no level-parameter basis");
}

}  // namespace

CNucSystem BuildSystemFromCNuc(CNuc &compound, const Config &configure, CNucSystemContent content) {
  const bool withParameters = (content == CNucSystemContent::Parameters);
  const bool allGroups = (content == CNucSystemContent::AllGroups);
  // Before TransformIn the level gammas are the file's partial widths (eV) and ANCs, not amplitudes.
  if (withParameters && (configure.paramMask & Config::TRANSFORM_PARAMETERS) && !compound.IsTransformedIn())
    throw std::runtime_error("BuildSystemFromCNuc: the compound nucleus has not been initialized "
                             "(CNuc::Initialize): its level parameters are still the input widths");

  CNucSystemMap map;

  std::vector<phys::Nucleus> nuclei;
  for (int i = 1; i <= compound.NumNuclei(); i++) {
    const CompoundNucleus *n = compound.GetNucleus(i);
    nuclei.push_back(phys::Nucleus{n->Z, n->A});
  }

  std::vector<phys::Partition> partitions;
  for (int a = 1; a <= compound.NumPairs(); a++) {
    const PPair *pair = compound.GetPair(a);
    const std::string what = "pair " + std::to_string(pair->GetPairKey());
    phys::Partition p;
    p.kind = PartitionKindOf(pair->GetPType(), pair->GetPairKey());
    p.nucleus = pair->GetNucleus() - 1;
    for (int i = 0; i < 2; i++) {
      p.Z[i] = pair->GetZ(i + 1);
      p.mass_u[i] = pair->GetM(i + 1);
      p.two_j[i] = TwiceSpin(pair->GetJ(i + 1), what + " particle " + std::to_string(i + 1) + " spin");
      p.parity[i] = pair->GetPi(i + 1);
      p.g_factor[i] = pair->GetG(i + 1);
    }
    p.excitation_MeV[0] = pair->GetLightExE();
    p.excitation_MeV[1] = pair->GetHeavyExE();
    p.separation_MeV = pair->GetSepE();
    partitions.push_back(p);
    map.pairOfPartition.push_back(a);
  }

  std::vector<phys::Channel> channels;
  std::vector<phys::JpiBlock> blocks;
  for (int j = 1; j <= compound.NumJGroups(); j++) {
    JGroup *group = compound.GetJGroup(j);
    if (!group->IsInRMatrix() && !allGroups) continue;  // a group made only to hold an external-capture final state
    const std::string what = "J-group " + std::to_string(j);

    std::vector<int> blockChannels, groupChannels;
    std::vector<double> boundary;
    for (int ch = 1; ch <= group->NumChannels(); ch++) {
      AChannel *channel = group->GetChannel(ch);
      phys::Channel c;
      c.partition = channel->GetPairNum() - 1;
      c.kind = ChannelKindOf(channel->GetRadType());
      c.l = channel->GetL();
      c.two_s = TwiceSpin(channel->GetS(), what + " channel " + std::to_string(ch) + " channel spin");
      c.radius_fm = compound.GetPair(channel->GetPairNum())->GetChRad();
      int index = -1;
      for (size_t k = 0; k < channels.size() && index < 0; k++)
        if (channels[k].partition == c.partition && channels[k].kind == c.kind && channels[k].l == c.l &&
            channels[k].two_s == c.two_s) index = (int)k;
      if (index < 0) {
        index = (int)channels.size();
        channels.push_back(c);
      }
      blockChannels.push_back(index);
      groupChannels.push_back(ch);
      boundary.push_back(withParameters ? channel->GetBoundaryCondition() : 0.0);
    }

    std::vector<double> energies, gammas;
    std::vector<int> groupLevels;
    for (int la = 1; la <= group->NumLevels(); la++) {
      ALevel *level = group->GetLevel(la);
      if (!level->IsInRMatrix() && !allGroups) continue;
      groupLevels.push_back(la);
      if (!withParameters) continue;
      energies.push_back(level->GetE());
      for (int ch = 1; ch <= group->NumChannels(); ch++) gammas.push_back(level->GetGamma(ch));
    }

    blocks.push_back(phys::JpiBlock(TwiceSpin(group->GetJ(), what + " J"), group->GetPi(), group->GetNucleus() - 1,
                                    std::move(blockChannels), std::move(boundary), std::move(energies),
                                    std::move(gammas)));
    map.groupOfBlock.push_back(j);
    map.levelOfBlockLevel.push_back(groupLevels);
    map.channelOfBlockChannel.push_back(groupChannels);
  }

  return CNucSystem{phys::System(BasisOf(configure), std::move(nuclei), std::move(partitions), std::move(channels),
                                 std::move(blocks)),
                    std::move(map)};
}
