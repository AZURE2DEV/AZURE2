// The reaction pathways of an entrance partition (Pathways.h; docs/upgrade/cnuc_split.md section 10).  The loops and
// their order are those of CNuc::SortPathways before step 6a, over System indices instead of CNuc numbers.
#include "azure2/physics/Pathways.h"

#include <cmath>
#include <cstdlib>
#include <set>
#include <tuple>

namespace azure2 {
namespace physics {

namespace {

// Index of the exit entry of partition alpha', appended when absent.
std::size_t ExitOf(std::vector<ExitPathways> &out, std::size_t partition) {
  for (std::size_t i = 0; i < out.size(); i++)
    if (out[i].partition == partition) return i;
  out.push_back(ExitPathways{partition, {}});
  return out.size() - 1;
}

const ExitPathways *FindExit(const std::vector<ExitPathways> &out, std::size_t partition, std::size_t &index) {
  for (std::size_t i = 0; i < out.size(); i++)
    if (out[i].partition == partition) {
      index = i;
      return &out[i];
    }
  return nullptr;
}

// Index of the group (s, s'[, s2']), appended when absent; a group without s2' matches on (s, s') alone.
std::size_t GroupOf(ExitPathways &exit, int two_s, int two_sp, bool has_sp2, int two_sp2) {
  for (std::size_t i = 0; i < exit.groups.size(); i++) {
    const SpinGroup &g = exit.groups[i];
    if (g.two_s == two_s && g.two_sp == two_sp && (!has_sp2 || g.two_sp2 == two_sp2)) return i;
  }
  SpinGroup g;
  g.two_s = two_s;
  g.two_sp = two_sp;
  g.has_sp2 = has_sp2;
  g.two_sp2 = has_sp2 ? two_sp2 : 0;
  exit.groups.push_back(g);
  return exit.groups.size() - 1;
}

// The internal pathway (b, k, k') of a group, appended when absent.
InternalPathway &PathwayOf(SpinGroup &group, std::size_t b, std::size_t k, std::size_t kp) {
  for (InternalPathway &p : group.internal)
    if (p.block == b && p.k_in == k && p.k_out == kp) return p;
  group.internal.push_back(InternalPathway{b, k, kp, false, 0.0});
  return group.internal.back();
}

}  // namespace

std::vector<ExitPathways> reaction_pathways(const System &system, std::size_t alpha,
                                            const std::vector<FinalState> &final_states,
                                            const PathwayOptions &options) {
  std::vector<ExitPathways> out;
  const std::vector<Partition> &partitions = system.partitions();
  const std::vector<Channel> &channels = system.channels();
  const Partition &entrance = partitions.at(alpha);
  const bool betaEntrance = entrance.kind == PartitionKind::Beta;

  auto channelOf = [&](const JpiBlock &block, std::size_t k) -> const Channel & { return channels[block.channel(k)]; };
  // (2J + 1) / ((2 I_1' + 1)(2 I_2' + 1)) of exit partition alpha', as the engine has always defined it.
  auto statSpinFactor = [&](const JpiBlock &block, std::size_t exit) {
    const Partition &p = partitions[exit];
    return (block.two_J() + 1.) * (1.0 / (p.two_j[0] + 1.0) / (p.two_j[1] + 1.0));
  };

  for (std::size_t ir = 0; ir < partitions.size(); ir++) {
    const Partition &exitPartition = partitions[ir];
    if (exitPartition.kind == PartitionKind::Beta) continue;
    if (betaEntrance) {
      for (int l = 0; l < 2; l++)
        for (std::size_t b = 0; b < system.blocks().size(); b++) {
          const JpiBlock &block = system.block(b);
          for (std::size_t k = 0; k < block.num_channels(); k++) {
            if ((std::size_t)channelOf(block, k).partition != alpha) continue;
            for (std::size_t kp = 0; kp < block.num_channels(); kp++) {
              if ((std::size_t)channelOf(block, kp).partition != ir || channelOf(block, k).l != l) continue;
              ExitPathways &exit = out[ExitOf(out, ir)];
              PathwayOf(exit.groups[GroupOf(exit, 2 * l, 0, false, 0)], b, k, kp);
            }
          }
        }
    } else if (exitPartition.kind == PartitionKind::Particle && ir == alpha) {
      for (int ts = std::abs(entrance.two_j[0] - entrance.two_j[1]); ts <= entrance.two_j[0] + entrance.two_j[1]; ts += 2)
        for (int tsp = std::abs(exitPartition.two_j[0] - exitPartition.two_j[1]);
             tsp <= exitPartition.two_j[0] + exitPartition.two_j[1]; tsp += 2)
          for (std::size_t b = 0; b < system.blocks().size(); b++) {
            const JpiBlock &block = system.block(b);
            for (std::size_t k = 0; k < block.num_channels(); k++) {
              if ((std::size_t)channelOf(block, k).partition != alpha) continue;
              for (std::size_t kp = 0; kp < block.num_channels(); kp++) {
                if ((std::size_t)channelOf(block, kp).partition != ir || channelOf(block, k).two_s != ts ||
                    channelOf(block, kp).two_s != tsp)
                  continue;
                ExitPathways &exit = out[ExitOf(out, ir)];
                InternalPathway &p = PathwayOf(exit.groups[GroupOf(exit, ts, tsp, false, 0)], b, k, kp);
                p.has_stat_spin_factor = true;
                p.stat_spin_factor = statSpinFactor(block, ir);
              }
            }
          }
    } else if (exitPartition.kind == PartitionKind::Particle) {
      // Unobserved primary, observed secondary: one group per (s, s', s2'); a group with s' != s2' holds the
      // pathways of both exit channel spins, a group with s' = s2' those of s' alone.
      const int tsMin = std::abs(exitPartition.two_j[0] - exitPartition.two_j[1]);
      const int tsMax = exitPartition.two_j[0] + exitPartition.two_j[1];
      for (int ts = std::abs(entrance.two_j[0] - entrance.two_j[1]); ts <= entrance.two_j[0] + entrance.two_j[1]; ts += 2)
        for (int tsp = tsMin; tsp <= tsMax; tsp += 2)
          for (int tsp2 = tsMin; tsp2 <= tsMax; tsp2 += 2)
            for (std::size_t b = 0; b < system.blocks().size(); b++) {
              const JpiBlock &block = system.block(b);
              for (std::size_t k = 0; k < block.num_channels(); k++) {
                if ((std::size_t)channelOf(block, k).partition != alpha) continue;
                for (std::size_t kp = 0; kp < block.num_channels(); kp++) {
                  if ((std::size_t)channelOf(block, kp).partition != ir || channelOf(block, k).two_s != ts ||
                      (channelOf(block, kp).two_s != tsp && channelOf(block, kp).two_s != tsp2))
                    continue;
                  ExitPathways &exit = out[ExitOf(out, ir)];
                  InternalPathway &p = PathwayOf(exit.groups[GroupOf(exit, ts, tsp, true, tsp2)], b, k, kp);
                  p.has_stat_spin_factor = true;
                  p.stat_spin_factor = statSpinFactor(block, ir);
                }
              }
            }
    } else if (exitPartition.kind == PartitionKind::Photon && !options.rmatrix_capture) {
      for (int ts = std::abs(entrance.two_j[0] - entrance.two_j[1]); ts <= entrance.two_j[0] + entrance.two_j[1]; ts += 2)
        for (std::size_t b = 0; b < system.blocks().size(); b++) {
          const JpiBlock &block = system.block(b);
          for (std::size_t k = 0; k < block.num_channels(); k++) {
            if ((std::size_t)channelOf(block, k).partition != alpha) continue;
            for (std::size_t kp = 0; kp < block.num_channels(); kp++) {
              if ((std::size_t)channelOf(block, kp).partition != ir || channelOf(block, k).two_s != ts) continue;
              ExitPathways &exit = out[ExitOf(out, ir)];
              InternalPathway &p = PathwayOf(exit.groups[GroupOf(exit, ts, 0, false, 0)], b, k, kp);
              p.has_stat_spin_factor = true;
              p.stat_spin_factor = statSpinFactor(block, ir);
            }
          }
        }
    }
  }

  if (betaEntrance) return out;

  // External capture.
  for (std::size_t f = 0; f < final_states.size(); f++) {
    const FinalState &fs = final_states[f];
    std::size_t photonExit = 0;
    if (!FindExit(out, (std::size_t)fs.photon_partition, photonExit)) continue;
    for (std::size_t g = 0; g < out[photonExit].groups.size(); g++) {
      const int two_s = out[photonExit].groups[g].two_s;
      for (std::size_t cf = 0; cf < fs.channels.size(); cf++) {
        const FinalStateChannel &finalChannel = fs.channels[cf];
        if (!finalChannel.particle) continue;
        std::size_t ccExit = 0;
        if (!FindExit(out, (std::size_t)finalChannel.partition, ccExit)) continue;
        // A UPOS decay holds, for each (s, s'), one group per s2' with identical pathways: count each pathway once.
        std::set<std::tuple<std::size_t, std::size_t, std::size_t>> seen;
        for (std::size_t kp = 0; kp < out[ccExit].groups.size(); kp++) {
          if (out[ccExit].groups[kp].two_s != two_s) continue;
          for (std::size_t mp = 0; mp < out[ccExit].groups[kp].internal.size(); mp++) {
            const InternalPathway cc = out[ccExit].groups[kp].internal[mp];
            if (!seen.insert(std::make_tuple(cc.block, cc.k_in, cc.k_out)).second) continue;
            const JpiBlock &block = system.block(cc.block);
            const Channel &ch = channelOf(block, cc.k_in);
            const Channel &chp = channelOf(block, cc.k_out);
            for (int L = 1; L <= 2; L++) {
              const Multipole X = (block.parity() * fs.parity == ((L % 2) ? -1 : 1)) ? Multipole::E : Multipole::M;
              if (!((X == Multipole::E && L == 1) && (fs.multipoles & kCaptureE1)) &&
                  !((X == Multipole::M && L == 1) && (fs.multipoles & kCaptureM1)) &&
                  !((X == Multipole::E && L == 2) && (fs.multipoles & kCaptureE2)))
                continue;
              if (std::abs(block.two_J() - 2 * L) > fs.two_J || fs.two_J > block.two_J() + 2 * L) continue;
              if (!couples(X, L, chp.l, chp.two_s / 2.0, finalChannel.l, finalChannel.two_s / 2.0, fs.two_J / 2.0))
                continue;
              ExternalPathway e;
              e.X = X;
              e.L = L;
              e.l = ch.l;
              e.two_J = block.two_J();
              e.final_state = f;
              e.final_channel = cf;
              e.cc_exit = ccExit;
              e.cc_group = kp;
              e.cc_pathway = mp;
              if (cc.k_in == cc.k_out) out[photonExit].groups[g].external.push_back(e);
              const ChannelKind photonKind = (X == Multipole::E) ? ChannelKind::Electric : ChannelKind::Magnetic;
              for (std::size_t k = 0; k < block.num_channels(); k++) {
                const Channel &c = channelOf(block, k);
                if (c.kind == photonKind && c.l == L && c.partition == fs.photon_partition) {
                  e.photon_channel = (int)k;
                  break;
                }
              }
              e.channel_capture = true;
              out[photonExit].groups[g].external.push_back(e);
            }
          }
        }
      }
    }
  }
  return out;
}

}  // namespace physics
}  // namespace azure2
