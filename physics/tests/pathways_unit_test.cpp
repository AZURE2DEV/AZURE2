// reaction_pathways (Pathways.h; docs/upgrade/cnuc_split.md section 10 (a)) on synthetic Systems, against properties
// the enumeration must have, checked without its loops:
//  1. every internal pathway joins an entrance channel of alpha and an exit channel of the entry's partition in one
//     block, with the group's channel spins (s; s' or, with s2', s' or s2'; beta: l and 0);
//  2. completeness: every pair of block channels (entrance alpha, exit alpha') appears exactly once per group it
//     belongs to, and the number of internal pathways equals a direct count over the block channels;
//  3. the statistical spin factor is (2J + 1) / ((2I1' + 1)(2I2' + 1)) of the exit;
//  4. every external pathway satisfies the parity, mask, J and coupling rules, and refers to an internal pathway with
//     the final channel's partition as exit and the group's entrance spin; direct capture only for k_in = k_out;
//  5. the R-matrix-capture option removes the photon exit and every external pathway.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "azure2/physics/Pathways.h"

using namespace azure2::physics;

namespace {

int failures = 0;
long checks = 0;

void check(const std::string &what, bool ok) {
  checks++;
  if (!ok) {
    std::printf("  FAIL  %s\n", what.c_str());
    failures++;
  }
}

Partition MakePartition(PartitionKind kind, int tj0, int tj1) {
  Partition p;
  p.kind = kind;
  p.Z = {{kind == PartitionKind::Particle ? 1 : 0, 6}};
  p.mass_u = {{1.007825, 12.0}};
  p.two_j = {{tj0, tj1}};
  p.separation_MeV = 1.0;
  return p;
}

// p + 12C-like entrance with spin (2I = 1, 3), an inelastic-like particle exit (1, 1), a photon exit.
struct Model {
  std::vector<Partition> partitions;
  std::vector<Channel> channels;
  std::vector<JpiBlock> blocks;
};

Model MakeModel() {
  Model m;
  m.partitions = {MakePartition(PartitionKind::Particle, 1, 3), MakePartition(PartitionKind::Particle, 1, 1),
                  MakePartition(PartitionKind::Photon, 2, 1)};
  auto add = [&](int partition, ChannelKind kind, int l, int two_s) {
    for (std::size_t c = 0; c < m.channels.size(); c++)
      if (m.channels[c].partition == partition && m.channels[c].kind == kind && m.channels[c].l == l &&
          m.channels[c].two_s == two_s)
        return (int)c;
    Channel c;
    c.partition = partition;
    c.kind = kind;
    c.l = l;
    c.two_s = two_s;
    c.radius_fm = 4.0;
    m.channels.push_back(c);
    return (int)m.channels.size() - 1;
  };
  // Blocks 2J^pi = 1+, 3+, 3-, 5-: all (l, s) of partitions 0 and 1 with l <= 3, and photon channels E1/M1/E2 to a
  // 1/2- final state (photon partition 2).
  const int twoJs[4] = {1, 3, 3, 5};
  const int parities[4] = {1, 1, -1, -1};
  for (int b = 0; b < 4; b++) {
    std::vector<int> ch;
    for (int a = 0; a < 2; a++)
      for (int ts = std::abs(m.partitions[a].two_j[0] - m.partitions[a].two_j[1]);
           ts <= m.partitions[a].two_j[0] + m.partitions[a].two_j[1]; ts += 2)
        for (int l = 0; l <= 3; l++)
          if ((l % 2 == 0 ? 1 : -1) == parities[b] && std::abs(2 * l - ts) <= twoJs[b] && twoJs[b] <= 2 * l + ts)
            ch.push_back(add(a, ChannelKind::Particle, l, ts));
    for (int L = 1; L <= 2; L++)
      if (std::abs(twoJs[b] - 1) <= 2 * L && 2 * L <= twoJs[b] + 1) {
        const bool electric = parities[b] * -1 == (L % 2 ? -1 : 1);
        ch.push_back(add(2, electric ? ChannelKind::Electric : ChannelKind::Magnetic, L, 1));
      }
    m.blocks.emplace_back(twoJs[b], parities[b], 0, ch, std::vector<double>(ch.size(), 0.0), std::vector<double>{},
                          std::vector<double>{});
  }
  return m;
}

}  // namespace

int main() {
  Model m = MakeModel();
  System system(LevelBasis::Brune, {Nucleus{7, 13}}, m.partitions, m.channels, m.blocks);
  // A 1/2- final state with particle channels of partitions 0 and 1 (all its l <= 3, s), multipoles E1, M1, E2.
  FinalState fs;
  fs.two_J = 1;
  fs.parity = -1;
  fs.photon_partition = 2;
  fs.multipoles = kCaptureE1 | kCaptureM1 | kCaptureE2;
  for (int a = 0; a < 2; a++)
    for (int ts = std::abs(m.partitions[a].two_j[0] - m.partitions[a].two_j[1]);
         ts <= m.partitions[a].two_j[0] + m.partitions[a].two_j[1]; ts += 2)
      for (int l = 0; l <= 3; l++)
        if (l % 2 == 1 && std::abs(2 * l - ts) <= 1 && 1 <= 2 * l + ts) fs.channels.push_back({a, true, l, ts});
  fs.channels.push_back({2, false, 1, 1});  // a non-particle channel, ignored

  for (int rmc = 0; rmc <= 1; rmc++) {
    PathwayOptions options;
    options.rmatrix_capture = rmc;
    const std::vector<ExitPathways> exits = reaction_pathways(system, 0, {fs}, options);
    const std::string tag = rmc ? "rmc: " : "";
    long internal = 0, external = 0;
    for (std::size_t e = 0; e < exits.size(); e++) {
      const ExitPathways &x = exits[e];
      const Partition &exitP = system.partitions()[x.partition];
      check(tag + "no beta exit", exitP.kind != PartitionKind::Beta);
      if (rmc) check(tag + "no photon exit", exitP.kind != PartitionKind::Photon);
      for (std::size_t g = 0; g < x.groups.size(); g++) {
        const SpinGroup &grp = x.groups[g];
        check(tag + "s2' only for a particle exit other than the entrance",
              grp.has_sp2 == (exitP.kind == PartitionKind::Particle && x.partition != 0));
        for (const InternalPathway &p : grp.internal) {
          internal++;
          const JpiBlock &block = system.block(p.block);
          const Channel &a = system.channels()[block.channel(p.k_in)];
          const Channel &b = system.channels()[block.channel(p.k_out)];
          check(tag + "1. entrance partition", a.partition == 0);
          check(tag + "1. exit partition", (std::size_t)b.partition == x.partition);
          check(tag + "1. entrance spin", a.two_s == grp.two_s);
          if (exitP.kind == PartitionKind::Photon)
            check(tag + "1. photon group s' = 0", grp.two_sp == 0);
          else
            check(tag + "1. exit spin", b.two_s == grp.two_sp || (grp.has_sp2 && b.two_s == grp.two_sp2));
          const double stat = (block.two_J() + 1.0) / ((exitP.two_j[0] + 1.0) * (exitP.two_j[1] + 1.0));
          check(tag + "3. statistical factor", p.has_stat_spin_factor && std::fabs(p.stat_spin_factor - stat) < 1e-15);
          int seen = 0;
          for (const InternalPathway &q : grp.internal) seen += (q.block == p.block && q.k_in == p.k_in && q.k_out == p.k_out);
          check(tag + "2. once per group", seen == 1);
        }
        for (const ExternalPathway &xp : grp.external) {
          external++;
          check(tag + "4. external only in the photon exit", exitP.kind == PartitionKind::Photon);
          check(tag + "4. final channel is a particle channel", fs.channels[xp.final_channel].particle);
          const ExitPathways &cx = exits[xp.cc_exit];
          check(tag + "4. channel capture through the final channel's partition",
                (int)cx.partition == fs.channels[xp.final_channel].partition);
          check(tag + "4. same entrance spin", cx.groups[xp.cc_group].two_s == grp.two_s);
          const InternalPathway &cc = cx.groups[xp.cc_group].internal[xp.cc_pathway];
          const JpiBlock &block = system.block(cc.block);
          check(tag + "4. compound J", xp.two_J == block.two_J());
          const bool electric = block.parity() * fs.parity == (xp.L % 2 ? -1 : 1);
          check(tag + "4. radiation type from the parities", (xp.X == Multipole::E) == electric);
          check(tag + "4. E1, M1 or E2", (xp.X == Multipole::E && xp.L <= 2) || (xp.X == Multipole::M && xp.L == 1));
          check(tag + "4. J triangle", std::abs(block.two_J() - 2 * xp.L) <= fs.two_J && fs.two_J <= block.two_J() + 2 * xp.L);
          const Channel &out = system.channels()[block.channel(cc.k_out)];
          const FinalStateChannel &fc = fs.channels[xp.final_channel];
          check(tag + "4. coupling rule", couples(xp.X, xp.L, out.l, out.two_s / 2.0, fc.l, fc.two_s / 2.0, fs.two_J / 2.0));
          check(tag + "4. entrance l", xp.l == system.channels()[block.channel(cc.k_in)].l);
          if (!xp.channel_capture) check(tag + "4. direct capture only for an elastic pathway", cc.k_in == cc.k_out);
          if (xp.channel_capture && xp.photon_channel >= 0) {
            const Channel &pc = system.channels()[block.channel(xp.photon_channel)];
            check(tag + "4. photon channel", pc.partition == 2 && pc.l == xp.L &&
                                                 pc.kind == (xp.X == Multipole::E ? ChannelKind::Electric : ChannelKind::Magnetic));
          }
        }
      }
    }
    // 2. Direct count: for each block, entrance channels of 0 times exit channels of each non-beta partition, once
    // per group the pair belongs to: elastic and photon one group; the particle exit one group per (s', s2') with s'
    // or s2' equal to the exit channel's spin.
    long expected = 0;
    for (std::size_t b = 0; b < system.blocks().size(); b++) {
      const JpiBlock &block = system.block(b);
      for (std::size_t k = 0; k < block.num_channels(); k++) {
        if (system.channels()[block.channel(k)].partition != 0) continue;
        for (std::size_t kp = 0; kp < block.num_channels(); kp++) {
          const Channel &c = system.channels()[block.channel(kp)];
          if (c.partition == 0) expected++;
          if (c.partition == 2 && !rmc) expected++;
          if (c.partition == 1) {
            const int nS = 2;  // 2s' in {0, 2} for spins 1/2, 1/2
            expected += 2 * nS - 1;  // (s', s2') pairs with s' or s2' = the channel's spin
          }
        }
      }
    }
    check(tag + "2. internal pathway count " + std::to_string(internal) + " = " + std::to_string(expected),
          internal == expected);
    if (rmc) check(tag + "5. no external pathways", external == 0);
    else check(tag + "external pathways found", external > 0);
    std::printf("  %s%ld exits, %ld internal, %ld external pathways\n", tag.c_str(), (long)exits.size(), internal, external);
  }
  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
