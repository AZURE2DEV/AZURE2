// Unit test of azure2::physics::System and JpiBlock, using physics/ alone.
//
// Checks the storage conventions (0-based indices, row-major gamma[lambda * n_c + k], doubled spins), that the
// objects are values (a copy is independent and equal), and that each structural error the constructors guard
// against is refused with std::invalid_argument rather than stored.  An error in the gamma layout, or a constructor
// that accepts a channel of another nucleus or a block with too few amplitudes, fails here.
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "azure2/physics/System.h"

using namespace azure2::physics;

namespace {

int failures = 0;

void check(const std::string &name, bool ok) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", name.c_str());
  if (!ok) failures++;
}

void refuses(const std::string &name, const std::function<void()> &build) {
  bool threw = false;
  try {
    build();
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  check("refuses " + name, threw);
}

// p + 12C -> 13N in miniature: one particle partition, one photon partition, two blocks.
std::vector<Nucleus> Nuclei() { return {Nucleus{7, 13}}; }

std::vector<Partition> Partitions() {
  Partition p;
  p.kind = PartitionKind::Particle;
  p.nucleus = 0;
  p.Z = {{1, 6}};
  p.mass_u = {{1.00782503223, 12.0}};
  p.two_j = {{1, 0}};
  p.parity = {{1, 1}};
  p.g_factor = {{5.5856946893, 0.0}};
  p.excitation_MeV = {{0.0, 0.0}};
  p.separation_MeV = 1.94349;
  Partition g = p;
  g.kind = PartitionKind::Photon;
  g.Z = {{0, 7}};
  g.mass_u = {{0.0, 13.005738609}};
  g.two_j = {{2, 1}};
  g.parity = {{1, -1}};
  g.separation_MeV = 0.0;
  return {p, g};
}

std::vector<Channel> Channels() {
  return {Channel{0, ChannelKind::Particle, 0, 1, 3.4},
          Channel{1, ChannelKind::Electric, 1, 1, 3.4},
          Channel{0, ChannelKind::Particle, 1, 1, 3.4}};
}

// 1/2+ : two levels, channels 0 (p, l = 0) and 1 (E1 to the ground state); 3/2- : one level, channel 2.
std::vector<JpiBlock> Blocks() {
  return {JpiBlock(1, 1, 0, {0, 1}, {-0.5, -0.5}, {2.365, 8.0}, {1.1, 0.01, 2.2, 0.02}),
          JpiBlock(3, -1, 0, {2}, {-0.9}, {3.5}, {0.7})};
}

}  // namespace

int main() {
  static_assert(std::is_copy_constructible<System>::value && std::is_copy_assignable<System>::value,
                "System is a value type");
  static_assert(std::is_copy_constructible<JpiBlock>::value, "JpiBlock is a value type");

  const System sys(LevelBasis::Brune, Nuclei(), Partitions(), Channels(), Blocks());
  std::printf("structure\n");
  check("counts", sys.num_nuclei() == 1 && sys.num_partitions() == 2 && sys.num_channels() == 3 &&
                      sys.num_blocks() == 2 && sys.num_levels() == 3);
  check("basis", sys.basis() == LevelBasis::Brune);
  const JpiBlock &b0 = sys.block(0);
  check("block 0: 2J = 1, parity +", b0.two_J() == 1 && b0.parity() == 1 && b0.nucleus() == 0);
  check("block 0: two levels, two channels", b0.num_levels() == 2 && b0.num_channels() == 2);
  // row-major: gamma(lambda, k) = gamma()[lambda * n_c + k]
  check("gamma(0,0), (0,1), (1,0), (1,1) row-major",
        b0.gamma(0, 0) == 1.1 && b0.gamma(0, 1) == 0.01 && b0.gamma(1, 0) == 2.2 && b0.gamma(1, 1) == 0.02);
  check("E and B by index", b0.E(1) == 8.0 && b0.B(1) == -0.5 && sys.block(1).B(0) == -0.9);
  check("block channel list indexes System channels", sys.channel(b0.channel(1)).kind == ChannelKind::Electric);
  check("find_channel in block", b0.find_channel(1) == 1 && b0.find_channel(2) == -1);
  check("find_channel in System", sys.find_channel(0, ChannelKind::Particle, 1, 1) == 2 &&
                                      sys.find_channel(1, ChannelKind::Magnetic, 1, 1) == -1);
  check("hard-sphere block allowed", JpiBlock(5, 1, 0, {2}, {0.0}, {}, {}).num_levels() == 0);

  std::printf("value semantics\n");
  System copy = sys;
  check("copy equal", copy.block(0).gamma() == sys.block(0).gamma() && copy.partition(0).mass_u == sys.partition(0).mass_u);
  check("copy has its own storage", &copy.block(0).gamma()[0] != &sys.block(0).gamma()[0]);

  std::printf("refusals\n");
  refuses("parity 0", [] { JpiBlock(1, 0, 0, {}, {}, {}, {}); });
  refuses("negative 2J", [] { JpiBlock(-1, 1, 0, {}, {}, {}, {}); });
  refuses("B count != channel count", [] { JpiBlock(1, 1, 0, {0, 1}, {0.0}, {1.0}, {1.0, 1.0}); });
  refuses("gamma size != levels x channels", [] { JpiBlock(1, 1, 0, {0, 1}, {0., 0.}, {1.0}, {1.0}); });
  refuses("channel listed twice in a block", [] { JpiBlock(1, 1, 0, {0, 0}, {0., 0.}, {}, {}); });
  refuses("block channel out of range", [] {
    System(LevelBasis::Formal, Nuclei(), Partitions(), Channels(), {JpiBlock(1, 1, 0, {3}, {0.}, {}, {})});
  });
  refuses("block nucleus out of range", [] {
    System(LevelBasis::Formal, Nuclei(), Partitions(), Channels(), {JpiBlock(1, 1, 1, {}, {}, {}, {})});
  });
  refuses("two blocks with the same J^pi and nucleus", [] {
    System(LevelBasis::Formal, Nuclei(), Partitions(), Channels(),
           {JpiBlock(1, 1, 0, {}, {}, {}, {}), JpiBlock(1, 1, 0, {}, {}, {}, {})});
  });
  refuses("channel partition out of range", [] {
    std::vector<Channel> c = Channels();
    c[0].partition = 2;
    System(LevelBasis::Formal, Nuclei(), Partitions(), c, {});
  });
  refuses("particle channel on a photon partition", [] {
    std::vector<Channel> c = Channels();
    c[1].kind = ChannelKind::Particle;
    System(LevelBasis::Formal, Nuclei(), Partitions(), c, {});
  });
  refuses("duplicate channel", [] {
    std::vector<Channel> c = Channels();
    c[2].l = 0;
    System(LevelBasis::Formal, Nuclei(), Partitions(), c, {});
  });
  refuses("partition nucleus out of range", [] {
    std::vector<Partition> p = Partitions();
    p[1].nucleus = 1;
    System(LevelBasis::Formal, Nuclei(), p, Channels(), {});
  });
  refuses("block channel of another nucleus", [] {
    std::vector<Nucleus> n = Nuclei();
    n.push_back(Nucleus{8, 14});
    System(LevelBasis::Formal, n, Partitions(), Channels(), {JpiBlock(1, 1, 1, {0}, {0.}, {}, {})});
  });

  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
