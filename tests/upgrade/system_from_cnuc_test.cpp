// The physics/ System built from today's CNuc (src/CNucSystem.cpp), on the upgrade test projects T1, T2 and T3.
//
// For each project, in each level-parameter basis (Brune, the default; Park; formal), the compound nucleus is loaded
// as the engine loads it (Config::ReadConfigFile, CNuc::Fill, CNuc::Initialize, with external capture on as in
// pyazr) and the System is built from it.  Checked:
//
//  1. every number of the System against the CNuc accessor it was taken from, with exact (==) equality: the nuclei
//     (Z, A); per partition Z, M, J, pi, g, e_light, e_heavy, S and the kind and nucleus; per channel the pair, the
//     radiation type, l, s and the radius; per block 2J, pi, nucleus, B_c, E_lambda and gamma_{lambda c};
//  2. completeness: every J^pi group in the R-matrix is one block and every R-matrix level of it one block level,
//     every group channel one block channel, and the System's channels are exactly the distinct (pair, kind, l, s)
//     of those groups;
//  3. the structure against the .azr files themselves, read off the <levels> lines by hand (the tables below): the
//     pair keys in order of first appearance, the number of channels, and per block 2J, pi, the channels as
//     (pair key, kind, l, 2s) and, in the Brune and Park bases (where the level energies are the observed ones), the
//     level energies E_x of the file;
//  4. that the adapter refuses a CNuc whose level parameters are still the input widths (Fill without Initialize).
//
// A wrong accessor, an off-by-one between CNuc's 1-based and the System's 0-based indices, a transposed gamma
// layout, a channel attached to the wrong block or a level dropped fails (1) to (3).
//
//   system_from_cnuc_test <T1.azr> <T2.azr> <T3.azr>      (ctest: physics_system_from_cnuc)
#include <cstdio>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "AZUREParams.h"
#include "CNuc.h"
#include "CNucSystem.h"
#include "Config.h"
#include "ECAmplitudeCache.h"

namespace phys = azure2::physics;

namespace {

int failures = 0;
long compared = 0;  // numbers compared with ==

void check(const std::string &name, bool ok, const std::string &detail = "") {
  if (!ok) {
    std::printf("  FAIL  %s%s%s\n", name.c_str(), detail.empty() ? "" : "  -- ", detail.c_str());
    failures++;
  }
}

template <class A, class B>
void same(const std::string &name, const A &system, const B &cnuc) {
  compared++;
  if (!(system == cnuc)) {
    std::ostringstream s;
    s.precision(17);
    s << "System " << system << ", CNuc " << cnuc;
    check(name, false, s.str());
  }
}

char KindLetter(phys::ChannelKind k) {
  switch (k) {
    case phys::ChannelKind::Particle: return 'P';
    case phys::ChannelKind::Electric: return 'E';
    case phys::ChannelKind::Magnetic: return 'M';
    case phys::ChannelKind::Fermi: return 'F';
    case phys::ChannelKind::GamowTeller: return 'G';
  }
  return '?';
}

int PairType(phys::PartitionKind k) {
  switch (k) {
    case phys::PartitionKind::Particle: return 0;
    case phys::PartitionKind::Photon: return 10;
    case phys::PartitionKind::Beta: return 20;
  }
  return -1;
}

// ---- the structure of each project, read off its <levels> lines (columns: io/include/azure2/io/LevelRows.h; 2s and 2l in the
// file).  Photon channels: E when pi_level * pi_final = (-1)^L, else M (AChannel's rule).
struct ExpChannel {
  int key;
  char kind;
  int l;
  int two_s;
};
struct ExpBlock {
  int two_J;
  int parity;
  std::vector<ExpChannel> channels;
  std::vector<double> energies;  // E_x (MeV) as written in the file
};
struct ExpProject {
  std::string name;
  std::vector<int> pairKeys;  // in order of first appearance
  std::vector<int> pairTypes;
  int numChannels;
  std::vector<ExpBlock> blocks;
};

const std::vector<ExpProject> kExpected = {
    {"T1",
     {1, 2, 3},
     {0, 10, 10},
     8,
     {{1, -1, {{1, 'P', 1, 1}}, {0.0}},
      {3, -1, {{1, 'P', 1, 1}, {2, 'M', 1, 1}, {2, 'E', 2, 1}, {3, 'M', 1, 3}}, {1.0, 3.502, 15.0}},
      {1, 1, {{1, 'P', 0, 1}, {2, 'E', 1, 1}, {3, 'E', 1, 3}}, {2.365, 15.0}},
      {5, 1, {{1, 'P', 2, 1}, {3, 'E', 1, 3}}, {3.547}}}},
    {"T2",
     {1, 2, 3},
     {0, 0, 0},
     6,
     {{2, -1, {{1, 'P', 0, 2}, {1, 'P', 2, 2}, {2, 'P', 1, 0}, {3, 'P', 1, 4}}, {12.4258, 13.105, 17.0}},
      {6, -1, {{1, 'P', 2, 2}, {1, 'P', 4, 2}, {2, 'P', 3, 0}, {3, 'P', 1, 4}}, {13.0734, 13.2652}}}},
    {"T3",
     {1, 4, 2, 3},
     {0, 20, 0, 0},
     11,
     {{4, 1, {{1, 'P', 2, 0}, {4, 'G', 1, 2}, {2, 'P', 1, 4}, {3, 'P', 1, 4}}, {3.03, 22.0}},
      {8, 1, {{1, 'P', 4, 0}}, {11.35}},
      {2, 1, {{2, 'P', 1, 2}, {2, 'P', 1, 4}}, {17.64}},
      {4, -1, {{2, 'P', 0, 4}, {3, 'P', 0, 4}}, {19.0}},
      {0, 1, {{1, 'P', 0, 0}}, {25.0}},
      {2, -1, {{2, 'P', 0, 2}, {3, 'P', 0, 2}}, {22.0}}}},
};

enum class Basis { Brune, Park, Formal };
const char *BasisName(Basis b) { return b == Basis::Brune ? "Brune" : b == Basis::Park ? "Park" : "formal"; }

struct Loaded {
  std::ostringstream sink;
  Config cfg;
  CNuc compound;
  LevelState state;  // the levels, as the run holds them
  LevelInputs inputs;
  Loaded() : cfg(sink) {}
};

// As AZUREAPI::Initialize loads the compound nucleus (no data: nothing of the System depends on them).
bool Load(Loaded &m, const std::string &azr, Basis basis, bool initialize) {
  m.cfg.configfile = azr;
  if (basis == Basis::Park) m.cfg.paramMask |= (Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  if (basis == Basis::Formal) m.cfg.paramMask &= ~(Config::USE_PARK_FORMALISM | Config::USE_BRUNE_FORMALISM);
  m.cfg.paramMask |= Config::USE_EXTERNAL_CAPTURE;
  if (m.cfg.ReadConfigFile() < 0) return false;
  m.cfg.ResolveMatrixRoute();
  InitializeECAmplitudeCache(m.cfg);
  if (m.compound.Fill(m.cfg, m.state, m.inputs) == -1) return false;
  if (initialize) m.compound.Initialize(m.cfg, m.state, m.inputs);
  return true;
}

void CheckProject(const ExpProject &exp, const std::string &azr, Basis basis) {
  const std::string tag = exp.name + " (" + BasisName(basis) + ")";
  Loaded m;
  if (!Load(m, azr, basis, true)) {
    check(tag + ": load", false, m.sink.str());
    return;
  }
  CNuc &cn = m.compound;
  const long before = compared;
  const int failuresBefore = failures;
  CNucSystem built = BuildSystemFromCNuc(cn, m.inputs, m.cfg);
  const phys::System &sys = built.system;
  const CNucSystemMap &map = built.map;

  // basis
  const int b = AZUREParams::BasisForConfig(m.cfg);
  check(tag + ": basis", (b == AZUREParams::kBasisBrune && sys.basis() == phys::LevelBasis::Brune) ||
                             (b == AZUREParams::kBasisPark && sys.basis() == phys::LevelBasis::Park) ||
                             (b == AZUREParams::kBasisStandard && sys.basis() == phys::LevelBasis::Formal));
  check(tag + ": basis is the one asked for", (basis == Basis::Brune && sys.basis() == phys::LevelBasis::Brune) ||
                                                 (basis == Basis::Park && sys.basis() == phys::LevelBasis::Park) ||
                                                 (basis == Basis::Formal && sys.basis() == phys::LevelBasis::Formal));

  // (1) nuclei
  same(tag + ": number of nuclei", (int)sys.num_nuclei(), cn.NumNuclei());
  for (int i = 0; i < (int)sys.num_nuclei() && i < cn.NumNuclei(); i++) {
    same(tag + ": nucleus Z", sys.nucleus(i).Z, cn.GetNucleus(i + 1)->Z);
    same(tag + ": nucleus A", sys.nucleus(i).A, cn.GetNucleus(i + 1)->A);
  }

  // (1) partitions
  same(tag + ": number of partitions", (int)sys.num_partitions(), cn.NumPairs());
  for (int a = 0; a < (int)sys.num_partitions() && a < cn.NumPairs(); a++) {
    const phys::Partition &p = sys.partition(a);
    const PPair *pair = cn.GetPair(map.pairOfPartition[a]);
    same(tag + ": partition -> pair number", map.pairOfPartition[a], a + 1);
    const std::string w = tag + ": pair " + std::to_string(pair->GetPairKey()) + " ";
    same(w + "type", PairType(p.kind), pair->GetPType());
    same(w + "nucleus", p.nucleus + 1, pair->GetNucleus());
    for (int i = 0; i < 2; i++) {
      same(w + "Z", p.Z[i], pair->GetZ(i + 1));
      same(w + "M", p.mass_u[i], pair->GetM(i + 1));
      same(w + "2J", (double)p.two_j[i], 2.0 * pair->GetJ(i + 1));
      same(w + "pi", p.parity[i], pair->GetPi(i + 1));
      same(w + "g", p.g_factor[i], pair->GetG(i + 1));
    }
    same(w + "e_light", p.excitation_MeV[0], pair->GetLightExE());
    same(w + "e_heavy", p.excitation_MeV[1], pair->GetHeavyExE());
    same(w + "e_light + e_heavy", p.excitation_MeV[0] + p.excitation_MeV[1], pair->GetExE());
    same(w + "S", p.separation_MeV, pair->GetSepE());
  }

  // (1), (2) blocks, their channels and levels
  std::vector<int> rmGroups;
  for (int j = 1; j <= cn.NumJGroups(); j++)
    if (cn.GetJGroup(j)->IsInRMatrix()) rmGroups.push_back(j);
  same(tag + ": blocks = J-groups in the R-matrix", (int)sys.num_blocks(), (int)rmGroups.size());
  same(tag + ": map blocks", (int)map.groupOfBlock.size(), (int)sys.num_blocks());
  std::set<std::tuple<int, char, int, double> > distinct;  // (pair number, radiation type, l, s) of those groups
  for (int bl = 0; bl < (int)sys.num_blocks() && bl < (int)rmGroups.size(); bl++) {
    const phys::JpiBlock &blk = sys.block(bl);
    same(tag + ": block -> J-group", map.groupOfBlock[bl], rmGroups[bl]);
    const JGroup *g = cn.GetJGroup(map.groupOfBlock[bl]);
    const std::string w = tag + ": block " + std::to_string(bl) + " ";
    same(w + "2J", (double)blk.two_J(), 2.0 * g->GetJ());
    same(w + "pi", blk.parity(), g->GetPi());
    same(w + "nucleus", blk.nucleus() + 1, g->GetNucleus());
    same(w + "channels", (int)blk.num_channels(), g->NumChannels());
    for (int k = 0; k < (int)blk.num_channels() && k < g->NumChannels(); k++) {
      same(w + "block channel -> group channel", map.channelOfBlockChannel[bl][k], k + 1);
      const AChannel *ch = g->GetChannel(map.channelOfBlockChannel[bl][k]);
      const phys::Channel &c = sys.channel(blk.channel(k));
      distinct.insert(std::make_tuple(ch->GetPairNum(), ch->GetRadType(), ch->GetL(), ch->GetS()));
      same(w + "channel pair", c.partition + 1, ch->GetPairNum());
      same(w + "channel kind", KindLetter(c.kind), ch->GetRadType());
      same(w + "channel l", c.l, ch->GetL());
      same(w + "channel 2s", (double)c.two_s, 2.0 * ch->GetS());
      same(w + "channel radius", c.radius_fm, cn.GetPair(ch->GetPairNum())->GetChRad());
      same(w + "B_c", blk.B(k), ch->GetBoundaryCondition());
    }
    std::vector<int> rmLevels;
    for (int la = 1; la <= g->NumLevels(); la++)
      if (g->GetLevel(la)->IsInRMatrix()) rmLevels.push_back(la);
    same(w + "levels = R-matrix levels of the group", (int)blk.num_levels(), (int)rmLevels.size());
    same(w + "gamma size", blk.gamma().size(), blk.num_levels() * blk.num_channels());
    for (int lam = 0; lam < (int)blk.num_levels() && lam < (int)rmLevels.size(); lam++) {
      same(w + "block level -> group level", map.levelOfBlockLevel[bl][lam], rmLevels[lam]);
      const LevelInputs &in = m.inputs;
      const int j = map.groupOfBlock[bl], la = rmLevels[lam];
      same(w + "E", blk.E(lam), in.E(j, la));
      for (int k = 0; k < (int)blk.num_channels(); k++) {
        same(w + "gamma", blk.gamma(lam, k), in.Gamma(j, la, k + 1));
        same(w + "gamma row-major", blk.gamma()[lam * blk.num_channels() + k], in.Gamma(j, la, k + 1));
      }
    }
  }
  same(tag + ": System channels = distinct (pair, type, l, s) of the blocks", (int)sys.num_channels(), (int)distinct.size());

  // (3) structure against the file
  same(tag + ": partitions (file)", (int)sys.num_partitions(), (int)exp.pairKeys.size());
  for (int a = 0; a < (int)sys.num_partitions() && a < (int)exp.pairKeys.size(); a++) {
    same(tag + ": pair key order (file)", cn.GetPair(map.pairOfPartition[a])->GetPairKey(), exp.pairKeys[a]);
    same(tag + ": pair type (file)", PairType(sys.partition(a).kind), exp.pairTypes[a]);
  }
  same(tag + ": channels (file)", (int)sys.num_channels(), exp.numChannels);
  same(tag + ": blocks (file)", (int)sys.num_blocks(), (int)exp.blocks.size());
  for (int bl = 0; bl < (int)sys.num_blocks() && bl < (int)exp.blocks.size(); bl++) {
    const phys::JpiBlock &blk = sys.block(bl);
    const ExpBlock &e = exp.blocks[bl];
    const std::string w = tag + ": block " + std::to_string(bl) + " (file) ";
    same(w + "2J", blk.two_J(), e.two_J);
    same(w + "pi", blk.parity(), e.parity);
    same(w + "channels", (int)blk.num_channels(), (int)e.channels.size());
    for (int k = 0; k < (int)blk.num_channels() && k < (int)e.channels.size(); k++) {
      const phys::Channel &c = sys.channel(blk.channel(k));
      same(w + "channel pair key", cn.GetPair(map.pairOfPartition[c.partition])->GetPairKey(), e.channels[k].key);
      same(w + "channel kind", KindLetter(c.kind), e.channels[k].kind);
      same(w + "channel l", c.l, e.channels[k].l);
      same(w + "channel 2s", c.two_s, e.channels[k].two_s);
    }
    same(w + "levels", (int)blk.num_levels(), (int)e.energies.size());
    if (basis != Basis::Formal)
      for (int lam = 0; lam < (int)blk.num_levels() && lam < (int)e.energies.size(); lam++)
        same(w + "E_x", blk.E(lam), e.energies[lam]);
  }

  std::printf("  %-14s %2zu partitions, %2zu channels, %zu blocks, %2zu levels: %ld numbers compared, %s\n",
              tag.c_str(), sys.num_partitions(), sys.num_channels(), sys.num_blocks(), sys.num_levels(),
              compared - before, failures == failuresBefore ? "all equal" : "MISMATCH");
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 4) {
    std::fprintf(stderr, "usage: %s T1.azr T2.azr T3.azr\n", argv[0]);
    return 2;
  }
  std::printf("System from CNuc, every number against the CNuc accessor (exact) and the structure against the file\n");
  for (int t = 0; t < 3; t++)
    for (Basis basis : {Basis::Brune, Basis::Park, Basis::Formal}) {
      try {
        CheckProject(kExpected[t], argv[t + 1], basis);
      } catch (const std::exception &e) {
        check(kExpected[t].name + " (" + BasisName(basis) + ")", false, std::string("threw: ") + e.what());
      }
    }

  std::printf("refusal\n");
  {
    Loaded m;
    bool threw = false;
    if (Load(m, argv[1], Basis::Brune, false)) {
      try {
        BuildSystemFromCNuc(m.compound, m.inputs, m.cfg);
      } catch (const std::runtime_error &) {
        threw = true;
      }
    }
    check("T1 before CNuc::Initialize (widths still in eV) is refused", threw);
    if (threw) std::printf("  ok    T1 before CNuc::Initialize (widths still in eV) is refused\n");
  }

  std::printf("%s: %ld numbers compared, %d failure(s)\n", failures ? "FAIL" : "PASS", compared, failures);
  return failures ? 1 : 0;
}
