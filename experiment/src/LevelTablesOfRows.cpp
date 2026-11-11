// The tables of a <levels> section made from its rows (ProjectReader.h, LevelTablesOfRows).
#include "azure2/experiment/ProjectReader.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>

#include "azure2/io/LevelGrouping.h"
#include "azure2/io/LevelRows.h"

namespace azure2 {
namespace experiment {

PairInfo PairInfoOfRow(const io::LevelRow &r) {
  PairInfo p;
  p.key = r.ir();
  p.kind = r.pType();
  if (p.kind != 10) {
    p.lightZ = r.z1();
    p.lightA = (int)std::lround(r.m1());
  }
  p.heavyZ = r.z2();
  p.heavyA = (int)std::lround(r.m2());
  p.lightExcitation = r.e1();
  p.heavyExcitation = r.e2();
  p.threshold = r.sepE() + (r.e1() + r.e2());  // PPair: GetSepE() + GetExE()
  p.lightMass = r.m1();
  p.heavyMass = r.m2();
  p.separationEnergy = r.sepE();
  p.radius = r.chRad();
  p.twoJLight = (int)std::lround(2.0 * r.j1());
  p.twoJHeavy = (int)std::lround(2.0 * r.j2());
  p.identical = io::IdenticalParticles(r);
  p.heavyParity = r.pi2();
  return p;
}

namespace {

// The compound nucleus (Z, A) a row's pair forms (PPair's constructor).
std::pair<int, int> CompoundOf(const io::LevelRow &r) {
  if (r.pType() == 20) return {r.z2() - r.z1(), (int)std::floor(r.m2() + 0.5)};
  if (r.pType() == 10) return {r.z2(), (int)std::floor(r.m2() + 0.5)};
  return {r.z1() + r.z2(), (int)std::floor(r.m1() + r.m2() + 0.5)};
}

struct LevelLines {
  std::map<int, std::pair<double, bool>> widths;  ///< group channel -> (width as written (eV), it is an amplitude)
};

}  // namespace

LevelTables LevelTablesOfRows(const std::vector<io::LevelRow> &rows) {
  LevelTables t;
  for (const io::LevelRow &r : rows)
    if (!PairByKey(r.ir(), t.pairs)) t.pairs.push_back(PairInfoOfRow(r));
  // External capture: the mask of a pair's first active row that has one, the pairs in the order of their first active
  // rows (CNuc::Fill numbers its pairs so, and fills its external-capture map from the same rows)
  std::vector<int> activeOrder;
  for (const io::LevelRow &r : rows) {
    if (r.isActive() != 1) continue;
    if (std::find(activeOrder.begin(), activeOrder.end(), r.ir()) == activeOrder.end()) activeOrder.push_back(r.ir());
    PairInfo &p = *std::find_if(t.pairs.begin(), t.pairs.end(), [&](const PairInfo &x) { return x.key == r.ir(); });
    if (p.externalCapture == 0 && r.ecMultMask() != 0) p.externalCapture = (unsigned int)r.ecMultMask();
  }
  for (int key : activeOrder)
    if (PairByKey(key, t.pairs)->externalCapture != 0) t.externalCapturePairs.push_back(key);

  // Groups (J, parity, compound nucleus) and their channels, in the order of the active rows, as CNuc::Fill.
  std::vector<std::tuple<double, int, std::pair<int, int>>> groups;
  std::vector<std::vector<LevelChannelInfo>> channels;
  std::vector<std::vector<LevelInfo>> levels;
  std::vector<std::vector<LevelLines>> lines;
  io::LevelGrouper grouper;
  for (const io::LevelRow &r : rows) {
    t.rowSlots.emplace_back();
    if (r.isActive() != 1) continue;
    const auto id = std::make_tuple(r.levelJ(), r.levelPi(), CompoundOf(r));
    int g = 0;
    while (g < (int)groups.size() && groups[g] != id) g++;
    if (g == (int)groups.size()) {
      groups.push_back(id);
      channels.emplace_back();
      levels.emplace_back();
      lines.emplace_back();
    }
    const LevelChannelInfo ch{r.ir(), r.l(), (int)std::lround(2.0 * r.s())};
    int c = 0;
    while (c < (int)channels[g].size() &&
           !(channels[g][c].pairKey == ch.pairKey && channels[g][c].l == ch.l && channels[g][c].twoS == ch.twoS))
      c++;
    if (c == (int)channels[g].size()) channels[g].push_back(ch);
    const io::LevelGrouper::Slot slot = grouper.Assign(g + 1, r.levelE(), r.levelID());
    if (slot.level > (int)levels[g].size()) {
      LevelInfo info;
      info.twoJ = (int)std::lround(2.0 * r.levelJ());
      info.parity = r.levelPi();
      info.energy = r.levelE();
      info.group = g + 1;
      info.level = slot.level;
      info.nucleusZ = std::get<2>(id).first;
      info.nucleusA = std::get<2>(id).second;
      levels[g].push_back(info);
      lines[g].emplace_back();
    }
    LevelInfo &level = levels[g][slot.level - 1];
    if (std::find(level.ids.begin(), level.ids.end(), r.levelID()) == level.ids.end()) level.ids.push_back(r.levelID());
    lines[g][slot.level - 1].widths[c + 1] = {r.gamma(), r.gammaIsRWA() == 1};
    t.rowSlots.back() = LevelRowSlot{g + 1, slot.level, c + 1};
  }
  for (std::size_t g = 0; g < groups.size(); g++)
    for (std::size_t la = 0; la < levels[g].size(); la++) {
      LevelInfo info = levels[g][la];
      info.channels = channels[g];
      // The total width as written (eV): open particle channels and photon channels given as widths; a static
      // moment (a photon channel to the level itself) is no width (src/ModelData.cpp reads it the same way).
      info.width = 0.0;
      for (const auto &w : lines[g][la].widths) {
        if (w.second.second) continue;
        const PairInfo &p = *PairByKey(channels[g][w.first - 1].pairKey, t.pairs);
        const bool open = p.kind == 0 && info.energy > p.threshold;
        const bool photon = p.kind == 10 && std::fabs(info.energy - p.heavyExcitation) >= 1e-3;
        if (open || photon) info.width += std::fabs(w.second.first) * 1e-6;
      }
      t.levels.push_back(info);
    }
  return t;
}

}  // namespace experiment
}  // namespace azure2
