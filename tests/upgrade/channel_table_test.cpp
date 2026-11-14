// The per-stage channel table (src/ChannelTable.cpp, docs/upgrade/channel_table.md) on the charged and neutral pairs
// of T1, T2 and T3 and the spin-orbit pair of tests/spin_orbit_potential (keys with both j = l -+ 1/2): rows for
// l = 0..4 at energies below, at round-off from and above each threshold.
//  1. a table built from the keys in two orders, with duplicates, on 1 and on 8 threads is bit-identical, and each
//     row is bit-identical to ChannelTable::Compute at its key (rows depend on their key only);
//  2. Find has no row for an energy one ulp from a key (no tolerance);
//  3. a table built from a previous one equals a table built from scratch, and copies the rows the previous holds;
//  4. a previous table of a changed system (another channel radius) is not reused.
// It fails if a row is computed by another formula than Compute, if the build depends on the schedule, or if a key
// matches a neighbouring energy, as the deleted memo (CoulFuncCache, 1e-12 MeV) did.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "CNuc.h"
#include "ChannelTable.h"
#include "Config.h"
#include "PPair.h"

namespace {

int failures = 0;

void check(const std::string &name, bool ok) {
  std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", name.c_str());
  if (!ok) failures++;
}

bool SameBits(const ChannelTable::Row &a, const ChannelTable::Row &b) { return std::memcmp(&a, &b, sizeof a) == 0; }

bool SameTable(const ChannelTable &a, const ChannelTable &b, const std::vector<ChannelTable::Key> &keys) {
  if (a.Size() != b.Size()) return false;
  for (const auto &k : keys) {
    const ChannelTable::Row *x = a.Find(k), *y = b.Find(k);
    if (!x || !y || !SameBits(*x, *y)) return false;
  }
  return true;
}

void CheckProject(const std::string &azr) {
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = azr;
  CNuc compound;
  if (cfg.ReadConfigFile() < 0 || compound.Fill(cfg) == -1) {
    check(azr + ": load", false);
    return;
  }
  compound.Initialize(cfg);
  std::vector<ChannelTable::Key> keys;
  for (int p = 1; p <= compound.NumPairs(); p++) {
    if (compound.GetPair(p)->GetPType() != 0) continue;
    for (int l = 0; l <= 4; l++)
      for (double e : {-2.0, -0.3, -1e-6, 2e-10, 1e-4, 0.01, 0.1, 0.35, 1.0, 2.5, 7.0})
        for (double J : {l + 0.5, std::fabs(l - 0.5)})  // both j of a spin-1/2 particle: distinct keys for a spin-orbit pair
          keys.push_back({p, l, compound.GetPair(p)->SpinOrbitTwoJ(J), e});
  }
  std::printf("%s: %zu keys\n", azr.c_str(), keys.size());
  std::vector<ChannelTable::Key> shuffled = keys;
  shuffled.insert(shuffled.end(), keys.begin(), keys.begin() + keys.size() / 2);
  std::mt19937 rng(7);
  std::shuffle(shuffled.begin(), shuffled.end(), rng);
#ifdef _OPENMP
  omp_set_num_threads(1);
#endif
  auto one = ChannelTable::Build(&compound, cfg, keys);
#ifdef _OPENMP
  omp_set_num_threads(8);
#endif
  auto eight = ChannelTable::Build(&compound, cfg, shuffled);
  check("1 thread in order == 8 threads shuffled with duplicates (bits)", SameTable(*one, *eight, keys));
  bool computeSame = true;
  for (const auto &k : keys)
    computeSame = computeSame && SameBits(*one->Find(k),
                                          ChannelTable::Compute(compound.GetPair(k.pair), k.l, k.two_j, k.energy, cfg));
  check("every row == ChannelTable::Compute at its key (bits)", computeSame);
  bool noNeighbour = true;
  for (const auto &k : keys)
    noNeighbour = noNeighbour && !one->Find({k.pair, k.l, k.two_j, std::nextafter(k.energy, 1e9)}) &&
                  !one->Find({k.pair, k.l, k.two_j, std::nextafter(k.energy, -1e9)});
  check("no row for E_c one ulp from a key", noNeighbour);
  // a pair with a spin-orbit term: rows of l >= 1 differ between j = l - 1/2 and l + 1/2 above threshold, and no row
  // is found under the other j's key (the key holds 2J)
  for (int p = 1; p <= compound.NumPairs(); p++) {
    const PPair *pair = compound.GetPair(p);
    if (pair->SpinOrbitTwoJ(0.5) == azure2::physics::kNoSpinOrbit) continue;
    const ChannelTable::Row *a = one->Find({p, 1, 1, 1.0}), *b = one->Find({p, 1, 3, 1.0});
    check("spin-orbit pair " + std::to_string(p) + ": p1/2 and p3/2 rows differ", a && b && a->P != b->P && a->S != b->S);
    check("spin-orbit pair " + std::to_string(p) + ": no row under kNoSpinOrbit",
          !one->Find({p, 1, azure2::physics::kNoSpinOrbit, 1.0}));
  }
  std::vector<ChannelTable::Key> more = keys;
  for (const auto &k : keys) more.push_back({k.pair, k.l, k.two_j, k.energy * 1.5 + 0.01});
  auto fromPrevious = ChannelTable::Build(&compound, cfg, more, one.get());
  auto fromScratch = ChannelTable::Build(&compound, cfg, more);
  check("built from a previous table == built from scratch (bits)", SameTable(*fromPrevious, *fromScratch, more));
  check("built from a previous table computes only the new rows",
        fromPrevious->Computed() == fromPrevious->Size() - one->Size());
  // The same project with the first key's pair at another channel radius (the structure is fixed after set-up: a
  // radius is changed by filling the compound again with the override, as the API's SetRadius does).
  CNuc moved;
  moved.Fill(cfg, {{keys.front().pair, compound.GetPair(keys.front().pair)->GetChRad() + 0.5}});
  moved.Initialize(cfg);
  auto changed = ChannelTable::Build(&moved, cfg, keys, one.get());
  check("a previous table of another channel radius is not reused", changed->Computed() == changed->Size());
}

}  // namespace

int main(int argc, char **argv) {
  for (int i = 1; i < argc; i++) CheckProject(argv[i]);
  std::printf(failures ? "FAILED (%d)\n" : "PASSED\n", failures);
  return failures ? 1 : 0;
}
