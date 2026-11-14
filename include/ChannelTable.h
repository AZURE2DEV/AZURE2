#ifndef CHANNELTABLE_H
#define CHANNELTABLE_H

#include <complex>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

class CNuc;
class PPair;
class Config;

/*!
 * The per-stage table of channel functions (docs/upgrade/channel_table.md): for each key (pair alpha, l, 2J, channel
 * energy E_c in MeV c.m.; 2J = PPair::SpinOrbitTwoJ, kNoSpinOrbit for every channel of a pair without a spin-orbit term) the penetrability P, the shift function S, the hard-sphere phase factor exp(-i phi) and the
 * Coulomb phase factor exp(i(sigma_l - sigma_0)), all dimensionless.  E_c is compared exactly (no tolerance).  The
 * table is built once, before the threads that read it start, and is never written afterwards; a new set of
 * energies gets a new table (Build with the previous one), never a modified one.
 */
class ChannelTable {
 public:
  struct Key {
    int pair;       ///< the engine's pair number alpha (1-based)
    int l;          ///< orbital angular momentum
    int two_j;      ///< 2J of the channel (PPair::SpinOrbitTwoJ): kNoSpinOrbit unless the pair has a spin-orbit term
    double energy;  ///< E_c = (E_x - S_alpha) - e_alpha, MeV c.m.
    bool operator<(const Key &o) const {
      if (pair != o.pair) return pair < o.pair;
      if (l != o.l) return l < o.l;
      if (two_j != o.two_j) return two_j < o.two_j;
      return energy < o.energy;
    }
    bool operator==(const Key &o) const { return pair == o.pair && l == o.l && two_j == o.two_j && energy == o.energy; }
  };
  struct Row {
    double P;                        ///< penetrability (0 for a closed channel)
    double S;                        ///< shift function
    std::complex<double> expHSP;     ///< exp(-i phi), 1 for a closed channel
    std::complex<double> expCP;      ///< exp(i(sigma_l - sigma_0)), 1 for a closed channel
  };

  /// The one formula of a row: an open channel from CoulFunc, a closed one (ChannelFunc::IsClosed) from
  /// ChannelFunc::Shift with P = 0.  Not for |E_c| <= ScatteringLength::kThresholdTolerance (Unit's own case).
  static Row Compute(const PPair *pair, int l, int two_j, double energy, const Config &configure);

  /// The table of the given keys (any order, duplicates allowed).  Rows of keys held by `previous` are copied from
  /// it when it was built for the same pairs (charges, masses, radii, thresholds, Coulomb route); the others are
  /// computed in parallel, each row by one thread.  Errors of the Coulomb routines are rethrown here.
  static std::shared_ptr<const ChannelTable> Build(const CNuc *cnuc, const Config &configure, std::vector<Key> keys,
                                                   const ChannelTable *previous = nullptr);

  /// The table of rows already computed (a free-shift refill records the rows its points used): sorted, duplicates
  /// removed (a row is a pure function of its key, so duplicates are equal).
  static std::shared_ptr<const ChannelTable> FromRows(const CNuc *cnuc, const Config &configure,
                                                      std::vector<std::pair<Key, Row>> rows);

  /// True when the table was built for the present pairs (charges, masses, radii, thresholds, Coulomb route).
  bool SameSystem(const CNuc *cnuc, const Config &configure) const;

  /// The row of a key, or nullptr when the table does not hold it.
  const Row *Find(const Key &key) const;

  std::size_t Size() const { return keys_.size(); }
  std::size_t Computed() const { return computed_; }  ///< rows computed (not copied) by Build

 private:
  std::vector<Key> keys_;  // sorted, unique
  std::vector<Row> rows_;
  std::vector<std::vector<double>> pairSignatures_;
  std::size_t computed_ = 0;
};

/// Where a point takes its rows from (Unit::CalcEDependentValues): the tables, searched in order by exact key, then
/// ChannelTable::Compute once per (pair, l, 2J, E_c) of the point; `used`, when set, receives every row the point used.
struct ChannelRowSource {
  std::vector<const ChannelTable *> tables;
  std::vector<std::pair<ChannelTable::Key, ChannelTable::Row>> *used = nullptr;
  const ChannelTable::Row *Find(const ChannelTable::Key &key) const {
    for (const ChannelTable *t : tables)
      if (const ChannelTable::Row *r = t->Find(key)) return r;
    return nullptr;
  }
};

#endif
