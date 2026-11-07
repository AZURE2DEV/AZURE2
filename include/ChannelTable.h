#ifndef CHANNELTABLE_H
#define CHANNELTABLE_H

#include <complex>
#include <cstddef>
#include <memory>
#include <vector>

class CNuc;
class PPair;
class Config;

/*!
 * The per-stage table of channel functions (docs/upgrade/channel_table.md): for each key (pair alpha, l, channel
 * energy E_c in MeV c.m.) the penetrability P, the shift function S, the hard-sphere phase factor exp(-i phi) and the
 * Coulomb phase factor exp(i(sigma_l - sigma_0)), all dimensionless.  E_c is compared exactly (no tolerance).  The
 * table is built once, before the threads that read it start, and is never written afterwards; a new set of
 * energies gets a new table (Build with the previous one), never a modified one.
 */
class ChannelTable {
 public:
  struct Key {
    int pair;       ///< the engine's pair number alpha (1-based)
    int l;          ///< orbital angular momentum
    double energy;  ///< E_c = (E_x - S_alpha) - e_alpha, MeV c.m.
    bool operator<(const Key &o) const {
      if (pair != o.pair) return pair < o.pair;
      if (l != o.l) return l < o.l;
      return energy < o.energy;
    }
    bool operator==(const Key &o) const { return pair == o.pair && l == o.l && energy == o.energy; }
  };
  struct Row {
    double P;                        ///< penetrability (0 for a closed channel)
    double S;                        ///< shift function
    std::complex<double> expHSP;     ///< exp(-i phi), 1 for a closed channel
    std::complex<double> expCP;      ///< exp(i(sigma_l - sigma_0)), 1 for a closed channel
  };

  /// The one formula of a row: an open channel from CoulFunc, a closed one (ChannelFunc::IsClosed) from
  /// ChannelFunc::Shift with P = 0.  Not for |E_c| <= ScatteringLength::kThresholdTolerance (EPoint's own case).
  static Row Compute(PPair *pair, int l, double energy, const Config &configure);

  /// The table of the given keys (any order, duplicates allowed).  Rows of keys held by `previous` are copied from
  /// it when it was built for the same pairs (charges, masses, radii, thresholds, Coulomb route); the others are
  /// computed in parallel, each row by one thread.  Errors of the Coulomb routines are rethrown here.
  static std::shared_ptr<const ChannelTable> Build(CNuc *cnuc, const Config &configure, std::vector<Key> keys,
                                                   const ChannelTable *previous = nullptr);

  /// The row of a key, or nullptr when the table does not hold it.
  const Row *Find(int pair, int l, double energy) const;

  std::size_t Size() const { return keys_.size(); }
  std::size_t Computed() const { return computed_; }  ///< rows computed (not copied) by Build

 private:
  std::vector<Key> keys_;  // sorted, unique
  std::vector<Row> rows_;
  std::vector<std::vector<double>> pairSignatures_;
  std::size_t computed_ = 0;
};

#endif
