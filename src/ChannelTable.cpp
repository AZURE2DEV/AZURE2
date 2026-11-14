#include "ChannelTable.h"

#include <algorithm>
#include <cmath>
#include <exception>

#include "CompoundSystem.h"
#include "ChannelFunc.h"
#include "Config.h"
#include "CoulFunc.h"
#include "Kinematics.h"
#include "PPair.h"

#ifdef _OPENMP
#include <omp.h>
#endif

ChannelTable::Row ChannelTable::Compute(const PPair *pair, int l, int two_j, double energy, const Config &configure) {
  const bool gsl = !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC);
  ChannelFunc channelFunc(pair, gsl);
  if (channelFunc.IsClosed(energy)) {
    // At and below threshold, and just above it where the Coulomb barrier makes P = 0 to double precision: the
    // shift function continued through E = 0 (ChannelFunc).  The Coulomb-wave routines are not usable there.
    return Row{0.0, channelFunc.Shift(l, two_j, energy), 1.0, 1.0};
  }
  CoulFunc coulomb(pair, gsl);
  const double radius = pair->GetChRad();
  Row row;
  row.P = coulomb.Penetrability(l, two_j, radius, energy);
  row.S = coulomb.PEShift(l, two_j, radius, energy);
  // exp(i(sigma_l - sigma_0)) = prod_{l'=1}^{l} (l' + i eta) / |l' + i eta|  (conventions.rst, conv-phases)
  const double eta = Kinematics::Channel(pair).Eta(energy);
  std::complex<double> expCP(1.0, 0.0);
  for (int ll = 1; ll <= l; ll++)
    expCP *= std::complex<double>((double)ll / sqrt(pow(eta, 2.0) + pow((double)ll, 2.0)),
                                  eta / sqrt(pow(eta, 2.0) + pow((double)ll, 2.0)));
  row.expCP = expCP;
  row.expHSP = coulomb.HardSpherePhase(l, two_j, radius, energy);
  return row;
}

namespace {

// Everything a row depends on besides (l, 2J, E_c), per pair: a previous table is reused only when these agree.
std::vector<std::vector<double>> PairSignatures(const CompoundSystem *compound, const Config &configure) {
  std::vector<std::vector<double>> s;
  for (int p = 1; p <= compound->NumPairs(); p++) {
    const PPair *pair = compound->GetPair(p);
    CoulFunc coulomb(pair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
    s.push_back({(double)pair->GetPType(), (double)pair->GetZ(1), (double)pair->GetZ(2), pair->GetM(1), pair->GetM(2),
                 pair->GetRedMass(), pair->GetChRad(), pair->GetSepE(), pair->GetExE(),
                 coulomb.UsesTable() ? 1.0 : 0.0, (configure.paramMask & Config::USE_GSL_COULOMB_FUNC) ? 1.0 : 0.0,
                 pair->FunctionSettings().relativistic ? 1.0 : 0.0});
    // the pair's hybrid potential, by value
    const std::vector<double> potential = azure2::physics::PotentialKey(pair->GetHybridPotential());
    s.back().insert(s.back().end(), potential.begin(), potential.end());
  }
  return s;
}

}  // namespace

std::shared_ptr<const ChannelTable> ChannelTable::Build(const CompoundSystem *compound, const Config &configure, std::vector<Key> keys,
                                                        const ChannelTable *previous) {
  auto table = std::make_shared<ChannelTable>();
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  table->keys_ = std::move(keys);
  table->rows_.resize(table->keys_.size());
  table->pairSignatures_ = PairSignatures(compound, configure);
  if (previous && previous->pairSignatures_ != table->pairSignatures_) previous = nullptr;

  std::vector<int> todo;
  for (int i = 0; i < (int)table->keys_.size(); i++) {
    const Key &k = table->keys_[i];
    const Row *old = previous ? previous->Find(k) : nullptr;
    if (old)
      table->rows_[i] = *old;
    else
      todo.push_back(i);
  }
  table->computed_ = todo.size();

  // Each row is written by the one thread that holds its index; a row is a pure function of its key, so the table
  // does not depend on the number of threads or the schedule.
  const int n = (int)todo.size();
  std::exception_ptr error;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) if (n > 32 && !omp_in_parallel())
#endif
  for (int t = 0; t < n; t++) {
    try {
      const Key &k = table->keys_[todo[t]];
      table->rows_[todo[t]] = Compute(compound->GetPair(k.pair), k.l, k.two_j, k.energy, configure);
    } catch (...) {
#ifdef _OPENMP
#pragma omp critical(channel_table_error)
#endif
      if (!error) error = std::current_exception();
    }
  }
  if (error) std::rethrow_exception(error);
  return table;
}

std::shared_ptr<const ChannelTable> ChannelTable::FromRows(const CompoundSystem *compound, const Config &configure,
                                                           std::vector<std::pair<Key, Row>> rows) {
  auto table = std::make_shared<ChannelTable>();
  std::sort(rows.begin(), rows.end(), [](const std::pair<Key, Row> &a, const std::pair<Key, Row> &b) { return a.first < b.first; });
  rows.erase(std::unique(rows.begin(), rows.end(), [](const std::pair<Key, Row> &a, const std::pair<Key, Row> &b) { return a.first == b.first; }),
             rows.end());
  table->keys_.reserve(rows.size());
  table->rows_.reserve(rows.size());
  for (const auto &kr : rows) {
    table->keys_.push_back(kr.first);
    table->rows_.push_back(kr.second);
  }
  table->pairSignatures_ = PairSignatures(compound, configure);
  return table;
}

bool ChannelTable::SameSystem(const CompoundSystem *compound, const Config &configure) const {
  return pairSignatures_ == PairSignatures(compound, configure);
}

const ChannelTable::Row *ChannelTable::Find(const Key &k) const {
  auto it = std::lower_bound(keys_.begin(), keys_.end(), k);
  if (it == keys_.end() || !(*it == k)) return nullptr;
  return &rows_[it - keys_.begin()];
}
