#include "ChannelTable.h"

#include <algorithm>
#include <cmath>
#include <exception>

#include "CNuc.h"
#include "ChannelFunc.h"
#include "Config.h"
#include "CoulFunc.h"
#include "Kinematics.h"
#include "PPair.h"

#ifdef _OPENMP
#include <omp.h>
#endif

ChannelTable::Row ChannelTable::Compute(PPair *pair, int l, double energy, const Config &configure) {
  const bool gsl = !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC);
  ChannelFunc channelFunc(pair, gsl);
  if (channelFunc.IsClosed(energy)) {
    // At and below threshold, and just above it where the Coulomb barrier makes P = 0 to double precision: the
    // shift function continued through E = 0 (ChannelFunc).  The Coulomb-wave routines are not usable there.
    return Row{0.0, channelFunc.Shift(l, energy), 1.0, 1.0};
  }
  CoulFunc coulomb(pair, gsl);
  const double radius = pair->GetChRad();
  Row row;
  row.P = coulomb.Penetrability(l, radius, energy);
  row.S = coulomb.PEShift(l, radius, energy);
  // exp(i(sigma_l - sigma_0)) = prod_{l'=1}^{l} (l' + i eta) / |l' + i eta|  (conventions.rst, conv-phases)
  const double eta = Kinematics::Channel(pair).Eta(energy);
  std::complex<double> expCP(1.0, 0.0);
  for (int ll = 1; ll <= l; ll++)
    expCP *= std::complex<double>((double)ll / sqrt(pow(eta, 2.0) + pow((double)ll, 2.0)),
                                  eta / sqrt(pow(eta, 2.0) + pow((double)ll, 2.0)));
  row.expCP = expCP;
  row.expHSP = coulomb.HardSpherePhase(l, radius, energy);
  return row;
}

namespace {

// Everything a row depends on besides (l, E_c), per pair: a previous table is reused only when these agree.
std::vector<std::vector<double>> PairSignatures(CNuc *cnuc, const Config &configure) {
  std::vector<std::vector<double>> s;
  for (int p = 1; p <= cnuc->NumPairs(); p++) {
    PPair *pair = cnuc->GetPair(p);
    CoulFunc coulomb(pair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
    s.push_back({(double)pair->GetPType(), (double)pair->GetZ(1), (double)pair->GetZ(2), pair->GetM(1), pair->GetM(2),
                 pair->GetRedMass(), pair->GetChRad(), pair->GetSepE(), pair->GetExE(), (double)coulomb.hybridTag(),
                 coulomb.UsesTable() ? 1.0 : 0.0, (configure.paramMask & Config::USE_GSL_COULOMB_FUNC) ? 1.0 : 0.0,
                 Kinematics::Relativistic() ? 1.0 : 0.0});
  }
  return s;
}

}  // namespace

std::shared_ptr<const ChannelTable> ChannelTable::Build(CNuc *cnuc, const Config &configure, std::vector<Key> keys,
                                                        const ChannelTable *previous) {
  auto table = std::make_shared<ChannelTable>();
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  table->keys_ = std::move(keys);
  table->rows_.resize(table->keys_.size());
  table->pairSignatures_ = PairSignatures(cnuc, configure);
  if (previous && previous->pairSignatures_ != table->pairSignatures_) previous = nullptr;

  std::vector<int> todo;
  for (int i = 0; i < (int)table->keys_.size(); i++) {
    const Key &k = table->keys_[i];
    const Row *old = previous ? previous->Find(k.pair, k.l, k.energy) : nullptr;
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
      table->rows_[todo[t]] = Compute(cnuc->GetPair(k.pair), k.l, k.energy, configure);
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

std::shared_ptr<const ChannelTable> ChannelTable::FromRows(CNuc *cnuc, const Config &configure,
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
  table->pairSignatures_ = PairSignatures(cnuc, configure);
  return table;
}

bool ChannelTable::SameSystem(CNuc *cnuc, const Config &configure) const {
  return pairSignatures_ == PairSignatures(cnuc, configure);
}

const ChannelTable::Row *ChannelTable::Find(int pair, int l, double energy) const {
  const Key k{pair, l, energy};
  auto it = std::lower_bound(keys_.begin(), keys_.end(), k);
  if (it == keys_.end() || !(*it == k)) return nullptr;
  return &rows_[it - keys_.begin()];
}
