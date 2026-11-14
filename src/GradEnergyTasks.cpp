// The analytic gradient over the energy tasks of the point pass (GradEnergyTasks.h): the forward pass, the
// cotangents of the units, one level-matrix adjoint per compound energy.
#include "GradEnergyTasks.h"

#include <algorithm>
#include <cstdlib>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "CompoundSystem.h"
#include "Config.h"
#include "JGroup.h"
#include "ModelPass.h"
#include "Profile.h"
#include "Results.h"
#include "Segment.h"
#include "THMMatrixFunc.h"
#include "ThmCotangents.h"
#include "Unit.h"
#include "WorkingCopy.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace azure2 {

namespace {

// The scalar gradient sums the tasks in this many interleaved chunks (task k in chunk k mod kChunks), each in task
// order, then the chunks in order: the same additions in the same order at any thread count.
constexpr int kChunks = 64;
// The row passes hold at most this many doubles of task contributions at once (32 MB) before adding them to the rows.
constexpr long long kRowBlockDoubles = 4 << 20;

/// AZURE_GRAD_DENSE (the run's Config::environment): the dense level-matrix adjoint.
bool DenseAdjointRequested(const Config &config) { return config.environment.denseAdjoint; }

void AddCotangent(TMatrixCotangent &dst, const TMatrixCotangent &src) {
  if (dst.t.empty() && dst.ec.empty()) {
    dst = src;
    return;
  }
  for (size_t k = 0; k < src.t.size() && k < dst.t.size(); k++)
    for (size_t m = 0; m < src.t[k].size() && m < dst.t[k].size(); m++) dst.t[k][m] += src.t[k][m];
  for (size_t k = 0; k < src.ec.size() && k < dst.ec.size(); k++)
    for (size_t m = 0; m < src.ec[k].size() && m < dst.ec[k].size(); m++) dst.ec[k][m] += src.ec[k][m];
  dst.screen.insert(dst.screen.end(), src.screen.begin(), src.screen.end());
}

// Adds the bilinear terms `add` (per J-group) into `merged`, by channel pair.
void MergeTerms(std::vector<std::vector<BilinearTerm>> &merged, const std::vector<std::vector<BilinearTerm>> &add) {
  if (merged.size() < add.size()) merged.resize(add.size());
  for (std::size_t j = 0; j < add.size(); j++)
    for (const BilinearTerm &t : add[j]) {
      bool found = false;
      for (BilinearTerm &m : merged[j])
        if (m.ch == t.ch && m.chp == t.chp) {
          m.beta += t.beta;
          found = true;
          break;
        }
      if (!found) merged[j].push_back(t);
    }
}

// The pass: its points, their models, and the units' cotangents (by unit: cot[cotFirst[u] .. cotFirst[u + 1])).
struct Pass {
  std::vector<const Segment *> segments;
  std::vector<GradPoint> points;
  std::vector<double> models;
  std::vector<UnitCotangent> cot;
  std::vector<int> cotFirst;
  std::vector<std::vector<ParamCotangent>> params;  // per point
};

// The points of the data segments (the THM ones only with `hoes`, the context of their adjoint), the forward pass over
// them, their models.
void Forward(const CompoundSystem *compound, const LevelState &state, WorkingCopy &data, Results &results, const Config &config, Pass &pass,
             const ThmAdjointContext *hoes) {
  AZUREProfile::Timer timer(config.profile, AZUREProfile::kGradForward);
  for (int s = 1; s <= data.Data().NumSegments(); s++) {
    const Segment *seg = data.Data().GetSegment(s);
    if (!seg || (seg->IsTHM() && !hoes)) continue;
    pass.segments.push_back(seg);
    for (int p = 0; p < seg->NumPoints(); p++)
      if (Unit *pt = seg->GetPoint(data, p + 1)) pass.points.push_back({seg, s, p, pt});
  }
  // The records of every unit at the parameters filled now: the point pass, or (no OpenMP, a failed task) each
  // point calculated alone, as an evaluation does.
  if (data.CalculatePointsInParallel(pass.segments, compound, state, config, results) != 0)
    for (const GradPoint &p : pass.points)
      EvaluateModelPoint(data, *p.segment, p.pointIdx, compound, state, config, results, LeafValues::Calculate);
  const int n = (int)pass.points.size();
  pass.models.assign(n, 0.0);
#pragma omp parallel for schedule(dynamic, 16)
  for (int k = 0; k < n; k++) {
    const GradPoint &p = pass.points[k];
    pass.models[k] = EvaluateModelPoint(data, *p.segment, p.pointIdx, compound, state, config, results, LeafValues::Ready);
  }
}

// Each point's cotangents for fitBar[k] (0: none), then the units' cotangents ordered by unit (each unit's in point
// order).  False, with bailSeg, where a point is outside the analytic path (the first such point in data order).
bool Cotangents(const CompoundSystem *compound, const LevelState &state, WorkingCopy &data, Results &results, const Config &config,
                const std::vector<double> &fitBar, Pass &pass, int &bailSeg) {
  const int n = (int)pass.points.size();
  std::vector<std::vector<UnitCotangent>> perPoint(n);
  pass.params.assign(n, {});
  std::vector<char> failed(n, 0);
#pragma omp parallel
  {
    GradAccum scratch;
    scratch.Init(compound, state);
#pragma omp for schedule(dynamic, 4)
    for (int k = 0; k < n; k++) {
      if (fitBar[k] == 0.0) continue;
      scratch.Zero();
      if (!PointCotangents(data, pass.points[k], compound, state, config, results, fitBar[k], k, scratch, perPoint[k])) {
        failed[k] = 1;
        continue;
      }
      for (std::size_t i = 0; i < scratch.abundance.size(); i++)
        if (scratch.abundance[i] != 0.0) pass.params[k].push_back({'a', (int)i, scratch.abundance[i]});
      for (std::size_t i = 0; i < scratch.doppler.size(); i++)
        if (scratch.doppler[i] != 0.0) pass.params[k].push_back({'d', (int)i, scratch.doppler[i]});
      for (std::size_t i = 0; i < scratch.thickness.size(); i++)
        if (scratch.thickness[i] != 0.0) pass.params[k].push_back({'t', (int)i, scratch.thickness[i]});
    }
  }
  for (int k = 0; k < n; k++)
    if (failed[k]) {
      bailSeg = pass.points[k].segNum;
      return false;
    }
  pass.cotFirst.assign(data.NumSites() + 1, 0);
  for (const auto &list : perPoint)
    for (const UnitCotangent &c : list) pass.cotFirst[c.unit + 1]++;
  for (int u = 0; u < data.NumSites(); u++) pass.cotFirst[u + 1] += pass.cotFirst[u];
  pass.cot.resize(pass.cotFirst.back());
  std::vector<int> next(pass.cotFirst.begin(), pass.cotFirst.end() - 1);
  for (const auto &list : perPoint)
    for (const UnitCotangent &c : list) pass.cot[next[c.unit]++] = c;
  return true;
}

// Whether unit u reads a T-matrix (every kind but kHoes, whose HOES formula reads the level matrix's elements).
bool ReadsTMatrix(WorkingCopy &data, UnitRef u) { return data.SiteAt(u).Observable() != SiteObservable::kHoes; }

bool HasCotangents(const Pass &pass, const EnergyTaskList &list, const EnergyTask &task) {
  for (int i = task.begin; i < task.begin + task.size; i++)
    if (pass.cotFirst[list.units[i]] < pass.cotFirst[list.units[i] + 1]) return true;
  return false;
}

// The first unit of `task` that reads a T-matrix and has cotangents, or -1.
int FirstTMatrixReader(WorkingCopy &data, const Pass &pass, const EnergyTaskList &list, const EnergyTask &task) {
  for (int i = task.begin; i < task.begin + task.size; i++)
    if (ReadsTMatrix(data, list.units[i]) && pass.cotFirst[list.units[i]] < pass.cotFirst[list.units[i] + 1])
      return list.units[i];
  return -1;
}

// The point whose row reads unit u's observable for cotangent k: the reader of a slot holding Legendre coefficients
// (the slot has no angle), the unit itself otherwise.
const Site &RowsOf(WorkingCopy &data, UnitRef u, const UnitCotangent &k) {
  const Site &site = data.SiteAt(u);
  return (k.rows && site.AngleFree()) ? *k.rows : site;
}
int RowSlotOf(WorkingCopy &data, UnitRef u, const UnitCotangent &k) {
  return (k.rows && data.SiteAt(u).AngleFree()) ? k.rowSlot : 0;
}

// The cotangents of one T-matrix group with the reads of one slot through one row set merged (RowSet.h): the readers
// of a row set see the slot under the same row, so the slot's adjoint (AMatrixFunc::PointTBar, linear in the
// cotangent) is called once with their cotangents summed, in their order.  Other cotangents stay one each.
struct SharedRead {
  int cotangent;  // the first cotangent of the read (its reader's row, components)
  UnitRef unit;
  double value;
};

void MergeSharedReads(WorkingCopy &data, const Pass &pass, const std::vector<std::pair<int, UnitRef>> &order,
                      std::vector<SharedRead> &reads) {
  reads.clear();
  // (unit, row set, position in it, component, A_y cross section) -> its entry in `reads`
  std::map<std::tuple<UnitRef, const RowSet *, int, int, bool>, std::size_t> at;
  for (const std::pair<int, UnitRef> &o : order) {
    const UnitCotangent &k = pass.cot[o.first];
    const RowSet *set = (k.rows && data.SiteAt(o.second).AngleFree()) ? k.rows->RowSetOf() : nullptr;
    if (!set) {
      reads.push_back({o.first, o.second, k.value});
      continue;
    }
    const auto key = std::make_tuple(o.second, set, k.rows->RowSetIndex(k.rowSlot), k.xsComponent, k.ayCrossSection);
    const auto found = at.find(key);
    if (found == at.end()) {
      at.emplace(key, reads.size());
      reads.push_back({o.first, o.second, k.value});
    } else {
      reads[found->second].value += k.value;
    }
  }
}

// The solver of a task on this thread, filled and factored at the task's first unit.
AMatrixFunc &SolveTask(const CompoundSystem *compound, const LevelState &state, const Config &config, Results &results, Site *first) {
  AMatrixFunc &amf = results.scratch.ForThisThread().GradientSolver(compound, state, config);
  AZUREProfile::Timer timer(config.profile, AZUREProfile::kGradForward);
  amf.ClearMatrices();
  amf.FillMatrices(first);
  amf.InvertMatrices();
  return amf;
}

// The cotangents of the units [begin, end) of `task` (a T-matrix group), by their readers' plan tasks (stable: each
// reader's in unit and point order) -- the order of the sum with one slot per point (UnitCotangent::reader).
void GroupCotangents(const Pass &pass, const EnergyTaskList &list, const EnergyTask &task, int begin, int end,
                     std::vector<std::pair<int, UnitRef>> &order) {
  order.clear();
  for (int i = begin; i < end; i++) {
    const UnitRef u = list.units[task.begin + i];
    for (int c = pass.cotFirst[u]; c < pass.cotFirst[u + 1]; c++) order.push_back({c, u});
  }
  std::stable_sort(order.begin(), order.end(), [&pass](const std::pair<int, UnitRef> &a, const std::pair<int, UnitRef> &b) {
    return pass.cot[a.first].reader < pass.cot[b.first].reader;
  });
}

// Calls group(begin, end) for each T-matrix group of `task` with cotangents, its T-matrix built first (from its first
// unit, the group's provider; directly for the first built, through RecalculateTMatrix after, as the point pass).
template <class Group>
bool ForEachGroup(AMatrixFunc &amf, WorkingCopy &data, const Pass &pass, const EnergyTaskList &list,
                  const EnergyTask &task, Group group) {
  const UnitRef *u = list.units.data() + task.begin;
  const int *g = list.tGroups.data() + task.begin;
  bool built = false;
  for (int i = 0; i < task.size;) {
    int end = i, provide = -1;  // the group's first unit that reads a T-matrix
    bool any = false;
    while (end < task.size && g[end] == g[i]) {
      if (ReadsTMatrix(data, u[end])) {
        if (provide < 0) provide = end;
        any = any || pass.cotFirst[u[end]] < pass.cotFirst[u[end] + 1];
      }
      end++;
    }
    if (any) {
      Site *provider = &data.SiteAt(u[provide]);
      {
        AZUREProfile::Timer timer(amf.configure().profile, AZUREProfile::kGradForward);
        if (!built)
          amf.CalculateTMatrix(provider);
        else
          amf.RecalculateTMatrix(provider);
      }
      built = true;
      if (!group(i, end, provider)) return false;
    }
    i = end;
  }
  return true;
}

// The THM sites of one task (kHoes) into `acc`: per site the cotangents of all its readers summed (in their order), and
// the HOES adjoint taken once per site on a level-space fill at the task's energy (THMMatrixFunc::THMPointAdjoint).
bool HoesTask(const CompoundSystem *compound, const LevelState &state, const Config &config, WorkingCopy &data, Results &results, const Pass &pass,
              const EnergyTaskList &list, const EnergyTask &task, const ThmAdjointContext &context, GradAccum &acc) {
  THMMatrixFunc &thm = results.scratch.ForThisThread().GradientThmSolver(compound, state, config);
  bool filled = false;
  for (int i = task.begin; i < task.begin + task.size; i++) {
    const UnitRef u = list.units[i];
    if (ReadsTMatrix(data, u) || pass.cotFirst[u] == pass.cotFirst[u + 1]) continue;
    double sum = 0.0;
    for (int c = pass.cotFirst[u]; c < pass.cotFirst[u + 1]; c++) sum += pass.cot[c].value;
    Site *site = &data.SiteAt(u);
    if (!filled) {
      AZUREProfile::Timer timer(config.profile, AZUREProfile::kGradForward);
      thm.Reset(compound, state, config);
      thm.ClearMatrices();
      thm.FillMatrices(site);
      thm.InvertMatrices();
      filled = true;
    }
    AZUREProfile::Timer timer(config.profile, AZUREProfile::kGradBackprop);
    if (!thm.THMPointAdjoint(site, sum, acc, context)) return false;
  }
  return true;
}

// The scalar adjoint of one task into `acc`: per T-matrix group the T-matrix cotangents of all its units summed;
// the bilinear terms of the decays without external capture merged over the task and back-propagated once per
// J-group, a decay with external capture on its own.
bool ScalarTask(const CompoundSystem *compound, const LevelState &state, const Config &config, WorkingCopy &data, Results &results, const Pass &pass,
                const EnergyTaskList &list, const EnergyTask &task, const vector_matrix_r *shiftDeriv,
                const ThmAdjointContext *hoes, GradAccum &acc) {
  if (!HasCotangents(pass, list, task)) return true;
  if (hoes && !HoesTask(compound, state, config, data, results, pass, list, task, *hoes, acc)) return false;
  if (FirstTMatrixReader(data, pass, list, task) < 0) return true;
  // filled with the task's first unit that reads a T-matrix, as the point pass (SolveAtOneEnergy)
  int firstReader = task.begin;
  while (!ReadsTMatrix(data, list.units[firstReader])) firstReader++;
  Site *first = &data.SiteAt(list.units[firstReader]);
  AMatrixFunc &amf = SolveTask(compound, state, config, results, first);
  std::vector<std::vector<BilinearTerm>> merged(compound->NumJGroups()), terms;
  TMatrixCotangent bar;
  const bool dense = DenseAdjointRequested(config);
  std::vector<std::pair<int, UnitRef>> order;
  std::vector<SharedRead> reads;
  auto group = [&](int begin, int end, Site *provider) {
    TMatrixCotangent sum;
    GroupCotangents(pass, list, task, begin, end, order);
    MergeSharedReads(data, pass, order, reads);
    for (const SharedRead &r : reads) {
      const UnitCotangent &k = pass.cot[r.cotangent];
      if (!ReadsTMatrix(data, r.unit)) continue;  // a THM site: HoesTask
      if (!amf.PointTBar(&data.SiteAt(r.unit), results, r.value, bar, k.xsComponent, k.ayCrossSection, RowsOf(data, r.unit, k), RowSlotOf(data, r.unit, k))) return false;
      AddCotangent(sum, bar);
    }
    if (dense || !amf.CollectBilinearTerms(provider, sum, terms)) {
      amf.TBarAdjoint(provider, sum, acc, shiftDeriv);
      return true;
    }
    for (const auto &s : sum.screen) acc.AddScreening(s.first, s.second);
    MergeTerms(merged, terms);
    return true;
  };
  if (!ForEachGroup(amf, data, pass, list, task, group)) return false;
  AZUREProfile::Timer timer(config.profile, AZUREProfile::kGradBackprop);
  for (int j = 1; j <= (int)merged.size(); j++)
    if (!merged[j - 1].empty()) amf.BilinearAdjoint(first, j, merged[j - 1], acc, shiftDeriv);
  return true;
}

// The rows of one task: for every data row with cotangents at this energy, its share of d(row)/d(level parameters)
// into rowAcc[slot] (rowsOf[slot] its row).  Per row the T-matrix cotangents of a group are summed, the bilinear
// terms merged over the task; then either one back-propagation per row, or, where it costs less, the derivatives of
// the task's bilinear forms by free parameter once (the basis: AMatrixFunc::LevelBasis or BilinearAdjoint with the
// cotangents 1 and i) and each row as a combination of them.
bool RowTask(const CompoundSystem *compound, const LevelState &state, const Config &config, WorkingCopy &data, Results &results, const Pass &pass,
             const EnergyTaskList &list, const EnergyTask &task, const vector_matrix_r *shiftDeriv,
             const AMatrixFunc::FreeMask *mask, std::vector<int> &slotOfRow, std::vector<int> &rowsOf,
             std::vector<GradAccum> &rowAcc, GradAccum &scratch) {
  rowsOf.clear();
  if (!HasCotangents(pass, list, task)) return true;
  Site *first = &data.SiteAt(list.units[task.begin]);
  AMatrixFunc &amf = SolveTask(compound, state, config, results, first);
  const bool dense = DenseAdjointRequested(config);
  std::vector<std::vector<std::vector<BilinearTerm>>> rowTerms;
  std::vector<TMatrixCotangent> sums;
  std::vector<int> groupSlots;
  std::vector<std::vector<BilinearTerm>> terms;
  TMatrixCotangent bar;
  auto slotOf = [&](int row) {
    if (slotOfRow[row] < 0) {
      slotOfRow[row] = (int)rowsOf.size();
      rowsOf.push_back(row);
      if (rowAcc.size() < rowsOf.size()) {
        rowAcc.emplace_back();
        rowAcc.back().Init(compound, state);
        rowAcc.back().mask = mask;
      }
      rowAcc[rowsOf.size() - 1].Zero();
      rowTerms.emplace_back(compound->NumJGroups());
    }
    return slotOfRow[row];
  };
  std::vector<std::pair<int, UnitRef>> order;
  auto group = [&](int begin, int end, Site *provider) {
    groupSlots.clear();
    GroupCotangents(pass, list, task, begin, end, order);
    for (const std::pair<int, UnitRef> &o : order) {
      const UnitCotangent &k = pass.cot[o.first];
      const int slot = slotOf(k.row);
      if (sums.size() < rowsOf.size()) sums.resize(rowsOf.size());
      if (std::find(groupSlots.begin(), groupSlots.end(), slot) == groupSlots.end()) {
        groupSlots.push_back(slot);
        sums[slot] = TMatrixCotangent();
      }
      if (!amf.PointTBar(&data.SiteAt(o.second), results, k.value, bar, k.xsComponent, k.ayCrossSection, RowsOf(data, o.second, k), RowSlotOf(data, o.second, k))) return false;
      AddCotangent(sums[slot], bar);
    }
    for (int slot : groupSlots) {
      if (dense || !amf.CollectBilinearTerms(provider, sums[slot], terms)) {
        amf.TBarAdjoint(provider, sums[slot], rowAcc[slot], shiftDeriv);
        continue;
      }
      for (const auto &s : sums[slot].screen) rowAcc[slot].AddScreening(s.first, s.second);
      MergeTerms(rowTerms[slot], terms);
    }
    return true;
  };
  const bool ok = ForEachGroup(amf, data, pass, list, task, group);
  for (int row : rowsOf) slotOfRow[row] = -1;
  if (!ok) return false;

  // Cost of the two routes, in units of level-matrix elements (as for the rows of one energy before).
  const int g = (int)rowsOf.size();
  struct Key {
    int j, ch, chp;
  };
  std::vector<Key> keys;
  std::vector<std::vector<int>> keyOf(g);
  double costMembers = 0.0, costBasis = 0.0;
  std::vector<double> freeCount;
  if (mask) {
    freeCount.assign(mask->e.size(), 0.0);
    for (size_t j = 0; j < mask->e.size(); j++)
      for (size_t la = 0; la < mask->e[j].size(); la++) {
        freeCount[j] += mask->e[j][la];
        for (char f : mask->g[j][la]) freeCount[j] += f;
      }
  }
  for (int q = 0; q < g; q++)
    for (int j = 1; j <= (int)rowTerms[q].size(); j++) {
      if (rowTerms[q][j - 1].empty()) continue;
      const JGroup *jg = compound->GetJGroup(j);
      const double n2 = (double)jg->NumLevels() * jg->NumLevels();
      costMembers += n2 * (rowTerms[q][j - 1].size() + jg->NumChannels());
      for (const BilinearTerm &t : rowTerms[q][j - 1]) {
        int k = 0;
        while (k < (int)keys.size() && !(keys[k].j == j && keys[k].ch == t.ch && keys[k].chp == t.chp)) k++;
        if (k == (int)keys.size()) {
          keys.push_back({j, t.ch, t.chp});
          costBasis += (mask ? 2.0 * freeCount[j - 1] * jg->NumLevels() : 2.0 * n2 * (1.0 + jg->NumChannels()));
        }
        keyOf[q].push_back(k);
        costBasis += 2.0 * jg->NumLevels() * (1.0 + jg->NumChannels());
      }
    }
  AZUREProfile::Timer timer(config.profile, AZUREProfile::kGradBackprop);
  if (!(costBasis < costMembers)) {
    for (int q = 0; q < g; q++)
      for (int j = 1; j <= (int)rowTerms[q].size(); j++)
        if (!rowTerms[q][j - 1].empty()) amf.BilinearAdjoint(first, j, rowTerms[q][j - 1], rowAcc[q], shiftDeriv);
    return true;
  }
  // The basis: for each channel pair, the J-group block (E then gamma, level by level) of the adjoint with cotangent
  // 1 and with cotangent i.
  std::vector<std::vector<double>> basisRe(keys.size()), basisIm(keys.size());
  for (size_t k = 0; k < keys.size(); k++) {
    const int j = keys[k].j;
    const JGroup *jg = compound->GetJGroup(j);
    const int n = jg->NumLevels(), nch = jg->NumChannels();
    for (int part = 0; part < 2; part++) {
      std::vector<double> &block = part ? basisIm[k] : basisRe[k];
      const BilinearTerm term{keys[k].ch, keys[k].chp, part ? complex(0.0, 1.0) : complex(1.0, 0.0)};
      if (mask && amf.LevelBasis(first, j, term, *mask, shiftDeriv, block)) continue;
      for (int la = 0; la < n; la++) {
        scratch.e[j - 1][la] = 0.0;
        std::fill(scratch.gamma[j - 1][la].begin(), scratch.gamma[j - 1][la].end(), 0.0);
      }
      amf.BilinearAdjoint(first, j, std::vector<BilinearTerm>{term}, scratch, shiftDeriv);
      block.assign((size_t)n * (1 + nch), 0.0);
      for (int la = 0; la < n; la++) {
        block[la] = scratch.e[j - 1][la];
        for (int d = 0; d < nch; d++) block[n + (size_t)la * nch + d] = scratch.gamma[j - 1][la][d];
      }
    }
  }
  AZUREProfile::Timer rowsTimer(config.profile, AZUREProfile::kGradRows);
  for (int q = 0; q < g; q++) {
    GradAccum &acc = rowAcc[q];
    size_t t = 0;
    for (int j = 1; j <= (int)rowTerms[q].size(); j++)
      for (const BilinearTerm &term : rowTerms[q][j - 1]) {
        const int k = keyOf[q][t++];
        const double a = std::real(term.beta), b = std::imag(term.beta);
        const std::vector<double> &re = basisRe[k], &im = basisIm[k];
        const JGroup *jg = compound->GetJGroup(j);
        const int n = jg->NumLevels(), nch = jg->NumChannels();
        for (int la = 0; la < n; la++) {
          acc.e[j - 1][la] += a * re[la] + b * im[la];
          double *gr = &acc.gamma[j - 1][la][0];
          const double *r0 = &re[n + (size_t)la * nch], *i0 = &im[n + (size_t)la * nch];
          for (int d = 0; d < nch; d++) gr[d] += a * r0[d] + b * i0[d];
        }
      }
  }
  return true;
}

}  // namespace

bool EnergyTaskGradient(const CompoundSystem *compound, const LevelState &state, WorkingCopy *data, Results &results, const Config &config,
                        const vector_matrix_r *shiftDeriv, const FitBarFn &fitBarFn, GradAccum &accum, int &bailSeg,
                        const ThmAdjointContext *hoes) {
  bailSeg = 0;
  Pass pass;
  Forward(compound, state, *data, results, config, pass, hoes);
  // The objective's cotangent of each model, in data order (fitBarFn may accumulate the norm gradient); a THM point's
  // through its profile (ThmCotangents.h).
  const int n = (int)pass.points.size();
  std::vector<double> fitBar(n);
  for (int k = 0; k < n; k++) {
    const GradPoint &p = pass.points[k];
    fitBar[k] = fitBarFn(p.segment, p.segNum, p.pointIdx, pass.models[k]);
  }
  if (hoes) ThmModelCotangents(*data, results, pass.points, pass.models, fitBar);
  if (!Cotangents(compound, state, *data, results, config, fitBar, pass, bailSeg)) return false;
  const EnergyTaskList &list = data->EnergyTasksOf(pass.segments, config);
  const int chunks = std::max(1, std::min(kChunks, list.size()));
  std::vector<GradAccum> chunkAcc(chunks);
  std::vector<int> failedTask(chunks, -1);
#pragma omp parallel for schedule(dynamic, 1)
  for (int c = 0; c < chunks; c++) {
    chunkAcc[c].Init(compound, state);
    // A scalar pass: the Brune and Park level-matrix adjoints are summed in the formal parameters and resolved once
    // by the caller (AMatrixFunc::ResolveFormalCotangents).
    chunkAcc[c].formal = true;
    for (int t = c; t < list.size(); t += chunks)
      if (!ScalarTask(compound, state, config, *data, results, pass, list, list[t], shiftDeriv, hoes, chunkAcc[c])) {
        failedTask[c] = t;
        break;
      }
  }
  for (int c = 0; c < chunks; c++)
    if (failedTask[c] >= 0) {
      const EnergyTask &task = list[failedTask[c]];
      for (int i = task.begin; i < task.begin + task.size; i++)
        if (pass.cotFirst[list.units[i]] < pass.cotFirst[list.units[i] + 1]) {
          bailSeg = pass.points[pass.cot[pass.cotFirst[list.units[i]]].row].segNum;
          break;
        }
      return false;
    }
  for (const GradAccum &c : chunkAcc) accum.Add(c);
  for (const auto &list : pass.params)
    for (const ParamCotangent &q : list) {
      if (q.kind == 'a') accum.AddAbundance(q.index + 1, q.value);
      else if (q.kind == 'd') accum.AddDoppler(q.index + 1, q.value);
      else accum.AddThickness(q.index + 1, q.value);
    }
  return true;
}

bool EnergyTaskRows(const CompoundSystem *compound, const LevelState &state, WorkingCopy *data, Results &results, const Config &config,
                    const vector_matrix_r *shiftDeriv, const ParamIndexMap &pmap, const AMatrixFunc::FreeMask *mask,
                    const std::function<double(const GradPoint &, double)> &rowFitBar, std::vector<GradPoint> &points,
                    std::vector<double> &models, std::vector<double> &rows, int &bailSeg, int &numTasks) {
  bailSeg = 0;
  numTasks = 0;
  Pass pass;
  Forward(compound, state, *data, results, config, pass, nullptr);
  const int n = (int)pass.points.size();
  std::vector<double> fitBar(n);
  for (int k = 0; k < n; k++) fitBar[k] = rowFitBar(pass.points[k], pass.models[k]);
  if (!Cotangents(compound, state, *data, results, config, fitBar, pass, bailSeg)) return false;
  const EnergyTaskList &list = data->EnergyTasksOf(pass.segments, config);
  numTasks = list.size();
  const int nCols = pmap.NumPacked();
  rows.assign((size_t)n * nCols, 0.0);
  if (nCols == 0) {  // no free parameter: the models only
    points.swap(pass.points);
    models.swap(pass.models);
    return true;
  }

  // The rows each task touches (for the size of its contributions).
  std::vector<int> rowsOfTask(list.size(), 0), seen(n, -1);
  for (int t = 0; t < list.size(); t++)
    for (int i = list[t].begin; i < list[t].begin + list[t].size; i++)
      for (int c = pass.cotFirst[list.units[i]]; c < pass.cotFirst[list.units[i] + 1]; c++)
        if (seen[pass.cot[c].row] != t) seen[pass.cot[c].row] = t, rowsOfTask[t]++;

  // Blocks of tasks whose contributions fit kRowBlockDoubles: the tasks of a block in parallel, each into its own
  // place, then added to the rows in task order.
  std::vector<double> buffer;
  std::vector<int> bufferRows;
  bool ok = true;
  for (int begin = 0; begin < list.size() && ok;) {
    int end = begin;
    long long size = 0;
    std::vector<long long> offset;
    while (end < list.size() && (end == begin || size + (long long)rowsOfTask[end] * nCols <= kRowBlockDoubles)) {
      offset.push_back(size / nCols);
      size += (long long)rowsOfTask[end] * nCols;
      end++;
    }
    buffer.assign((size_t)size, 0.0);
    bufferRows.assign((size_t)(size / nCols), -1);
    std::vector<char> failed(end - begin, 0);
#pragma omp parallel
    {
      std::vector<int> slotOfRow(n, -1), rowsOf;
      std::vector<GradAccum> rowAcc;
      GradAccum scratch;
      scratch.Init(compound, state);
      vector_r fullRow;
#pragma omp for schedule(dynamic, 1)
      for (int t = begin; t < end; t++) {
        if (!RowTask(compound, state, config, *data, results, pass, list, list[t], shiftDeriv, mask, slotOfRow, rowsOf,
                     rowAcc, scratch)) {
          failed[t - begin] = 1;
          continue;
        }
        for (std::size_t q = 0; q < rowsOf.size(); q++) {
          const long long r = offset[t - begin] + (long long)q;
          bufferRows[r] = rowsOf[q];
          fullRow.assign(pmap.NumFull(), 0.0);
          rowAcc[q].Scatter(pmap, fullRow);
          double *out = buffer.data() + r * nCols;
          for (int f = 0; f < pmap.NumFull(); f++) {
            const int packed = pmap.FullToPacked(f);
            if (packed >= 0 && packed < nCols) out[packed] = fullRow[f];
          }
        }
      }
    }
    for (int t = begin; t < end && ok; t++)
      if (failed[t - begin]) {
        ok = false;
        const EnergyTask &task = list[t];
        for (int i = task.begin; i < task.begin + task.size && !bailSeg; i++)
          if (pass.cotFirst[list.units[i]] < pass.cotFirst[list.units[i] + 1])
            bailSeg = pass.points[pass.cot[pass.cotFirst[list.units[i]]].row].segNum;
      }
    for (std::size_t r = 0; r < bufferRows.size() && ok; r++) {
      if (bufferRows[r] < 0) continue;
      double *row = rows.data() + (size_t)bufferRows[r] * nCols;
      const double *add = buffer.data() + r * nCols;
      for (int c = 0; c < nCols; c++) row[c] += add[c];
    }
    begin = end;
  }
  if (!ok) return false;
  // The parameter parts of each row (abundances, Doppler temperatures, thicknesses).
  for (int k = 0; k < n; k++)
    for (const ParamCotangent &q : pass.params[k]) {
      const int full = q.kind == 'a' ? pmap.AbundanceIndex(q.index + 1)
                       : q.kind == 'd' ? pmap.DopplerIndex(q.index + 1) : pmap.ThicknessIndex(q.index + 1);
      if (full < 0 || full >= pmap.NumFull()) continue;
      const int packed = pmap.FullToPacked(full);
      if (packed >= 0 && packed < nCols) rows[(size_t)k * nCols + packed] += q.value;
    }
  points.swap(pass.points);
  models.swap(pass.models);
  return true;
}

}  // namespace azure2
