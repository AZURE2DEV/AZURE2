// The analytic gradient over the energy tasks of the point pass (GradEnergyTasks.h): the forward pass, the
// cotangents of the units, one level-matrix adjoint per compound energy.
#include "GradEnergyTasks.h"

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

#include "CNuc.h"
#include "Config.h"
#include "JGroup.h"
#include "ModelPass.h"
#include "Profile.h"
#include "Results.h"
#include "Segment.h"
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

bool DenseAdjointRequested() {
  static const bool dense = [] {
    const char *e = std::getenv("AZURE_GRAD_DENSE");
    return e && *e && std::string(e) != "0";
  }();
  return dense;
}

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

// The points of the data segments but the THM ones, the forward pass over them, their models.
void Forward(CNuc *compound, WorkingCopy &data, Results &results, const Config &config, Pass &pass) {
  AZUREProfile::Timer timer(AZUREProfile::kGradForward);
  for (int s = 1; s <= data.Data().NumSegments(); s++) {
    const Segment *seg = data.Data().GetSegment(s);
    if (!seg || seg->IsTHM()) continue;
    pass.segments.push_back(seg);
    for (int p = 0; p < seg->NumPoints(); p++)
      if (Unit *pt = seg->GetPoint(data, p + 1)) pass.points.push_back({seg, s, p, pt});
  }
  // The records of every unit at the parameters filled now: the point pass, or (no OpenMP, a failed task) each
  // point calculated alone, as an evaluation does.
  if (data.CalculatePointsInParallel(pass.segments, compound, config, results) != 0)
    for (const GradPoint &p : pass.points)
      EvaluateModelPoint(data, *p.segment, p.pointIdx, compound, config, results, LeafValues::Calculate);
  const int n = (int)pass.points.size();
  pass.models.assign(n, 0.0);
#pragma omp parallel for schedule(dynamic, 16)
  for (int k = 0; k < n; k++) {
    const GradPoint &p = pass.points[k];
    pass.models[k] = EvaluateModelPoint(data, *p.segment, p.pointIdx, compound, config, results, LeafValues::Ready);
  }
}

// Each point's cotangents for fitBar[k] (0: none), then the units' cotangents ordered by unit (each unit's in point
// order).  False, with bailSeg, where a point is outside the analytic path (the first such point in data order).
bool Cotangents(CNuc *compound, WorkingCopy &data, Results &results, const Config &config,
                const std::vector<double> &fitBar, Pass &pass, int &bailSeg) {
  const int n = (int)pass.points.size();
  std::vector<std::vector<UnitCotangent>> perPoint(n);
  pass.params.assign(n, {});
  std::vector<char> failed(n, 0);
#pragma omp parallel
  {
    GradAccum scratch;
    scratch.Init(compound);
#pragma omp for schedule(dynamic, 4)
    for (int k = 0; k < n; k++) {
      if (fitBar[k] == 0.0) continue;
      scratch.Zero();
      if (!PointCotangents(data, pass.points[k], compound, config, results, fitBar[k], k, scratch, perPoint[k])) {
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
  pass.cotFirst.assign(data.NumUnits() + 1, 0);
  for (const auto &list : perPoint)
    for (const UnitCotangent &c : list) pass.cotFirst[c.unit + 1]++;
  for (int u = 0; u < data.NumUnits(); u++) pass.cotFirst[u + 1] += pass.cotFirst[u];
  pass.cot.resize(pass.cotFirst.back());
  std::vector<int> next(pass.cotFirst.begin(), pass.cotFirst.end() - 1);
  for (const auto &list : perPoint)
    for (const UnitCotangent &c : list) pass.cot[next[c.unit]++] = c;
  return true;
}

bool HasCotangents(const Pass &pass, const EnergyTaskList &list, const EnergyTask &task) {
  for (int i = task.begin; i < task.begin + task.size; i++)
    if (pass.cotFirst[list.units[i]] < pass.cotFirst[list.units[i] + 1]) return true;
  return false;
}

// The solver of a task on this thread, filled and factored at the task's first unit.
AMatrixFunc &SolveTask(CNuc *compound, const Config &config, Unit *first) {
  static thread_local AMatrixFunc amf(compound, config);
  amf.Reset(compound, config);
  AZUREProfile::Timer timer(AZUREProfile::kGradForward);
  amf.ClearMatrices();
  amf.FillMatrices(first);
  amf.InvertMatrices();
  return amf;
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
    int end = i;
    bool any = false;
    while (end < task.size && g[end] == g[i]) {
      any = any || pass.cotFirst[u[end]] < pass.cotFirst[u[end] + 1];
      end++;
    }
    if (any) {
      Unit *provider = &data.UnitAt(u[i]);
      {
        AZUREProfile::Timer timer(AZUREProfile::kGradForward);
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

// The scalar adjoint of one task into `acc`: per T-matrix group the T-matrix cotangents of all its units summed;
// the bilinear terms of the decays without external capture merged over the task and back-propagated once per
// J-group, a decay with external capture on its own.
bool ScalarTask(CNuc *compound, const Config &config, WorkingCopy &data, Results &results, const Pass &pass,
                const EnergyTaskList &list, const EnergyTask &task, const vector_matrix_r *shiftDeriv, GradAccum &acc) {
  if (!HasCotangents(pass, list, task)) return true;
  Unit *first = &data.UnitAt(list.units[task.begin]);
  AMatrixFunc &amf = SolveTask(compound, config, first);
  std::vector<std::vector<BilinearTerm>> merged(compound->NumJGroups()), terms;
  TMatrixCotangent bar;
  const bool dense = DenseAdjointRequested();
  auto group = [&](int begin, int end, Unit *provider) {
    TMatrixCotangent sum;
    for (int i = begin; i < end; i++) {
      const UnitRef u = list.units[task.begin + i];
      for (int c = pass.cotFirst[u]; c < pass.cotFirst[u + 1]; c++) {
        const UnitCotangent &k = pass.cot[c];
        if (!amf.PointTBar(&data.UnitAt(u), results, k.value, bar, k.xsComponent, k.ayCrossSection)) return false;
        AddCotangent(sum, bar);
      }
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
  AZUREProfile::Timer timer(AZUREProfile::kGradBackprop);
  for (int j = 1; j <= (int)merged.size(); j++)
    if (!merged[j - 1].empty()) amf.BilinearAdjoint(first, j, merged[j - 1], acc, shiftDeriv);
  return true;
}

// The rows of one task: for every data row with cotangents at this energy, its share of d(row)/d(level parameters)
// into rowAcc[slot] (rowsOf[slot] its row).  Per row the T-matrix cotangents of a group are summed, the bilinear
// terms merged over the task; then either one back-propagation per row, or, where it costs less, the derivatives of
// the task's bilinear forms by free parameter once (the basis: AMatrixFunc::LevelBasis or BilinearAdjoint with the
// cotangents 1 and i) and each row as a combination of them.
bool RowTask(CNuc *compound, const Config &config, WorkingCopy &data, Results &results, const Pass &pass,
             const EnergyTaskList &list, const EnergyTask &task, const vector_matrix_r *shiftDeriv,
             const AMatrixFunc::FreeMask *mask, std::vector<int> &slotOfRow, std::vector<int> &rowsOf,
             std::vector<GradAccum> &rowAcc, GradAccum &scratch) {
  rowsOf.clear();
  if (!HasCotangents(pass, list, task)) return true;
  Unit *first = &data.UnitAt(list.units[task.begin]);
  AMatrixFunc &amf = SolveTask(compound, config, first);
  const bool dense = DenseAdjointRequested();
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
        rowAcc.back().Init(compound);
        rowAcc.back().mask = mask;
      }
      rowAcc[rowsOf.size() - 1].Zero();
      rowTerms.emplace_back(compound->NumJGroups());
    }
    return slotOfRow[row];
  };
  auto group = [&](int begin, int end, Unit *provider) {
    groupSlots.clear();
    for (int i = begin; i < end; i++) {
      const UnitRef u = list.units[task.begin + i];
      for (int c = pass.cotFirst[u]; c < pass.cotFirst[u + 1]; c++) {
        const UnitCotangent &k = pass.cot[c];
        const int slot = slotOf(k.row);
        if (sums.size() < rowsOf.size()) sums.resize(rowsOf.size());
        if (std::find(groupSlots.begin(), groupSlots.end(), slot) == groupSlots.end()) {
          groupSlots.push_back(slot);
          sums[slot] = TMatrixCotangent();
        }
        if (!amf.PointTBar(&data.UnitAt(u), results, k.value, bar, k.xsComponent, k.ayCrossSection)) return false;
        AddCotangent(sums[slot], bar);
      }
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
      JGroup *jg = compound->GetJGroup(j);
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
  AZUREProfile::Timer timer(AZUREProfile::kGradBackprop);
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
    JGroup *jg = compound->GetJGroup(j);
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
  AZUREProfile::Timer rowsTimer(AZUREProfile::kGradRows);
  for (int q = 0; q < g; q++) {
    GradAccum &acc = rowAcc[q];
    size_t t = 0;
    for (int j = 1; j <= (int)rowTerms[q].size(); j++)
      for (const BilinearTerm &term : rowTerms[q][j - 1]) {
        const int k = keyOf[q][t++];
        const double a = std::real(term.beta), b = std::imag(term.beta);
        const std::vector<double> &re = basisRe[k], &im = basisIm[k];
        JGroup *jg = compound->GetJGroup(j);
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

bool EnergyTaskGradient(CNuc *compound, WorkingCopy *data, Results &results, const Config &config,
                        const vector_matrix_r *shiftDeriv, const FitBarFn &fitBarFn, GradAccum &accum, int &bailSeg) {
  bailSeg = 0;
  Pass pass;
  Forward(compound, *data, results, config, pass);
  // The objective's cotangent of each model, in data order (fitBarFn may accumulate the norm gradient).
  const int n = (int)pass.points.size();
  std::vector<double> fitBar(n);
  for (int k = 0; k < n; k++) {
    const GradPoint &p = pass.points[k];
    fitBar[k] = fitBarFn(p.segment, p.segNum, p.pointIdx, pass.models[k]);
  }
  if (!Cotangents(compound, *data, results, config, fitBar, pass, bailSeg)) return false;
  const EnergyTaskList &list = data->EnergyTasksOf(pass.segments, compound, config);
  const int chunks = std::max(1, std::min(kChunks, list.size()));
  std::vector<GradAccum> chunkAcc(chunks);
  std::vector<int> failedTask(chunks, -1);
#pragma omp parallel for schedule(dynamic, 1)
  for (int c = 0; c < chunks; c++) {
    chunkAcc[c].Init(compound);
    // A scalar pass: the Brune and Park level-matrix adjoints are summed in the formal parameters and resolved once
    // by the caller (AMatrixFunc::ResolveFormalCotangents).
    chunkAcc[c].formal = true;
    for (int t = c; t < list.size(); t += chunks)
      if (!ScalarTask(compound, config, *data, results, pass, list, list[t], shiftDeriv, chunkAcc[c])) {
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

bool EnergyTaskRows(CNuc *compound, WorkingCopy *data, Results &results, const Config &config,
                    const vector_matrix_r *shiftDeriv, const ParamIndexMap &pmap, const AMatrixFunc::FreeMask *mask,
                    const std::function<double(const GradPoint &, double)> &rowFitBar, std::vector<GradPoint> &points,
                    std::vector<double> &models, std::vector<double> &rows, int &bailSeg, int &numTasks) {
  bailSeg = 0;
  numTasks = 0;
  Pass pass;
  Forward(compound, *data, results, config, pass);
  const int n = (int)pass.points.size();
  std::vector<double> fitBar(n);
  for (int k = 0; k < n; k++) fitBar[k] = rowFitBar(pass.points[k], pass.models[k]);
  if (!Cotangents(compound, *data, results, config, fitBar, pass, bailSeg)) return false;
  const EnergyTaskList &list = data->EnergyTasksOf(pass.segments, compound, config);
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
      scratch.Init(compound);
      vector_r fullRow;
#pragma omp for schedule(dynamic, 1)
      for (int t = begin; t < end; t++) {
        if (!RowTask(compound, config, *data, results, pass, list, list[t], shiftDeriv, mask, slotOfRow, rowsOf,
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
