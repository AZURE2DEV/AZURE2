// The calculation plan (CalculationPlan.h; docs/upgrade/data_model.md step S5).
#include "CalculationPlan.h"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AMatrixFunc.h"
#include "CNuc.h"
#include "Config.h"
#include "NodeTable.h"
#include "Profile.h"
#include "TargetEffect.h"
#include "WorkingCopy.h"
#include "Unit.h"
#include "Segment.h"
#include "Results.h"

namespace azure2 {

namespace {

std::string Where(WorkingCopy &data, PointRef r) {
  std::ostringstream s;
  s << "unit " << r;
  if (r >= 0 && r < data.NumUnits())
    s << " (key " << data.UnitAt(r).GetSegmentKey() << ", E_lab " << data.UnitAt(r).GetLabEnergy() << " MeV)";
  return s.str();
}

}  // namespace

Unit *CalculationPlan::Point(WorkingCopy &data, int task) const { return data.PlanPoint(points_[task]); }

void CalculationPlan::Build(WorkingCopy &data) {
  points_ = data.PlanUnits();
  energy_tasks_ = EnergyTaskList(), energy_work_.clear(), energy_units_.clear(), energy_keys_.clear(), energy_sources_.clear();
  for (int t = 0; t < NumTasks(); t++) data.PlanPoint(points_[t])->SetPlanTask(t);
  const std::string why = Check(data);
  if (!why.empty()) throw std::logic_error("WorkingCopy: the calculation plan does not give every point one writer: " + why);
}

void CalculationPlan::AddRecordRange(const Unit &u, std::vector<std::pair<std::size_t, std::size_t>> &ranges) {
  if (u.ResultSlot() == Unit::kNoRecord) return;
  ranges.push_back({u.ResultSlot(), u.ResultSlot() + Results::RecordSize(u.ResultAngularSlots(), 0)});
}

void CalculationPlan::RunChecked(
    Results &results, int numSteps,
    const std::function<void(int, std::vector<std::pair<std::size_t, std::size_t>> &)> &run,
    const std::function<std::string(int)> &name) {
  // The slots are compared after every batch of steps: one step per batch up to 2,000 steps, more in a larger pass
  // (natural_neon has 56,703 energy tasks), so the check stays linear in the slots times 2,000.  Inside a batch a
  // changed slot is attributed to the step whose records hold it; one step's write into another's record of the same
  // batch is then not seen (the records of two steps never overlap: that is checked).
  const int batch = std::max(1, (numSteps + 1999) / 2000);
  std::vector<int> owner(results.slot.size(), -1);   // the step whose records hold the slot (latest such step)
  std::vector<int> writer(results.slot.size(), -1);  // the step that wrote it so far (-1: none)
  std::vector<double> before;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  for (int first = 0; first < numSteps; first += batch) {
    const int last = std::min(numSteps, first + batch);
    before = results.slot;
    for (int k = first; k < last; k++) {
      ranges.clear();
      run(k, ranges);
      for (const auto &r : ranges)
        for (std::size_t i = r.first; i < r.second && i < owner.size(); i++) {
          if (owner[i] >= first && owner[i] != k)
            throw std::logic_error("AZURE_PLAN_CHECK: the records of " + name(owner[i]) + " and " + name(k) +
                                   " overlap at slot " + std::to_string(i));
          owner[i] = k;
        }
    }
    const std::string steps = batch == 1 ? name(first) : "one of " + name(first) + " .. " + name(last - 1);
    for (std::size_t i = 0; i < results.slot.size() && i < before.size(); i++) {
      if (std::memcmp(&before[i], &results.slot[i], sizeof(double)) == 0) continue;
      const int k = owner[i];
      if (k < first) throw std::logic_error("AZURE_PLAN_CHECK: " + steps + " writes slot " + std::to_string(i) + ", outside its records");
      if (writer[i] >= 0 && writer[i] != k)
        throw std::logic_error("AZURE_PLAN_CHECK: " + name(k) + " writes slot " + std::to_string(i) + ", which " +
                               name(writer[i]) + " wrote");
      writer[i] = k;
    }
  }
}

std::string CalculationPlan::Check(WorkingCopy &data) const {
  // writer[unit]: the task that calculates the point (-1: none yet; -2: a unit that is no point of a term).
  const std::vector<PointRef> points = data.PlanUnits();
  std::vector<int> writer(data.NumUnits(), -2);
  for (PointRef u : points) writer[u] = -1;
  for (int t = 0; t < NumTasks(); t++) {
    const PointRef r = points_[t];
    if (r < 0 || r >= (int)writer.size() || writer[r] == -2)
      return "task " + std::to_string(t) + " names a unit that is no point of the data set";
    if (writer[r] >= 0)
      return Where(data, r) + " is calculated by two tasks (" + std::to_string(writer[r]) + " and " + std::to_string(t) + ")";
    writer[r] = t;
    if (data.PlanPoint(r)->PlanTask() != t)
      return Where(data, r) + " names task " + std::to_string(data.PlanPoint(r)->PlanTask()) + ", the plan " +
             std::to_string(t);
  }
  for (PointRef u : points)
    if (writer[u] < 0) return Where(data, u) + " is calculated by no task";
  return std::string();
}

bool HasSubPointFold(const WorkingCopy &data, const Unit &point) {
  return point.IsTargetEffect() && data.Data().GetTargetEffect(point.GetTargetEffectNum())->IsSubPointEffect();
}

void CalculationPlan::TaskUnits(WorkingCopy &data, int task, std::vector<UnitRef> &units) const {
  const PointRef point = points_[task];
  const Unit &p = data.UnitAt(point);
  if (!HasSubPointFold(data, p)) {
    units.push_back(point);
    return;
  }
  for (int i = 0; i < p.NumSubPoints(); i++) units.push_back(p.SubPointFirst() + i);
}

long long ECAmplitudeSource(const Unit &unit, const Config &configure) {
  if (!(configure.paramMask & Config::USE_EXTERNAL_CAPTURE)) return 0;
  if (configure.useECTable) {
    const std::uint64_t held = unit.ECAmplitudesHeldAt(unit.GetCMEnergy());
    return held ? (long long)(held | 0x8000000000000000ULL) : 1;  // a held hash is negative, never 1 or 2 + key
  }
  return 2 + ((long long)unit.GetSegmentKey() - (long long)INT32_MIN);
}

namespace {

std::uint64_t Bits(double x) {
  std::uint64_t b;
  std::memcpy(&b, &x, sizeof b);
  return b;
}

// Sorts `v` by `less` (a strict total order): pieces of it in parallel, then merged pairwise.
template <class T, class Less>
void SortInPieces(std::vector<T> &v, Less less) {
  const int n = (int)v.size();
  int pieces = 1;
#ifdef _OPENMP
  if (n > 4096) pieces = std::min(8, omp_get_max_threads());
#endif
  std::vector<int> bound(pieces + 1);
  for (int p = 0; p <= pieces; p++) bound[p] = (int)((long long)n * p / pieces);
#pragma omp parallel for schedule(static) if (pieces > 1)
  for (int p = 0; p < pieces; p++) std::sort(v.begin() + bound[p], v.begin() + bound[p + 1], less);
  for (int width = 1; width < pieces; width *= 2) {
#pragma omp parallel for schedule(static) if (pieces > 2 * width)
    for (int p = 0; p < pieces - width; p += 2 * width)
      std::inplace_merge(v.begin() + bound[p], v.begin() + bound[p + width],
                         v.begin() + bound[std::min(pieces, p + 2 * width)], less);
  }
}

char RouteOf(bool thm, const Config &configure) {
  if (thm) return 'T';
  return (configure.paramMask & Config::USE_AMATRIX) ? 'A' : 'R';
}

}  // namespace

const EnergyTaskList &CalculationPlan::EnergyTasks(WorkingCopy &data, const Config &configure,
                                                  const std::vector<int> &work) {
  // The units of the work in plan order, with their plan tasks, their nodes (the slots of the working copy: the node
  // tables give the nucleus, the route and E_x), and their amplitude sources now.
  const NodeTable *stage = data.Nodes();
  const NodeSlots *slots = data.Slots();
  if (!stage || !slots) throw std::logic_error("CalculationPlan::EnergyTasks: the working copy has no node table");
  const NodeOverflow &overflow = data.Overflow();
  const int numStage = stage->Size();
  std::vector<UnitRef> units;
  std::vector<int> owners, nodes;
  units.reserve(energy_units_.size()), owners.reserve(energy_units_.size()), nodes.reserve(energy_units_.size());
  for (int t : work) {
    TaskUnits(data, t, units);
    owners.resize(units.size(), t);
    nodes.insert(nodes.end(), slots->node.begin() + slots->begin[t], slots->node.begin() + slots->begin[t + 1]);
    if (nodes.size() != units.size())
      throw std::logic_error("CalculationPlan::EnergyTasks: the slots of plan task " + std::to_string(t) +
                             " do not match its units");
  }
  const int n = (int)units.size();
  auto nodeList = [&](int node) -> const NodeList & {
    return node < numStage ? static_cast<const NodeList &>(*stage) : overflow;
  };
  auto at = [numStage](int node) { return node < numStage ? node : node - numStage; };
  std::vector<double> keys(n);
  std::vector<long long> sources(n);
  // Each unit's key is its own (reads only), so they are formed in parallel: a fit that moves an energy shift
  // rebuilds the tasks at every evaluation, and this serial section is what the other threads wait for.
#pragma omp parallel for schedule(static) if (n > 4096)
  for (int k = 0; k < n; k++) {
    keys[k] = nodeList(nodes[k]).energy[at(nodes[k])];
    sources[k] = ECAmplitudeSource(data.UnitAt(units[k]), configure);
  }
  const bool same = work == energy_work_ && units == energy_units_ && sources == energy_sources_ &&
                    keys.size() == energy_keys_.size() &&
                    (keys.empty() || std::memcmp(keys.data(), energy_keys_.data(), keys.size() * sizeof(double)) == 0);
  if (same) return energy_tasks_;
  AZUREProfile::Count(AZUREProfile::kEnergyTaskBuilds);

  struct Key {
    int nucleus;
    char route;
    double energy;       // E_x, MeV
    std::uint64_t bits;  // of E_x
    int entrance, exit;
    long long source;
    int ordinal;
  };
  std::vector<Key> order(n);
#pragma omp parallel for schedule(static) if (n > 4096)
  for (int k = 0; k < n; k++) {
    const Unit &u = data.UnitAt(units[k]);
    const NodeList &list = nodeList(nodes[k]);
    const int i = at(nodes[k]);
    order[k] = {list.nucleus[i], RouteOf(list.thm[i] != 0, configure), keys[k], Bits(keys[k]), u.GetEntranceKey(),
                u.GetExitKey(), sources[k], k};
  }
  // By task key (nucleus, route, E_x, bits), then by T-matrix group (entrance, exit, amplitude source), then in plan
  // order.  The order is total (the ordinal last), so sorting pieces in parallel and merging them gives the one serial
  // result.
  auto before = [](const Key &a, const Key &b) {
    if (a.nucleus != b.nucleus) return a.nucleus < b.nucleus;
    if (a.route != b.route) return a.route < b.route;
    if (a.energy != b.energy) return a.energy < b.energy;
    if (a.bits != b.bits) return a.bits < b.bits;
    if (a.entrance != b.entrance) return a.entrance < b.entrance;
    if (a.exit != b.exit) return a.exit < b.exit;
    if (a.source != b.source) return a.source < b.source;
    return a.ordinal < b.ordinal;
  };
  SortInPieces(order, before);
  // The tasks in energy order, each a run of `order` (its first index in `start`).
  EnergyTaskList list;
  std::vector<int> start, group(n);
  int head = 0;  // the first index (in `order`) of the current T-matrix group
  for (int k = 0; k < n; k++) {
    const Key &o = order[k];
    if (k == 0 || o.nucleus != order[k - 1].nucleus || o.route != order[k - 1].route || o.bits != order[k - 1].bits) {
      list.tasks.emplace_back();
      EnergyTask &t = list.tasks.back();
      t.nucleus = o.nucleus, t.route = o.route, t.energy = o.energy, t.node = nodes[o.ordinal];
      start.push_back(k);
    }
    EnergyTask &task = list.tasks.back();
    if (task.size == 0 || o.entrance != order[k - 1].entrance || o.exit != order[k - 1].exit ||
        o.source != order[k - 1].source) {
      task.numTMatrices++;
      head = k;
    } else if (o.source < 0) {
      // Amplitudes held: the hash agrees with the group's first unit; so must the amplitudes, bit for bit.
      const Unit &self = data.UnitAt(units[o.ordinal]);
      if (!self.SameECAmplitudesAs(data.UnitAt(units[order[head].ordinal])))
        throw std::logic_error("CalculationPlan: two sets of external-capture amplitudes with one hash at E_x " +
                               std::to_string(o.energy) + " MeV");
    }
    group[k] = task.numTMatrices - 1;
    task.size++;
  }
  // Largest first (one fill and factorisation, one T-matrix per group, one observable per unit); ties in energy order.
  std::vector<int> rank(list.tasks.size());
  for (std::size_t t = 0; t < rank.size(); t++) rank[t] = (int)t;
  std::stable_sort(rank.begin(), rank.end(), [&list](int a, int b) {
    return 1 + list.tasks[a].numTMatrices + list.tasks[a].size > 1 + list.tasks[b].numTMatrices + list.tasks[b].size;
  });
  std::vector<EnergyTask> sorted(rank.size());
  list.units.resize(n), list.tGroups.resize(n), list.owners.resize(n);
  int next = 0;
  for (std::size_t r = 0; r < rank.size(); r++) {
    sorted[r] = list.tasks[rank[r]];
    sorted[r].begin = next;
    for (int i = 0; i < sorted[r].size; i++, next++) {
      const int k = start[rank[r]] + i, o = order[k].ordinal;
      list.units[next] = units[o], list.tGroups[next] = group[k];
      list.owners[next] = owners[o];
    }
  }
  list.tasks.swap(sorted);
  energy_tasks_ = std::move(list);
  energy_work_ = work;
  energy_units_.swap(units);
  energy_keys_.swap(keys);
  energy_sources_.swap(sources);
  return energy_tasks_;
}

bool CalculationPlan::RunTimeCheckRequested() {
  static const bool on = [] {
    const char *v = std::getenv("AZURE_PLAN_CHECK");
    return v && *v && std::strcmp(v, "0") != 0;
  }();
  return on;
}

}  // namespace azure2
