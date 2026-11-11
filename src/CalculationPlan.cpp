// The calculation plan (CalculationPlan.h; docs/upgrade/data_model.md step S5).
#include "CalculationPlan.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdexcept>

#include "AMatrixFunc.h"
#include "CNuc.h"
#include "Config.h"
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

// The slots a point and its sub-points hold apart from the natural terms (the model pass's).
void AddUnitRanges(WorkingCopy &data, Unit &p, std::vector<std::pair<std::size_t, std::size_t>> &ranges) {
  CalculationPlan::AddRecordRange(p, ranges);
  for (Unit &sub : p.GetSubPoints(data)) CalculationPlan::AddRecordRange(sub, ranges);
}

}  // namespace

Unit *MappedPoints::operator[](int i) const { return data_->PlanPoint((*refs_)[i]); }

Unit *CalculationPlan::Host(WorkingCopy &data, int task) const { return data.PlanPoint(tasks_[task].host); }

void CalculationPlan::Build(WorkingCopy &data, const std::vector<std::pair<PointRef, PointRef>> &mapping) {
  tasks_.clear();
  energy_tasks_.clear(), energy_work_.clear(), energy_units_.clear(), energy_keys_.clear();
  for (PointRef u : data.PlanUnits()) {
    Unit *point = data.PlanPoint(u);
    point->SetPlanTask(-1);
    if (point->IsMapped()) continue;
    tasks_.push_back({u, {}});
    point->SetPlanTask((int)tasks_.size() - 1);
  }
  for (const std::pair<PointRef, PointRef> &m : mapping) {
    const int task = data.PlanPoint(m.first)->PlanTask();
    if (task < 0) throw std::logic_error("WorkingCopy: the calculation plan maps " + Where(data, m.second) + " onto " +
                                         Where(data, m.first) + ", which is not the host of a task");
    tasks_[task].mapped.push_back(m.second);
    data.PlanPoint(m.second)->SetPlanTask(task);
  }
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
  auto claim = [&](PointRef r, int task, std::string &why) {
    if (r < 0 || r >= (int)writer.size() || writer[r] == -2) {
      why = "task " + std::to_string(task) + " names a unit that is no point of the data set";
      return false;
    }
    int &w = writer[r];
    if (w >= 0) {
      why = Where(data, r) + " is calculated by two tasks (" + std::to_string(w) + " and " + std::to_string(task) + ")";
      return false;
    }
    w = task;
    if (data.PlanPoint(r)->PlanTask() != task) {
      why = Where(data, r) + " names task " + std::to_string(data.PlanPoint(r)->PlanTask()) + ", the plan " +
            std::to_string(task);
      return false;
    }
    return true;
  };
  std::string why;
  for (int t = 0; t < NumTasks(); t++) {
    const CalculationTask &task = tasks_[t];
    if (!claim(task.host, t, why)) return why;
    Unit *host = data.PlanPoint(task.host);
    if (host->IsMapped()) return Where(data, task.host) + " is a host and mapped onto another point";
    for (PointRef m : task.mapped) {
      if (!claim(m, t, why)) return why;
      Unit *mp = data.PlanPoint(m);
      if (!mp->IsMapped()) return Where(data, m) + " is calculated by a host but not mapped";
      if (host->NumSubPoints() > 0 && mp->NumSubPoints() != host->NumSubPoints())
        return Where(data, m) + " has " + std::to_string(mp->NumSubPoints()) + " sub-points, its host " +
               std::to_string(host->NumSubPoints()) + ": the host would write past them";
    }
  }
  for (PointRef u : points)
    if (writer[u] < 0) return Where(data, u) + " is calculated by no task";
  return std::string();
}

void CalculationPlan::WriteRanges(WorkingCopy &data, int task, std::vector<std::pair<std::size_t, std::size_t>> &ranges) const {
  ranges.clear();
  AddUnitRanges(data, *data.PlanPoint(tasks_[task].host), ranges);
  for (PointRef m : tasks_[task].mapped) AddUnitRanges(data, *data.PlanPoint(m), ranges);
}

bool HasSubPointFold(const WorkingCopy &data, const Unit &point) {
  return point.IsTargetEffect() && data.Data().GetTargetEffect(point.GetTargetEffectNum())->IsSubPointEffect();
}

void CalculationPlan::TaskUnits(WorkingCopy &data, int task, std::vector<UnitRef> &units) const {
  auto add = [&](PointRef r) {
    const Unit &p = data.UnitAt(r);
    if (!HasSubPointFold(data, p)) {
      units.push_back(r);
      return;
    }
    for (int i = 0; i < p.NumSubPoints(); i++) units.push_back(p.SubPointFirst() + i);
  };
  add(tasks_[task].host);
  for (PointRef m : tasks_[task].mapped) add(m);
}

namespace {

std::uint64_t Bits(double x) {
  std::uint64_t b;
  std::memcpy(&b, &x, sizeof b);
  return b;
}

char RouteOf(const Unit &u, const Config &configure) {
  if (u.IsTHM()) return 'T';
  return (configure.paramMask & Config::USE_AMATRIX) ? 'A' : 'R';
}

}  // namespace

const std::vector<EnergyTask> &CalculationPlan::EnergyTasks(WorkingCopy &data, CNuc *cnuc, const Config &configure,
                                                            const std::vector<int> &work) {
  // The units of the work in plan order, with the plan task of each, and their E_x now.
  std::vector<UnitRef> units;
  std::vector<int> owners;
  for (int t : work) {
    TaskUnits(data, t, units);
    owners.resize(units.size(), t);
  }
  std::vector<double> keys(units.size());
  for (std::size_t k = 0; k < units.size(); k++) keys[k] = AMatrixFunc::InEnergy(cnuc, &data.UnitAt(units[k]));
  const bool same = work == energy_work_ && units == energy_units_ && keys.size() == energy_keys_.size() &&
                    (keys.empty() || std::memcmp(keys.data(), energy_keys_.data(), keys.size() * sizeof(double)) == 0);
  if (same) return energy_tasks_;

  struct Key {
    int nucleus;
    char route;
    double energy;
    std::uint64_t bits;
    int entrance, exit;
    std::size_t ordinal;
  };
  std::vector<Key> order(units.size());
  for (std::size_t k = 0; k < units.size(); k++) {
    Unit &u = data.UnitAt(units[k]);
    order[k] = {cnuc->GetPair(cnuc->GetPairNumFromKey(u.GetEntranceKey()))->GetNucleus(), RouteOf(u, configure), keys[k],
                Bits(keys[k]), u.GetEntranceKey(), u.GetExitKey(), k};
  }
  // By task key (nucleus, route, E_x, bits), then by (entrance, exit) pair, then in plan order.
  std::sort(order.begin(), order.end(), [](const Key &a, const Key &b) {
    if (a.nucleus != b.nucleus) return a.nucleus < b.nucleus;
    if (a.route != b.route) return a.route < b.route;
    if (a.energy != b.energy) return a.energy < b.energy;
    if (a.bits != b.bits) return a.bits < b.bits;
    if (a.entrance != b.entrance) return a.entrance < b.entrance;
    if (a.exit != b.exit) return a.exit < b.exit;
    return a.ordinal < b.ordinal;
  });
  std::vector<EnergyTask> tasks;
  for (std::size_t k = 0; k < order.size(); k++) {
    const Key &o = order[k];
    if (k == 0 || o.nucleus != order[k - 1].nucleus || o.route != order[k - 1].route || o.bits != order[k - 1].bits) {
      tasks.emplace_back();
      tasks.back().nucleus = o.nucleus, tasks.back().route = o.route, tasks.back().energy = o.energy;
    }
    EnergyTask &task = tasks.back();
    if (task.units.empty() || o.entrance != order[k - 1].entrance || o.exit != order[k - 1].exit) task.numPairs++;
    task.units.push_back(units[o.ordinal]);
    task.owners.push_back(owners[o.ordinal]);
  }
  // Largest first (one fill and factorisation, one T-matrix per pair, one observable per unit); ties in energy order.
  std::stable_sort(tasks.begin(), tasks.end(), [](const EnergyTask &a, const EnergyTask &b) {
    return 1 + a.numPairs + a.units.size() > 1 + b.numPairs + b.units.size();
  });
  energy_tasks_.swap(tasks);
  energy_work_ = work;
  energy_units_.swap(units);
  energy_keys_.swap(keys);
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
