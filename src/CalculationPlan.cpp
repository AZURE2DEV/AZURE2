// The calculation plan (CalculationPlan.h; docs/upgrade/data_model.md step S5).
#include "CalculationPlan.h"

#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdexcept>

#include "WorkingCopy.h"
#include "Unit.h"
#include "Segment.h"
#include "Results.h"

namespace azure2 {

namespace {

std::string Where(EData &data, PointRef r) {
  std::ostringstream s;
  s << "unit " << r;
  if (r >= 0 && r < data.NumUnits())
    s << " (key " << data.Unit(r).GetSegmentKey() << ", E_lab " << data.Unit(r).GetLabEnergy() << " MeV)";
  return s.str();
}

// The slots a point and its sub-points hold apart from the natural terms (the model pass's).
void AddUnitRanges(EData &data, EPoint &p, std::vector<std::pair<std::size_t, std::size_t>> &ranges) {
  auto add = [&ranges](EPoint &u) {
    if (u.ResultSlot() == EPoint::kNoRecord) return;
    const std::size_t n = Results::RecordSize(u.ResultAngularSlots(), 0);
    ranges.push_back({u.ResultSlot(), u.ResultSlot() + n});
  };
  add(p);
  for (EPoint &sub : p.GetSubPoints(data)) add(sub);
}

}  // namespace

EPoint *MappedPoints::operator[](int i) const { return data_->PlanPoint((*refs_)[i]); }

void CalculationPlan::Build(EData &data, const std::vector<std::pair<PointRef, PointRef>> &mapping) {
  tasks_.clear();
  for (PointRef u : data.PlanUnits()) {
    EPoint *point = data.PlanPoint(u);
    point->SetPlanTask(-1);
    if (point->IsMapped()) continue;
    tasks_.push_back({u, {}});
    point->SetPlanTask((int)tasks_.size() - 1);
  }
  for (const std::pair<PointRef, PointRef> &m : mapping) {
    const int task = data.PlanPoint(m.first)->PlanTask();
    if (task < 0) throw std::logic_error("EData: the calculation plan maps " + Where(data, m.second) + " onto " +
                                         Where(data, m.first) + ", which is not the host of a task");
    tasks_[task].mapped.push_back(m.second);
    data.PlanPoint(m.second)->SetPlanTask(task);
  }
  const std::string why = Check(data);
  if (!why.empty()) throw std::logic_error("EData: the calculation plan does not give every point one writer: " + why);
}

namespace {

// One task's writes: `before` the slots before it ran, `writer` the task that wrote each slot so far (-1: none).
void CheckTaskWrites(EData &data, int task, const std::vector<double> &before, const Results &results,
                     std::vector<int> &writer) {
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  if (task >= 0) data.Plan().WriteRanges(data, task, ranges);
  if (writer.size() != results.slot.size()) writer.assign(results.slot.size(), -1);
  for (std::size_t i = 0; i < results.slot.size() && i < before.size(); i++) {
    if (std::memcmp(&before[i], &results.slot[i], sizeof(double)) == 0) continue;
    bool own = false;
    for (const auto &r : ranges) own = own || (i >= r.first && i < r.second);
    std::string why;
    if (!own) why = "writes slot " + std::to_string(i) + ", outside its records";
    else if (writer[i] >= 0 && writer[i] != task)
      why = "writes slot " + std::to_string(i) + ", which task " + std::to_string(writer[i]) + " wrote";
    if (!why.empty())
      throw std::logic_error("AZURE_PLAN_CHECK: task " + std::to_string(task) +
                             (task >= 0 ? " (host " + Where(data, data.Plan().Task(task).host) + ")" : std::string()) +
                             " " + why);
    writer[i] = task;
  }
}

}  // namespace

void CalculationPlan::RunChecked(EData &data, const std::vector<EPoint *> &hosts, Results &results,
                                 const std::function<void(EPoint *)> &run) {
  std::vector<int> writer;
  std::vector<double> before;
  for (EPoint *host : hosts) {
    before = results.slot;
    run(host);
    CheckTaskWrites(data, host->PlanTask(), before, results, writer);
  }
}

std::string CalculationPlan::Check(EData &data) const {
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
    EPoint *host = data.PlanPoint(task.host);
    if (host->IsMapped()) return Where(data, task.host) + " is a host and mapped onto another point";
    for (PointRef m : task.mapped) {
      if (!claim(m, t, why)) return why;
      EPoint *mp = data.PlanPoint(m);
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

void CalculationPlan::WriteRanges(EData &data, int task, std::vector<std::pair<std::size_t, std::size_t>> &ranges) const {
  ranges.clear();
  AddUnitRanges(data, *data.PlanPoint(tasks_[task].host), ranges);
  for (PointRef m : tasks_[task].mapped) AddUnitRanges(data, *data.PlanPoint(m), ranges);
}

bool CalculationPlan::RunTimeCheckRequested() {
  static const bool on = [] {
    const char *v = std::getenv("AZURE_PLAN_CHECK");
    return v && *v && std::strcmp(v, "0") != 0;
  }();
  return on;
}

}  // namespace azure2
