// The calculation plan gives every point one writer (include/CalculationPlan.h, docs/upgrade/data_model.md S5), on
// every project given: the set-up of a calculation (Config, CNuc::Fill, WorkingCopy::Fill, the Initialize steps), which
// builds the plan with the mapping (WorkingCopy::MapData) and checks it there -- a plan that breaks the assertion throws.
// For each project that sets up, also: the plan checks again here (CalculationPlan::Check), every point of the data
// set is a task's host or a point the plan maps onto one, a host's mapped points are those whose Unit::PlanTask is
// the host's, and a copy of the data set (WorkingCopy::Clone) holds the same plan, valid in the copy, whose hosts are the
// copy's own points; and the energy tasks of the point pass (CalculationPlan::EnergyTasks, CheckEnergyTasks below):
// every unit in one task, one nucleus, route and E_x per task.  No evaluation is run.  A project that does not set up (missing data files: the projects whose
// data their check.sh writes) is listed and counted, not a failure.
//
//   plan_test <dir/project.azr> ...   (each project is read in its own directory: data paths are relative)
#include <cstdio>
#include <cstring>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unistd.h>

#include "AMatrixFunc.h"
#include "AZUREParams.h"
#include "CNuc.h"
#include "CalculationPlan.h"
#include "Config.h"
#include "PPair.h"
#include "WorkingCopy.h"
#include "Unit.h"
#include "Segment.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace {

int failures = 0, setUp = 0, notSetUp = 0;
long tasks = 0, mapped = 0;

void Fail(const std::string &project, const std::string &why) {
  std::printf("  FAIL  %s: %s\n", project.c_str(), why.c_str());
  failures++;
}

// What the check of this test adds to CalculationPlan::Check: the points' task indices against the plan, on `data`.
std::string CrossCheck(WorkingCopy &data) {
  const azure2::CalculationPlan &plan = data.Plan();
  for (int t = 0; t < plan.NumTasks(); t++) {
    const azure2::CalculationTask &task = plan.Task(t);
    Unit *host = data.PlanPoint(task.host);
    const azure2::MappedPoints m = host->PlanMapped(&data);
    if (m.size() != (int)task.mapped.size()) return "task " + std::to_string(t) + ": its host does not see its mapped points";
    for (int k = 0; k < m.size(); k++)
      if (m[k] != data.PlanPoint(task.mapped[k]) || m[k]->PlanTask() != t || m[k]->PlanMapped(&data).size() != 0)
        return "task " + std::to_string(t) + ": mapped point " + std::to_string(k) + " is not the plan's";
  }
  return plan.Check(data);
}

// The energy tasks of every plan task (CalculationPlan::EnergyTasks): every unit of the plan tasks -- a point, or each
// sub-point of a point with a sub-point fold -- in exactly one energy task, with the plan task of its point; inside a
// task one nucleus, one route and the same bits of E_x (AMatrixFunc::InEnergy), the units by (entrance, exit) pair
// (and external-capture amplitude source: the T-matrix groups) and numTMatrices their count; the tasks by decreasing cost.  Counts the solves of a pass (energy tasks) against the
// units.
long energyTasks = 0, energyUnits = 0;
std::string CheckEnergyTasks(WorkingCopy &data, CNuc &compound, const Config &cfg) {
  azure2::CalculationPlan &plan = const_cast<azure2::CalculationPlan &>(data.Plan());
  std::vector<int> work(plan.NumTasks());
  for (int t = 0; t < plan.NumTasks(); t++) work[t] = t;
  const azure2::EnergyTaskList &list = plan.EnergyTasks(data, &compound, cfg, work);
  const std::vector<azure2::EnergyTask> &tasks = list.tasks;
  std::vector<int> owner(data.NumUnits(), -2);  // the plan task of each unit of the pass (-2: no unit of it)
  long expected = 0;
  for (int t = 0; t < plan.NumTasks(); t++) {
    std::vector<azure2::UnitRef> units;
    plan.TaskUnits(data, t, units);
    for (azure2::UnitRef u : units) {
      if (owner[u] != -2) return "unit " + std::to_string(u) + " is a unit of two plan tasks";
      owner[u] = t;
    }
    expected += (long)units.size();
  }
  long seen = 0;
  double lastCost = 1e300;
  for (std::size_t k = 0; k < tasks.size(); k++) {
    const azure2::EnergyTask &task = tasks[k];
    const std::string where = "energy task " + std::to_string(k);
    if (task.size <= 0 || task.begin < 0 || task.begin + task.size > (int)list.units.size() ||
        list.owners.size() != list.units.size() || list.tGroups.size() != list.units.size() ||
        list.providers.size() != list.units.size())
      return where + " is empty, outside the arrays, or has no owners";
    if (task.begin != (int)seen) return where + " does not follow the task before it in the arrays";
    const azure2::UnitRef *taskUnits = list.units.data() + task.begin, *taskProviders = list.providers.data() + task.begin;
    const int *taskOwners = list.owners.data() + task.begin, *taskGroups = list.tGroups.data() + task.begin;
    const double cost = 1.0 + task.numTMatrices + task.size;
    if (cost > lastCost) return where + " is out of the order of cost";
    lastCost = cost;
    int pairs = 0;
    for (int i = 0; i < task.size; i++) {
      const azure2::UnitRef r = taskUnits[i];
      if (r < 0 || r >= data.NumUnits() || owner[r] < 0) return where + " names a unit of no plan task (or twice)";
      if (owner[r] != taskOwners[i]) return where + ": unit " + std::to_string(r) + " has another owner";
      owner[r] = -1;
      Unit &u = data.UnitAt(r);
      const double e = AMatrixFunc::InEnergy(&compound, &u);
      if (std::memcmp(&e, &task.energy, sizeof e) != 0) return where + ": a unit at another E_x";
      if (compound.GetPair(compound.GetPairNumFromKey(u.GetEntranceKey()))->GetNucleus() != task.nucleus)
        return where + ": a unit of another nucleus";
      if ((u.IsTHM() ? 'T' : 'A') != task.route && task.route != 'R') return where + ": a unit of another route";
      Unit &prov = data.UnitAt(taskProviders[i]);
      const std::tuple<int, int, long long> b{u.GetEntranceKey(), u.GetExitKey(), azure2::ECAmplitudeSource(prov, cfg)};
      if (prov.GetEntranceKey() != u.GetEntranceKey() || prov.GetExitKey() != u.GetExitKey())
        return where + ": a provider of another pair";
      if (i > 0) {
        Unit &p = data.UnitAt(taskUnits[i - 1]);
        const std::tuple<int, int, long long> a{p.GetEntranceKey(), p.GetExitKey(),
                                                azure2::ECAmplitudeSource(data.UnitAt(taskProviders[i - 1]), cfg)};
        if (b < a) return where + ": units not ordered by T-matrix group";
        if ((b != a) != (taskGroups[i] != taskGroups[i - 1])) return where + ": T-matrix groups wrong";
        if (b != a) pairs++;
      } else
        pairs = 1;
      // a group's first unit is its own provider
      if ((i == 0 || taskGroups[i] != taskGroups[i - 1]) && taskProviders[i] != taskUnits[i])
        return where + ": a T-matrix group whose first unit is not its own provider";
    }
    if (pairs != task.numTMatrices)
      return where + ": numTMatrices " + std::to_string(task.numTMatrices) + ", groups " + std::to_string(pairs);
    seen += (long)task.size;
  }
  if (seen != expected) return std::to_string(expected - seen) + " units of the plan tasks in no energy task";
  energyTasks += (long)tasks.size();
  energyUnits += seen;
  return std::string();
}

void RunProject(const std::string &path) {
  const size_t slash = path.find_last_of('/');
  const std::string dir = slash == std::string::npos ? "." : path.substr(0, slash);
  const std::string file = slash == std::string::npos ? path : path.substr(slash + 1);
  if (chdir(dir.c_str()) != 0) return Fail(path, "no directory");
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = file;
  g_config = &cfg;
  CNuc compound;
  WorkingCopy data;
  AZUREParams params;
  try {
    if (cfg.ReadConfigFile() < 0 || compound.Fill(cfg) == -1 || data.Fill(cfg, &compound) == -1) {
      notSetUp++;
      std::printf("  --    %s: does not set up\n", path.c_str());
      return;
    }
    compound.Initialize(cfg);
    compound.FillMnParams(params, &cfg);
    data.FillMnParams(params);
    if (data.Initialize(&compound, cfg) == -1) {
      notSetUp++;
      std::printf("  --    %s: does not initialize\n", path.c_str());
      return;
    }
  } catch (const std::logic_error &e) {
    return Fail(path, e.what());  // the plan's assertion (or another logic error of the set-up)
  } catch (const std::exception &e) {
    notSetUp++;
    std::printf("  --    %s: set-up stops (%s)\n", path.c_str(), e.what());
    return;
  }
  setUp++;
  std::string why = CrossCheck(data);
  if (!why.empty()) return Fail(path, why);
  why = CheckEnergyTasks(data, compound, cfg);
  if (!why.empty()) return Fail(path, why);
  std::unique_ptr<WorkingCopy> copy(data.Clone());
  if (copy->Plan().NumTasks() != data.Plan().NumTasks()) return Fail(path, "the copy has another plan");
  why = CrossCheck(*copy);
  if (!why.empty()) return Fail(path, "in the copy: " + why);
  for (int t = 0; t < copy->Plan().NumTasks(); t++)
    if (copy->PlanPoint(copy->Plan().Task(t).host) == data.PlanPoint(data.Plan().Task(t).host))
      return Fail(path, "the copy's task " + std::to_string(t) + " is hosted by a point of the original");
  tasks += data.Plan().NumTasks();
  for (int t = 0; t < data.Plan().NumTasks(); t++) mapped += (long)data.Plan().Task(t).mapped.size();
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: plan_test <dir/project.azr> ...\n");
    return 2;
  }
  char start[4096];
  if (!getcwd(start, sizeof(start))) return 2;
  for (int k = 1; k < argc; k++) {
    if (chdir(start) != 0) return 2;
    RunProject(argv[k]);
  }
  std::printf("plan: %d projects set up (%ld tasks, %ld mapped points; %ld energy tasks for %ld units), %d do not set "
              "up; %s (%d failure%s)\n", setUp, tasks, mapped, energyTasks, energyUnits, notSetUp,
              failures ? "FAIL" : "one writer per point", failures, failures == 1 ? "" : "s");
  return failures || setUp == 0 ? 1 : 0;
}
