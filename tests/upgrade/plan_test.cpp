// The calculation plan gives every point one writer (include/CalculationPlan.h, docs/upgrade/data_model.md S5), on
// every project given: the set-up of a calculation (Config, CNuc::Fill, EData::Fill, the Initialize steps), which
// builds the plan with the mapping (EData::MapData) and checks it there -- a plan that breaks the assertion throws.
// For each project that sets up, also: the plan checks again here (CalculationPlan::Check), every point of the data
// set is a task's host or a point the plan maps onto one, a host's mapped points are those whose EPoint::PlanTask is
// the host's, and a copy of the data set (EData::Clone) holds the same plan, valid in the copy, whose hosts are the
// copy's own points.  No evaluation is run.  A project that does not set up (missing data files: the projects whose
// data their check.sh writes) is listed and counted, not a failure.
//
//   plan_test <dir/project.azr> ...   (each project is read in its own directory: data paths are relative)
#include <cstdio>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

#include "AZUREParams.h"
#include "CNuc.h"
#include "CalculationPlan.h"
#include "Config.h"
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
std::string CrossCheck(EData &data) {
  const azure2::CalculationPlan &plan = data.Plan();
  for (int t = 0; t < plan.NumTasks(); t++) {
    const azure2::CalculationTask &task = plan.Task(t);
    EPoint *host = data.PlanPoint(task.host);
    const azure2::MappedPoints m = host->PlanMapped(&data);
    if (m.size() != (int)task.mapped.size()) return "task " + std::to_string(t) + ": its host does not see its mapped points";
    for (int k = 0; k < m.size(); k++)
      if (m[k] != data.PlanPoint(task.mapped[k]) || m[k]->PlanTask() != t || m[k]->PlanMapped(&data).size() != 0)
        return "task " + std::to_string(t) + ": mapped point " + std::to_string(k) + " is not the plan's";
  }
  return plan.Check(data);
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
  EData data;
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
  std::unique_ptr<EData> copy(data.Clone());
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
  std::printf("plan: %d projects set up (%ld tasks, %ld mapped points), %d do not set up; %s (%d failure%s)\n", setUp,
              tasks, mapped, notSetUp, failures ? "FAIL" : "one writer per point", failures, failures == 1 ? "" : "s");
  return failures || setUp == 0 ? 1 : 0;
}
