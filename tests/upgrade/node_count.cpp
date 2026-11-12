// The node count of a folded project (docs/upgrade/energy_grid.md section 2.2): how many distinct compound energies
// the units of the point pass lie on today, and how many nodes one shared mesh per compound nucleus at today's node
// density would have.  A measuring tool for the shared mesh (branch 7, part 2), not a test: it prints, it asserts
// nothing.  Each project is set up as a calculation is (Config, CNuc::Fill, WorkingCopy::Fill, the Initialize steps; no
// evaluation).
//
//   node_count <dir/project.azr> ...   (each project is read in its own directory: data paths are relative)
//
// Per project, one line: the points, the points with a sub-point fold, the units of the pass (a point without a fold,
// or each sub-point of one), the distinct (nucleus, E_x bits) of those units ("today": what the energy tasks solve), and
// the estimate of a shared mesh,
//     sum over nuclei of  [ integral over the union of the folded points' windows of max_g 1/h_g(E_x) dE_x
//                           + 2 per folded point (the kernel's edges)
//                           + the distinct E_x of the points without a fold ],
// where g runs over the sub-point grids of the folded points whose window covers E_x and h_g is the spacing of g's
// nodes there (section 1.2's definition, on the grids the set-up placed).  E_x = AMatrixFunc::InEnergy, MeV.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unistd.h>
#include <vector>

#include "AMatrixFunc.h"
#include "AZUREParams.h"
#include "CNuc.h"
#include "CalculationPlan.h"
#include "Config.h"
#include "PPair.h"
#include "Segment.h"
#include "Unit.h"
#include "WorkingCopy.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace {

std::uint64_t Bits(double x) {
  std::uint64_t b;
  std::memcpy(&b, &x, sizeof b);
  return b;
}

// The integral of max over the intervals covering E of their density (1/width), over the union of the intervals:
// a sweep over the interval ends with a multiset of the densities alive.
double MeshNodes(std::vector<std::tuple<double, double, double>> &intervals) {  // (low, high, density)
  std::vector<std::tuple<double, int, double>> events;  // (E, +1 open / -1 close, density)
  for (const auto &iv : intervals) {
    if (!(std::get<1>(iv) > std::get<0>(iv))) continue;
    events.emplace_back(std::get<0>(iv), +1, std::get<2>(iv));
    events.emplace_back(std::get<1>(iv), -1, std::get<2>(iv));
  }
  std::sort(events.begin(), events.end());
  std::multiset<double> alive;
  double total = 0.0, last = 0.0;
  for (const auto &ev : events) {
    const double e = std::get<0>(ev);
    if (!alive.empty()) total += (e - last) * *alive.rbegin();
    last = e;
    if (std::get<1>(ev) > 0)
      alive.insert(std::get<2>(ev));
    else
      alive.erase(alive.find(std::get<2>(ev)));
  }
  return total;
}

void RunProject(const std::string &path) {
  const size_t slash = path.find_last_of('/');
  const std::string dir = slash == std::string::npos ? "." : path.substr(0, slash);
  const std::string file = slash == std::string::npos ? path : path.substr(slash + 1);
  if (chdir(dir.c_str()) != 0) return (void)std::printf("%s\tno directory\n", path.c_str());
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = file;
  g_config = &cfg;
  CNuc compound;
  WorkingCopy data;
  AZUREParams params;
  try {
    if (cfg.ReadConfigFile() < 0 || compound.Fill(cfg) == -1 || data.Fill(cfg, &compound) == -1)
      return (void)std::printf("%s\tdoes not set up\n", path.c_str());
    compound.Initialize(cfg);
    compound.FillMnParams(params, &cfg);
    data.FillMnParams(params);
    if (data.Initialize(&compound, cfg) == -1) return (void)std::printf("%s\tdoes not initialize\n", path.c_str());
  } catch (const std::exception &e) {
    return (void)std::printf("%s\tset-up stops (%s)\n", path.c_str(), e.what());
  }
  const azure2::CalculationPlan &plan = data.Plan();
  auto nucleusOf = [&](const Site &u) { return compound.GetPair(compound.GetPairNumFromKey(u.GetEntranceKey()))->GetNucleus(); };
  std::set<std::pair<int, std::uint64_t>> today, unfolded;
  std::map<int, std::vector<std::tuple<double, double, double>>> intervals;  // by nucleus
  long units = 0, folded = 0, edges = 0;
  for (int t = 0; t < plan.NumTasks(); t++) {
    Unit &p = *plan.Point(data, t);
    std::vector<azure2::UnitRef> refs;
    plan.TaskUnits(data, t, refs);
    for (azure2::UnitRef r : refs) {
      Site &u = data.SiteAt(r);
      today.insert({nucleusOf(u), Bits(AMatrixFunc::InEnergy(&compound, &u))});
    }
    units += (long)refs.size();
    if (!azure2::HasSubPointFold(data, p)) {
      unfolded.insert({nucleusOf(p), Bits(AMatrixFunc::InEnergy(&compound, &p))});
      continue;
    }
    folded++;
    edges += 2;
    std::vector<double> ex;
    for (azure2::UnitRef r : refs) ex.push_back(AMatrixFunc::InEnergy(&compound, &data.SiteAt(r)));
    std::sort(ex.begin(), ex.end());
    for (size_t k = 0; k + 1 < ex.size(); k++)
      if (ex[k + 1] > ex[k]) intervals[nucleusOf(p)].emplace_back(ex[k], ex[k + 1], 1.0 / (ex[k + 1] - ex[k]));
  }
  if (folded == 0) return;
  double mesh = (double)edges + (double)unfolded.size();
  for (auto &kv : intervals) mesh += MeshNodes(kv.second);
  std::printf("%s\t%d\t%ld\t%ld\t%zu\t%.0f\t%.2f\n", path.c_str(), plan.NumTasks(), folded, units, today.size(), mesh,
              mesh / (double)today.size());
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: node_count <dir/project.azr> ...\n");
    return 2;
  }
  char start[4096];
  if (!getcwd(start, sizeof(start))) return 2;
  std::printf("project\tpoints\tfolded points\tunits\tdistinct E_x today\tshared mesh estimate\tmesh / today\n");
  for (int k = 1; k < argc; k++) {
    if (chdir(start) != 0) return 2;
    RunProject(argv[k]);
    std::fflush(stdout);
  }
  return 0;
}
