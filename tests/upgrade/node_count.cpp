// The node count of a folded project (docs/upgrade/energy_grid.md section 2.2): how many distinct compound energies
// the units of the point pass lie on today, and how many nodes one shared mesh per compound nucleus at today's node
// density would have.  A measuring tool for the shared mesh (branch 7, part 2), not a test: it prints, it asserts
// nothing.  Each project is set up as a calculation is (Config, CNuc::Fill, WorkingCopy::Fill, the Initialize steps; no
// evaluation).
//
//   node_count [--plan] [--detail] <dir/project.azr> ...   (each project is read in its own directory: data paths are relative)
//
// Per project, one line: the points, the points with a sub-point fold, the units of the pass (a point without a fold,
// or each sub-point of one), the distinct (nucleus, E_x bits) of those units ("today": what the energy tasks solve), and
// the estimate of a shared mesh,
//     sum over nuclei of  [ integral over the union of the folded points' windows of max_g 1/h_g(E_x) dE_x
//                           + 2 per folded point (the kernel's edges)
//                           + the distinct E_x of the points without a fold ],
// where g runs over the sub-point grids of the folded points whose window covers E_x and h_g is the spacing of g's
// nodes there (section 1.2's definition, on the grids the set-up placed).  E_x = AMatrixFunc::InEnergy, MeV.
//
// With --plan, the shared mesh the planner (include/MeshPlanner.h) places on the same supports and anchors, and more
// columns (section 2.8): per folded point the support [min, max] of its sub-points' E_x with the settings of its grid
// (Unit::SubGridSpec: smooth step, points per width, core widths; a uniform grid without a spec has no lattices), per
// nucleus the anchors of its levels quantized in E_x (AdaptiveIntegrationGrid::Anchors with the frame shift -(S + E_x)
// of the entrance pair of its first folded point, before the widths are transformed, as the set-up saw them); the planned nodes plus the exact energies of the points without a fold; the
// over-coverage (planned nodes with no node of today's grids, of any point, within the mesh's own spacing there); the
// nodes outside every support (the outer nodes of the edge intervals); and the cost of the ladder alone, the integral
// of max 1/L over the union of the supports divided by that of max 1/S (S today's smooth step, L its power of two);
// the supports in which the mesh has fewer nodes (inside, plus the two outer ones) than the point's own grid today,
// and the smallest such ratio.  --detail (implies --plan) prints the supports, the anchors and the supports with fewer
// nodes than today to stderr.
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
#include "ProjectSetup.h"
#include "AdaptiveIntegrationGrid.h"
#include "AZUREParams.h"
#include "CNuc.h"
#include "CalculationPlan.h"
#include "DataSet.h"
#include "MeshPlanner.h"
#include "Config.h"
#include "PPair.h"
#include "Segment.h"
#include "Unit.h"
#include "WorkingCopy.h"

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

bool g_plan = false;
bool g_detail = false;  // --detail: the supports, anchors and short supports on stderr

/// The planner's columns for one project: supports and anchors as the set-up has them (see the head of the file).
/// The quantized anchors in E_x per entrance pair key, as the set-up's grids saw the levels: taken after
/// WorkingCopy::Fill (which places the grids) and before CNuc::Initialize transforms the widths.
std::map<int, std::vector<azure2::MeshAnchor>> AnchorsByPair(CNuc &compound, const Config &cfg) {
  std::map<int, std::vector<azure2::MeshAnchor>> out;
  for (int i = 1; i <= compound.NumPairs(); i++) {
    const PPair *pair = compound.GetPair(i);
    AdaptiveIntegrationGrid::GridConfig gc;
    gc.entranceKey = pair->GetPairKey();
    gc.inputWidthsArePhysical = (cfg.paramMask & Config::TRANSFORM_PARAMETERS) && !compound.IsTransformedIn();
    std::vector<azure2::MeshAnchor> &anchors = out[pair->GetPairKey()];
    for (const auto &r : AdaptiveIntegrationGrid(gc).Anchors(&compound, -(pair->GetSepE() + pair->GetExE())))
      anchors.push_back({r.energy, r.totalWidth, r.particleWidth, r.id});
  }
  return out;
}

void PlanColumns(CNuc &compound, WorkingCopy &data, const std::map<int, std::vector<azure2::MeshAnchor>> &anchorsByPair,
                 const std::set<std::pair<int, std::uint64_t>> &today,
                 const std::set<std::pair<int, std::uint64_t>> &unfolded) {
  const azure2::CalculationPlan &plan = data.Plan();
  auto nucleusOf = [&](const Unit &u) { return compound.GetPair(compound.GetPairNumFromKey(u.GetEntranceKey()))->GetNucleus(); };
  std::map<int, std::vector<azure2::MeshSupport>> supports;
  std::map<int, std::vector<std::vector<double>>> grids;  // per nucleus, per support: today's nodes (E_x, ascending)
  std::map<int, int> entranceOf;  // a pair key of each nucleus
  std::map<int, std::vector<std::tuple<double, double, double>>> smoothDensity, ladderDensity;
  for (int t = 0; t < plan.NumTasks(); t++) {
    Unit &p = *plan.Point(data, t);
    if (!azure2::HasSubPointFold(data, p)) continue;
    std::vector<azure2::UnitRef> refs;
    plan.TaskUnits(data, t, refs);
    std::vector<double> grid;
    for (azure2::UnitRef r : refs) grid.push_back(AMatrixFunc::InEnergy(&compound, &data.UnitAt(r)));
    std::sort(grid.begin(), grid.end());
    const double lo = grid.front(), hi = grid.back();
    const int nucleus = nucleusOf(p);
    grids[nucleus].push_back(grid);
    entranceOf.emplace(nucleus, p.GetEntranceKey());
    const Unit::SubGridSpec &spec = data.Data().GridSpecOf(p);
    azure2::MeshSupport s;
    s.low = lo;
    s.high = hi;
    s.dataSet = p.GetEntranceKey();
    if (spec.valid) {
      s.smoothStep = spec.baseEnergyStep;
      s.settings.pointsPerWidth = spec.pointsPerWidth;
      s.settings.coreWidths = spec.resonanceWidthMultiplier;
    } else {
      s.smoothStep = refs.size() > 1 ? (hi - lo) / (double)(refs.size() - 1) : 0.0;
      s.settings.pointsPerWidth = 0.0;  // a uniform grid: no lattices
      s.settings.coreWidths = 0.0;
    }
    supports[nucleus].push_back(s);
    if (g_detail)
      std::fprintf(stderr, "support %zu [%.9f, %.9f] valid %d smooth %.4g p %g m %g ladder %.4g nodes %zu\n",
                   supports[nucleus].size() - 1, lo, hi, (int)spec.valid, s.smoothStep, s.settings.pointsPerWidth,
                   s.settings.coreWidths, azure2::MeshLadderStep(s), grid.size());
    const double ladder = azure2::MeshLadderStep(s);
    if (ladder > 0.0) {
      const double range = hi - lo;
      double base = s.smoothStep;
      if (!(base > 0.0) || base > range) base = range;
      const double smooth = range / std::max(1.0, std::ceil(range / base - 1.0e-9));
      smoothDensity[nucleus].emplace_back(lo, hi, 1.0 / smooth);
      ladderDensity[nucleus].emplace_back(lo, hi, 1.0 / ladder);
    }
  }
  azure2::MeshTailConstants tails;
  tails.tailRatio = AdaptiveIntegrationGrid::tailRatio;
  tails.convolutionCoreWidths = AdaptiveIntegrationGrid::convolutionCoreWidths;
  long planned = 0, over = 0, outside = 0, anchorsUsed = 0;
  double smoothNodes = 0.0, ladderNodes = 0.0, smallest = 1e300;
  long fewer = 0;
  for (auto &kv : supports) {
    const int nucleus = kv.first;
    const std::vector<azure2::MeshAnchor> &anchors = anchorsByPair.at(entranceOf[nucleus]);
    anchorsUsed += (long)anchors.size();
    if (g_detail)
      for (const auto &a : anchors) std::fprintf(stderr, "anchor nucleus %d E_x %.9f Gamma %.4g id %d\n", nucleus, a.energy, a.particleWidth, a.id);
    const azure2::MeshPlan mesh = azure2::PlanMesh(anchors, kv.second, tails);
    // The nucleus's mesh: the planned nodes and the exact energies of its unfolded points.
    std::set<double> all(mesh.nodes.begin(), mesh.nodes.end());
    for (const auto &u : unfolded) {
      if (u.first != nucleus) continue;
      double e;
      std::memcpy(&e, &u.second, sizeof e);
      all.insert(e);
    }
    planned += (long)all.size();
    std::vector<double> mine;
    for (const auto &u : today) {
      if (u.first != nucleus) continue;
      double e;
      std::memcpy(&e, &u.second, sizeof e);
      mine.push_back(e);
    }
    std::sort(mine.begin(), mine.end());
    std::vector<std::pair<double, double>> covered;
    for (const azure2::MeshSupport &s : kv.second) covered.emplace_back(s.low, s.high);
    std::sort(covered.begin(), covered.end());
    const std::vector<double> &x = mesh.nodes;
    for (size_t j = 0; j < x.size(); j++) {
      double spacing = 1e300;
      if (j > 0) spacing = std::min(spacing, x[j] - x[j - 1]);
      if (j + 1 < x.size()) spacing = std::min(spacing, x[j + 1] - x[j]);
      const auto it = std::lower_bound(mine.begin(), mine.end(), x[j] - spacing);
      if (it == mine.end() || *it > x[j] + spacing) over++;
      bool in = false;
      for (const auto &c : covered) {
        if (c.first > x[j]) break;
        if (c.second >= x[j]) { in = true; break; }
      }
      if (!in) outside++;
    }
    // Per support: the mesh nodes inside it against the nodes of the point's own grid today.
    for (size_t i = 0; i < kv.second.size(); i++) {
      const std::vector<double> &own = grids[nucleus][i];
      const azure2::MeshSupport &s = kv.second[i];
      if (own.size() < 2) continue;
      const long inside = (long)(std::upper_bound(x.begin(), x.end(), s.high) - std::lower_bound(x.begin(), x.end(), s.low));
      const double ratio = (double)(inside + 2) / (double)own.size();  // + the two outer nodes the quadrature uses
      if (ratio < 1.0) fewer++;
      smallest = std::min(smallest, ratio);
      if (g_detail && ratio < 1.0)
        std::fprintf(stderr, "fewer: support %zu [%.9f, %.9f] mesh %ld + 2, own %zu\n", i, s.low, s.high, inside, own.size());
    }
    smoothNodes += MeshNodes(smoothDensity[nucleus]);
    ladderNodes += MeshNodes(ladderDensity[nucleus]);
  }
  std::printf("\t%ld\t%.2f\t%ld\t%ld\t%ld\t%.2f\t%ld\t%.2f", planned, (double)planned / (double)today.size(), over, outside,
              anchorsUsed, smoothNodes > 0.0 ? ladderNodes / smoothNodes : 0.0, fewer, smallest);
}

void RunProject(const std::string &path) {
  const size_t slash = path.find_last_of('/');
  const std::string dir = slash == std::string::npos ? "." : path.substr(0, slash);
  const std::string file = slash == std::string::npos ? path : path.substr(slash + 1);
  if (chdir(dir.c_str()) != 0) return (void)std::printf("%s\tno directory\n", path.c_str());
  std::ostringstream sink;
  Config cfg(sink);
  cfg.configfile = file;
  CNuc compound;
  WorkingCopy data;
  AZUREParams params;
  std::map<int, std::vector<azure2::MeshAnchor>> anchorsByPair;
  try {
    if (cfg.ReadConfigFile() < 0 || compound.Fill(cfg) == -1 || data.Fill(cfg, &compound) == -1)
      return (void)std::printf("%s\tdoes not set up\n", path.c_str());
    if (g_plan) anchorsByPair = AnchorsByPair(compound, cfg);
    RequestDataEntrances(compound, cfg);
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
  std::printf("%s\t%d\t%ld\t%ld\t%zu\t%.0f\t%.2f", path.c_str(), plan.NumTasks(), folded, units, today.size(), mesh,
              mesh / (double)today.size());
  if (g_plan) PlanColumns(compound, data, anchorsByPair, today, unfolded);
  std::printf("\n");
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: node_count [--plan] [--detail] <dir/project.azr> ...\n");
    return 2;
  }
  char start[4096];
  if (!getcwd(start, sizeof(start))) return 2;
  int first = 1;
  for (; first < argc && argv[first][0] == '-'; first++) {
    const std::string option = argv[first];
    if (option == "--plan") g_plan = true;
    if (option == "--detail") g_detail = g_plan = true;
  }
  std::printf("project\tpoints\tfolded points\tunits\tdistinct E_x today\tshared mesh estimate\tmesh / today");
  if (g_plan)
    std::printf("\tplanned mesh\tplanned / today\tover-coverage\toutside supports\tanchors\tladder cost"
                "\tsupports with fewer nodes than today\tsmallest mesh/own nodes in a support");
  std::printf("\n");
  for (int k = first; k < argc; k++) {
    if (chdir(start) != 0) return 2;
    RunProject(argv[k]);
    std::fflush(stdout);
  }
  return 0;
}
