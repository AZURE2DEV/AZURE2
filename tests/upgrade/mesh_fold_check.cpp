// The quadrature of a folded point on today's private grid and on the planned shared mesh, against a converged
// reference (docs/upgrade/energy_grid.md §2.8).  A measuring tool, no engine: a Breit-Wigner cross section (one level,
// or two interfering levels as a sum of amplitudes) folded with a Gaussian of width sigma (support E0 +- 5 sigma, a
// convolution's settings: 2 core widths) and with a rectangle [E0 - sigma, E0] (a thick-target-like kernel: 20 core
// widths), by the rule of §2.4 (azure2::experiment::GridQuadrature: sigma linear between nodes, two-point
// Gauss-Legendre per interval, clipped at the support), on
//   (a) today's grid of each point, AdaptiveIntegrationGrid::GenerateGrid(high, low, anchors), and
//   (b) the mesh PlanMesh places for all the points of the case together (one data set),
// against the same rule on the union of both, each interval cut in 64.  Per case: the largest relative error over 13
// points E0 = E_r + k w / 2 (k = -6 ... 6, w = max(Gamma, sigma)), and the nodes of (a) (distinct over the points)
// and (b).
//
//   mesh_fold_check                 prints the table (tab-separated)
//   mesh_fold_check --with-edges    the same with each support's edges added to its mesh nodes (a diagnosis: the
//                                   mesh clips the edge intervals, today's grids have the edges as nodes)
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <set>
#include <string>
#include <vector>

#include "AdaptiveIntegrationGrid.h"
#include "Config.h"
#include "MeshPlanner.h"

namespace {

struct Level {
  double energy, width;  // MeV
  std::complex<double> phase;
};

double CrossSection(const std::vector<Level> &levels, double e) {
  std::complex<double> amplitude = 0.0;
  for (const Level &l : levels) amplitude += l.phase * (0.5 * l.width) / std::complex<double>(e - l.energy, 0.5 * l.width);
  return std::norm(amplitude);
}

/// The rule of §2.4 on nodes in descending order: num / den of the kernel K over [lo, hi].
double Fold(const std::vector<double> &nodes, const std::vector<Level> &levels, const std::function<double(double)> &kernel,
            double lo, double hi) {
  const double x1 = -1.0 / std::sqrt(3.0), x2 = 1.0 / std::sqrt(3.0);
  double num = 0.0, den = 0.0;
  for (size_t j = 0; j + 1 < nodes.size(); j++) {
    const double ea = nodes[j], eb = nodes[j + 1];
    if (eb >= hi || ea <= lo) continue;
    const double a = std::min(ea, hi), b = std::max(eb, lo);
    const double mid = 0.5 * (a + b), half = 0.5 * (a - b);
    if (half < 1.0e-12) continue;
    const double sa = CrossSection(levels, ea), sb = CrossSection(levels, eb);
    for (double x : {x1, x2}) {
      const double e = mid + half * x;
      const double frac = (ea - e) / (ea - eb);
      num += half * kernel(e) * ((1.0 - frac) * sa + frac * sb);
      den += half * kernel(e);
    }
  }
  return num / den;
}

std::vector<double> Refined(std::vector<double> nodes, int parts) {  // ascending in, descending out
  std::sort(nodes.begin(), nodes.end());
  nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
  std::vector<double> out;
  for (size_t j = 0; j + 1 < nodes.size(); j++)
    for (int k = 0; k < parts; k++) out.push_back(nodes[j] + (nodes[j + 1] - nodes[j]) * k / parts);
  out.push_back(nodes.back());
  std::reverse(out.begin(), out.end());
  return out;
}

struct Row {
  double today = 0.0, mesh = 0.0;
  size_t todayNodes = 0, meshNodes = 0;
};

bool g_withEdges = false;

Row RunCase(const std::vector<Level> &levels, double sigma, bool gaussian) {
  const double coreWidths = gaussian ? AdaptiveIntegrationGrid::convolutionCoreWidths : 20.0;
  const int numPoints = gaussian ? 100 : 50;
  std::vector<AdaptiveIntegrationGrid::ResonanceInfo> info;
  std::vector<azure2::MeshAnchor> anchors;
  for (size_t i = 0; i < levels.size(); i++) {
    AdaptiveIntegrationGrid::ResonanceInfo r{levels[i].energy, levels[i].width, levels[i].width, 1001 + (int)i};
    r = AdaptiveIntegrationGrid::Quantize(r);
    info.push_back(r);
    anchors.push_back({r.energy, r.totalWidth, r.particleWidth, r.id});
  }
  double wmax = sigma;
  for (const Level &l : levels) wmax = std::max(wmax, l.width);
  std::vector<double> centres;
  for (int k = -6; k <= 6; k++) centres.push_back(levels[0].energy + 0.5 * k * wmax);
  std::vector<azure2::MeshSupport> supports;
  for (double e0 : centres) {
    azure2::MeshSupport s;
    s.low = gaussian ? e0 - 5.0 * sigma : e0 - sigma;
    s.high = gaussian ? e0 + 5.0 * sigma : e0;
    s.smoothStep = (s.high - s.low) / numPoints;
    s.settings.pointsPerWidth = 50.0;
    s.settings.coreWidths = coreWidths;
    supports.push_back(s);
  }
  azure2::MeshTailConstants tails;
  tails.tailRatio = AdaptiveIntegrationGrid::tailRatio;
  tails.convolutionCoreWidths = AdaptiveIntegrationGrid::convolutionCoreWidths;
  const azure2::MeshPlan plan = azure2::PlanMesh(anchors, supports, tails);
  Row row;
  row.meshNodes = plan.nodes.size();
  std::set<double> todayAll;
  for (size_t i = 0; i < supports.size(); i++) {
    const azure2::MeshSupport &s = supports[i];
    const double e0 = centres[i];
    std::function<double(double)> kernel;
    if (gaussian)
      kernel = [&](double e) { return std::exp(-0.5 * (e - e0) * (e - e0) / (sigma * sigma)); };
    else
      kernel = [](double) { return 1.0; };
    AdaptiveIntegrationGrid::GridConfig gc;
    gc.maxPoints = numPoints;
    gc.baseEnergyStep = s.smoothStep;
    gc.resonanceWidthMultiplier = coreWidths;
    gc.pointsPerWidth = 50.0;
    const std::vector<double> today = AdaptiveIntegrationGrid(gc).GenerateGrid(s.high, s.low, info);
    todayAll.insert(today.begin(), today.end());
    std::vector<double> mesh(plan.nodes.begin() + (long)plan.ranges[i].first, plan.nodes.begin() + (long)plan.ranges[i].last + 1);
    if (g_withEdges) {
      mesh.push_back(s.low);
      mesh.push_back(s.high);
      std::sort(mesh.begin(), mesh.end());
      mesh.erase(std::unique(mesh.begin(), mesh.end()), mesh.end());
    }
    std::vector<double> both = mesh;
    for (double e : today) both.push_back(e);
    for (double &e : both) e = std::min(std::max(e, s.low), s.high);
    const double reference = Fold(Refined(both, 64), levels, kernel, s.low, s.high);
    std::reverse(mesh.begin(), mesh.end());
    const double a = Fold(today, levels, kernel, s.low, s.high), b = Fold(mesh, levels, kernel, s.low, s.high);
    row.today = std::max(row.today, std::fabs(a - reference) / std::fabs(reference));
    row.mesh = std::max(row.mesh, std::fabs(b - reference) / std::fabs(reference));
  }
  row.todayNodes = todayAll.size();
  return row;
}

}  // namespace

int main(int argc, char **argv) {
  g_withEdges = argc > 1 && std::string(argv[1]) == "--with-edges";
  const double sigma = 0.01;  // MeV
  std::printf("kernel\tlevels\tGamma/sigma\terror today\terror mesh\tnodes today\tnodes mesh\tmesh <= today\n");
  for (int kernelKind = 0; kernelKind < 2; kernelKind++) {
    const bool gaussian = kernelKind == 0;
    for (int two = 0; two < 2; two++) {
      for (double ratio : {0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0}) {
        const double width = ratio * sigma;
        std::vector<Level> levels = {{1.0, width, 1.0}};
        // A second level 3 widths above, twice as wide, interfering with a phase of 2 radians.
        if (two) levels.push_back({1.0 + 3.0 * width, 2.0 * width, std::polar(0.7, 2.0)});
        const Row r = RunCase(levels, sigma, gaussian);
        std::printf("%s\t%s\t%g\t%.2e\t%.2e\t%zu\t%zu\t%s\n", gaussian ? "Gaussian" : "rectangle", two ? "two" : "one", ratio,
                    r.today, r.mesh, r.todayNodes, r.meshNodes, r.mesh <= r.today * (1.0 + 1.0e-6) + 1.0e-13 ? "yes" : "NO");
      }
    }
  }
  return 0;
}
