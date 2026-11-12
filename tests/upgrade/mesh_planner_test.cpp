// The planner of the shared mesh (include/MeshPlanner.h, docs/upgrade/energy_grid.md §2.8), tested alone against
// today's private grids (AdaptiveIntegrationGrid::GenerateGrid on the same anchors).  Properties:
//   1. never coarser than today: every mesh interval inside a support is no longer than today's step D(E) of that
//      support anywhere on it (D as AdaptiveIntegrationGrid.cpp computes it: the smooth step of the window, the cores
//      and the tails), and today's lattice nodes are mesh nodes, bit for bit, wherever the anchor's lattice is at least
//      as fine as the mesh's smooth step;
//   2. shared: overlapping supports read the same nodes; N copies of a support plan the nodes of one;
//   3. a function of anchors, supports and settings: permutations change no bit; disjoint groups plan alike together
//      and apart;
//   4. local: nodes outside MeshChangedRegion do not move when an anchor changes;
//   5. deterministic: the nodes of 1 and 3 (core E_r + k pitch, tail recurrence) are today's bits.
// Cases: an isolated narrow level, two overlapping levels of different widths, a level at a support's edge, a level
// outside the support but in reach, no level, a support of zero width, and 200 random cases (fixed seed).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "AdaptiveIntegrationGrid.h"
#include "Config.h"
#include "MeshPlanner.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

namespace {

using azure2::MeshAnchor;
using azure2::MeshPlan;
using azure2::MeshSupport;

int failures = 0;
int checks = 0;
void Expect(bool ok, const std::string &what) {
  checks++;
  if (!ok) {
    failures++;
    if (failures <= 20) std::printf("FAIL: %s\n", what.c_str());
  }
}

azure2::MeshTailConstants Tails() {
  azure2::MeshTailConstants t;
  t.tailRatio = AdaptiveIntegrationGrid::tailRatio;
  t.convolutionCoreWidths = AdaptiveIntegrationGrid::convolutionCoreWidths;
  return t;
}

MeshAnchor Anchor(double energy, double width, int id) {
  AdaptiveIntegrationGrid::ResonanceInfo r{energy, width, width, id};
  r = AdaptiveIntegrationGrid::Quantize(r);
  return MeshAnchor{r.energy, r.totalWidth, r.particleWidth, r.id};
}

MeshSupport Support(double low, double high, int numPoints, double p, double m, int dataSet = 0) {
  MeshSupport s;
  s.low = low;
  s.high = high;
  s.smoothStep = (high - low) / numPoints;
  s.settings.pointsPerWidth = p;
  s.settings.coreWidths = m;
  s.dataSet = dataSet;
  return s;
}

AdaptiveIntegrationGrid::GridConfig TodayConfig(const MeshSupport &s) {
  AdaptiveIntegrationGrid::GridConfig gc;
  gc.maxPoints = 1000;
  gc.baseEnergyStep = s.smoothStep;
  gc.resonanceWidthMultiplier = s.settings.coreWidths;
  gc.pointsPerWidth = s.settings.pointsPerWidth;
  return gc;
}

/// Today's grid of a support, ascending.
std::vector<double> TodayGrid(const MeshSupport &s, const std::vector<MeshAnchor> &anchors) {
  std::vector<AdaptiveIntegrationGrid::ResonanceInfo> info;
  for (const MeshAnchor &a : anchors) info.push_back({a.energy, a.totalWidth, a.particleWidth, a.id});
  std::vector<double> g = AdaptiveIntegrationGrid(TodayConfig(s)).GenerateGrid(s.high, s.low, info);
  std::reverse(g.begin(), g.end());
  return g;
}

/// Today's step of a support (the rule of AdaptiveIntegrationGrid::GenerateGrid): the smooth step of the window and,
/// per anchor, its core pitch out to m widths and its tail beyond, each capped at the configured base step.
struct TodayStep {
  double smooth = 0.0, base = 0.0, p = 0.0, m = 0.0, tail = 0.0;
  std::vector<MeshAnchor> anchors;
  TodayStep(const MeshSupport &s, const std::vector<MeshAnchor> &a) : anchors(a) {
    const double range = s.high - s.low;
    base = s.smoothStep;
    if (!(base > 0.0) || base > range) base = range;
    int n = (int)std::ceil(range / base - 1.0e-9);
    if (n < 1) n = 1;
    smooth = range / n;
    p = s.settings.pointsPerWidth;
    m = s.settings.coreWidths;
    tail = AdaptiveIntegrationGrid::tailRatio / p *
           std::sqrt(std::max(AdaptiveIntegrationGrid::convolutionCoreWidths, m) /
                     AdaptiveIntegrationGrid::convolutionCoreWidths);
  }
  /// `tolerance` is today's relative widening of the core (1e-9), which only decides on which side a node lying on
  /// the core's edge falls; between nodes the core ends at m widths.
  double Local(const MeshAnchor &r, double e, double tolerance) const {
    const double pitch = r.particleWidth / p;
    const double d = std::fabs(e - r.energy);
    const double step = (d <= std::max(0.0, m) * r.particleWidth * (1.0 + tolerance)) ? pitch : std::max(pitch, tail * d);
    return std::min(step, base);
  }
  /// The smallest step on [lo, hi] (inside an interval between nodes): each anchor's step grows with the distance, so
  /// its minimum is at the point of the stretch nearest to it.
  double MinOn(double lo, double hi) const {
    double s = smooth;
    for (const MeshAnchor &r : anchors) {
      if (!(r.particleWidth > 0.0)) continue;
      s = std::min(s, Local(r, std::min(std::max(r.energy, lo), hi), 0.0));
    }
    return s;
  }
  /// The step at a node, as GenerateGrid decides it.
  double At(double e) const {
    double s = smooth;
    for (const MeshAnchor &r : anchors)
      if (r.particleWidth > 0.0) s = std::min(s, Local(r, e, 1.0e-9));
    return s;
  }
};

/// Property 1 for one support of a plan.  A positive `finestLadder` (the finest smooth step of the plan's supports)
/// also checks that today's lattice nodes finer than it are mesh nodes (only when the plan has a single class of
/// settings, whose lattices are then exactly today's).
struct Coarseness {
  double worstRatio = 0.0;     // largest (mesh interval / today's real interval) at an E
  double worstAt = 0.0;
  double coarserLength = 0.0;  // length of the support where the mesh interval is longer than today's real one
  double length = 0.0;
};

void CheckNeverCoarser(const MeshPlan &plan, size_t index, const MeshSupport &s, const std::vector<MeshAnchor> &anchors,
                       double finestLadder, const std::string &label, Coarseness *stats = nullptr) {
  const azure2::MeshNodeRange r = plan.ranges[index];
  const std::vector<double> &x = plan.nodes;
  if (s.high - s.low < azure2::MeshZeroWidth()) {
    Expect(r.first == r.last && x[r.first] == s.high, label + ": a zero-width support is its own exact node");
    return;
  }
  Expect(x[r.first] <= s.low && x[r.last] >= s.high, label + ": the range encloses the support");
  Expect(r.first + 1 < x.size() && x[r.first + 1] > s.low && x[r.last - 1] < s.high,
         label + ": the range is the tightest enclosing one");
  const TodayStep today(s, anchors);
  for (size_t j = r.first; j < r.last; j++) {
    const double lo = std::max(x[j], s.low), hi = std::min(x[j + 1], s.high);
    if (!(hi > lo)) continue;
    const double len = hi - lo, inset = 1.0e-9 * len;
    const double allowed = today.MinOn(lo + inset, hi - inset);
    if (!(len <= allowed * (1.0 + 1.0e-9) + 1.0e-14)) {
      char buf[256];
      std::snprintf(buf, sizeof buf, ": mesh interval [%.12g, %.12g] of %.4g MeV, today's step there %.4g MeV", lo, hi,
                    len, allowed);
      Expect(false, label + buf);
    } else {
      checks++;
    }
  }
  const std::vector<double> real = TodayGrid(s, anchors);
  // Today's real intervals are at most about twice its step (where two lattices cross, the last node of one and the
  // first of the other); a check that TodayStep is today's rule.
  for (size_t k = 0; k + 1 < real.size(); k++) {
    const double mid = 0.5 * (real[k] + real[k + 1]);
    Expect(real[k + 1] - real[k] <= 2.5 * today.At(mid) * (1.0 + 1.0e-9) + 1.0e-14,
           label + ": today's real grid is within its rule");
  }
  if (finestLadder > 0.0) {
    // Today's smooth nodes, bit for bit (GenerateGrid's formula), are the only nodes the mesh need not have.
    const double range = s.high - s.low;
    int n = (int)std::ceil(range / std::min(today.base, range) - 1.0e-9);
    if (n < 1) n = 1;
    std::set<double> smooth;
    for (int i = 0; i <= n; i++) smooth.insert(s.high - range * i / n);
    for (size_t k = 1; k + 1 < real.size(); k++) {
      const double e = real[k];
      if (smooth.count(e)) continue;
      if (!(today.At(e) <= finestLadder)) continue;
      Expect(std::binary_search(x.begin(), x.end(), e), label + ": today's lattice node is a mesh node, bit for bit");
    }
  }
  if (stats) {
    // The comparison with today's real grid, point by point (measured, not asserted: where today keeps its smooth
    // lattice beside a tail, or two lattices cross, its intervals are shorter than its own step): over the pieces
    // between the nodes of both grids, the length where the mesh interval is the longer, and the largest ratio.
    std::vector<double> cuts(real.begin(), real.end());
    for (size_t j = r.first; j <= r.last; j++) cuts.push_back(std::min(std::max(x[j], s.low), s.high));
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    for (size_t q = 0; q + 1 < cuts.size(); q++) {
      const double e = 0.5 * (cuts[q] + cuts[q + 1]), piece = cuts[q + 1] - cuts[q];
      const size_t k = (size_t)(std::upper_bound(real.begin(), real.end(), e) - real.begin());
      const size_t j = (size_t)(std::upper_bound(x.begin(), x.end(), e) - x.begin());
      if (k == 0 || k >= real.size() || j == 0 || j >= x.size()) continue;
      const double mine = std::min(x[j], s.high) - std::max(x[j - 1], s.low), theirs = real[k] - real[k - 1];
      stats->length += piece;
      if (mine > theirs * (1.0 + 1.0e-9)) {
        stats->coarserLength += piece;
        if (mine / theirs > stats->worstRatio) {
          stats->worstRatio = mine / theirs;
          stats->worstAt = e;
        }
      }
    }
  }
}

bool SameBits(const std::vector<double> &a, const std::vector<double> &b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++)
    if (!(a[i] == b[i])) return false;
  return true;
}

std::vector<double> Inside(const std::vector<double> &nodes, double lo, double hi) {
  std::vector<double> out;
  for (double e : nodes)
    if (e >= lo && e <= hi) out.push_back(e);
  return out;
}

void CheckCase(const std::string &label, const std::vector<MeshAnchor> &anchors, const std::vector<MeshSupport> &supports,
               bool print) {
  const MeshPlan plan = azure2::PlanMesh(anchors, supports, Tails());
  Expect(plan.ranges.size() == supports.size(), label + ": a range per support");
  Expect(std::is_sorted(plan.nodes.begin(), plan.nodes.end()) &&
             std::adjacent_find(plan.nodes.begin(), plan.nodes.end()) == plan.nodes.end(),
         label + ": nodes ascending and distinct");
  std::set<std::pair<double, double>> classes;
  for (const MeshSupport &s : supports) classes.insert({s.settings.pointsPerWidth, s.settings.coreWidths});
  double finestLadder = 1e300;
  for (const MeshSupport &s : supports)
    if (azure2::MeshLadderStep(s) > 0.0) finestLadder = std::min(finestLadder, azure2::MeshLadderStep(s));
  if (classes.size() != 1) finestLadder = 0.0;
  for (size_t i = 0; i < supports.size(); i++) {
    Coarseness stats;
    CheckNeverCoarser(plan, i, supports[i], anchors, finestLadder, label + " support " + std::to_string(i),
                      &stats);
    if (print)
      std::printf("  %-28s support %zu: mesh nodes in range %zu, today %zu; mesh longer than today's real interval on "
                  "%.3g of the support (largest ratio %.3g at %.9g)\n",
                  label.c_str(), i, plan.ranges[i].last - plan.ranges[i].first + 1, TodayGrid(supports[i], anchors).size(),
                  stats.length > 0 ? stats.coarserLength / stats.length : 0.0, stats.worstRatio, stats.worstAt);
  }
}

void NamedCases() {
  const double p = 50.0, m = 20.0;
  // An isolated narrow level.
  CheckCase("isolated narrow level", {Anchor(1.0, 1.0e-5, 1001)}, {Support(0.95, 1.05, 100, p, m)}, true);
  // Two overlapping levels of different widths.
  CheckCase("two overlapping levels", {Anchor(1.0, 1.0e-5, 1001), Anchor(1.0008, 2.0e-4, 2001)},
            {Support(0.95, 1.05, 100, p, m)}, true);
  // A level at the edge of a support.
  CheckCase("level at the edge", {Anchor(1.05, 1.0e-5, 1001)}, {Support(0.95, 1.05, 100, p, m)}, true);
  // A level outside the support, in reach of its tail.
  CheckCase("level outside, in reach", {Anchor(1.053, 1.0e-5, 1001)}, {Support(0.95, 1.05, 100, p, m)}, true);
  // No level.
  CheckCase("no level", {}, {Support(0.95, 1.05, 100, p, m)}, true);
  // A Gaussian-like core of 2 widths.
  CheckCase("convolution core", {Anchor(1.0, 3.3e-5, 1001)}, {Support(0.98, 1.02, 200, p, 2.0)}, true);
  // A support of zero width beside an ordinary one.
  {
    const std::vector<MeshSupport> s = {Support(1.0, 1.0, 100, p, m), Support(0.95, 1.05, 100, p, m)};
    CheckCase("zero-width support", {Anchor(1.0, 1.0e-5, 1001)}, s, true);
  }
}

std::vector<MeshAnchor> RandomAnchors(std::mt19937 &rng) {
  std::uniform_int_distribution<int> count(0, 5);
  std::uniform_real_distribution<double> energy(0.9, 1.1), logWidth(std::log(1.0e-7), std::log(1.0e-2));
  std::vector<MeshAnchor> a;
  const int n = count(rng);
  for (int i = 0; i < n; i++) a.push_back(Anchor(energy(rng), std::exp(logWidth(rng)), 1001 + i));
  return a;
}

std::vector<MeshSupport> RandomSupports(std::mt19937 &rng, bool oneClass) {
  std::uniform_int_distribution<int> count(1, 4), pick(0, 2), points(0, 3);
  std::uniform_real_distribution<double> centre(0.9, 1.1), logWidth(std::log(1.0e-4), std::log(0.1));
  const double ps[3] = {25.0, 50.0, 100.0}, ms[3] = {2.0, 5.0, 20.0};
  const int ns[4] = {50, 100, 200, 1000};
  const int fixedP = pick(rng), fixedM = pick(rng);
  std::vector<MeshSupport> s;
  const int n = count(rng);
  for (int i = 0; i < n; i++) {
    const double c = centre(rng), w = std::exp(logWidth(rng));
    const double p = ps[oneClass ? fixedP : pick(rng)], m = ms[oneClass ? fixedM : pick(rng)];
    s.push_back(Support(c - 0.5 * w, c + 0.5 * w, ns[points(rng)], p, m, i));
  }
  return s;
}

void RandomCases() {
  std::mt19937 rng(20261112u);
  for (int c = 0; c < 200; c++) {
    const std::vector<MeshAnchor> anchors = RandomAnchors(rng);
    const std::vector<MeshSupport> supports = RandomSupports(rng, c % 2 == 0);
    CheckCase("random " + std::to_string(c), anchors, supports, false);
  }
}

void Shared() {
  const double p = 50.0, m = 20.0;
  const std::vector<MeshAnchor> anchors = {Anchor(1.0, 1.0e-5, 1001), Anchor(1.02, 4.0e-4, 2001)};
  // Two data sets with the same settings, overlapping: the same nodes over the overlap.
  const std::vector<MeshSupport> two = {Support(0.95, 1.01, 100, p, m, 0), Support(0.99, 1.05, 80, p, m, 1)};
  const MeshPlan plan = azure2::PlanMesh(anchors, two, Tails());
  auto nodesOf = [&](size_t i, double lo, double hi) {
    std::vector<double> range(plan.nodes.begin() + (long)plan.ranges[i].first,
                              plan.nodes.begin() + (long)plan.ranges[i].last + 1);
    return Inside(range, lo, hi);
  };
  Expect(SameBits(nodesOf(0, 0.99, 1.01), nodesOf(1, 0.99, 1.01)), "shared: two data sets read the same overlap nodes");
  // N identical supports plan the nodes of one.
  const std::vector<MeshSupport> one = {Support(0.95, 1.05, 100, p, m, 0)};
  std::vector<MeshSupport> many;
  for (int k = 0; k < 7; k++) many.push_back(Support(0.95, 1.05, 100, p, m, k));
  const MeshPlan a = azure2::PlanMesh(anchors, one, Tails()), b = azure2::PlanMesh(anchors, many, Tails());
  Expect(SameBits(a.nodes, b.nodes), "shared: N identical supports have the nodes of one");
  bool sameRanges = true;
  for (const azure2::MeshNodeRange &r : b.ranges)
    sameRanges = sameRanges && r.first == a.ranges[0].first && r.last == a.ranges[0].last;
  Expect(sameRanges, "shared: N identical supports have the range of one");
}

void Permutations() {
  std::mt19937 rng(7u);
  for (int c = 0; c < 40; c++) {
    std::vector<MeshAnchor> anchors = RandomAnchors(rng);
    std::vector<MeshSupport> supports = RandomSupports(rng, false);
    const MeshPlan base = azure2::PlanMesh(anchors, supports, Tails());
    for (int k = 0; k < 3; k++) {
      std::vector<size_t> order(supports.size());
      for (size_t i = 0; i < order.size(); i++) order[i] = i;
      std::shuffle(order.begin(), order.end(), rng);
      std::vector<MeshSupport> s2;
      for (size_t i : order) {
        s2.push_back(supports[i]);
        s2.back().dataSet = (int)(order.size() - i);  // other data-set numbers too
      }
      std::vector<MeshAnchor> a2 = anchors;
      std::shuffle(a2.begin(), a2.end(), rng);
      const MeshPlan p2 = azure2::PlanMesh(a2, s2, Tails());
      Expect(SameBits(base.nodes, p2.nodes), "permutation " + std::to_string(c) + ": the same nodes, bit for bit");
      bool ranges = true;
      for (size_t i = 0; i < order.size(); i++)
        ranges = ranges && p2.ranges[i].first == base.ranges[order[i]].first &&
                 p2.ranges[i].last == base.ranges[order[i]].last;
      Expect(ranges, "permutation " + std::to_string(c) + ": the ranges follow their supports");
    }
  }
  // Two disjoint groups (each with its levels) planned together and apart.
  const double p = 50.0, m = 20.0;
  const std::vector<MeshAnchor> anchors = {Anchor(0.92, 1.0e-5, 1001), Anchor(0.93, 2.0e-4, 1002),
                                           Anchor(1.07, 3.0e-6, 2001)};
  const std::vector<MeshSupport> g1 = {Support(0.90, 0.94, 100, p, m), Support(0.91, 0.95, 60, p, 2.0)};
  const std::vector<MeshSupport> g2 = {Support(1.05, 1.10, 200, 25.0, m)};
  std::vector<MeshSupport> both = g1;
  both.insert(both.end(), g2.begin(), g2.end());
  const MeshPlan together = azure2::PlanMesh(anchors, both, Tails());
  std::vector<double> apart = azure2::PlanMesh(anchors, g1, Tails()).nodes;
  const std::vector<double> second = azure2::PlanMesh(anchors, g2, Tails()).nodes;
  apart.insert(apart.end(), second.begin(), second.end());
  std::sort(apart.begin(), apart.end());
  apart.erase(std::unique(apart.begin(), apart.end()), apart.end());
  Expect(SameBits(together.nodes, apart), "disjoint groups: the same nodes together and apart");
}

void Locality() {
  std::mt19937 rng(11u);
  int moved = 0;
  for (int c = 0; c < 100; c++) {
    std::vector<MeshAnchor> before = RandomAnchors(rng);
    if (before.empty()) continue;
    const std::vector<MeshSupport> supports = RandomSupports(rng, c % 2 == 0);
    std::vector<MeshAnchor> after = before;
    std::uniform_int_distribution<int> which(0, (int)before.size() - 1), how(0, 2), quanta(-6, 6);
    MeshAnchor &a = after[which(rng)];
    const int kind = how(rng);
    if (kind == 0) a.energy += 0.25 * a.particleWidth * quanta(rng);                  // moved by quarter widths
    if (kind == 1) a = Anchor(a.energy, a.particleWidth * 1.25 * 1.25, a.id);        // two notches wider
    if (kind == 2) a = Anchor(a.energy + 0.003, a.particleWidth / 1.25, a.id);       // moved far and narrower
    const MeshPlan p0 = azure2::PlanMesh(before, supports, Tails()), p1 = azure2::PlanMesh(after, supports, Tails());
    const auto region = azure2::MeshChangedRegion(before, after, supports, Tails());
    auto outside = [&](const std::vector<double> &nodes) {
      std::vector<double> out;
      for (double e : nodes) {
        bool in = false;
        for (const auto &r : region) in = in || (e >= r.first && e <= r.second);
        if (!in) out.push_back(e);
      }
      return out;
    };
    Expect(SameBits(outside(p0.nodes), outside(p1.nodes)), "locality " + std::to_string(c) +
                                                               ": the nodes outside the changed region are the same bits");
    if (!SameBits(p0.nodes, p1.nodes)) moved++;
  }
  std::printf("  locality: %d of the cases changed nodes inside their region\n", moved);
  // No change, no region.
  const std::vector<MeshAnchor> same = {Anchor(1.0, 1.0e-5, 1001)};
  Expect(azure2::MeshChangedRegion(same, same, {Support(0.95, 1.05, 100, 50.0, 20.0)}, Tails()).empty(),
         "locality: unchanged anchors change nothing");
}

void Ladder() {
  // The ladder step: a power of two, never above today's smooth step, above half of it.
  std::mt19937 rng(3u);
  std::uniform_real_distribution<double> logWidth(std::log(1.0e-5), std::log(1.0));
  std::uniform_int_distribution<int> points(1, 3000);
  for (int c = 0; c < 1000; c++) {
    const double w = std::exp(logWidth(rng));
    const MeshSupport s = Support(2.0, 2.0 + w, points(rng), 50.0, 20.0);
    const double l = azure2::MeshLadderStep(s), today = TodayStep(s, {}).smooth;
    int e = 0;
    Expect(std::frexp(l, &e) == 0.5 && l <= today && l > 0.5 * today, "ladder: a power of two in (S/2, S]");
  }
}

}  // namespace

int main() {
  std::printf("mesh planner: named cases\n");
  NamedCases();
  RandomCases();
  Shared();
  Permutations();
  Locality();
  Ladder();
  std::printf("mesh planner: %d checks, %d failures\n", checks, failures);
  return failures == 0 ? 0 : 1;
}
