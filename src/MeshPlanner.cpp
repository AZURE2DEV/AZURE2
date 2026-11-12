#include "MeshPlanner.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <tuple>

namespace azure2 {

namespace {

using Stretch = std::pair<double, double>;  // [low, high], MeV

/// Sorted, disjoint stretches (touching ones merged).
std::vector<Stretch> Merge(std::vector<Stretch> s) {
  std::sort(s.begin(), s.end());
  std::vector<Stretch> merged;
  for (const Stretch &x : s) {
    if (!merged.empty() && x.first <= merged.back().second)
      merged.back().second = std::max(merged.back().second, x.second);
    else
      merged.push_back(x);
  }
  return merged;
}

/// The parts of `all` not covered by `part` (both merged; part inside all).
std::vector<Stretch> Minus(const std::vector<Stretch> &all, const std::vector<Stretch> &part) {
  std::vector<Stretch> out;
  size_t j = 0;
  for (const Stretch &a : all) {
    double from = a.first;
    while (j < part.size() && part[j].second <= from) j++;
    size_t k = j;
    while (k < part.size() && part[k].first < a.second) {
      if (part[k].first > from) out.emplace_back(from, part[k].first);
      from = std::max(from, part[k].second);
      k++;
    }
    if (from < a.second) out.emplace_back(from, a.second);
  }
  return out;
}

/// Whether [lo, hi] shares a stretch of positive length with the merged list s.
bool Meets(const std::vector<Stretch> &s, double lo, double hi) {
  auto it = std::upper_bound(s.begin(), s.end(), lo, [](double e, const Stretch &x) { return e < x.second; });
  return it != s.end() && std::min(hi, it->second) > std::max(lo, it->first);
}

/// The settings of a class of supports, compared bit for bit.
using ClassKey = std::pair<double, double>;  // (points per width, core widths)

/// What a support contributes once its settings and step are worked out.
struct SupportFacts {
  bool exact = false;  // narrower than MeshZeroWidth: one exact node at `high`
  double ladder = 0.0;
  ClassKey key;
};

/// The offsets of an anchor's lattice from its energy, with the length of each interval: the core k * pitch for
/// k = 0 ... nCore, then the tail recurrence while its step is below `longest`.  offsets[0] = 0; step[k] is the
/// length of [offsets[k], offsets[k + 1]].
struct Offsets {
  std::vector<double> offset, step;
  double coreEdge = 0.0;
};

Offsets AnchorOffsets(double pitch, double coreWidths, double pointsPerWidth, double tail, double longest) {
  Offsets o;
  const long nCore = std::max(0L, (long)std::ceil(std::max(0.0, coreWidths) * pointsPerWidth - 1.0e-9));
  for (long k = 0; k <= nCore; k++) o.offset.push_back(k * pitch);
  for (long k = 0; k < nCore; k++) o.step.push_back(pitch);
  o.coreEdge = o.offset.back();
  double at = o.coreEdge;
  // The tail: at most a few hundred steps (log(longest / pitch) / tail); the bound only guards a degenerate input.
  for (int guard = 0; guard < 10000000; guard++) {
    const double step = std::max(pitch, tail * at);
    if (!(step < longest)) break;
    o.step.push_back(step);
    at = at + step;
    o.offset.push_back(at);
  }
  return o;
}

/// One lattice: an anchor's for one class of settings, or the ladder of one smooth step.
struct Lattice {
  bool ladder = false;
  double center = 0.0, pitch = 0.0, tail = 0.0, coreEdge = 0.0, reach = 0.0;  // an anchor's
  double step = 0.0;                                                          // a ladder's
  const std::vector<Stretch> *covered = nullptr;  // where the lattice is active
  const std::vector<Stretch> *gaps = nullptr;     // the parts of the supports' union where it is not

  /// Its step at E: +inf where it has no node (outside an anchor's reach); the ladder's is checked by coverage.
  double StepAt(double e) const {
    if (ladder) return step;
    const double d = std::fabs(e - center);
    if (d > reach) return std::numeric_limits<double>::infinity();
    if (d <= coreEdge) return pitch;
    return std::max(pitch, tail * d);
  }
};

/// One interval of a lattice, a candidate for the mesh.
struct Piece {
  double lo, hi, length;
  int lattice;
};

/// The classes of settings and ladder steps of the supports, with their coverage.
struct Coverage {
  std::vector<SupportFacts> facts;
  std::vector<Stretch> all;  // the union of the supports (zero-width ones excluded)
  std::map<ClassKey, std::vector<Stretch>> classCovered, classGaps;
  std::map<ClassKey, double> classLongest;  // the largest ladder step of the class's supports
  std::map<double, std::vector<Stretch>> ladderCovered, ladderGaps;
  double longest = 0.0;  // the largest ladder step of all
};

Coverage Cover(const std::vector<MeshSupport> &supports) {
  Coverage c;
  std::map<ClassKey, std::vector<Stretch>> byClass;
  std::map<double, std::vector<Stretch>> byLadder;
  for (const MeshSupport &s : supports) {
    SupportFacts f;
    f.key = ClassKey(s.settings.pointsPerWidth, s.settings.coreWidths);
    f.ladder = MeshLadderStep(s);
    f.exact = !(f.ladder > 0.0);
    c.facts.push_back(f);
    if (f.exact) continue;
    const Stretch x(std::min(s.low, s.high), std::max(s.low, s.high));
    c.all.push_back(x);
    byClass[f.key].push_back(x);
    byLadder[f.ladder].push_back(x);
    double &longest = c.classLongest[f.key];
    longest = std::max(longest, f.ladder);
    c.longest = std::max(c.longest, f.ladder);
  }
  c.all = Merge(c.all);
  for (auto &kv : byClass) {
    c.classCovered[kv.first] = Merge(kv.second);
    c.classGaps[kv.first] = Minus(c.all, c.classCovered[kv.first]);
  }
  for (auto &kv : byLadder) {
    c.ladderCovered[kv.first] = Merge(kv.second);
    c.ladderGaps[kv.first] = Minus(c.all, c.ladderCovered[kv.first]);
  }
  return c;
}

/// Bitwise order of anchors for the comparison of two anchor lists (the lattice depends on energy and width only).
bool AnchorLess(const MeshAnchor &a, const MeshAnchor &b) {
  return std::make_tuple(a.energy, a.particleWidth) < std::make_tuple(b.energy, b.particleWidth);
}

}  // namespace

double MeshLadderStep(const MeshSupport &support) {
  const double range = std::fabs(support.high - support.low);
  if (!(range >= MeshZeroWidth())) return 0.0;
  // Today's smooth step (AdaptiveIntegrationGrid::GenerateGrid): the window in ceil(window / step) equal parts.
  double base = support.smoothStep;
  if (!(base > 0.0) || base > range) base = range;
  int n = (int)std::ceil(range / base - 1.0e-9);
  if (n < 1) n = 1;
  const double step = range / n;
  int exponent = 0;
  std::frexp(step, &exponent);  // step = f 2^exponent, f in [0.5, 1)
  return std::ldexp(1.0, exponent - 1);
}

double MeshTailFactor(const MeshSettings &settings, const MeshTailConstants &tails) {
  const double core = std::max(tails.convolutionCoreWidths, settings.coreWidths);
  return tails.tailRatio / settings.pointsPerWidth * std::sqrt(core / tails.convolutionCoreWidths);
}

MeshPlan PlanMesh(const std::vector<MeshAnchor> &anchors, const std::vector<MeshSupport> &supports,
                  const MeshTailConstants &tails) {
  const Coverage cover = Cover(supports);
  std::vector<Lattice> lattices;
  std::vector<Piece> pieces;
  std::vector<double> nodes;

  // The ladders: n L for the n whose intervals meet a support of that step.
  for (const auto &kv : cover.ladderCovered) {
    Lattice x;
    x.ladder = true;
    x.step = kv.first;
    x.covered = &kv.second;
    x.gaps = &cover.ladderGaps.at(kv.first);
    const int index = (int)lattices.size();
    lattices.push_back(x);
    double last = -std::numeric_limits<double>::infinity();
    for (const Stretch &s : kv.second) {
      double n = std::max(std::floor(s.first / x.step), std::floor(last / x.step));
      const double nEnd = std::ceil(s.second / x.step);
      for (; n < nEnd; n += 1.0) pieces.push_back({n * x.step, (n + 1.0) * x.step, x.step, index});
      last = nEnd * x.step;
    }
  }

  // The anchors' lattices, one per class of settings whose ladder is coarser than the anchor's pitch.
  for (const MeshAnchor &a : anchors) {
    if (!(a.particleWidth > 0.0)) continue;
    for (const auto &kv : cover.classCovered) {
      const ClassKey &key = kv.first;
      const MeshSettings settings{key.first, key.second};
      if (!(settings.pointsPerWidth > 0.0)) continue;
      const double longest = cover.classLongest.at(key);
      Lattice x;
      x.center = a.energy;
      x.pitch = a.particleWidth / settings.pointsPerWidth;
      if (!(x.pitch < longest)) continue;
      x.tail = MeshTailFactor(settings, tails);
      const Offsets o = AnchorOffsets(x.pitch, settings.coreWidths, settings.pointsPerWidth, x.tail, longest);
      x.coreEdge = o.coreEdge;
      x.reach = o.offset.back();
      if (!Meets(kv.second, x.center - x.reach, x.center + x.reach)) continue;
      x.covered = &kv.second;
      x.gaps = &cover.classGaps.at(key);
      const int index = (int)lattices.size();
      lattices.push_back(x);
      std::vector<double> above(o.offset.size()), below(o.offset.size());
      for (size_t k = 0; k < o.offset.size(); k++) {
        above[k] = x.center + o.offset[k];
        below[k] = x.center - o.offset[k];
      }
      for (size_t k = 0; k < o.step.size(); k++) {
        pieces.push_back({above[k], above[k + 1], o.step[k], index});
        pieces.push_back({below[k + 1], below[k], o.step[k], index});
      }
    }
  }

  // Keep a piece unless its lattice is inactive on it or a single other lattice, active over all of the piece's part
  // inside the supports, is finer than its length everywhere on it.  A sweep in E finds the anchors' lattices whose
  // reach spans a piece.
  std::vector<int> byStart;
  for (int i = 0; i < (int)lattices.size(); i++)
    if (!lattices[i].ladder) byStart.push_back(i);
  std::sort(byStart.begin(), byStart.end(), [&](int i, int j) {
    return lattices[i].center - lattices[i].reach < lattices[j].center - lattices[j].reach;
  });
  std::vector<int> ladders;
  for (int i = 0; i < (int)lattices.size(); i++)
    if (lattices[i].ladder) ladders.push_back(i);
  std::sort(pieces.begin(), pieces.end(),
            [](const Piece &a, const Piece &b) { return std::tie(a.lo, a.hi, a.length, a.lattice) <
                                                        std::tie(b.lo, b.hi, b.length, b.lattice); });
  std::vector<int> alive;
  size_t next = 0;
  const double margin = 1.0 - 1.0e-12;  // a lattice dominates only if finer beyond rounding
  for (const Piece &p : pieces) {
    while (next < byStart.size() && lattices[byStart[next]].center - lattices[byStart[next]].reach <= p.lo)
      alive.push_back(byStart[next++]);
    alive.erase(std::remove_if(alive.begin(), alive.end(),
                               [&](int i) { return lattices[i].center + lattices[i].reach < p.lo; }),
                alive.end());
    const Lattice &own = lattices[p.lattice];
    if (!(p.hi > p.lo) || !Meets(*own.covered, p.lo, p.hi)) continue;
    bool dominated = false;
    auto check = [&](int y) {
      if (dominated || y == p.lattice) return;
      const Lattice &other = lattices[y];
      const double finest = std::max(other.StepAt(p.lo), other.StepAt(p.hi));  // the step is quasi-convex in E
      if (!(finest < p.length * margin)) return;
      if (Meets(*other.gaps, p.lo, p.hi)) return;
      dominated = true;
    };
    for (int y : ladders) check(y);
    for (int y : alive)
      if (lattices[y].center + lattices[y].reach >= p.hi) check(y);
    if (dominated) continue;
    nodes.push_back(p.lo);
    nodes.push_back(p.hi);
  }

  // Supports narrower than MeshZeroWidth: their exact high edge.
  for (size_t i = 0; i < supports.size(); i++)
    if (cover.facts[i].exact) nodes.push_back(std::max(supports[i].low, supports[i].high));

  std::sort(nodes.begin(), nodes.end());
  nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());

  MeshPlan plan;
  plan.nodes = std::move(nodes);
  for (size_t i = 0; i < supports.size(); i++) {
    const double low = std::min(supports[i].low, supports[i].high), high = std::max(supports[i].low, supports[i].high);
    MeshNodeRange r;
    if (cover.facts[i].exact) {
      r.first = r.last = (size_t)(std::lower_bound(plan.nodes.begin(), plan.nodes.end(), high) - plan.nodes.begin());
    } else {
      const auto above = std::upper_bound(plan.nodes.begin(), plan.nodes.end(), low);
      r.first = above == plan.nodes.begin() ? 0 : (size_t)(above - plan.nodes.begin()) - 1;
      r.last = (size_t)(std::lower_bound(plan.nodes.begin(), plan.nodes.end(), high) - plan.nodes.begin());
      if (r.last >= plan.nodes.size() && !plan.nodes.empty()) r.last = plan.nodes.size() - 1;
    }
    plan.ranges.push_back(r);
  }
  return plan;
}

std::vector<std::pair<double, double>> MeshChangedRegion(const std::vector<MeshAnchor> &before,
                                                         const std::vector<MeshAnchor> &after,
                                                         const std::vector<MeshSupport> &supports,
                                                         const MeshTailConstants &tails) {
  std::vector<MeshAnchor> a = before, b = after, changed;
  std::sort(a.begin(), a.end(), AnchorLess);
  std::sort(b.begin(), b.end(), AnchorLess);
  std::set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(changed), AnchorLess);
  const Coverage cover = Cover(supports);
  std::vector<Stretch> region;
  for (const MeshAnchor &r : changed) {
    if (!(r.particleWidth > 0.0)) continue;
    for (const auto &kv : cover.classLongest) {
      const MeshSettings settings{kv.first.first, kv.first.second};
      if (!(settings.pointsPerWidth > 0.0)) continue;
      const double pitch = r.particleWidth / settings.pointsPerWidth;
      if (!(pitch < kv.second)) continue;
      const Offsets o = AnchorOffsets(pitch, settings.coreWidths, settings.pointsPerWidth,
                                      MeshTailFactor(settings, tails), kv.second);
      const double reach = o.offset.back() + cover.longest;
      region.emplace_back(r.energy - reach, r.energy + reach);
    }
  }
  return Merge(region);
}

}  // namespace azure2
