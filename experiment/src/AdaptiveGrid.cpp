// The refinement of an adaptive extrapolation grid (AdaptiveGrid.h).
#include "azure2/experiment/AdaptiveGrid.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace azure2 {
namespace experiment {

namespace {

struct Split {
  double deviation;  ///< the largest relative deviation over the components
  double a, b, m;
  std::vector<double> value;  ///< the model at m
};

// The largest |value| of each component over the grid and the batch.
std::vector<double> Floors(const std::map<double, std::vector<double>> &grid, const std::vector<std::vector<double>> &batch,
                           std::size_t n) {
  std::vector<double> top(n, 0.0);
  auto add = [&](const std::vector<double> &v) {
    for (std::size_t c = 0; c < n; c++)
      if (std::isfinite(v[c])) top[c] = std::max(top[c], std::fabs(v[c]));
  };
  for (const auto &g : grid) add(g.second);
  for (const std::vector<double> &v : batch) add(v);
  for (double &t : top) t *= kFloorFraction;
  return top;
}

std::vector<std::vector<double>> Evaluate(const GridBatchModel &model, const std::vector<double> &x, std::size_t &n) {
  std::vector<std::vector<double>> v = model(x);
  if (v.size() != x.size()) throw std::runtime_error("RefineGrid: the model returned a batch of another length");
  for (const std::vector<double> &r : v) {
    if (n == 0) n = r.size();
    if (r.size() != n || n == 0) throw std::runtime_error("RefineGrid: the model's values have different lengths");
  }
  return v;
}

}  // namespace

AdaptiveGridResult RefineGrid(std::vector<double> start, double tolerance, long maxPoints, const GridBatchModel &model) {
  AdaptiveGridResult out;
  std::sort(start.begin(), start.end());
  start.erase(std::unique(start.begin(), start.end()), start.end());
  if (start.empty()) return out;
  std::size_t n = 0;  // components of the model's value
  std::map<double, std::vector<double>> grid;
  {
    std::vector<std::vector<double>> v = Evaluate(model, start, n);
    out.evaluations += (long)start.size();
    for (std::size_t i = 0; i < start.size(); i++) grid.emplace(start[i], std::move(v[i]));
  }
  std::vector<std::pair<double, double>> test;  // the intervals of this sweep, ascending
  for (auto it = grid.begin(); std::next(it) != grid.end(); ++it) test.push_back({it->first, std::next(it)->first});
  while (true) {
    std::vector<std::pair<double, double>> intervals;
    std::vector<double> mids;
    for (const auto &ab : test) {
      const double a = ab.first, b = ab.second, m = 0.5 * (a + b);
      if (!(m > a && m < b) || b - a <= 1e-12 * std::max(std::fabs(a), std::fabs(b))) continue;
      intervals.push_back(ab);
      mids.push_back(m);
    }
    if (mids.empty()) break;
    out.sweeps++;
    const std::vector<std::vector<double>> fm = Evaluate(model, mids, n);
    out.evaluations += (long)mids.size();
    const std::vector<double> floors = Floors(grid, fm, n);
    std::vector<Split> splits;
    for (std::size_t i = 0; i < mids.size(); i++) {
      const std::vector<double> &fa = grid.at(intervals[i].first), &fb = grid.at(intervals[i].second);
      double worst = 0.0;
      for (std::size_t c = 0; c < n; c++) {
        if (!std::isfinite(fa[c]) || !std::isfinite(fb[c]) || !std::isfinite(fm[i][c])) continue;
        const double deviation = std::fabs(fm[i][c] - 0.5 * (fa[c] + fb[c]));
        const double scale = std::max(std::fabs(fm[i][c]), floors[c]);
        const double relative = scale > 0.0 ? deviation / scale : (deviation > 0.0 ? HUGE_VAL : 0.0);
        worst = std::max(worst, relative);
      }
      if (worst > tolerance) splits.push_back({worst, intervals[i].first, intervals[i].second, mids[i], fm[i]});
    }
    if (splits.empty()) break;
    const long room = std::max(0L, maxPoints - (long)grid.size());
    if ((long)splits.size() > room) {
      // The largest deviations first; the order of equal ones by energy (the splits are in ascending order).
      std::stable_sort(splits.begin(), splits.end(), [](const Split &p, const Split &q) { return p.deviation > q.deviation; });
      out.budgetReached = true;
      out.largestDeviation = splits[room].deviation;
      out.largestDeviationAt = splits[room].m;
      splits.resize(room);
    }
    test.clear();
    for (Split &s : splits) {
      grid.emplace(s.m, std::move(s.value));
      test.push_back({s.a, s.m});
      test.push_back({s.m, s.b});
    }
    std::sort(test.begin(), test.end());
    if (out.budgetReached) break;
  }
  for (auto &g : grid) {
    out.x.push_back(g.first);
    out.values.push_back(std::move(g.second));
  }
  return out;
}

}  // namespace experiment
}  // namespace azure2
