// The calculation units of a working copy (EData::Units, docs/upgrade/data_structure.md section 1.4): sub-points
// given to points during the set-up or a grid re-placement, and the layout of the one unit vector.
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "EData.h"

void EData::AddSubPoint(EPoint &point, EPoint sub) {
  const EPoint *base = units_.data();
  if (&point < base || &point >= base + units_.size())
    throw std::logic_error("EData::AddSubPoint: the point is not a unit of this data set");
  point.MakeSubPoint(sub);
  pendingSubPoints_.emplace_back(static_cast<int>(&point - base), std::move(sub));
}

void EData::CommitSubPoints() {
  // each point's sub-points in the order given, the points in the order of their first sub-point
  SubPointGrids grids(1);
  std::vector<int> slot(units_.size(), -1);
  for (auto &p : pendingSubPoints_) {
    if (slot[p.first] < 0) slot[p.first] = static_cast<int>(grids[0].size()), grids[0].push_back({p.first, {}});
    grids[0][slot[p.first]].second.push_back(std::move(p.second));
  }
  pendingSubPoints_.clear();
  PlaceSubPointGrids(grids);
}

void EData::PlaceSubPointGrids(SubPointGrids &grids) {
  for (auto &item : grids)
    for (auto &grid : item) {
      const int first = NumUnits();
      for (EPoint &s : grid.second) units_.push_back(std::move(s));
      units_[grid.first].SetSubPointRange(first, static_cast<int>(grid.second.size()));
    }
  CompactUnits();
}

void EData::CompactUnits() {
  std::vector<EPoint> next;
  next.reserve(units_.size());
  // the points of every term, in data order
  for (ESegment &segment : segments_)
    for (int t = 0; t < segment.NumTerms(); t++) {
      const int first = static_cast<int>(next.size());
      for (EPoint &p : segment.TermPoints(*this, t)) next.push_back(std::move(p));
      segment.SetTermRange(t, first);
    }
  // then their sub-points, in the same order
  const int numPoints = static_cast<int>(next.size());
  for (int i = 0; i < numPoints; i++) {
    const int from = next[i].SubPointFirst(), count = next[i].NumSubPoints();
    const int first = static_cast<int>(next.size());
    for (int j = 0; j < count; j++) next.push_back(std::move(units_[from + j]));
    next[i].SetSubPointRange(first, count);
  }
  units_.swap(next);
}

std::vector<azure2::PointRef> EData::PlanUnits() const {
  std::vector<azure2::PointRef> units;
  for (const ESegment &segment : segments_)
    for (int p = 0; p < segment.NumPoints(); p++) units.push_back(segment.TermFirstUnit(0) + p);
  for (const LeafRef &leaf : leafRefs_) {
    const ESegment &segment = segments_[leaf.first];
    for (int p = 0; p < segment.NumTermPoints(leaf.second); p++) units.push_back(segment.TermFirstUnit(leaf.second) + p);
  }
  return units;
}
