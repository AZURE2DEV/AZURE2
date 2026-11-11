// The calculation units of a working copy (EData::Units, docs/upgrade/data_structure.md section 1.4): sub-points
// given to points during the set-up or a grid re-placement, and the layout of the one unit vector.
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "WorkingCopy.h"

void EData::AddSubPoint(EPoint &point, EPoint sub) {
  const EPoint *base = units_.data();
  if (&point < base || &point >= base + units_.size())
    throw std::logic_error("EData::AddSubPoint: the point is not a unit of this data set");
  point.MakeSubPoint(sub);
  const std::size_t unit = static_cast<std::size_t>(&point - base);
  if (pendingSubPoints_.size() <= unit) pendingSubPoints_.resize(units_.size());
  pendingSubPoints_[unit].push_back(std::move(sub));
}

void EData::CommitSubPoints() {
  std::vector<std::vector<EPoint> *> fresh(units_.size(), nullptr);
  for (std::size_t u = 0; u < pendingSubPoints_.size(); u++)
    if (!pendingSubPoints_[u].empty()) fresh[u] = &pendingSubPoints_[u];
  LayOutUnits(&fresh);
  std::vector<std::vector<EPoint>>().swap(pendingSubPoints_);
}

void EData::PlaceSubPointGrids(SubPointGrids &grids) {
  // A grid with as many nodes as the one it replaces takes its place; any other is a new range, and the units are
  // laid out again.
  std::vector<std::vector<EPoint> *> fresh(units_.size(), nullptr);
  bool layOut = false;
  for (auto &item : grids)
    for (auto &grid : item) {
      EPoint &point = units_[grid.first];
      if ((int)grid.second.size() == point.NumSubPoints()) {
        for (int j = 0; j < point.NumSubPoints(); j++)
          units_[point.SubPointFirst() + j] = std::move(grid.second[j]);
      } else {
        fresh[grid.first] = &grid.second;
        layOut = true;
      }
    }
  if (layOut) LayOutUnits(&fresh);
}

void EData::CompactUnits() { LayOutUnits(nullptr); }

// The one layout of the units: the points of every term in data order, then their sub-points in the same order --
// a point's new sub-points from `fresh` (by its present unit) where it has some, else its present ones -- each unit
// moved once into a vector of exactly the size needed.
void EData::LayOutUnits(std::vector<std::vector<EPoint> *> *fresh) {
  auto newSubs = [&](int unit) { return fresh && (*fresh)[unit] ? (*fresh)[unit] : nullptr; };
  std::size_t size = 0;
  for (const ESegment &segment : set_->segments_)
    for (int t = 0; t < segment.NumTerms(); t++)
      for (int p = 0; p < segment.NumTermPoints(t); p++) {
        const int unit = segment.TermFirstUnit(*this, t) + p;
        size += 1 + (newSubs(unit) ? newSubs(unit)->size() : (std::size_t)units_[unit].NumSubPoints());
      }
  std::vector<EPoint> next;
  next.reserve(size);
  std::vector<int> from;  // the present unit of each point placed
  for (const ESegment &segment : set_->segments_)
    for (int t = 0; t < segment.NumTerms(); t++) {
      const int first = static_cast<int>(next.size());
      for (int p = 0; p < segment.NumTermPoints(t); p++) {
        from.push_back(segment.TermFirstUnit(*this, t) + p);
        next.push_back(std::move(units_[from.back()]));
      }
      State(segment).terms[t].first = first;
    }
  for (std::size_t i = 0; i < from.size(); i++) {
    const int first = static_cast<int>(next.size());
    if (std::vector<EPoint> *subs = newSubs(from[i])) {
      for (EPoint &sub : *subs) next.push_back(std::move(sub));
      next[i].SetSubPointRange(first, static_cast<int>(subs->size()));
    } else {
      const int old = next[i].SubPointFirst(), count = next[i].NumSubPoints();
      for (int j = 0; j < count; j++) next.push_back(std::move(units_[old + j]));
      next[i].SetSubPointRange(first, count);
    }
  }
  units_.swap(next);
}

std::vector<azure2::PointRef> EData::PlanUnits() const {
  std::vector<azure2::PointRef> units;
  for (const ESegment &segment : set_->segments_)
    for (int p = 0; p < segment.NumPoints(); p++) units.push_back(segment.TermFirstUnit(*this, 0) + p);
  for (const LeafRef &leaf : set_->leafRefs_) {
    const ESegment &segment = set_->segments_[leaf.first];
    for (int p = 0; p < segment.NumTermPoints(leaf.second); p++) units.push_back(segment.TermFirstUnit(*this, leaf.second) + p);
  }
  return units;
}
