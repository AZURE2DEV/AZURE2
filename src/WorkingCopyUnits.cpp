// The points and slots of a working copy (WorkingCopy::UnitAt, SlotAt; docs/upgrade/data_structure.md section 1.4):
// slots given to points during the set-up or a grid re-placement, and the layout of the two vectors.
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "WorkingCopy.h"

void WorkingCopy::AddSubPoint(Unit &point, Site sub) {
  const Unit *base = units_.data();
  if (&point < base || &point >= base + units_.size())
    throw std::logic_error("WorkingCopy::AddSubPoint: the point is not a unit of this data set");
  point.MakeSubPoint(sub);
  const std::size_t unit = static_cast<std::size_t>(&point - base);
  if (pendingSubPoints_.size() <= unit) pendingSubPoints_.resize(units_.size());
  pendingSubPoints_[unit].push_back(std::move(sub));
}

void WorkingCopy::CommitSubPoints() {
  std::vector<std::vector<Site> *> fresh(units_.size(), nullptr);
  for (std::size_t u = 0; u < pendingSubPoints_.size(); u++)
    if (!pendingSubPoints_[u].empty()) fresh[u] = &pendingSubPoints_[u];
  LayOutUnits(&fresh);
  std::vector<std::vector<Site>>().swap(pendingSubPoints_);
}

void WorkingCopy::PlaceSubPointGrids(SubPointGrids &grids) {
  // A grid with as many nodes as the one it replaces takes its place; any other is a new range, and the slots are
  // laid out again.
  std::vector<std::vector<Site> *> fresh(units_.size(), nullptr);
  bool layOut = false;
  for (auto &item : grids)
    for (auto &grid : item) {
      Unit &point = units_[grid.first];
      if ((int)grid.second.size() == point.NumSubPoints()) {
        for (int j = 0; j < point.NumSubPoints(); j++) slots_[point.SubPointFirst() + j] = std::move(grid.second[j]);
      } else {
        fresh[grid.first] = &grid.second;
        layOut = true;
      }
    }
  if (layOut) LayOutUnits(&fresh);
}

void WorkingCopy::CompactUnits() { LayOutUnits(nullptr); }

// The one layout: the points of every term in data order in the units, and their slots in the same order in the
// slots -- a point's new slots from `fresh` (by its present unit) where it has some, else its present ones -- each
// moved once into a vector of exactly the size needed.
void WorkingCopy::LayOutUnits(std::vector<std::vector<Site> *> *fresh) {
  auto newSubs = [&](int unit) { return fresh && (*fresh)[unit] ? (*fresh)[unit] : nullptr; };
  std::size_t points = 0, slots = 0;
  for (const Segment &segment : set_->segments_)
    for (int t = 0; t < segment.NumTerms(); t++)
      for (int p = 0; p < segment.NumTermPoints(t); p++) {
        const int unit = segment.TermFirstUnit(*this, t) + p;
        points++;
        slots += newSubs(unit) ? newSubs(unit)->size() : (std::size_t)units_[unit].NumSubPoints();
      }
  std::vector<Unit> next;
  next.reserve(points);
  std::vector<int> from;  // the present unit of each point placed
  for (const Segment &segment : set_->segments_)
    for (int t = 0; t < segment.NumTerms(); t++) {
      const int first = static_cast<int>(next.size());
      for (int p = 0; p < segment.NumTermPoints(t); p++) {
        from.push_back(segment.TermFirstUnit(*this, t) + p);
        next.push_back(std::move(units_[from.back()]));
      }
      State(segment).terms[t].first = first;
    }
  std::vector<Site> nextSlots;
  nextSlots.reserve(slots);
  for (std::size_t i = 0; i < from.size(); i++) {
    const int first = static_cast<int>(nextSlots.size());
    if (std::vector<Site> *subs = newSubs(from[i])) {
      for (Site &sub : *subs) nextSlots.push_back(std::move(sub));
      next[i].SetSubPointRange(first, static_cast<int>(subs->size()));
    } else {
      const int old = next[i].SubPointFirst(), count = next[i].NumSubPoints();
      for (int j = 0; j < count; j++) nextSlots.push_back(std::move(slots_[old + j]));
      next[i].SetSubPointRange(first, count);
    }
  }
  units_.swap(next);
  slots_.swap(nextSlots);
}

std::vector<azure2::PointRef> WorkingCopy::PlanUnits() const {
  std::vector<azure2::PointRef> units;
  for (const Segment &segment : set_->segments_)
    for (int p = 0; p < segment.NumPoints(); p++) units.push_back(segment.TermFirstUnit(*this, 0) + p);
  for (const LeafRef &leaf : set_->leafRefs_) {
    const Segment &segment = set_->segments_[leaf.first];
    for (int p = 0; p < segment.NumTermPoints(leaf.second); p++) units.push_back(segment.TermFirstUnit(*this, leaf.second) + p);
  }
  return units;
}
