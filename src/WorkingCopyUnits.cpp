// The points and slots of a working copy (WorkingCopy::UnitAt, SlotAt; docs/upgrade/data_structure.md section 1.4):
// slots given to points during the set-up or a grid re-placement, and the layout of the two vectors.
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cstdint>
#include <cstring>
#include <iostream>

#include "MeshSlots.h"
#include "WorkingCopy.h"

void WorkingCopy::AddSubPoint(Unit &point, Site sub) {
  const Unit *base = units_.data();
  if (&point < base || &point >= base + units_.size())
    throw std::logic_error("WorkingCopy::AddSubPoint: the point is not a unit of this data set");
  point.MakeSubPoint(sub);
  const std::size_t unit = static_cast<std::size_t>(&point - base);
  if (pendingOwnBlock_.size() <= unit) pendingOwnBlock_.resize(units_.size(), -1);
  if (pendingOwnBlock_[unit] < 0) {
    pendingOwnBlock_[unit] = static_cast<int>(pendingBlocks_.size());
    pendingBlocks_.emplace_back();
    pendingBlocks_.back().members.push_back({static_cast<int>(unit), 0, 0});
  }
  azure2::SlotBlock &block = pendingBlocks_[pendingOwnBlock_[unit]];
  block.sites.push_back(std::move(sub));
  block.members[0].count++;
}

void WorkingCopy::AddSlotGroup(azure2::SlotBlock block) {
  for (const azure2::SlotMember &m : block.members)
    if (m.unit < 0 || m.unit >= NumUnits() || m.first < 0 || m.count < 0 || m.first + m.count > (int)block.sites.size())
      throw std::logic_error("WorkingCopy::AddSlotGroup: a member is no unit of this data set or lies outside the block");
  pendingBlocks_.push_back(std::move(block));
}

void WorkingCopy::CommitSubPoints() {
  LayOutUnits(&pendingBlocks_);
  std::vector<azure2::SlotBlock>().swap(pendingBlocks_);
  std::vector<int>().swap(pendingOwnBlock_);
}

void WorkingCopy::PlaceSubPointGrids(SubPointGrids &grids) {
  // A block whose members keep their places in the block they read now, with as many sites, takes that block's place;
  // any other is laid out anew with the units.
  std::vector<azure2::SlotBlock> fresh;
  for (auto &item : grids)
    for (azure2::SlotBlock &block : item) {
      bool inPlace = !block.members.empty();
      int first = -1;
      for (const azure2::SlotMember &m : block.members) {
        const Unit &point = units_[m.unit];
        if (!inPlace) break;
        if (point.NumSubPoints() == 0 || point.NumSubPoints() != m.count) {
          inPlace = false;
          break;
        }
        const int b = slotBlock_[point.SubPointFirst()];
        if (first < 0) first = b;
        inPlace = b == first && point.SubPointFirst() - b == m.first;
      }
      if (inPlace) {
        const int size = (int)block.sites.size();
        inPlace = first + size <= NumSlots() && slotBlock_[first + size - 1] == first &&
                  (first + size == NumSlots() || slotBlock_[first + size] != first);
      }
      if (inPlace) {
        for (std::size_t j = 0; j < block.sites.size(); j++) slots_[first + j] = std::move(block.sites[j]);
      } else {
        fresh.push_back(std::move(block));
      }
    }
  if (!fresh.empty()) LayOutUnits(&fresh);
}

void WorkingCopy::CompactUnits() { LayOutUnits(nullptr); }

void WorkingCopy::SeparateSlotsByAngle(const std::vector<int> &units) {
  auto same = [](double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; };
  bool any = false;
  for (int u : units) {
    Unit &point = units_[u];
    bool differ = false;
    for (int i = 1; i <= point.NumSubPoints() && !differ; i++) {
      const Site &slot = *point.GetSubPoint(*this, i);
      differ = !same(slot.GetLabAngle(), point.GetLabAngle()) || !same(slot.GetCMAngle(), point.GetCMAngle());
    }
    if (!differ) continue;
    for (int i = 1; i <= point.NumSubPoints(); i++) {
      Site copy = *point.GetSubPoint(*this, i);
      copy.SetLabAngle(point.GetLabAngle());
      copy.SetCMAngle(point.GetCMAngle());
      AddSubPoint(point, std::move(copy));
    }
    any = true;
  }
  if (any) CommitSubPoints();
}

void WorkingCopy::SetSlotRange(int unit, int first, int count) {
  units_[unit].SetSubPointRange(first, count);
  for (int k = first; k < first + count; k++) {
    if (slotWriter_[k] < 0) slotWriter_[k] = unit;
    slotReaders_[k]++;
  }
}

// The one layout: the points of every term in data order in the units, and their blocks in the same order in the
// slots -- a point's new block from `fresh` where it is a member of one, else the block it reads now -- each block
// moved once into a vector of exactly the size needed (a block no point reads is dropped).  Then the writer of each
// slot: the first point, in this order, whose range holds it.
void WorkingCopy::LayOutUnits(std::vector<azure2::SlotBlock> *fresh) {
  std::vector<std::pair<int, int>> member(units_.size(), {-1, -1});  // per present unit: (fresh block, member)
  if (fresh)
    for (std::size_t b = 0; b < fresh->size(); b++)
      for (std::size_t m = 0; m < (*fresh)[b].members.size(); m++) member[(*fresh)[b].members[m].unit] = {(int)b, (int)m};
  std::size_t points = 0;
  for (const Segment &segment : set_->segments_)
    for (int t = 0; t < segment.NumTerms(); t++) points += segment.NumTermPoints(t);
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
  // The slots kept: every fresh block with a member placed, every present block a point outside them reads.
  std::size_t slots = 0;
  {
    std::vector<char> freshKept(fresh ? fresh->size() : 0, 0), blockKept(slots_.size(), 0);
    for (std::size_t i = 0; i < from.size(); i++) {
      const int b = member[from[i]].first;
      if (b >= 0) {
        if (!freshKept[b]) freshKept[b] = 1, slots += (*fresh)[b].sites.size();
        continue;
      }
      const Unit &point = next[i];
      for (int k = point.SubPointFirst(); k < point.SubPointFirst() + point.NumSubPoints(); k++)
        if (!blockKept[slotBlock_[k]]) {
          blockKept[slotBlock_[k]] = 1;
          for (int j = slotBlock_[k]; j < (int)slots_.size() && slotBlock_[j] == slotBlock_[k]; j++) slots++;
        }
    }
  }
  std::vector<Site> nextSlots;
  std::vector<int> nextBlock;
  nextSlots.reserve(slots), nextBlock.reserve(slots);
  std::vector<int> freshAt(fresh ? fresh->size() : 0, -1);  // the new first slot of each fresh block placed
  std::vector<int> movedTo(slots_.size(), -1);                // the new index of each present slot moved
  for (std::size_t i = 0; i < from.size(); i++) {
    Unit &point = next[i];
    const std::pair<int, int> m = member[from[i]];
    if (m.first >= 0) {
      azure2::SlotBlock &block = (*fresh)[m.first];
      if (freshAt[m.first] < 0) {
        freshAt[m.first] = static_cast<int>(nextSlots.size());
        for (Site &site : block.sites) nextSlots.push_back(std::move(site)), nextBlock.push_back(freshAt[m.first]);
      }
      point.SetSubPointRange(freshAt[m.first] + block.members[m.second].first, block.members[m.second].count);
      continue;
    }
    const int old = point.SubPointFirst(), count = point.NumSubPoints();
    if (count == 0) {
      point.SetSubPointRange(static_cast<int>(nextSlots.size()), 0);
      continue;
    }
    for (int k = old; k < old + count; k++) {
      if (movedTo[k] >= 0) continue;
      const int b = slotBlock_[k], at = static_cast<int>(nextSlots.size());
      for (int j = b; j < (int)slots_.size() && slotBlock_[j] == b; j++) {
        movedTo[j] = static_cast<int>(nextSlots.size());
        nextSlots.push_back(std::move(slots_[j]));
        nextBlock.push_back(at);
      }
    }
    for (int k = old; k < old + count; k++)
      if (movedTo[k] != movedTo[old] + (k - old))
        throw std::logic_error("WorkingCopy::LayOutUnits: the slots of a point do not stay contiguous");
    point.SetSubPointRange(movedTo[old], count);
  }
  units_.swap(next);
  slots_.swap(nextSlots);
  slotBlock_.swap(nextBlock);
  slotWriter_.assign(slots_.size(), -1);
  slotReaders_.assign(slots_.size(), 0);
  for (int u = 0; u < (int)units_.size(); u++) {
    const Unit &point = units_[u];
    for (int k = point.SubPointFirst(); k < point.SubPointFirst() + point.NumSubPoints(); k++) {
      if (slotWriter_[k] < 0) slotWriter_[k] = u;
      slotReaders_[k]++;
    }
  }
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

long WorkingCopy::CheckSlotGroups() const {
  long reads = 0, outside = 0, writers = 0, otherGroup = 0, offPlan = 0;
  std::vector<std::string> reports;
  auto report = [&reports](const std::string &r) {
    if (reports.size() < 10) reports.push_back(r);
  };
  auto bits = [](double x) {
    std::uint64_t b;
    std::memcpy(&b, &x, sizeof b);
    return b;
  };
  // What the points that read one block have in common (SlotGroup.h), compared with the block's first reader.
  auto key = [&](const Unit &p) {
    const Unit::SubGridSpec &spec = Data().GridSpecOf(p);
    return std::vector<std::uint64_t>{(std::uint64_t)p.GetSegmentKey(), bits(p.GetCMAngle()),
                                      (std::uint64_t)(p.IsCoefficient() ? p.GetCoefficientOrder() : -1),
                                      (std::uint64_t)p.GetTargetEffectNum(), bits(p.GetScreeningUe()),
                                      (std::uint64_t)p.GetScreeningParam(), (std::uint64_t)spec.mesh,
                                      bits(spec.frameShift), bits(spec.stoppingConversion)};
  };
  std::vector<int> firstReader(slots_.size(), -1);  // per block (its first slot): its first reader
  // The plans of the shared meshes for the present anchors and supports.
  const std::vector<DataSet::NucleusMesh> &meshes = Data().Meshes();
  std::vector<azure2::MeshPlan> plans(meshes.size());
  for (std::size_t m = 0; m < meshes.size() && m < meshAnchors_.size(); m++)
    plans[m] = azure2::PlanMesh(meshAnchors_[m], meshSupports_[m], azure2::MeshTails());
  for (int u = 0; u < NumUnits(); u++) {
    const Unit &p = units_[u];
    const int f = p.SubPointFirst(), c = p.NumSubPoints();
    if (c == 0) continue;
    reads += c;
    if (f < 0 || f + c > NumSlots()) {
      outside++;
      report("unit " + std::to_string(u) + ": slots [" + std::to_string(f) + ", " + std::to_string(f + c) +
             ") beyond the " + std::to_string(NumSlots()) + " slots");
      continue;
    }
    const int b = slotBlock_[f];
    bool inBlock = true;
    for (int k = f; k < f + c; k++) inBlock = inBlock && slotBlock_[k] == b;
    if (!inBlock) {
      outside++;
      report("unit " + std::to_string(u) + " (segment key " + std::to_string(p.GetSegmentKey()) +
             "): its slots leave the block at slot " + std::to_string(b));
    }
    for (int k = f; k < f + c; k++) {
      const int w = slotWriter_[k];
      const bool readsIt = w >= 0 && w < NumUnits() && units_[w].SubPointFirst() <= k &&
                           k < units_[w].SubPointFirst() + units_[w].NumSubPoints();
      if (!readsIt || w > u) {
        writers++;
        report("slot " + std::to_string(k) + ": its writer, unit " + std::to_string(w) + ", is no earlier reader of it");
      }
      // the sites of a block: one term and one angle
      const Site &slot = slots_[k], &head = slots_[slotBlock_[k]];
      if (slot.GetSegmentKey() != head.GetSegmentKey() || bits(slot.GetCMAngle()) != bits(head.GetCMAngle()) ||
          bits(slot.GetLabAngle()) != bits(head.GetLabAngle())) {
        otherGroup++;
        report("slot " + std::to_string(k) + " has another term or angle than its block's first slot");
      }
    }
    if (firstReader[b] < 0) {
      firstReader[b] = u;
    } else if (key(units_[firstReader[b]]) != key(p)) {
      otherGroup++;
      report("unit " + std::to_string(u) + " reads the block of unit " + std::to_string(firstReader[b]) +
             ", another group");
    }
    const Unit::SubGridSpec &spec = Data().GridSpecOf(p);
    if (spec.mesh >= 0 && (std::size_t)spec.mesh < plans.size()) {
      const std::vector<double> energies = azure2::MeshSlotEnergies(plans[spec.mesh], spec);
      bool same = (int)energies.size() == c;
      for (int k = 0; same && k < c; k++) same = slots_[f + k].GetOriginalEnergy() == energies[k];
      if (!same) {
        offPlan++;
        report("unit " + std::to_string(u) + " (segment key " + std::to_string(p.GetSegmentKey()) + "): " +
               std::to_string(c) + " slots, not the " + std::to_string(energies.size()) + " nodes of its mesh range");
      }
    }
  }
  long blocks = 0;
  for (int k = 0; k < NumSlots(); k++) blocks += slotBlock_[k] == k;
  std::cerr << "AZURE_NODE_CHECK: slot groups: " << blocks << " blocks, " << NumSlots() << " slots, " << reads
            << " reads; " << outside << " ranges outside a block, " << writers << " slots with a wrong writer, "
            << otherGroup << " reads of another group, " << offPlan << " ranges off their mesh plan" << std::endl;
  for (const std::string &r : reports) std::cerr << "AZURE_NODE_CHECK:   " << r << std::endl;
  return outside + writers + otherGroup + offPlan;
}
