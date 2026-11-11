// The node table of a fit stage and its check (include/NodeTable.h; docs/upgrade/energy_grid.md section 2.6).
#include "NodeTable.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <tuple>

#include "AMatrixFunc.h"
#include "CNuc.h"
#include "CalculationPlan.h"
#include "Config.h"
#include "PPair.h"
#include "Unit.h"
#include "WorkingCopy.h"

namespace azure2 {

namespace {

std::uint64_t Bits(double x) {
  std::uint64_t b;
  std::memcpy(&b, &x, sizeof b);
  return b;
}

int NucleusOf(CNuc *cnuc, const Unit &u) { return cnuc->GetPair(cnuc->GetPairNumFromKey(u.GetEntranceKey()))->GetNucleus(); }

template <class V>
bool SameBits(const V &a, const V &b) {
  return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(a[0])) == 0);
}

template <class M>
bool SameRows(const M &a, const M &b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); i++)
    if (!SameBits(a[i], b[i])) return false;
  return true;
}

// The channel rows of two units at one node, bit for bit.
bool SameTables(const EDependentTables &a, const EDependentTables &b) {
  return SameBits(a.rowStart, b.rowStart) && SameBits(a.lo_elements, b.lo_elements) &&
         SameBits(a.penetrabilities, b.penetrabilities) && SameBits(a.coulombphase, b.coulombphase) &&
         SameBits(a.hardspherephase, b.hardspherephase) && SameRows(a.thm_jl, b.thm_jl) &&
         SameRows(a.thm_rhodjl, b.thm_rhodjl) && SameRows(a.thm_coul, b.thm_coul);
}

std::size_t TableBytes(const EDependentTables &t) {
  std::size_t n = sizeof t + t.rowStart.capacity() * sizeof(int) +
                  (t.lo_elements.capacity() + t.coulombphase.capacity() + t.hardspherephase.capacity()) * sizeof(complex) +
                  t.penetrabilities.capacity() * sizeof(double);
  for (const auto &r : t.thm_jl) n += r.capacity() * sizeof(double);
  for (const auto &r : t.thm_rhodjl) n += r.capacity() * sizeof(double);
  for (const auto &r : t.thm_coul) n += r.capacity() * sizeof(complex);
  return n;
}

NodeKey KeyOf(CNuc *cnuc, Unit &u) {
  return NodeKey(NucleusOf(cnuc, u), u.IsTHM() ? 1 : 0, Bits(AMatrixFunc::InEnergy(cnuc, &u)));
}

}  // namespace

NodeKey NodeList::Key(int n) const { return NodeKey(nucleus[n], thm[n], Bits(energy[n])); }

void NodeList::Append(const NodeKey &key, double ex, SharedTable<EDependentTables> nodeRows) {
  nucleus.push_back(std::get<0>(key));
  thm.push_back(std::get<1>(key));
  energy.push_back(ex);
  rows.push_back(std::move(nodeRows));
}

void NodeList::Clear() {
  nucleus.clear();
  energy.clear();
  rows.clear();
  thm.clear();
}

std::size_t NodeList::Bytes() const {
  std::size_t n = nucleus.capacity() * sizeof(int) + energy.capacity() * sizeof(double) +
                  rows.capacity() * sizeof(rows[0]) + thm.capacity();
  for (const auto &r : rows)
    if (r) n += TableBytes(*r);
  return n;
}

int NodeTable::Find(const NodeKey &key) const {
  // nodes are sorted by their key
  int lo = 0, hi = Size();
  while (lo < hi) {
    const int mid = (lo + hi) / 2;
    if (Key(mid) < key)
      lo = mid + 1;
    else
      hi = mid;
  }
  return (lo < Size() && Key(lo) == key) ? lo : -1;
}

std::size_t NodeTable::Bytes() const {
  return sizeof *this + NodeList::Bytes() +
         (slots.begin.capacity() + slots.node.capacity()) * sizeof(int) + slots.overflow.capacity();
}

std::size_t NodeOverflow::KeyHash::operator()(const NodeKey &k) const {
  const std::uint64_t h = std::get<2>(k) * 0x9E3779B97F4A7C15ULL ^
                          ((std::uint64_t)(unsigned)std::get<0>(k) << 1 | (std::uint64_t)std::get<1>(k));
  return (std::size_t)(h ^ (h >> 29));
}

int NodeOverflow::Place(const NodeKey &key, double ex, const SharedTable<EDependentTables> &nodeRows) {
  const auto it = index.emplace(key, Size());
  if (it.second) Append(key, ex, nodeRows);
  return it.first->second;
}

NodeState BuildNodeTable(WorkingCopy &data, CNuc *cnuc, const Config &) {
  const CalculationPlan &plan = data.Plan();
  // Every unit of the pass with its key, in plan order; the first unit at a key gives the node its rows.
  std::vector<std::pair<NodeKey, UnitRef>> all;
  std::vector<UnitRef> units;
  auto table = std::make_shared<NodeTable>();
  for (int t = 0; t < plan.NumTasks(); t++) {
    units.clear();
    plan.TaskUnits(data, t, units);
    for (UnitRef r : units) all.push_back({KeyOf(cnuc, data.UnitAt(r)), r});
    table->slots.begin.push_back((int)all.size());
  }
  std::vector<int> order(all.size());
  for (std::size_t k = 0; k < all.size(); k++) order[k] = (int)k;
  std::stable_sort(order.begin(), order.end(), [&all](int a, int b) { return all[a].first < all[b].first; });
  // The slots, in plan order: the node of each unit.
  table->slots.node.assign(all.size(), -1);
  table->slots.overflow.assign(plan.NumTasks(), 0);
  for (std::size_t k = 0; k < order.size(); k++) {
    const auto &u = all[order[k]];
    if (k == 0 || u.first != all[order[k - 1]].first) {
      Unit &unit = data.UnitAt(u.second);
      table->Append(u.first, AMatrixFunc::InEnergy(cnuc, &unit), unit.ChannelRowsShared());
    }
    table->slots.node[order[k]] = table->Size() - 1;
  }
  NodeState state;
  state.slots = std::shared_ptr<const NodeSlots>(table, &table->slots);
  state.table = std::move(table);
  return state;
}

bool NodeCheckRequested() {
  static const bool on = [] {
    const char *v = std::getenv("AZURE_NODE_CHECK");
    return v && *v && std::strcmp(v, "0") != 0;
  }();
  return on;
}

}  // namespace azure2

void WorkingCopy::PlaceOverflowNodes(CNuc *cnuc, const std::vector<int> &changedTasks, bool rowsReleased) {
  const azure2::NodeTable *stage = nodes_.table.get();
  if (!stage || !cnuc || (changedTasks.empty() && !rowsReleased)) return;
  const azure2::NodeSlots &old = *nodes_.slots;
  std::vector<char> changed(plan_.NumTasks(), 0);
  for (int t : changedTasks) changed[t] = 1;
  // Every slot of a task that moved, or that was on the overflow table (the table is placed again from nothing, so it
  // holds only what the copy's units need now, in plan order), is looked up; every other slot is the old one.
  auto slots = std::make_shared<azure2::NodeSlots>();
  slots->node.reserve(old.node.size());
  slots->begin.reserve(old.begin.size());
  slots->overflow.assign(plan_.NumTasks(), 0);
  nodes_.overflow.Clear();
  std::vector<azure2::UnitRef> refs;
  for (int t = 0; t < plan_.NumTasks(); t++) {
    if (!changed[t] && !old.overflow[t]) {
      slots->node.insert(slots->node.end(), old.node.begin() + old.begin[t], old.node.begin() + old.begin[t + 1]);
    } else {
      refs.clear();
      plan_.TaskUnits(*this, t, refs);
      // the units of a task share their point's entrance pair, and so its compound nucleus
      const int nucleus = refs.empty() ? 0 : azure2::NucleusOf(cnuc, units_[refs[0]]);
      for (azure2::UnitRef r : refs) {
        Unit &u = units_[r];
        const double ex = AMatrixFunc::InEnergy(cnuc, &u);
        const azure2::NodeKey key(nucleus, u.IsTHM() ? 1 : 0, azure2::Bits(ex));
        int n = stage->Find(key);
        if (n < 0) {
          n = stage->Size() + nodes_.overflow.Place(key, ex, u.ChannelRowsShared());
          slots->overflow[t] = 1;
        }
        slots->node.push_back(n);
      }
    }
    slots->begin.push_back((int)slots->node.size());
  }
  nodes_.slots = std::move(slots);
}

void WorkingCopy::ReleaseOverflowRows() {
  // The rows of the overflow table are holders of the units' own; a unit about to be re-energised writes them in
  // place only when it is their sole holder (SharedTable copy-on-write), so they are let go first.
  for (auto &r : nodes_.overflow.rows) r = SharedTable<EDependentTables>();
}

long WorkingCopy::CheckNodeTable(CNuc *cnuc, const std::vector<int> &work) {
  const azure2::NodeTable *stage = nodes_.table.get();
  if (!stage || !nodes_.slots) throw std::logic_error("AZURE_NODE_CHECK: the working copy has no node table");
  const azure2::NodeSlots &slots = *nodes_.slots;
  long units = 0, onStage = 0, onOverflow = 0, offTables = 0, wrongNode = 0, rowsDiffer = 0;
  std::vector<std::string> reports;
  auto report = [&reports](const std::string &r) {
    if (reports.size() < 10) reports.push_back(r);
  };
  std::vector<azure2::UnitRef> refs;
  for (int t : work) {
    refs.clear();
    plan_.TaskUnits(*this, t, refs);
    const int b = slots.begin[t], e = slots.begin[t + 1];
    if (e - b != (int)refs.size()) {
      offTables += (long)refs.size();
      report("plan task " + std::to_string(t) + " has " + std::to_string(e - b) + " slots for " +
             std::to_string(refs.size()) + " units");
      continue;
    }
    for (int k = 0; k < e - b; k++) {
      Unit &u = units_[refs[k]];
      units++;
      const int n = slots.node[b + k];
      const bool inStage = n >= 0 && n < stage->Size();
      const bool inOverflow = n >= stage->Size() && n - stage->Size() < nodes_.overflow.Size();
      if (!inStage && !inOverflow) {
        offTables++;
        report("slot " + std::to_string(k) + " of plan task " + std::to_string(t) + " names no node (" +
               std::to_string(n) + ")");
        continue;
      }
      const azure2::NodeList &list = inStage ? static_cast<const azure2::NodeList &>(*stage) : nodes_.overflow;
      const int at = inStage ? n : n - stage->Size();
      if (list.Key(at) != azure2::KeyOf(cnuc, u)) {
        wrongNode++;
        report("slot " + std::to_string(k) + " of plan task " + std::to_string(t) + " names node " + std::to_string(n) +
               ", its unit " + std::to_string(refs[k]) + " is elsewhere");
        continue;
      }
      (inStage ? onStage : onOverflow)++;
      const EDependentTables *mine = u.ChannelRows();
      if (!mine || !list.rows[at] || !azure2::SameTables(*mine, *list.rows[at])) {
        rowsDiffer++;
        report("channel rows of unit " + std::to_string(refs[k]) + " (segment key " + std::to_string(u.GetSegmentKey()) +
               ") differ from node " + std::to_string(n));
      }
    }
  }
  std::cerr << "AZURE_NODE_CHECK: " << stage->Size() << " stage nodes (" << stage->Bytes() << " bytes), "
            << nodes_.overflow.Size() << " overflow nodes (" << nodes_.overflow.Bytes() << " bytes); " << units
            << " units of the pass: " << onStage << " on the stage table, " << onOverflow << " on the overflow table, "
            << offTables << " off the tables, " << wrongNode << " on a wrong node; channel rows " << rowsDiffer
            << " differ" << std::endl;
  for (const std::string &r : reports) std::cerr << "AZURE_NODE_CHECK:   " << r << std::endl;
  return offTables + wrongNode + rowsDiffer;
}
