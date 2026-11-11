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

// The key of a node: nucleus, route (THM or not: a THM unit's rows carry the vertex pieces), bits of E_x.
using Key = std::tuple<int, char, std::uint64_t>;

Key KeyOf(CNuc *cnuc, Unit &u) { return Key(NucleusOf(cnuc, u), u.IsTHM() ? 1 : 0, Bits(AMatrixFunc::InEnergy(cnuc, &u))); }

}  // namespace

int NodeTable::Find(int n, double e) const {
  const std::uint64_t b = Bits(e);
  // nodes are sorted by (nucleus, thm, bits); a classic node is looked for first
  for (char thmFlag = 0; thmFlag <= 1; thmFlag++) {
    int lo = 0, hi = Size();
    while (lo < hi) {
      const int mid = (lo + hi) / 2;
      const Key k(nucleus[mid], thm[mid], Bits(energy[mid]));
      if (k < Key(n, thmFlag, b))
        lo = mid + 1;
      else
        hi = mid;
    }
    if (lo < Size() && nucleus[lo] == n && thm[lo] == thmFlag && Bits(energy[lo]) == b) return lo;
  }
  return -1;
}

std::size_t NodeTable::Bytes() const {
  std::size_t n = sizeof *this + nucleus.capacity() * sizeof(int) + energy.capacity() * sizeof(double) +
                  rows.capacity() * sizeof(rows[0]) + thm.capacity() + (slotBegin.capacity() + slotNode.capacity()) * sizeof(int);
  for (const auto &r : rows)
    if (r) n += TableBytes(*r);
  return n;
}

std::shared_ptr<const NodeTable> BuildNodeTable(WorkingCopy &data, CNuc *cnuc, const Config &) {
  const CalculationPlan &plan = data.Plan();
  // Every unit of the pass with its key, in plan order; the first unit at a key gives the node its rows.
  std::vector<std::pair<Key, UnitRef>> all;
  std::vector<UnitRef> units;
  for (int t = 0; t < plan.NumTasks(); t++) {
    units.clear();
    plan.TaskUnits(data, t, units);
    for (UnitRef r : units) all.push_back({KeyOf(cnuc, data.UnitAt(r)), r});
  }
  std::stable_sort(all.begin(), all.end(),
                   [](const std::pair<Key, UnitRef> &a, const std::pair<Key, UnitRef> &b) { return a.first < b.first; });
  auto table = std::make_shared<NodeTable>();
  for (std::size_t k = 0; k < all.size(); k++) {
    if (k > 0 && all[k].first == all[k - 1].first) continue;
    const Unit &u = data.UnitAt(all[k].second);
    table->nucleus.push_back(std::get<0>(all[k].first));
    table->thm.push_back(std::get<1>(all[k].first));
    table->energy.push_back(AMatrixFunc::InEnergy(cnuc, &const_cast<Unit &>(u)));
    table->rows.push_back(u.ChannelRowsShared());
  }
  // The slots of the folded points.
  table->slotBegin.assign(1, 0);
  for (int t = 0; t < plan.NumTasks(); t++) {
    Unit &p = *plan.Point(data, t);
    if (HasSubPointFold(data, p))
      for (int i = 1; i <= p.NumSubPoints(); i++) {
        Unit &s = *p.GetSubPoint(data, i);
        table->slotNode.push_back(table->Find(NucleusOf(cnuc, s), AMatrixFunc::InEnergy(cnuc, &s)));
        if (table->slotNode.back() >= 0 && table->thm[table->slotNode.back()] != (s.IsTHM() ? 1 : 0))
          table->slotNode.back() = -1;
      }
    table->slotBegin.push_back((int)table->slotNode.size());
  }
  return table;
}

bool NodeCheckRequested() {
  static const bool on = [] {
    const char *v = std::getenv("AZURE_NODE_CHECK");
    return v && *v && std::strcmp(v, "0") != 0;
  }();
  return on;
}

}  // namespace azure2

long WorkingCopy::CheckNodeTable(CNuc *cnuc, const std::vector<int> &work) {
  const azure2::NodeTable *table = nodeTable_.get();
  if (!table) throw std::logic_error("AZURE_NODE_CHECK: the working copy has no node table");
  long units = 0, onTable = 0, offTable = 0, rowsDiffer = 0, slots = 0, slotsOnUnit = 0, slotsMoved = 0, slotsDiffer = 0;
  std::vector<std::string> reports;
  std::vector<azure2::UnitRef> refs;
  for (int t : work) {
    refs.clear();
    plan_.TaskUnits(*this, t, refs);
    for (azure2::UnitRef r : refs) {
      Unit &u = units_[r];
      units++;
      const azure2::Key key = azure2::KeyOf(cnuc, u);
      const int n = table->Find(std::get<0>(key), AMatrixFunc::InEnergy(cnuc, &u));
      if (n < 0 || table->thm[n] != std::get<1>(key)) {
        offTable++;  // re-placed or shifted inside the stage: part 2's overflow table
        continue;
      }
      onTable++;
      const EDependentTables *mine = u.ChannelRows();
      if (!mine || !table->rows[n] || !azure2::SameTables(*mine, *table->rows[n])) {
        rowsDiffer++;
        if (reports.size() < 10)
          reports.push_back("channel rows of unit " + std::to_string(r) + " (segment key " +
                            std::to_string(u.GetSegmentKey()) + ") differ from node " + std::to_string(n));
      }
    }
    // The slots of a folded point: as many as its sub-points; a slot on a node holds the node its sub-point is at,
    // unless the sub-point has moved off the stage's node (re-placed or shifted: counted, not compared).
    Unit &p = *plan_.Point(*this, t);
    if (!azure2::HasSubPointFold(*this, p)) continue;
    const int b = table->slotBegin[t], e = table->slotBegin[t + 1];
    if (e - b != p.NumSubPoints()) {
      slotsMoved += p.NumSubPoints();  // a re-placed grid with another node count
      continue;
    }
    for (int i = 0; i < e - b; i++) {
      slots++;
      const int n = table->slotNode[b + i];
      Unit &s = *p.GetSubPoint(*this, i + 1);
      const double ex = AMatrixFunc::InEnergy(cnuc, &s);
      if (n >= 0 && std::memcmp(&table->energy[n], &ex, sizeof ex) == 0) {
        slotsOnUnit++;
        continue;
      }
      const int now = table->Find(azure2::NucleusOf(cnuc, s), ex);
      if (n < 0 && now < 0) continue;  // off the table from the start, still off
      if (now != n) {
        slotsMoved++;  // the sub-point moved since the set-up
        continue;
      }
      slotsDiffer++;
      if (reports.size() < 10) reports.push_back("slot " + std::to_string(i) + " of plan task " + std::to_string(t) +
                                                 " names node " + std::to_string(n) + ", its sub-point is at node " +
                                                 std::to_string(now));
    }
  }
  std::cerr << "AZURE_NODE_CHECK: " << table->Size() << " nodes (" << table->Bytes() << " bytes); " << units
            << " units of the pass: " << onTable << " on a node of the table, " << offTable << " off it; channel rows "
            << rowsDiffer << " differ; " << slots << " slots: " << slotsOnUnit << " at their sub-point's node, "
            << slotsMoved << " moved, " << slotsDiffer << " differ" << std::endl;
  for (const std::string &r : reports) std::cerr << "AZURE_NODE_CHECK:   " << r << std::endl;
  return rowsDiffer + slotsDiffer;
}
