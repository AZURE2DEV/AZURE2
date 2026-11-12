#ifndef AZURE2_NODETABLE_H
#define AZURE2_NODETABLE_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "SharedTable.h"
#include "UnitTables.h"

class WorkingCopy;
class CNuc;
class Config;

namespace azure2 {

/// A node by (compound nucleus, route, bits of E_x): route 1 for a THM unit (its rows carry the vertex pieces), else 0.
using NodeKey = std::tuple<int, char, std::uint64_t>;

/*!
 * A list of nodes: compound energies at which the point pass solves the level matrix, each with the channel rows there
 * (L_o, sqrt(P), Coulomb and hard-sphere phases of every channel of every J group: they depend on E_x alone).  The
 * rows are a holder of the rows the first unit at the node computed, not a second copy.  E_x in MeV
 * (AMatrixFunc::InEnergy).
 */
struct NodeList {
  std::vector<int> nucleus;                         ///< compound nucleus of each node (1-based, PPair::GetNucleus)
  std::vector<double> energy;                       ///< E_x of each node, MeV
  std::vector<SharedTable<EDependentTables>> rows;  ///< the channel rows at the node (shared with a unit there)
  std::vector<char> thm;                            ///< 1: a THM node (its rows carry the THM vertex pieces as well)
  int Size() const { return (int)energy.size(); }
  NodeKey Key(int n) const;
  void Append(const NodeKey &key, double ex, SharedTable<EDependentTables> nodeRows);
  void Clear();
  /// Bytes held (the rows counted once per node, shared or not).
  std::size_t Bytes() const;
};

/*!
 * The slots of the plan tasks of a working copy: for task t, [begin[t], begin[t + 1]) of `node`, one slot per unit the
 * task calculates (CalculationPlan::TaskUnits: its point, or its sub-points in order when it carries a sub-point fold),
 * holding the node that unit is at.  A node index below the stage table's size is a node of the stage table, any
 * other one of the working copy's overflow table (NodeOverflow), at index - stage size.
 */
struct NodeSlots {
  std::vector<int> begin{0};
  std::vector<int> node;
  std::vector<char> overflow;  ///< per task: 1 when a slot of it is on the overflow table
};

/*!
 * The node table of a fit stage (docs/upgrade/energy_grid.md section 2.6): every compound energy at which the point
 * pass solves the level matrix -- a point without a sub-point fold, or a sub-point of one -- once per (compound nucleus,
 * route, bits of E_x), exact duplicates merged as the energy tasks group them (CalculationPlan::EnergyTasks), sorted by
 * that key.  Built once at the end of the set-up (WorkingCopy::Initialize) from the units the set-up placed, immutable,
 * and held by every copy of the working copy through a shared_ptr<const>; with it the slots of every plan task as the
 * set-up placed them (each on a node of this table).
 */
struct NodeTable : NodeList {
  NodeSlots slots;
  /// The node at `key` bit for bit, or -1.
  int Find(const NodeKey &key) const;
  /// Bytes held by the nodes and the slots.
  std::size_t Bytes() const;
};

/*!
 * The overflow table of a working copy (docs/upgrade/energy_grid.md section 2.6): the nodes of its units that the stage
 * table does not hold, placed in the copy's hand-over (WorkingCopy::PlaceOverflowNodes) and a function of the
 * parameter vector the copy was handed alone.  Two kinds of node live here:
 *  - nodes re-placed inside a stage: a level's quantized anchor moved and WorkingCopy::ReplaceSubPointGrids laid the
 *    point's sub-point grid out again (the rule kept bit for bit);
 *  - the node of a point without a fold in a segment with an energy shift: it moves with the shift (energy_grid.md
 *    section 2.5, Segment::UpdatePointEnergyWithShift).  A folded point with a shift keeps the nodes of the mesh and
 *    its kernel moves (ShiftedFold.h): its slots never land here for the shift.
 */
struct NodeOverflow : NodeList {
  /// A node at `key`: the last one placed when it is at `key` (the slots of a grid come in order), else a new one with
  /// `nodeRows`.  An exact duplicate further back is a node of its own; the energy tasks group sites by key, not by
  /// node, so it costs a row holder and no solve.  No look-up structure: the table is placed again at every hand-over
  /// that moves a site, and its vectors keep their capacity.
  int Place(const NodeKey &key, double ex, const SharedTable<EDependentTables> &nodeRows);
  void Clear() { NodeList::Clear(); }
};

/// What a working copy holds of the nodes: the stage table (shared by every copy), its slots (the stage table's, or the
/// copy's own once a unit moved off it), and its overflow table.
struct NodeState {
  std::shared_ptr<const NodeTable> table;
  std::shared_ptr<const NodeSlots> slots;
  NodeOverflow overflow;
};

/// The node state of `data` as its units are now (one pass over the plan tasks' units; see NodeTable): the stage table
/// with its slots, an empty overflow table.
NodeState BuildNodeTable(WorkingCopy &data, CNuc *cnuc, const Config &configure);

/// AZURE_NODE_CHECK=1 or AZURE_PLAN_CHECK=1: the node tables and slots of `data` against its units, bit for bit
/// (WorkingCopy::CheckNodeTable: every slot names a node of a table, the stage nodes are in key order, a slot's node
/// has the key -- nucleus, route, bits of E_x computed from the site's c.m. energy -- of its site, and its channel
/// rows are the site's).  At every hand-over and every point pass.
bool NodeCheckRequested();

}  // namespace azure2

#endif
