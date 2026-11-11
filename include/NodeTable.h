#ifndef AZURE2_NODETABLE_H
#define AZURE2_NODETABLE_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "SharedTable.h"
#include "UnitTables.h"

class WorkingCopy;
class CNuc;
class Config;

namespace azure2 {

/*!
 * The node table of a fit stage (docs/upgrade/energy_grid.md section 2.6; branch 7, part 1): every compound energy at
 * which the point pass solves the level matrix -- a point without a sub-point fold, or a sub-point of one -- once per
 * (compound nucleus, bits of E_x), exact duplicates merged as the energy tasks group them (CalculationPlan::EnergyTasks),
 * with the channel rows there (L_o, sqrt(P), Coulomb and hard-sphere phases of every channel of every J group: they
 * depend on E_x alone).  Built once at the end of the set-up (WorkingCopy::Initialize) from the units the set-up
 * placed, immutable, and held by every copy of the working copy through a shared_ptr<const>.
 *
 * Beside it, per plan task of a point with a sub-point fold, the point's slots: for each of its sub-points (in order) the
 * node it is at, as an index of the table, or -1 for a node the stage table does not hold (a sub-point re-placed or
 * moved by an energy shift inside the stage: what part 2 puts in the working copy's own overflow table).
 *
 * In this part the table and the slots stand beside the sub-point units, which still carry the calculation; the check
 * mode AZURE_NODE_CHECK=1 shows, at every point pass, that they say what the units say, bit for bit
 * (WorkingCopy::CheckNodeTable).  E_x in MeV (AMatrixFunc::InEnergy).
 */
struct NodeTable {
  std::vector<int> nucleus;                         ///< compound nucleus of each node (1-based, PPair::GetNucleus)
  std::vector<double> energy;                       ///< E_x of each node, MeV; sorted by (nucleus, bits)
  std::vector<SharedTable<EDependentTables>> rows;  ///< the channel rows at the node (shared with the unit they came from)
  std::vector<char> thm;                            ///< 1: a THM node (its rows carry the THM vertex pieces as well)
  /// The slots of the point of each plan task: [slotBegin[t], slotBegin[t + 1]) of slotNode (empty without a fold).
  std::vector<int> slotBegin;
  std::vector<int> slotNode;  ///< node index of each slot, -1 when no node of the table
  int Size() const { return (int)energy.size(); }
  /// The node at (nucleus, E_x) bit for bit, or -1.
  int Find(int nucleus, double energy) const;
  /// Bytes held (the rows counted once per node, shared or not).
  std::size_t Bytes() const;
};

/// The node table of `data` as its units are now (one pass over the plan tasks' units; see NodeTable).
std::shared_ptr<const NodeTable> BuildNodeTable(WorkingCopy &data, CNuc *cnuc, const Config &configure);

/// AZURE_NODE_CHECK=1: the node table of `data` against its units, bit for bit (see WorkingCopy::CheckNodeTable).
bool NodeCheckRequested();

}  // namespace azure2

#endif
