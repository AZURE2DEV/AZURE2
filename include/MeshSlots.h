#ifndef AZURE2_MESHSLOTS_H
#define AZURE2_MESHSLOTS_H

// The slots of folded points on the shared energy mesh of their compound nucleus (docs/upgrade/energy_grid.md
// sections 2.2, 2.3, 2.8): which points take them, the supports and anchors handed to the planner (MeshPlanner.h), and
// the c.m. energies of a point's slots.  Used by WorkingCopy::SetUpTargetEffects (one plan per nucleus) and by the
// hand-over, WorkingCopy::ReplaceSubPointGrids (a new plan when the quantized anchors change).

#include <utility>
#include <vector>

#include "DataSet.h"
#include "MeshPlanner.h"
#include "Segment.h"
#include "Unit.h"

class CNuc;
class Config;
class PPair;
class TargetEffect;
class WorkingCopy;

namespace azure2 {

/// A folded point waiting for its slots from the shared mesh of its nucleus (WorkingCopy::SetUpTargetEffects).
struct MeshPoint {
  Unit *point;
  TermRef term;
  int mesh;  ///< index in DataSet::meshes_
};

/// The pair whose separation and excitation turn a point's c.m. energy into E_x (AMatrixFunc::InEnergy): the exit
/// pair of a beta-delayed spectrum (entrance of type 20), the entrance pair otherwise.
const PPair *FramePair(const PPair *entrance, const PPair *exit);

/// Whether a segment's folded points take their slots from the shared mesh.  Not a segment with an energy shift (its
/// slots move rigidly with the shift until the shifted fold puts it on fixed nodes, branch 8), not THM (its kinematics
/// on the mesh are branch 10).
bool OnSharedMesh(const Segment &segment, const PPair *entrance, const PPair *exit);

/// Whether a support's lower edge lies at a physical floor (a thick target's 1 keV, the integration floor or a
/// threshold): a hard edge of the mesh (MeshSupport::lowHard).
bool EdgeAtFloor(const TargetEffect &effect, double low, double floor);

/// Records a folded point's support (E_x) in the mesh of its nucleus and the mesh's place in its grid record.
MeshPoint QueueMeshPoint(std::vector<DataSet::NucleusMesh> &meshes, Unit &point, const TermRef &term,
                         Unit::SubGridSpec &spec, const PPair *entrance, const PPair *exit, bool lowHard);

/// The quantized anchors of a nucleus's mesh, in E_x: the levels as each entrance pair of its folded points sees them
/// (AdaptiveIntegrationGrid::Anchors in that pair's frame), merged and sorted.  `inputWidthsArePhysical` as the grid
/// configuration's: at the set-up the input widths, in the hand-over the observed parameters of the current values.
std::vector<MeshAnchor> MeshAnchorsOf(CNuc *nucleus, const DataSet::NucleusMesh &mesh,
                                      bool inputWidthsArePhysical);

/// A point's slots on a mesh: the nodes of its support's range, as c.m. energies of its pair, descending.
std::vector<double> MeshSlotEnergies(const MeshPlan &plan, const Unit::SubGridSpec &spec);

/// Gives a folded point the slots at the c.m. energies of `grid` (descending), each with the stopping of its
/// effect's target at that energy.
void AddGridSlots(WorkingCopy &data, Unit &point, const TermRef &term, const std::vector<double> &grid,
                  const TargetEffect &effect, double conversion, const Config &configure);

/// The tail constants of AdaptiveIntegrationGrid, for the planner.
MeshTailConstants MeshTails();

/// Whether two anchor lists are the same, bit for bit (energy and particle width).
bool SameAnchors(const std::vector<MeshAnchor> &a, const std::vector<MeshAnchor> &b);

/// Whether [low, high] meets one of the stretches (MeshChangedRegion's).
bool MeetsStretches(const std::vector<std::pair<double, double>> &stretches, double low, double high);

}  // namespace azure2

#endif  // AZURE2_MESHSLOTS_H
