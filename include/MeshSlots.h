#ifndef AZURE2_MESHSLOTS_H
#define AZURE2_MESHSLOTS_H

// The slots of folded points on the shared energy mesh of their compound nucleus (docs/upgrade/energy_grid.md
// sections 2.2, 2.3, 2.8): which points take them, the supports and anchors handed to the planner (MeshPlanner.h), and
// the c.m. energies of a point's slots.  Used by WorkingCopy::SetUpTargetEffects (one plan per nucleus) and by the
// hand-over, WorkingCopy::ReplaceSubPointGrids (a new plan when the quantized anchors change).

#include <functional>
#include <utility>
#include <vector>

#include "DataSet.h"
#include "MeshPlanner.h"
#include "Segment.h"
#include "SlotGroup.h"
#include "Unit.h"

class CompoundSystem;
class LevelInputs;
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

/// Whether a segment's folded points take their slots from the shared mesh: every segment whose pairs give an E_x
/// frame.  A segment with an energy shift is folded there with the kernel moving on the fixed nodes (ShiftedFold.h).
///
/// THM data are folded points of the mesh like any other.  Their support is the kernel's (EffectSupport: the Gaussian
/// resolution or a beam profile, tails reaching kThmConvolutionRange sigma) at the point's c.m. energy E_0 of x + A,
/// clipped at the data kind's floor (Segment::IntegrationFloor): a HOES excitation function reaches below the
/// entrance threshold, so its floor in E_cm is minIntegrationEnergy - (S + E_x of the entrance pair), the compound
/// system's E_x = minIntegrationEnergy, and E_cm of its nodes may be negative.  The caller hands the planner the
/// support in E_x, [max(E_0 - r sigma, floor) + S + E_x, E_0 + r sigma + S + E_x], with lowHard when the floor clips
/// it: lowHard means what it means for a thick target -- a physical floor that is itself a node, below which no node is
/// evaluated for this support.  The planner works in E_x and needs nothing else.  Nodes: a THM site at E_x and a
/// direct site at the same E_x hold one table of channel rows (L_o, P, the phases; WorkingCopy::FillChannelRows), and
/// a THM site its vertex pieces besides (j_l, rho j_l', the Coulomb term, the window's table: Site::GetThmFormFactor).
/// The level matrix at E_x is the same for both; NodeKey's route flag (NodeTable.h) keeps the THM units apart from the
/// direct ones until the energy pass solves a node once for every reader (CalculationPlan.h).
bool OnSharedMesh(const Segment &segment, const PPair *entrance, const PPair *exit);

/// Whether a support's lower edge lies at a physical floor (a thick target's 1 keV, the integration floor or a
/// threshold): a hard edge of the mesh (MeshSupport::lowHard).
bool EdgeAtFloor(const TargetEffect &effect, double low, double floor);

/// Records a folded point's support (E_x) in the mesh of its nucleus and the mesh's place in its grid record; the
/// support of a point with an energy shift (spec.shifted) is the planned one at its start centre (ShiftedFold.h), its
/// window's growth at a shift given by `thicknessChange` (azure2::SupportThicknessChange).
MeshPoint QueueMeshPoint(std::vector<DataSet::NucleusMesh> &meshes, Unit &point, const TermRef &term,
                         Unit::SubGridSpec &spec, const PPair *entrance, const PPair *exit, bool lowHard,
                         const std::function<double(double)> &thicknessChange);

/// The support of a folded point (E_x) as its kernels give it at the nominal point, from its grid record.
MeshSupport NominalMeshSupport(const Unit::SubGridSpec &spec, int dataSet);

/// The quantized anchors of a nucleus's mesh, in E_x: the levels as each entrance pair of its folded points sees them
/// (AdaptiveIntegrationGrid::Anchors in that pair's frame), merged and sorted.  `inputWidthsArePhysical` as the grid
/// configuration's: at the set-up the input widths, in the hand-over the observed parameters of the current values.
std::vector<MeshAnchor> MeshAnchorsOf(const CompoundSystem *nucleus, const LevelInputs &inputs, const DataSet::NucleusMesh &mesh,
                                      bool inputWidthsArePhysical);

/// A point's slots on a mesh: the nodes of its support's range, as c.m. energies of its pair, descending.
std::vector<double> MeshSlotEnergies(const MeshPlan &plan, const Unit::SubGridSpec &spec);

/// The key of a folded point's sharing group on mesh `mesh` (SlotGroup.h).
SlotGroupKey MeshGroupKey(const WorkingCopy &data, const Unit &point, const TermRef &term, const Unit::SubGridSpec &spec);

/// The nodes of a sharing group on a plan: the union of its members' ranges (`supports`: each member's support in
/// the plan), descending as the slots are, and each member's range in it, (first, count).
std::vector<std::size_t> GroupNodes(const MeshPlan &plan, const std::vector<int> &supports,
                                    std::vector<std::pair<int, int>> &ranges);

/// Gives the points queued for mesh `mesh` (`pending`, those of other meshes skipped) their slots on `plan`: one
/// block per sharing group (SlotGroup.h), each slot with the stopping of its effect's target at its energy.
void AddMeshSlotGroups(WorkingCopy &data, const MeshPlan &plan, int mesh, const std::vector<MeshPoint> &pending,
                       const Config &configure);

/// The tail constants of AdaptiveIntegrationGrid, for the planner.
MeshTailConstants MeshTails();

/// Whether two anchor lists are the same, bit for bit (energy and particle width).
bool SameAnchors(const std::vector<MeshAnchor> &a, const std::vector<MeshAnchor> &b);

/// Whether [low, high] meets one of the stretches (MeshChangedRegion's).
bool MeetsStretches(const std::vector<std::pair<double, double>> &stretches, double low, double high);

}  // namespace azure2

#endif  // AZURE2_MESHSLOTS_H
