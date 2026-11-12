// The slots of folded points on the shared mesh of their nucleus (include/MeshSlots.h).
#include "MeshSlots.h"

#include <algorithm>
#include <tuple>

#include "AdaptiveIntegrationGrid.h"
#include "CNuc.h"
#include "Config.h"
#include "PPair.h"
#include "TargetEffect.h"
#include "WorkingCopy.h"

namespace azure2 {

azure2::MeshTailConstants MeshTails() {
  azure2::MeshTailConstants tails;
  tails.tailRatio = AdaptiveIntegrationGrid::tailRatio;
  tails.convolutionCoreWidths = AdaptiveIntegrationGrid::convolutionCoreWidths;
  return tails;
}

/// The pair whose separation and excitation turn a point's c.m. energy into E_x (AMatrixFunc::InEnergy): the exit
/// pair of a beta-delayed spectrum (entrance of type 20), the entrance pair otherwise.
const PPair *FramePair(const PPair *entrance, const PPair *exit) {
  return entrance->GetPType() == 20 ? exit : entrance;
}

/// Whether a segment's folded points take their slots from the shared mesh.  Not a segment with an energy shift (its
/// slots move rigidly with the shift until the shifted fold puts it on fixed nodes, branch 8), not THM (its kinematics
/// on the mesh are branch 10).
bool OnSharedMesh(const Segment &segment, const PPair *entrance, const PPair *exit) {
  if (segment.IsTHM()) return false;
  if (segment.IsVaryEnergyShift() || segment.StartEnergyShift() != 0.0 || segment.GetShiftLeaderPos() >= 0 ||
      segment.HasShiftFn() || segment.IsVaryEnergyShiftSqrt() || segment.StartEnergyShiftSqrt() != 0.0)
    return false;
  return FramePair(entrance, exit) != nullptr;
}

/// Whether a support's lower edge lies at a physical floor (a thick target's 1 keV, the integration floor or a
/// threshold): a hard edge of the mesh (MeshSupport::lowHard).
bool EdgeAtFloor(const TargetEffect &effect, double low, double floor) {
  return (effect.IsTargetIntegration() && effect.IsThickTarget()) || low <= floor;
}

/// Records a folded point's support (E_x) in the mesh of its nucleus and the mesh's place in its grid record.
MeshPoint QueueMeshPoint(std::vector<DataSet::NucleusMesh> &meshes, Unit &point, const TermRef &term,
                         Unit::SubGridSpec &spec, const PPair *entrance, const PPair *exit, bool lowHard) {
  const PPair *frame = FramePair(entrance, exit);
  const int nucleus = entrance->GetNucleus();
  spec.frameShift = frame->GetSepE() + frame->GetExE();
  size_t m = 0;
  while (m < meshes.size() && meshes[m].nucleus != nucleus) m++;
  if (m == meshes.size()) {
    meshes.emplace_back();
    meshes.back().nucleus = nucleus;
  }
  DataSet::NucleusMesh &mesh = meshes[m];
  const std::pair<int, double> frameKey(term.GetEntranceKey(), spec.frameShift);
  if (std::find(mesh.anchorFrames.begin(), mesh.anchorFrames.end(), frameKey) == mesh.anchorFrames.end())
    mesh.anchorFrames.push_back(frameKey);
  azure2::MeshSupport support;
  support.low = spec.endEnergy + spec.frameShift;
  support.high = spec.startEnergy + spec.frameShift;
  support.smoothStep = spec.baseEnergyStep;
  support.settings.pointsPerWidth = spec.pointsPerWidth;
  support.settings.coreWidths = spec.resonanceWidthMultiplier;
  support.lowHard = lowHard;
  support.dataSet = term.GetSegmentKey();
  spec.mesh = (int)m;
  spec.meshSupport = (int)mesh.supports.size();
  mesh.supports.push_back(support);
  return {&point, term, (int)m};
}

/// The quantized anchors of a nucleus's mesh, in E_x: the levels as each entrance pair of its folded points sees them
/// (AdaptiveIntegrationGrid::Anchors in that pair's frame), merged and sorted.  `inputWidthsArePhysical` as the grid
/// configuration's: at the set-up the input widths, in the hand-over the observed parameters of the current values.
std::vector<azure2::MeshAnchor> MeshAnchorsOf(CNuc *nucleus, const DataSet::NucleusMesh &mesh,
                                              bool inputWidthsArePhysical) {
  std::vector<azure2::MeshAnchor> anchors;
  for (const std::pair<int, double> &frame : mesh.anchorFrames) {
    AdaptiveIntegrationGrid::GridConfig gc;
    gc.entranceKey = frame.first;
    gc.inputWidthsArePhysical = inputWidthsArePhysical;
    for (const auto &r : AdaptiveIntegrationGrid(gc).Anchors(nucleus, -frame.second))
      anchors.push_back({r.energy, r.totalWidth, r.particleWidth, r.id});
  }
  auto key = [](const azure2::MeshAnchor &a) { return std::make_tuple(a.energy, a.particleWidth, a.totalWidth, a.id); };
  std::sort(anchors.begin(), anchors.end(),
            [&](const azure2::MeshAnchor &a, const azure2::MeshAnchor &b) { return key(a) < key(b); });
  anchors.erase(std::unique(anchors.begin(), anchors.end(),
                            [](const azure2::MeshAnchor &a, const azure2::MeshAnchor &b) {
                              return a.energy == b.energy && a.particleWidth == b.particleWidth;
                            }),
                anchors.end());
  return anchors;
}

bool SameAnchors(const std::vector<azure2::MeshAnchor> &a, const std::vector<azure2::MeshAnchor> &b) {
  if (a.size() != b.size()) return false;
  for (size_t k = 0; k < a.size(); k++)
    if (a[k].energy != b[k].energy || a[k].particleWidth != b[k].particleWidth) return false;
  return true;
}

bool MeetsStretches(const std::vector<std::pair<double, double>> &stretches, double low, double high) {
  for (const auto &x : stretches)
    if (x.first <= high && x.second >= low) return true;
  return false;
}

/// A point's slots on a mesh: the nodes of its support's range, as c.m. energies of its pair, descending.
std::vector<double> MeshSlotEnergies(const azure2::MeshPlan &plan, const Unit::SubGridSpec &spec) {
  const azure2::MeshNodeRange r = plan.ranges[(size_t)spec.meshSupport];
  std::vector<double> energies;
  for (size_t k = r.last + 1; k-- > r.first;) energies.push_back(plan.nodes[k] - spec.frameShift);
  return energies;
}

/// Gives a folded point the slots at the c.m. energies of `grid` (descending), each with the stopping of its
/// effect's target at that energy.
void AddGridSlots(WorkingCopy &data, Unit &point, const TermRef &term, const std::vector<double> &grid,
                  const TargetEffect &effect, double conversion, const Config &configure) {
  for (double e : grid) {
    Site slot(point.GetCMAngle(), e, term.segment, term.term);
    slot.BelongToRowOf(point);
    if (effect.IsTargetIntegration())
      slot.SetStoppingPower(conversion * effect.StoppingAt(e / conversion, configure,
                                                           data.EffectValues(term.GetTargetEffectNum()).stoppingScale));
    data.AddSubPoint(point, slot);
  }
}

}  // namespace azure2
