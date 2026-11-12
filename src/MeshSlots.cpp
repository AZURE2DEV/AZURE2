// The slots of folded points on the shared mesh of their nucleus (include/MeshSlots.h).
#include "MeshSlots.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <functional>
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

SlotGroupKey MeshGroupKey(const WorkingCopy &data, const Unit &point, const TermRef &term, const Unit::SubGridSpec &spec) {
  auto bits = [](double x) {
    std::uint64_t b;
    std::memcpy(&b, &x, sizeof b);
    return b;
  };
  SlotGroupKey key;
  key.segment = term.segment;
  key.term = term.term;
  key.mesh = spec.mesh;
  key.targetEffect = point.GetTargetEffectNum();
  key.legendreOrder = point.IsCoefficient() ? point.GetCoefficientOrder() : -1;
  key.angle = bits(point.GetCMAngle());
  key.conversion = bits(spec.stoppingConversion);
  key.frameShift = bits(spec.frameShift);
  key.screeningUe = bits(point.GetScreeningUe());
  key.screeningParam = point.GetScreeningParam();
  return key;
}

std::vector<std::size_t> GroupNodes(const MeshPlan &plan, const std::vector<int> &supports,
                                    std::vector<std::pair<int, int>> &ranges) {
  std::vector<std::size_t> nodes;
  for (int s : supports) {
    const MeshNodeRange r = plan.ranges[(std::size_t)s];
    for (std::size_t k = r.first; k <= r.last && k < plan.nodes.size(); k++) nodes.push_back(k);
  }
  std::sort(nodes.begin(), nodes.end(), std::greater<std::size_t>());
  nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
  ranges.clear();
  for (int s : supports) {
    const MeshNodeRange r = plan.ranges[(std::size_t)s];
    const int count = r.last >= r.first ? (int)(r.last - r.first + 1) : 0;
    const auto at = std::lower_bound(nodes.begin(), nodes.end(), r.last, std::greater<std::size_t>());
    ranges.push_back({count ? (int)(at - nodes.begin()) : 0, count});
  }
  return nodes;
}

void AddMeshSlotGroups(WorkingCopy &data, const MeshPlan &plan, int mesh, const std::vector<MeshPoint> &pending,
                       const Config &configure) {
  // The groups in the order of their first points, each with its points in the order queued.
  std::vector<SlotGroupKey> keys;
  std::vector<std::vector<const MeshPoint *>> groups;
  for (const MeshPoint &q : pending) {
    if (q.mesh != mesh) continue;
    const SlotGroupKey key = MeshGroupKey(data, *q.point, q.term, data.Data().GridSpecOf(*q.point));
    std::size_t g = 0;
    while (g < keys.size() && !(keys[g] == key)) g++;
    if (g == keys.size()) keys.push_back(key), groups.emplace_back();
    groups[g].push_back(&q);
  }
  for (const auto &group : groups) {
    const MeshPoint &lead = *group.front();
    const Unit::SubGridSpec &spec = data.Data().GridSpecOf(*lead.point);
    const TargetEffect &effect = *data.Data().GetTargetEffect(lead.term.GetTargetEffectNum());
    std::vector<int> supports;
    for (const MeshPoint *q : group) supports.push_back(data.Data().GridSpecOf(*q->point).meshSupport);
    std::vector<std::pair<int, int>> ranges;
    const std::vector<std::size_t> nodes = GroupNodes(plan, supports, ranges);
    SlotBlock block;
    block.sites.reserve(nodes.size());
    for (std::size_t k : nodes) {
      const double e = plan.nodes[k] - spec.frameShift;
      Site slot(lead.point->GetCMAngle(), e, lead.term.segment, lead.term.term);
      slot.BelongToRowOf(*lead.point);
      if (effect.IsTargetIntegration())
        slot.SetStoppingPower(spec.stoppingConversion *
                              effect.StoppingAt(e / spec.stoppingConversion, configure,
                                                data.EffectValues(lead.term.GetTargetEffectNum()).stoppingScale));
      lead.point->MakeSubPoint(slot);
      block.sites.push_back(std::move(slot));
    }
    for (std::size_t i = 0; i < group.size(); i++)
      if (ranges[i].second > 0)
        block.members.push_back({static_cast<int>(group[i]->point - &data.UnitAt(0)), ranges[i].first, ranges[i].second});
    if (!block.members.empty()) data.AddSlotGroup(std::move(block));
  }
}

}  // namespace azure2
