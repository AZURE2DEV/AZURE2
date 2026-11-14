// The target effects of the data (lines of <targetInfo>): which points they apply to, the sub-point grid of each
// target-effect point (energy-loss integration, beam profile, UDR, Doppler, screening), and the re-placement of the
// adaptive grids for the current fit parameters.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "AdaptiveIntegrationGrid.h"
#include "CalculationPlan.h"
#include "MeshSlots.h"
#include "ShiftedFold.h"
#include "CNuc.h"
#include "Config.h"
#include "TargetLayer.h"
#include "WorkingCopy.h"
#include "RunDiagnostics.h"
#include "azure2/support/Environment.h"
#include "azure2/experiment/Kernels.h"

/*!
 * Reads the target effects input file and creates the TargetEffect objects
 * to be applied to the data.
 */

using azure2::EdgeAtFloor;
using azure2::MeetsStretches;
using azure2::MeshAnchorsOf;
using azure2::MeshPoint;
using azure2::MeshSlotEnergies;
using azure2::MeshTails;
using azure2::OnSharedMesh;
using azure2::QueueMeshPoint;
using azure2::SameAnchors;

namespace {
Unit::SubGridSpec MakeSubGridSpec(const AdaptiveIntegrationGrid::GridConfig &gc, double start, double end,
                                    const TargetEffect *effect, double conversion) {
  Unit::SubGridSpec spec;
  spec.valid = true;
  spec.startEnergy = start;
  spec.endEnergy = end;
  spec.maxPoints = gc.maxPoints;
  spec.entranceKey = gc.entranceKey;
  spec.baseEnergyStep = gc.baseEnergyStep;
  spec.resonanceWidthMultiplier = gc.resonanceWidthMultiplier;
  spec.pointsPerWidth = gc.pointsPerWidth;
  spec.targetIntegration = effect->IsTargetIntegration();
  spec.stoppingConversion = conversion;
  return spec;
}

/// What the fold of a point with an energy shift on the shared mesh needs of its layout (ShiftedFold.h): the c.m.
/// energy the support was laid out at, the integration floor, and for a data set whose shift is active its c.m. shift
/// at the start values and its reach.
void SetShiftLayout(Unit::SubGridSpec &spec, const std::vector<Segment> &segments, const Segment &segment,
                    const SegmentState &state, const Unit &point, double floor) {
  spec.centreEnergy = point.GetCMEnergy();
  spec.floorEnergy = floor;
  spec.shifted = azure2::HasEnergyShift(segment);
  if (!spec.shifted) return;
  spec.startShift = azure2::PointShift(segment, state, point.GetOriginalEnergy(), spec.stoppingConversion);
  spec.shiftReach = azure2::ShiftReach(segments, segment, state, point.GetOriginalEnergy(), spec.stoppingConversion,
                                       spec.startEnergy - spec.endEnergy);
}
}  // namespace

WorkingCopy::GridNucleus::GridNucleus() {}
WorkingCopy::GridNucleus::GridNucleus(const GridNucleus &) {}
WorkingCopy::GridNucleus &WorkingCopy::GridNucleus::operator=(const GridNucleus &) { return *this; }
WorkingCopy::GridNucleus::~GridNucleus() {}

namespace {
// Does this segment need external-capture amplitudes at its sub-points?  They are
// integrals computed once at initialization (CalculateECAmplitudes), so its grid
// cannot be re-placed.
bool SegmentNeedsEC(const TermRef &segment, CNuc *theCNuc, const Config &configure) {
  if (!(configure.paramMask & Config::USE_EXTERNAL_CAPTURE)) return false;
  const int ir = theCNuc->GetPairNumFromKey(segment.GetExitKey());
  for (int j = 1; j <= theCNuc->NumJGroups(); j++)
    for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
      const ALevel *level = theCNuc->GetJGroup(j)->GetLevel(la);
      if (level->IsECLevel() && level->GetECPairNum() == ir) return true;
    }
  return false;
}
}  // namespace

/*!
 * Re-places the adaptive sub-point grids (DeBoer 2026-10-19, Coulomb table stage 3).
 *
 * SetUpTargetEffects lays each grid out once, from the level energies and
 * widths in the input file.  A fit that moves a narrow level then integrates it
 * on a grid built for another position, and the objective depends on which file
 * the session started from (natural neon refit1: 25,773.96 from the start file,
 * 25,776.08 from the fitted one, same parameters).  Here the grid is laid out
 * again for the current parameters: the fit values are converted to observed
 * energies and widths, as a saved file would hold them (CNuc::TransformOut on a
 * private copy of the compound nucleus), and fed to the same generator with the
 * same window, so the grid of a parameter set is the one a session started from
 * its saved file would place.  Only points whose grid changed get new
 * sub-points; their Coulomb functions are table lookups.
 *
 * Energy shifts: a folded point on the shared mesh keeps its nodes; its support
 * moves with its shift, and when the shift leaves the reach its nodes were planned
 * for, the mesh is planned again for the moved support (ShiftedFold.h, question 2).
 * Every folded point is on a shared mesh (SetUpTargetEffects refuses one that has none).
 *
 * A grid is kept only if every node is bit for bit where the new one puts it;
 * otherwise the point gets the new grid, so the grid depends on the current
 * parameters alone (docs/upgrade/fit_threads.md).
 */
int WorkingCopy::ReplaceSubPointGrids(CNuc *theCNuc, const Config &configure, const std::vector<TermRef> &shiftActive,
                                      std::vector<int> &replacedTasks) {
  if (!theCNuc) return 0;
  // Without grid replacement (--no-grid-replacement) the anchors stay; the supports of shifted points still follow
  // their shifts (ShiftedFold.h, question 2), which their folds need.
  const bool replaceGrids = configure.replaceSubPointGrids;
  const auto t0 = std::chrono::steady_clock::now();
  // Per shared mesh: the supports whose centre moves now, with the grid record of their point.
  struct MovedSupport { const Unit::SubGridSpec *spec; double centre; const Unit *point; };
  std::vector<std::vector<MovedSupport>> movedSupports(Data().Meshes().size());

  // The points of a sharing group on a shared mesh (SlotGroup.h), by their units (UnitAt), which are re-placed
  // together -- of every term that reads the group's block.
  struct Item { TermRef seg; std::vector<int> group; };
  std::vector<Item> work;
  std::vector<std::pair<int, std::size_t>> blocks;  // (first slot of a block, its work item)
  auto collect = [&](const TermRef &seg) {
    bool any = false;
    for (int i = 1; i <= seg.NumPoints() && !any; i++) any = Data().GridSpecOf(*seg.GetPoint(*this, i)).valid;
    if (!any) return;
    const bool shifted = std::find(shiftActive.begin(), shiftActive.end(), seg) != shiftActive.end();
    // The centre of each shifted point's support on its mesh, a function of its c.m. shift alone.
    bool supportMoved = false;
    for (int i = 1; i <= seg.NumPoints(); i++) {
      const Unit *point = seg.GetPoint(*this, i);
      const Unit::SubGridSpec &spec = Data().GridSpecOf(*point);
      if (spec.mesh < 0 || (!shifted && MeshCentre(spec) == 0.0)) continue;
      const double delta = shifted ? azure2::PointShift(*seg.segment, State(*seg.segment), point->GetOriginalEnergy(),
                                                        spec.stoppingConversion)
                                   : 0.0;
      const double centre = azure2::ShiftCentre(delta, spec.startShift, spec.shiftReach);
      if (centre == MeshCentre(spec)) continue;
      movedSupports[(size_t)spec.mesh].push_back({&spec, centre, point});
      supportMoved = true;
    }
    // External capture: grids are not re-placed for a moved level (the set-up's integrals are kept), but the blocks of
    // a shifted point whose support moved are: their new sites take their capture amplitudes at their energies from the
    // tables, as a shifted site always has (Site::GetECAmplitudeWithShift).
    const bool ec = SegmentNeedsEC(seg, theCNuc, configure);
    if (ec && !supportMoved) {
      if (configure.FirstTime("EC grids not re-placed"))
        configure.outStream << "Note: the sub-point grids of segments with external capture are not re-placed "
                               "(their capture integrals are computed once)." << std::endl;
      return;
    }
    for (int i = 1; i <= seg.NumPoints(); i++) {
      Unit *point = seg.GetPoint(*this, i);
      if (!Data().GridSpecOf(*point).valid) continue;
      if (Data().GridSpecOf(*point).mesh < 0)
        throw std::logic_error("ReplaceSubPointGrids: a folded point of segment " + std::to_string(point->GetSegmentKey()) +
                               " is off the shared mesh (MeshSlots.h, OnSharedMesh)");
      // With the other points of the block it reads (a shifted point's slots do not move with its shift:
      // ShiftedFold.h).
      const int block = point->NumSubPoints() ? SlotBlockOf(point->SubPointFirst()) : -1;
      std::size_t at = 0;
      while (at < blocks.size() && (block < 0 || blocks[at].first != block)) at++;
      if (at == blocks.size()) {
        blocks.push_back({block, work.size()});
        work.push_back({seg, {}});
      }
      work[blocks[at].second].group.push_back(static_cast<int>(point - &this->UnitAt(0)));
    }
  };
  for (const Segment &seg : set_->segments_) collect({&seg, 0});
  for (const LeafRef &leaf : set_->leafRefs_) collect({&set_->segments_[leaf.first], leaf.second});
  if (work.empty()) return 0;

  // The observed parameters of the current fit values, on the private copy.
  CNuc *g = nullptr;
  if (replaceGrids) {
    if (!gridNuc_.nuc) gridNuc_.nuc.reset(theCNuc->Clone());
    g = gridNuc_.nuc.get();
    for (int j = 1; j <= theCNuc->NumJGroups(); j++)
      for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
        const ALevel *from = theCNuc->GetJGroup(j)->GetLevel(la);
        ALevel *to = g->MutableLevel(j, la);
        to->SetFitE(from->GetFitE());
        for (int ch = 1; ch <= theCNuc->GetJGroup(j)->NumChannels(); ch++) to->SetFitGamma(ch, from->GetFitGamma(ch));
      }
    g->FillCompoundFromParamsPhysical(g->GetTransformParams(g->TransformOut(configure)));
  }
  const auto t1 = std::chrono::steady_clock::now();

  // The shared meshes: the anchors of the current parameters and the supports of the shifted points at their present
  // centres; where they differ from the copy's, the new plan and the region where its nodes differ from the present
  // one's (for a moved support, the whole mesh: every group is compared node by node).
  const std::vector<DataSet::NucleusMesh> &meshes = Data().Meshes();
  std::vector<std::vector<azure2::MeshAnchor>> newAnchors(meshes.size());
  std::vector<azure2::MeshPlan> newPlans(meshes.size());
  std::vector<std::vector<std::pair<double, double>>> regions(meshes.size());
  std::vector<char> meshMoved(meshes.size(), 0);
  std::size_t supportsMoved = 0;
  for (size_t m = 0; m < meshes.size(); m++) {
    newAnchors[m] = replaceGrids ? MeshAnchorsOf(g, meshes[m], true) : meshAnchors_[m];
    const bool anchorsMoved = !SameAnchors(newAnchors[m], meshAnchors_[m]);
    if (!anchorsMoved && movedSupports[m].empty()) continue;
    meshMoved[m] = 1;
    for (const MovedSupport &moved : movedSupports[m]) {
      const std::size_t k = (std::size_t)moved.spec->meshSupport;
      const TargetEffect &effect = *Data().GetTargetEffect(moved.point->GetTargetEffectNum());
      meshSupports_[m][k] = azure2::PlannedShiftSupport(
          azure2::NominalMeshSupport(*moved.spec, meshes[m].supports[k].dataSet), *moved.spec, moved.centre,
          [&](double delta) { return azure2::SupportThicknessChange(*moved.point, *moved.spec, effect, configure, delta); });
      meshSupports_[m][k].label = meshes[m].supports[k].label;
      meshCentres_[m][k] = moved.centre;
      supportsMoved++;
    }
    newPlans[m] = azure2::PlanMesh(newAnchors[m], meshSupports_[m], MeshTails());
    if (movedSupports[m].empty())
      regions[m] = azure2::MeshChangedRegion(meshAnchors_[m], newAnchors[m], meshSupports_[m], MeshTails());
    else
      regions[m] = {{-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()}};
  }
  if (supportsMoved > 0) {
#pragma omp critical(azure_shift_replan)
    configure.outStream << "Note: the energy shifts moved " << supportsMoved << " folded point(s) beyond the reach their "
                        << "mesh nodes were planned for; the mesh was planned again for them (ShiftedFold.h)." << std::endl;
  }
  const bool planCheck = azure2::CalculationPlan::RunTimeCheckRequested();
  std::atomic<long> outsideRegionMoved{0};
  const auto t2 = std::chrono::steady_clock::now();

  std::atomic<int> changed{0};
  std::atomic<long> newSubPoints{0}, reusedSubPoints{0};
  // For AZURE_GRID_TIMING: per work item, -1 = a different number of nodes, else the
  // largest node move (MeV) of a grid that changed.
  std::vector<double> moved(work.size(), 0.0);
  SubPointGrids placed(work.size());  // per work item; the units never grow in the parallel loop
  // The sites of a grid (descending c.m. energies) for the points of a term (placed after the loop), each with the
  // angle, row and settings of `target`.  A node that stays where it was keeps its tables (Unit::AdoptEDependentValues
  // from `present`, the sites the points read now: the same energy, angle and pairs, so exactly what
  // CalcEDependentValues and CalcLegendreP would compute); they are usually still shared with the copy's source, so a
  // copy re-placing its grids adds tables only for the nodes that moved.
  auto laySites = [&](Unit *target, const TermRef &seg, const std::vector<double> &grid,
                      const UnitRange<Site> &present, long &reused) {
    const Unit::SubGridSpec &tspec = Data().GridSpecOf(*target);
    const TargetEffect *effect = Data().GetTargetEffect(seg.GetTargetEffectNum());
    const TargetEffect *qEffect = effect->IsQCoefficients() ? effect : nullptr;
    std::vector<Unit::EDependentShare> old;
    old.reserve(present.size());
    for (const Site &sp : present) old.push_back(sp.ShareEDependentValues());
    std::sort(old.begin(), old.end(),
              [](const Unit::EDependentShare &x, const Unit::EDependentShare &y) { return x.energy < y.energy; });
    std::vector<Site> fresh;
    fresh.reserve(grid.size());
    for (size_t k = 0; k < grid.size(); k++) {
      const double e = grid[k];
      Site sub(target->GetCMAngle(), e, seg.segment, seg.term);
      sub.BelongToRowOf(*target);
      if (tspec.targetIntegration)
        sub.SetStoppingPower(tspec.stoppingConversion *
                             effect->StoppingAt(e / tspec.stoppingConversion, configure,
                                                EffectValues(seg.GetTargetEffectNum()).stoppingScale));
      auto o = std::lower_bound(old.begin(), old.end(), sub.GetCMEnergy(),
                                [](const Unit::EDependentShare &x, double v) { return x.energy < v; });
      if (o != old.end() && sub.AdoptEDependentValues(*o)) {
        sub.CalcCoulombAmplitude(this, theCNuc);  // depends on the angle; cheap
        reused++;
      } else {
        sub.CalcLegendreP(this, configure.maxLOrder, theCNuc, qEffect);
        sub.CalcEDependentValues(this, theCNuc, configure);
        sub.CalcCoulombAmplitude(this, theCNuc);
      }
      target->MakeSubPoint(sub);
      fresh.push_back(std::move(sub));
    }
    return fresh;
  };
  // The points of a sharing group on a shared mesh: their slots move only if their nucleus's anchors changed and the
  // support (or the present slots) of one of them meets the region where the nodes move (MeshChangedRegion);
  // elsewhere the new plan has the same nodes bit for bit, which AZURE_PLAN_CHECK=1 checks.  The group's sites are
  // placed again together, on the union of its points' new ranges (SlotGroup.h).
  auto replaceGroup = [&](const Item &it, size_t q) {
    Unit *lead = &this->UnitAt(it.group.front());
    const Unit::SubGridSpec &spec = Data().GridSpecOf(*lead);
    const size_t m = (size_t)spec.mesh;
    if (!meshMoved[m]) return;
    bool inRegion = false;
    std::vector<int> supports;
    for (int u : it.group) {
      const Unit *point = &this->UnitAt(u);
      const Unit::SubGridSpec &pspec = Data().GridSpecOf(*point);
      const azure2::MeshSupport &support = meshSupports_[m][(size_t)pspec.meshSupport];
      double low = support.low, high = support.high;
      if (point->NumSubPoints() > 0) {
        low = std::min(low, point->GetSubPoint(*this, point->NumSubPoints())->GetOriginalEnergy() + pspec.frameShift);
        high = std::max(high, point->GetSubPoint(*this, 1)->GetOriginalEnergy() + pspec.frameShift);
      }
      inRegion = inRegion || MeetsStretches(regions[m], low, high);
      supports.push_back(pspec.meshSupport);
    }
    if (!inRegion && !planCheck) return;
    std::vector<std::pair<int, int>> ranges;
    const std::vector<std::size_t> nodes = azure2::GroupNodes(newPlans[m], supports, ranges);
    std::vector<double> grid(nodes.size());
    for (size_t k = 0; k < nodes.size(); k++) grid[k] = newPlans[m].nodes[nodes[k]] - spec.frameShift;
    // Whether each point's slots are at the new energies, bit for bit; the group's block as a whole.
    bool same = true;
    for (size_t j = 0; j < it.group.size(); j++) {
      const Unit *point = &this->UnitAt(it.group[j]);
      bool mine = point->NumSubPoints() == ranges[j].second;
      for (int k = 0; mine && k < ranges[j].second; k++)
        mine = point->GetSubPoint(*this, k + 1)->GetOriginalEnergy() == grid[ranges[j].first + k];
      if (!mine) {
        same = false;
        if (point->NumSubPoints() != ranges[j].second) moved[q] = -1.0;
        else if (moved[q] >= 0.0)
          for (int k = 0; k < ranges[j].second; k++)
            moved[q] = std::max(moved[q], std::fabs(point->GetSubPoint(*this, k + 1)->GetOriginalEnergy() -
                                                    grid[ranges[j].first + k]));
        if (inRegion) changed++;
      }
    }
    const int first = lead->NumSubPoints() ? SlotBlockOf(lead->SubPointFirst()) : NumSlots();
    int size = 0;
    while (first + size < NumSlots() && SlotBlockOf(first + size) == first) size++;
    if (!inRegion) {
      if (!same) outsideRegionMoved++;
      return;
    }
    if (same && size == (int)grid.size()) {
      for (int k = 0; same && k < size; k++) same = SlotAt(first + k).GetOriginalEnergy() == grid[k];
      if (same) return;
    }
    long reused = 0;
    azure2::SlotBlock block;
    block.sites = laySites(lead, it.seg, grid, Slots(first, size), reused);
    for (size_t j = 0; j < it.group.size(); j++)
      if (ranges[j].second > 0)
        block.members.push_back({it.group[j], ranges[j].first, ranges[j].second});
    placed[q].push_back(std::move(block));
    reusedSubPoints += reused;
    newSubPoints += (long)grid.size();
  };
#ifdef _OPENMP
  const bool serial = omp_in_parallel();
#else
  const bool serial = true;
#endif
  if (serial) {
    for (size_t q = 0; q < work.size(); q++) replaceGroup(work[q], q);
  } else {
#pragma omp parallel for schedule(dynamic, 4)
    for (int q = 0; q < (int)work.size(); q++)
      replaceGroup(work[q], (size_t)q);
  }
  if (changed > 0) {
    for (const auto &item : placed)
      for (const azure2::SlotBlock &block : item)
        for (const azure2::SlotMember &member : block.members) replacedTasks.push_back(this->UnitAt(member.unit).PlanTask());
    this->PlaceSubPointGrids(placed);
  }
  for (size_t m = 0; m < meshes.size(); m++)
    if (meshMoved[m]) meshAnchors_[m] = std::move(newAnchors[m]);
  if (outsideRegionMoved > 0)
    throw std::logic_error("AZURE_PLAN_CHECK: " + std::to_string(outsideRegionMoved.load()) +
                           " points on the shared mesh outside the changed region would get other nodes");

  if (azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kGridTiming)) {
    const int c = ++configure.diagnostics->gridPlacements;
    if (c <= 20 || c % 50 == 0) {
      auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
      int resized = 0;
      double maxMove = 0.0;
      for (double m : moved) {
        if (m < 0.0) resized++;
        else maxMove = std::max(maxMove, m);
      }
      std::cerr << "[grid] call " << c << ": " << work.size() << " grids, " << changed.load() << " re-placed ("
                << newSubPoints.load() << " sub-points; " << resized << " with a new node count, largest node move "
                << maxMove * 1.0e6 << " eV otherwise; " << reusedSubPoints.load()
                << " sub-points kept their tables); observed parameters " << ms(t0, t1) << " ms, mesh plans "
                << ms(t1, t2) << " ms (" << std::count(meshMoved.begin(), meshMoved.end(), 1) << " of "
                << meshes.size() << " moved), grids " << ms(t2, std::chrono::steady_clock::now()) << " ms"
                << std::endl;
    }
  }
  return changed.load();
}

int WorkingCopy::SetUpTargetEffects(const Config &configure, CNuc *compound,
                              const azure2::experiment::SpecList<azure2::experiment::TargetEffectSpec> &specs) {
  int activeLines = 0;
  for (const auto &item : specs.items) {
    const std::string &line = item.spec.source;
    {
      TargetEffect targetEffect(item.spec, configure);
      if (targetEffect.HasScreeningError()) {
        std::cout << "ERROR: <targetInt> line \"" << line << "\": a screening block is 'screening Ue_eV vary sigma_eV'"
                  << std::endl;
        return -1;
      }
      if (targetEffect.HasDopplerError()) {
        std::cout << "ERROR: <targetInt> line \"" << line << "\": " << targetEffect.GetDopplerErrorText() << std::endl;
        return -1;
      }
      if (targetEffect.HasKernelError()) {
        std::cout << "ERROR: <targetInt> line \"" << line << "\": " << targetEffect.GetKernelErrorText() << std::endl;
        return -1;
      }
      if (targetEffect.HasDensityFitError()) {
        std::cout << "ERROR: <targetInt> line \"" << line << "\": " << targetEffect.GetDensityFitErrorText() << std::endl;
        return -1;
      }
      if (item.readFailed) return -1;
      if (targetEffect.IsActive()) {
        targetEffect.SetLineNumber(++activeLines);
        // Every listed segment gets its own copy of the effect.  The effect's
        // energies (sigma, beam profile) are converted lab -> c.m. below with
        // the segment's own factor; one object shared by "1,2" was converted
        // once per listed segment, and both segments were then folded with
        // the twice-converted sigma (17 keV c.m. became 16.1 keV for 18O+p).
        // A copy per segment also lets segments of different entrance pairs
        // share one line.
        std::vector<int> segmentsList = targetEffect.GetSegmentsList();
        for (int i = 1; i <= segmentsList.size(); i++) {
          if (Draft().IsSegmentKey(segmentsList[i - 1])) {
            Segment *segment = Draft().GetSegmentFromKey(segmentsList[i - 1]);
            if (segment && segment->IsTHM() && targetEffect.IsTargetIntegration()) {
              // A THM energy is the x + A energy reconstructed from the
              // three-body kinematics, not a beam energy slowing down in the
              // target; its sub-points also go below the entrance threshold,
              // where a stopping power is evaluated at E < 0 (NaN).
              configure.outStream << "ERROR: <targetInt>: segment " << segmentsList[i - 1]
                                  << " is a THM segment: target integration does not apply to it "
                                     "(use the Gaussian convolution for the THM energy resolution)."
                                  << std::endl;
              return -1;
            }
            if (segment) {
              this->AddTargetEffect(targetEffect);
              segment->SetTargetEffectNum(Draft().NumTargetEffects());
              // If segment has components, set the target effect to these as
              // well -- again a copy each, since each is converted below.
              for (int t = 1; t < segment->NumTerms(); t++) {
                this->AddTargetEffect(targetEffect);
                segment->SetTermTargetEffectNum(t, Draft().NumTargetEffects());
              }
            }
          }
        }
      }
    }
  }
  if (specs.failed) return -1;
  if (this->SetUpDensityParams(configure) == -1) return -1;
  this->SetUpStoppingScaleParams(configure);
  // The working copy starts from each effect's values as the set-up leaves them (a density-group member takes its
  // group's density).
  for (int n = 1; n <= Draft().NumTargetEffects(); n++) {
    const TargetEffect &te = *Draft().GetTargetEffect(n);
    effectValues_[n - 1] = TargetEffectValues{te.GetScreeningUe(), te.GetDensity(), te.GetDopplerTemperature(), te.GetStoppingScaleNominal()};
  }
  // a stopping from the library (evaluation, table): its data and its value at every beam energy, before any is used
  if (!this->CheckTargetStoppings(configure)) return -1;
  Draft().meshes_.clear();
  std::vector<MeshPoint> meshPending;  // the folded points whose slots the shared mesh gives, in set-up order
  // The folded points of one term (the base segment, t = 0, or a component), each with the pairs, the lab<->c.m.
  // factor and the kinematics of the pair the term calculates: their target thickness, the support of their kernels,
  // and their place in the shared mesh of their nucleus, whose nodes are their slots (placed below).
  auto layOutTerm = [&](const TermRef &term, const PPair *entrance, const PPair *exit, double conversion) -> int {
    const Segment &segment = *term.segment;
    TargetEffect *targetEffect = Draft().GetTargetEffect(term.GetTargetEffectNum());
    const std::string what = std::string(term.term ? "a component of " : "") + "segment #" +
                             std::to_string(segment.GetSegmentKey());
    // the line's energies stay in the lab frame; every point converts them with its own pair's factor
    const int bohrZ1 = (entrance->GetPType() == 20) ? 0 : entrance->GetZ(1);
    const int bohrZ2 = entrance->GetZ(2);
    if (targetEffect->IsUdr()) {
      double eMin = 1.0e300, eMax = 0.0;
      for (UnitIterator point = term.GetPoints(*this).begin(); point < term.GetPoints(*this).end(); point++) {
        eMin = std::min(eMin, point->GetLabEnergy());
        eMax = std::max(eMax, point->GetLabEnergy());
      }
      if (eMax > 0.0 && !targetEffect->LoadUdrTable(eMin, eMax, configure.outStream, configure.cacheInputFiles)) return -1;
    }
    for (UnitIterator point = term.GetPoints(*this).begin(); point < term.GetPoints(*this).end(); point++) {
      point->SetTargetEffectNum(term.GetTargetEffectNum());
      point->SetTargetKinematics(conversion, bohrZ1, bohrZ2);
      if (!targetEffect->IsSubPointEffect()) continue;
      double targetThickness = 0.0;
      if (targetEffect->IsTargetIntegration()) {
        // an infinitely thick target: the point integrates from its energy down to 1 keV; otherwise the grid spans
        // the target at its layout density (a fitted density's largest expected value), and a layer yield's grid a
        // margin beyond it (its E_out comes from the grid's own quadrature)
        targetThickness = targetEffect->SupportThickness(point->GetLabEnergy(), point->GetCMEnergy(), conversion, configure);
        if (!targetEffect->IsThickTarget()) {
          if (!(targetThickness > 0.0) && targetEffect->IsLayerYield())
            return RefuseLayerAtBeamEnergy(*targetEffect, *point, segment.GetDataFile(), configure);
          if (!(targetThickness > 0.0)) {
            configure.outStream << "ERROR: target integration of " << what << " at E_lab = " << point->GetLabEnergy()
                                << " MeV: the stopping cross section times the density is " << targetThickness / conversion
                                << " (not positive).  A density of 0 means an infinitely thick target; otherwise the"
                                << " stopping power must be positive at every point." << std::endl;
            return -1;
          }
        }
        point->SetTargetLabEnergy(point->GetLabEnergy());
        point->SetTargetThickness((targetEffect->HasDensityFit() || targetEffect->IsLayerYield())
                                      ? conversion * targetEffect->TargetThickness(point->GetLabEnergy(), configure)
                                      : targetThickness);
      }
      // the term's own entrance pair (an advanced component can have another entrance channel)
      const double minSubEnergy = segment.IntegrationFloor(entrance);
      // What the kernels read beyond the effect's parameters: the point's UDR kernel, from the effect's table, and the
      // inverse kinematics of a beam profile's detailed-balance weight, Q = E_gamma - E_cm: the entrance pair's threshold
      // (S plus its total excitation) minus the energy of the final level the gamma pair feeds (its heavy excitation).
      azure2::experiment::SupportPoint support;
      if (targetEffect->IsUdr()) {
        const double mass = entrance->GetM(1) * uconv;
        DataSet::PointSetUp &setUp = Draft().SetUpFor(*point);
        setUp.udr = targetEffect->BuildUdrKernel(point->GetLabEnergy(), mass);
        setUp.udrMass = mass;
        support.udr = setUp.udr.get();
        support.udrMass = mass;
      }
      if (targetEffect->IsBeamProfile())
        point->SetPhotoKinematics(entrance->GetSepE() + entrance->GetExE() - exit->GetHeavyExE(),
                                  (entrance->GetM(1) + entrance->GetM(2)) * uconv);
      // The interval the nodes span: the support of the effect's kernels (azure2/experiment/Kernels.h).
      support.energy = point->GetCMEnergy();
      support.labEnergy = point->GetLabEnergy();
      support.cmFactor = conversion;
      support.z1 = bohrZ1;
      support.z2 = bohrZ2;
      support.targetThickness = targetThickness;
      support.minimumEnergy = minSubEnergy;
      support.thm = segment.IsTHM();
      support.binLow = point->GetBinLowCM();
      support.binHigh = point->GetBinHighCM();
      double startEnergy = 0.0, endEnergy = 0.0;
      azure2::experiment::EffectSupport(targetEffect->Kernels(), targetEffect->KernelParametersFor(configure), support,
                                        startEnergy, endEnergy);
      // AZURE_REFINE_QUADRATURE=k: k times the smooth intervals and the points per width
      const int numPoints = targetEffect->NumSubPoints() * configure.refineQuadrature;
      AdaptiveIntegrationGrid::GridConfig gridConfig;
      gridConfig.maxPoints = numPoints;
      // The term's own entrance pair: its resonances, on its energy axis, are what the mesh must resolve.
      gridConfig.entranceKey = term.GetEntranceKey();
      gridConfig.inputWidthsArePhysical = (configure.paramMask & Config::TRANSFORM_PARAMETERS) && !compound->IsTransformedIn();
      gridConfig.baseEnergyStep = (startEnergy - endEnergy) / numPoints;
      gridConfig.resonanceWidthMultiplier = targetEffect->GetResonanceWidthMultiplier();
      gridConfig.pointsPerWidth = targetEffect->GetPointsPerWidth() * configure.refineQuadrature;
      Unit::SubGridSpec &spec = Draft().SetUpFor(*point).grid;
      spec = MakeSubGridSpec(gridConfig, startEnergy, endEnergy, targetEffect, conversion);
      spec.supportThickness = targetThickness;
      SetShiftLayout(spec, Draft().GetSegments(), segment, State(segment), *point, minSubEnergy);
      // Every folded point is a support of the shared mesh; a pair without an E_x frame has no place on it.
      if (!OnSharedMesh(segment, entrance, exit)) {
        configure.outStream << "ERROR: " << what << ": its points are folded (a target, beam or detector effect), but its"
                            << " pairs give no excitation energy of the compound nucleus (a beta-delayed spectrum without an exit pair), so the fold has"
                            << " no nodes." << std::endl;
        return -1;
      }
      const Unit &queued = *point;
      meshPending.push_back(QueueMeshPoint(
          Draft().meshes_, *point, term, spec, entrance, exit, EdgeAtFloor(*targetEffect, endEnergy, minSubEnergy),
          [&](double delta) { return azure2::SupportThicknessChange(queued, spec, *targetEffect, configure, delta); }));
    }
    return 0;
  };
  for (SegmentIterator segment = Draft().GetSegments().begin(); segment < Draft().GetSegments().end(); segment++) {
    const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(segment->GetEntranceKey()));
    const PPair *exitPair = compound->GetPairFromKey(segment->GetExitKey());  // nullptr for total capture (-1)
    if (segment->IsTargetEffect()) {
      const TargetEffect *targetEffect = Draft().GetTargetEffect(segment->GetTargetEffectNum());
      if (targetEffect->IsDoppler() && entrancePair->GetPType() != 0) {
        configure.outStream << "ERROR: segment #" << segment->GetSegmentKey() << ": free-gas Doppler broadening"
                            << " (doppler) needs a particle entrance pair (a beam on target atoms in thermal motion)"
                            << std::endl;
        return -1;
      }
      // A beta-delayed spectrum has no beam: a target integration or a beam profile would fold the decay energy with
      // a beam's energy loss, screening would use the beta pair's charges.  A detector resolution (Gaussian, UDR) applies.
      if (entrancePair->GetPType() == 20 &&
          (targetEffect->IsTargetIntegration() || targetEffect->IsBeamProfile() || targetEffect->HasScreening())) {
        configure.outStream << "ERROR: segment #" << segment->GetSegmentKey() << ": a beta-delayed spectrum (entrance"
                            << " pair of type 20) takes a detector resolution only (Gaussian convolution, UDR), not a"
                            << " target integration, a beam profile or electron screening." << std::endl;
        return -1;
      }
      if (layOutTerm({&*segment, 0}, entrancePair, exitPair, Segment::CMEnergyFactor(entrancePair, exitPair)) != 0)
        return -1;
    }
    for (int t = 1; t < segment->NumTerms(); t++) {
      const TermRef component{&*segment, t};
      // The component's grid, kinematics and lab<->c.m. factor belong to the
      // pair the component calculates -- for a natural target another
      // compound nucleus with its own masses -- not to the base segment's.
      // The points were copied from the base segment with its c.m.
      // energies; convert them for this pair before the grid is laid out
      // (WorkingCopy::Initialize repeats this, idempotently, from the lab energy).
      const PPair *compEntrance = compound->GetPair(compound->GetPairNumFromKey(component.GetEntranceKey()));
      const PPair *compExit = compound->GetPairFromKey(component.GetExitKey());
      for (int i = 0; i < component.NumPoints(); i++) segment->ConvertRowEnergyToCM(*this, t, i, compEntrance, compExit);
      if (segment->IsTargetEffect() &&
          layOutTerm(component, compEntrance, compExit, Segment::CMEnergyFactor(compEntrance, compExit)) != 0)
        return -1;
    }
  }
  // The shared mesh of each compound nucleus: one plan from its supports and the anchors of the input levels, and
  // every point queued above takes the mesh's nodes of its range as its slots.
  meshAnchors_.assign(Draft().meshes_.size(), {});
  meshSupports_.assign(Draft().meshes_.size(), {});
  meshCentres_.assign(Draft().meshes_.size(), {});
  for (size_t m = 0; m < Draft().meshes_.size(); m++) {
    DataSet::NucleusMesh &mesh = Draft().meshes_[m];
    mesh.anchors = MeshAnchorsOf(compound, mesh,
                                 (configure.paramMask & Config::TRANSFORM_PARAMETERS) && !compound->IsTransformedIn());
    meshAnchors_[m] = mesh.anchors;
    meshSupports_[m] = mesh.supports;
    meshCentres_[m].assign(mesh.supports.size(), 0.0);
    const azure2::MeshPlan plan = azure2::PlanMesh(mesh.anchors, mesh.supports, MeshTails());
    azure2::AddMeshSlotGroups(*this, plan, (int)m, meshPending, configure);
  }
  for (const MeshPoint &q : meshPending) {  // the centres of the shifted supports: their c.m. shifts at the start
    const Unit::SubGridSpec &spec = Draft().GridSpecOf(*q.point);
    if (spec.shifted) meshCentres_[(size_t)spec.mesh][(size_t)spec.meshSupport] = spec.startShift;
  }
  this->CommitSubPoints();  // the sub-points given above take their places in the units
  if (!this->CheckSubPointStoppings(configure)) return -1;  // and at every energy of the sub-point grids
  if (!CheckLayerYields(*this, configure)) return -1;  // a yield on a layer: its stopping and E_out on the grid
  // One screening parameter per <targetInt> line, shared by the copies of its effect.
  Draft().screeningEffects_.clear();
  Draft().screeningCopies_.clear();
  for (int n = 1; n <= Draft().NumTargetEffects(); n++) {
    const TargetEffect *te = Draft().GetTargetEffect(n);
    if (!te->HasScreening()) continue;
    size_t k = 0;
    while (k < Draft().screeningEffects_.size() && Draft().screeningEffects_[k] != te->GetLineNumber()) k++;
    if (k == Draft().screeningEffects_.size()) {
      Draft().screeningEffects_.push_back(te->GetLineNumber());
      Draft().screeningCopies_.emplace_back();
    }
    Draft().screeningCopies_[k].push_back(n);
  }
  this->ApplyScreening();
  // One Doppler parameter per <targetInt> line with the fit form, shared by the copies of its effect (as screening).
  Draft().dopplerEffects_.clear();
  Draft().dopplerCopies_.clear();
  for (int n = 1; n <= Draft().NumTargetEffects(); n++) {
    const TargetEffect *te = Draft().GetTargetEffect(n);
    if (!te->HasDopplerFit()) continue;
    size_t k = 0;
    while (k < Draft().dopplerEffects_.size() && Draft().dopplerEffects_[k] != te->GetLineNumber()) k++;
    if (k == Draft().dopplerEffects_.size()) {
      Draft().dopplerEffects_.push_back(te->GetLineNumber());
      Draft().dopplerCopies_.emplace_back();
    }
    Draft().dopplerCopies_[k].push_back(n);
  }
  return 0;
}
