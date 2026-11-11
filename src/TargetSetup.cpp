// The target effects of the data (lines of <targetInfo>): which points they apply to, the sub-point grid of each
// target-effect point (energy-loss integration, beam profile, UDR, Doppler, screening), the energy windows of a
// restricted integration, and the re-placement of the adaptive grids for the current fit parameters.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "AdaptiveIntegrationGrid.h"
#include "CNuc.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "UdrKernel.h"

/*!
 * Reads the target effects input file and creates the TargetEffect objects
 * to be applied to the data.
 */

namespace {
/*!
 * The sub-points of a thin-target point outside the windows of an energy-restricted integration
 * (Unit::SetTargetFallbackOnly): one at E_bar = E - D for a fixed density; for a fitted one three, at E_bar of the
 * smallest density the fit may reach (0.01 N0, WorkingCopy::FillDensitiesFromParams), of the nominal N0 and of the layout
 * density, through which the kernel interpolates E_bar(N) -- exact at N0, smooth in N.
 */
static void AddTargetFallbackNodes(WorkingCopy &data, Unit *point, TargetEffect *te, double cmConversion,
                                   const Config &configure, const TermRef &segment) {
  const double e = point->GetCMEnergy();
  std::vector<double> nodes;
  if (te->HasDensityFit()) {
    const double n0 = te->GetDensityNominal(), eLab = point->GetLabEnergy();
    const double dMin = cmConversion * te->MeanEnergyLoss(eLab, configure, 0.01 * n0);
    const double d0 = cmConversion * te->MeanEnergyLoss(eLab, configure, n0);
    double dMax = cmConversion * te->MeanEnergyLoss(eLab, configure, te->GetLayoutDensity());
    if (!(dMax > d0)) dMax = 1.5 * d0;
    nodes = {e - std::min(dMin, 0.5 * d0), e - d0, e - dMax};
  } else {
    nodes = {e - point->GetTargetFallbackDepth()};
  }
  for (double sub : nodes) {
    Unit subPoint(point->GetCMAngle(), sub, segment.segment, segment.term);
    subPoint.BelongToRowOf(*point);
    subPoint.SetStoppingPower(cmConversion * te->GetStoppingPowerEq()->Evaluate(configure, sub / cmConversion));
    data.AddSubPoint(*point, subPoint);
  }
}

/*!
 * Refuses the combinations the sigma(E_bar) fallback of a windowed target integration does not cover.
 */
static bool CheckTargetWindows(TargetEffect *te, int segmentKey, std::ostream &out) {
  if (te->IsAutoWindows() && (!te->IsTargetIntegration() || te->IsThickTarget())) {
    out << "ERROR: segment #" << segmentKey << ": automatic energy windows (\"auto\") need a THIN-target integration: "
        << "they decide where the depth-averaged energy is good enough, and only a thin target has one." << std::endl;
    return false;
  }
  if (!te->IsTargetIntegration()) return true;
  const bool windows = te->HasWindows();
  if (te->IsThickTarget() && windows) {
    out << "ERROR: segment #" << segmentKey << ": energy windows on a THICK-target integration (density 0): a thick-"
        << "target yield has no single equivalent energy to evaluate outside the windows.  Give the density or "
        << "drop the windows." << std::endl;
    return false;
  }
  if (te->IsDoppler() && !te->IsThickTarget() && (windows || te->GetAutoTolerance() > 0.)) {
    out << "ERROR: segment #" << segmentKey << ": energy windows or an automatic tolerance on a target integration "
        << "with free-gas Doppler broadening (doppler) are not supported." << std::endl;
    return false;
  }
  return true;
}

/*!
 * The lab-energy windows of one segment (or component) of a target-effect line: the line's own ranges, plus,
 * with "auto" (TargetEffect::IsAutoWindows, 2026-11-06), windows placed from the level scheme of the
 * segment's entrance pair.
 *
 * Outside a window a thin-target point is sigma at its depth-averaged energy E_bar, whose error is about
 * |sigma''| Delta^2 / 24.  For a resonance of total width Gamma (a Lorentzian) that is Delta^2 / (3 Gamma^2) of
 * the peak, so a resonance is "narrow" -- it needs the integration -- when Gamma < Delta / sqrt(3 tol), with
 * Delta the target thickness at the resonance (lab, at the layout density) plus twice the beam spread.  A point
 * at E0 samples [E0 - Delta, E0]; its window reaches margin x max(Gamma, Delta) beyond the points whose target
 * holds the resonance: [E_R - w, E_R + Delta + w].  The widths are the observed total widths the sub-point
 * grids are laid out with (AdaptiveIntegrationGrid::GetResonances, from the parameters the session starts
 * from).  The windows are fixed for the session: a fit that moves a narrow level by more than w leaves it
 * outside; a new session from the fitted file places them again.
 */
static std::vector<std::pair<double, double> > TargetWindowsFor(WorkingCopy &data, TargetEffect *te, CNuc *compound, const TermRef &segment,
                                                                 double cmConversion, const Config &configure,
                                                                 int reportKey, bool report) {
  std::vector<std::pair<double, double> > windows = te->GetRanges();
  if (!te->IsAutoWindows() || !compound || segment.NumPoints() == 0) return windows;
  double eLo = 1.0e300, eHi = 0.0;
  for (UnitIterator point = segment.GetPoints(data).begin(); point < segment.GetPoints(data).end(); point++) {
    eLo = std::min(eLo, point->GetLabEnergy());
    eHi = std::max(eHi, point->GetLabEnergy());
  }
  AdaptiveIntegrationGrid::GridConfig gc;
  gc.entranceKey = segment.GetEntranceKey();
  gc.inputWidthsArePhysical = (configure.paramMask & Config::TRANSFORM_PARAMETERS) && !compound->IsTransformedIn();
  AdaptiveIntegrationGrid finder(gc);
  // every level up to well above the data: a window reaches below its resonance only by w
  const std::vector<AdaptiveIntegrationGrid::ResonanceInfo> resonances =
      finder.GetResonances(cmConversion * eHi + 10.0, 0.0, compound);
  const double tol = te->GetAutoWindowTolerance(), margin = te->GetAutoWindowMargin();
  std::vector<std::pair<double, double> > found;
  int narrow = 0;
  for (const auto &r : resonances) {
    const double eR = r.energy / cmConversion;
    if (!(eR > 0.0) || !(r.totalWidth > 0.0)) continue;
    // the target and the spread as the segment's own points see them: the stopping power is an input fit over
    // the data's energies and means nothing far outside them
    const double eAt = std::min(std::max(eR, eLo), eHi);
    double delta = te->TargetThickness(eAt, configure, te->GetLayoutDensity());
    if (!(delta > 0.0) || !(delta < eAt)) continue;
    if (te->IsGaussian()) delta += 2.0 * te->SigmaCM(cmConversion * eAt, configure, cmConversion) / cmConversion;
    const double gamma = r.totalWidth / cmConversion;
    if (!(gamma < delta / std::sqrt(3.0 * tol))) continue;
    const double w = margin * std::max(gamma, delta);
    const double lo = eR - w, hi = eR + delta + w;
    if (hi < eLo || lo > eHi) continue;
    narrow++;
    found.push_back(std::make_pair(lo, hi));
  }
  // Particle thresholds of the same compound nucleus: a target that straddles one sees the cusp of the channel
  // opening, which no single energy reproduces (13C+alpha ND2021r: 15 % at the 16O(6.05) n1 threshold).
  int thresholds = 0;
  PPair *entrance = compound->GetPairFromKey(segment.GetEntranceKey());
  if (entrance) {
    const double eEntrance = entrance->GetSepE() + entrance->GetExE();
    for (int i = 1; i <= compound->NumPairs(); i++) {
      PPair *pair = compound->GetPair(i);
      if (pair == entrance || pair->GetPType() != 0 || pair->GetNucleus() != entrance->GetNucleus()) continue;
      const double eTh = (pair->GetSepE() + pair->GetExE() - eEntrance) / cmConversion;
      if (!(eTh > 0.0)) continue;
      const double eAt = std::min(std::max(eTh, eLo), eHi);
      double delta = te->TargetThickness(eAt, configure, te->GetLayoutDensity());
      if (!(delta > 0.0) || !(delta < eAt)) continue;
      if (te->IsGaussian()) delta += 2.0 * te->SigmaCM(cmConversion * eAt, configure, cmConversion) / cmConversion;
      const double lo = eTh - margin * delta, hi = eTh + delta + margin * delta;
      if (hi < eLo || lo > eHi) continue;
      thresholds++;
      found.push_back(std::make_pair(lo, hi));
    }
  }
  windows.insert(windows.end(), found.begin(), found.end());
  std::sort(windows.begin(), windows.end());
  std::vector<std::pair<double, double> > merged;
  for (const auto &w : windows) {
    if (!merged.empty() && w.first <= merged.back().second)
      merged.back().second = std::max(merged.back().second, w.second);
    else
      merged.push_back(w);
  }
  int inside = 0;
  for (UnitIterator point = segment.GetPoints(data).begin(); point < segment.GetPoints(data).end(); point++)
    if (te->BlendWeight(point->GetLabEnergy(), merged) > 0.0) inside++;
  std::ostringstream note;
  note << "Note: segment #" << reportKey << ": automatic target-integration windows (tolerance " << tol
       << ", margin " << margin << "): " << narrow << " narrow level(s), " << thresholds << " threshold(s), " << inside << " of " << segment.NumPoints()
       << " points integrated";
  if (!merged.empty()) {
    note << "; windows (lab MeV)";
    for (size_t i = 0; i < merged.size(); i++)
      note << (i ? ", " : " ") << std::setprecision(6) << merged[i].first << "-" << merged[i].second;
  }
  if (report) configure.Note(note.str());
  return merged;
}

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
      ALevel *level = theCNuc->GetJGroup(j)->GetLevel(la);
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
 * Energy shifts: the sub-points are laid out at the unshifted energies and moved
 * with the data by the shift step; the resonances are offset by the point's c.m.
 * shift so the dense part of the grid lands on each resonance after the move.
 *
 * A grid is kept only if every node is bit for bit where the new one puts it;
 * otherwise the point gets the new grid, so the grid depends on the current
 * parameters alone (docs/upgrade/fit_threads.md).
 */
int WorkingCopy::ReplaceSubPointGrids(CNuc *theCNuc, const Config &configure, const std::vector<TermRef> &shiftActive,
                                std::vector<TermRef> &reshift) {
  if (!configure.replaceSubPointGrids || !configure.useAdaptiveGrid || !theCNuc) return 0;
  const auto t0 = std::chrono::steady_clock::now();

  struct Item { TermRef seg; int i; double offset; bool shifted; };
  std::vector<Item> work;
  struct Owner { TermRef seg; bool shifted; bool replaceable; };
  std::unordered_map<const Unit *, Owner> mappedSegment;  // the segment of each mapped point
  static std::atomic<bool> ecNoted{false};
  auto collect = [&](const TermRef &seg) {
    // Mapped points are registered first, also in segments whose grids are not
    // re-placed: a host with such a mapped point keeps its grid too.
    for (int i = 1; i <= seg.NumPoints(); i++)
      if (seg.GetPoint(*this, i)->IsMapped()) mappedSegment[seg.GetPoint(*this, i)] = {seg, false, false};
    bool any = false;
    for (int i = 1; i <= seg.NumPoints() && !any; i++) any = Data().GridSpecOf(*seg.GetPoint(*this, i)).valid;
    if (!any) return;
    if (SegmentNeedsEC(seg, theCNuc, configure)) {
      if (!ecNoted.exchange(true))
        configure.outStream << "Note: the sub-point grids of segments with external capture are not re-placed "
                               "(their capture integrals are computed once)." << std::endl;
      return;
    }
    const bool shifted = std::find(shiftActive.begin(), shiftActive.end(), seg) != shiftActive.end();
    PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(seg.GetEntranceKey()));
    PPair *exitPair = theCNuc->GetPairFromKey(seg.GetExitKey());  // nullptr for total capture (-1)
    // The c.m./lab factor the shift step uses for sub-points (Segment::UpdatePointEnergyWithShift).
    const double cmFactor = Segment::CMEnergyFactor(entrancePair, exitPair);
    for (int i = 1; i <= seg.NumPoints(); i++) {
      Unit *point = seg.GetPoint(*this, i);
      if (!Data().GridSpecOf(*point).valid) continue;
      // A mapped point keeps its host's grid: the mapping (MapData, MapLeaves) was made because the
      // two grids are identical, and the host's calculation fills the mapped point's sub-points node by node
      // (Unit::Calculate).  It is re-placed with its host below, never on its own.
      if (point->IsMapped()) {
        mappedSegment[point] = {seg, shifted, true};
        continue;
      }
      // the parent's shift, which is what moves every sub-point (Segment::UpdatePointEnergyWithShift)
      const double offset = shifted ? cmFactor * seg.segment->TotalEnergyShift(State(*seg.segment), point->GetOriginalEnergy()) : 0.0;
      work.push_back({seg, i, offset, shifted});
    }
  };
  for (const Segment &seg : set_->segments_) collect({&seg, 0});
  for (const LeafRef &leaf : set_->leafRefs_) collect({&set_->segments_[leaf.first], leaf.second});
  if (work.empty()) return 0;

  // The observed parameters of the current fit values, on the private copy.
  if (!gridNuc_.nuc) gridNuc_.nuc.reset(theCNuc->Clone());
  CNuc *g = gridNuc_.nuc.get();
  for (int j = 1; j <= theCNuc->NumJGroups(); j++)
    for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
      ALevel *from = theCNuc->GetJGroup(j)->GetLevel(la);
      ALevel *to = g->GetJGroup(j)->GetLevel(la);
      to->SetFitE(from->GetFitE());
      for (int ch = 1; ch <= theCNuc->GetJGroup(j)->NumChannels(); ch++) to->SetFitGamma(ch, from->GetFitGamma(ch));
    }
  g->TransformOut(configure);
  g->FillCompoundFromParamsPhysical(g->GetTransformParams(configure));
  const auto t1 = std::chrono::steady_clock::now();

  std::atomic<int> changed{0};
  std::atomic<long> newSubPoints{0}, reusedSubPoints{0}, sharedSubPoints{0};
  // For AZURE_GRID_TIMING: per work item, -1 = a different number of nodes, else the
  // largest node move (MeV) of a grid that changed.
  std::vector<double> moved(work.size(), 0.0);
  std::mutex reshiftMutex;
  std::vector<TermRef> mappedShifted;  // shifted segments of mapped points re-placed with their hosts
  SubPointGrids placed(work.size());  // per work item; the units never grow in the parallel loop
  auto replace = [&](const Item &it, size_t q) {
    Unit *point = it.seg.GetPoint(*this, it.i);
    const Unit::SubGridSpec &spec = Data().GridSpecOf(*point);
    AdaptiveIntegrationGrid::GridConfig gc;
    gc.maxPoints = spec.maxPoints;
    gc.entranceKey = spec.entranceKey;
    gc.inputWidthsArePhysical = true;
    gc.baseEnergyStep = spec.baseEnergyStep;
    gc.resonanceWidthMultiplier = spec.resonanceWidthMultiplier;
    gc.pointsPerWidth = spec.pointsPerWidth;
    gc.resonanceOffset = it.offset;
    const azure2::MappedPoints mapped = point->PlanMapped(this);
    for (int m = 0; m < mapped.size(); m++) {
      auto own = mappedSegment.find(mapped[m]);
      if (own == mappedSegment.end() || !own->second.replaceable) return;  // its mapped point keeps the grid
    }
    std::vector<double> grid = AdaptiveIntegrationGrid(gc).GenerateGrid(spec.startEnergy, spec.endEnergy, g);
    const UnitRange<Unit> subs = point->GetSubPoints(*this);
    // The grid kept is the one placed now, bit for bit.  A grid kept within a tolerance would depend on the
    // parameter vectors the working copy evaluated before, which Minuit's parallel gradient hands out by thread:
    // the objective at one p, and so a fit, would depend on the thread count (docs/upgrade/fit_threads.md).
    bool same = ((size_t)subs.size() == grid.size());
    for (size_t k = 0; same && k < grid.size(); k++) same = subs[k].GetOriginalEnergy() == grid[k];
    if (same) return;
    if ((size_t)subs.size() != grid.size())
      moved[q] = -1.0;
    else
      for (int k = 0; k < subs.size(); k++) moved[q] = std::max(moved[q], std::fabs(subs[k].GetOriginalEnergy() - grid[k]));
    // Lays a point's sub-points out on the grid (placed after the loop).  A node that stays where it was keeps its
    // sub-point's tables (Unit::AdoptEDependentValues: the same energy, angle and pairs, so exactly what
    // CalcEDependentValues and CalcLegendreP would compute); they are usually still shared with the copy's source, so
    // a copy re-placing its grids adds tables only for the nodes that moved.
    long reused = 0, fromOtherCopies = 0;
    auto layOut = [&](Unit *target, const TermRef &seg, bool shifted) {
      const Unit::SubGridSpec &tspec = Data().GridSpecOf(*target);
      const TargetEffect *effect = Data().GetTargetEffect(seg.GetTargetEffectNum());
      const TargetEffect *qEffect = effect->IsQCoefficients() ? effect : nullptr;
      const UnitRange<Unit> tsubs = target->GetSubPoints(*this);
      std::vector<Unit::EDependentShare> old;
      if (!shifted) {
        old.reserve(tsubs.size());
        for (const Unit &sp : tsubs) old.push_back(sp.ShareEDependentValues());
        std::sort(old.begin(), old.end(),
                  [](const Unit::EDependentShare &x, const Unit::EDependentShare &y) { return x.energy < y.energy; });
      }
      std::vector<Unit> fresh;
      for (size_t k = 0; k < grid.size(); k++) {
        const double e = grid[k];
        Unit sub(target->GetCMAngle(), e, seg.segment, seg.term);
        sub.BelongToRowOf(*target);
        if (tspec.targetIntegration)
          sub.SetStoppingPower(tspec.stoppingConversion *
                               effect->GetStoppingPowerEq()->Evaluate(configure, e / tspec.stoppingConversion));
        auto o = std::lower_bound(old.begin(), old.end(), sub.GetCMEnergy(),
                                  [](const Unit::EDependentShare &x, double v) { return x.energy < v; });
        Unit::EDependentShare cached;
        if (o != old.end() && sub.AdoptEDependentValues(*o)) {
          sub.CalcCoulombAmplitude(this, theCNuc);  // depends on the angle; cheap
          reused++;
        } else if (!shifted && nodeCache_->Find(seg.GetSegmentKey(), sub.GetCMEnergy(), sub.GetCMAngle(), cached) &&
                   sub.AdoptEDependentValues(cached)) {
          sub.CalcCoulombAmplitude(this, theCNuc);  // another copy placed this node: its tables
          fromOtherCopies++;
        } else {
          sub.CalcLegendreP(this, configure.maxLOrder, theCNuc, qEffect);
          // A shifted segment's sub-points are moved and energised by the shift step.
          if (!shifted) {
            sub.CalcEDependentValues(this, theCNuc, configure, true, nullptr);
            sub.CalcCoulombAmplitude(this, theCNuc);
          }
        }
        if (!shifted) nodeCache_->Store(seg.GetSegmentKey(), sub.ShareEDependentValues());
        target->MakeSubPoint(sub);
        fresh.push_back(std::move(sub));
      }
      placed[q].emplace_back(static_cast<int>(target - &this->UnitAt(0)), std::move(fresh));
    };
    layOut(point, it.seg, it.shifted);
    for (int m = 0; m < mapped.size(); m++) {
      Unit *mp = mapped[m];
      auto own = mappedSegment.find(mp);  // present and replaceable: checked above
      layOut(mp, own->second.seg, own->second.shifted);
      newSubPoints += (long)grid.size();
      if (own->second.shifted) {
        std::lock_guard<std::mutex> lock(reshiftMutex);
        mappedShifted.push_back(own->second.seg);
      }
    }
    reusedSubPoints += reused;
    sharedSubPoints += fromOtherCopies;
    changed++;
    newSubPoints += (long)grid.size();
  };
#ifdef _OPENMP
  const bool serial = omp_in_parallel();
#else
  const bool serial = true;
#endif
  if (serial) {
    for (size_t q = 0; q < work.size(); q++) replace(work[q], q);
  } else {
#pragma omp parallel for schedule(dynamic, 4)
    for (int q = 0; q < (int)work.size(); q++) replace(work[q], (size_t)q);
  }
  if (changed > 0) this->PlaceSubPointGrids(placed);
  if (changed > 0)
    for (const Item &it : work)
      if (it.shifted && std::find(reshift.begin(), reshift.end(), it.seg) == reshift.end()) reshift.push_back(it.seg);
  for (const TermRef &seg : mappedShifted)
    if (std::find(reshift.begin(), reshift.end(), seg) == reshift.end()) reshift.push_back(seg);

  // The invariant point mapping rests on (SharesCalculation): a host and its
  // mapped points have the same sub-point grid, which Unit::Calculate fills
  // node by node.  Before 2026-11-05 every grid was re-placed on its own and a
  // component's grid could come out with another node count than its host's:
  // the host then wrote past the end of the component's sub-points (heap
  // corruption, natural neon fits at 4 threads) or filled nodes at other
  // energies.  Fail loudly rather than calculate on mismatched grids.
  if (changed > 0)
    for (const Item &it : work) {
      Unit *host = it.seg.GetPoint(*this, it.i);
      const UnitRange<Unit> hs = host->GetSubPoints(*this);
      const azure2::MappedPoints mapped = host->PlanMapped(this);
      for (int m = 0; m < mapped.size(); m++) {
        const UnitRange<Unit> ms = mapped[m]->GetSubPoints(*this);
        bool same = ms.size() == hs.size();
        for (int k = 0; same && k < hs.size(); k++) same = ms[k].GetOriginalEnergy() == hs[k].GetOriginalEnergy();
        if (!same)
          throw std::logic_error("WorkingCopy::ReplaceSubPointGrids: segment " + std::to_string(it.seg.GetSegmentKey()) +
                                 " point " + std::to_string(it.i) + " and a point mapped onto it have different "
                                 "sub-point grids (" + std::to_string(hs.size()) + " and " +
                                 std::to_string(ms.size()) + " nodes)");
      }
    }
  if (std::getenv("AZURE_GRID_TIMING")) {
    static std::atomic<int> calls{0};
    const int c = ++calls;
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
                << " sub-points kept their tables, " << sharedSubPoints.load()
                << " took another copy's); observed parameters " << ms(t0, t1) << " ms, grids "
                << ms(t1, std::chrono::steady_clock::now()) << " ms" << std::endl;
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
  // The working copy starts from each effect's values as the set-up leaves them (a density-group member takes its
  // group's density).
  for (int n = 1; n <= Draft().NumTargetEffects(); n++) {
    const TargetEffect &te = *Draft().GetTargetEffect(n);
    effectValues_[n - 1] = TargetEffectValues{te.GetScreeningUe(), te.GetDensity(), te.GetDopplerTemperature()};
  }
  for (SegmentIterator segment = Draft().GetSegments().begin(); segment < Draft().GetSegments().end(); segment++) {
    PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(segment->GetEntranceKey()));
    PPair *exitPair = compound->GetPairFromKey(segment->GetExitKey());  // nullptr for total capture (-1)
    const double cmConversion = Segment::CMEnergyFactor(entrancePair, exitPair);

    if (segment->IsTargetEffect()) {
      TargetEffect *targetEffect = Draft().GetTargetEffect(segment->GetTargetEffectNum());
      // the line's energies stay in the lab frame; every point converts them with its own pair's factor
      const int bohrZ1 = (entrancePair->GetPType() == 20) ? 0 : entrancePair->GetZ(1);
      const int bohrZ2 = entrancePair->GetZ(2);
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
      if (targetEffect->IsUdr()) {
        double eMin = 1.0e300, eMax = 0.0;
        for (UnitIterator point = segment->GetPoints(*this).begin(); point < segment->GetPoints(*this).end(); point++) {
          eMin = std::min(eMin, point->GetLabEnergy());
          eMax = std::max(eMax, point->GetLabEnergy());
        }
        if (eMax > 0.0 && !targetEffect->LoadUdrTable(eMin, eMax, configure.outStream, configure.cacheInputFiles)) return -1;
      }

      if (!CheckTargetWindows(targetEffect, segment->GetSegmentKey(), configure.outStream)) return -1;
      const bool thinTarget = targetEffect->IsTargetIntegration() && !targetEffect->IsThickTarget();
      const std::vector<std::pair<double, double> > windows =
          TargetWindowsFor(*this, targetEffect, compound, {&*segment, 0}, cmConversion, configure, segment->GetSegmentKey(), true);
      for (UnitIterator point = segment->GetPoints(*this).begin(); point < segment->GetPoints(*this).end(); point++) {
        // An effect restricted to energy ranges leaves points outside them
        // untouched: they never receive the effect number, so the whole
        // downstream machinery treats them as ordinary points -- except under a
        // thin-target integration, whose data carry the beam energy: there a
        // point outside the windows is evaluated at its depth-averaged energy
        // E_bar (Unit::SetTargetFallbackOnly).
        double blendWeight = targetEffect->BlendWeight(point->GetLabEnergy(), windows);
        if (blendWeight <= 0. && !thinTarget) continue;
        point->SetTargetEffectNum(segment->GetTargetEffectNum());
        point->SetTargetBlendWeight(std::max(blendWeight, 0.0));
        point->SetTargetKinematics(cmConversion, bohrZ1, bohrZ2);

        if (targetEffect->IsSubPointEffect()) {
          double forwardDepth = 0.0;
          double backwardDepth = 0.0;

          double targetThickness = 0.0;
          if (targetEffect->IsTargetIntegration()) {
            if (targetEffect->IsThickTarget()) {
              // an infinitely thick target: the point integrates from its energy down to 1 keV
              targetThickness = std::max(point->GetCMEnergy() - 0.001, 0.0);
            } else {
              // the grid spans the target at its layout density (a fitted density's largest expected value)
              targetThickness = cmConversion * targetEffect->TargetThickness(point->GetLabEnergy(), configure,
                                                                             targetEffect->GetLayoutDensity());
              if (!(targetThickness > 0.0)) {
                configure.outStream << "ERROR: target integration of segment #" << segment->GetSegmentKey()
                                    << " at E_lab = " << point->GetLabEnergy() << " MeV: the stopping cross section times the"
                                    << " density is " << targetThickness / cmConversion << " (not positive).  A density of 0"
                                    << " means an infinitely thick target; otherwise the stopping power must be positive"
                                    << " at every point." << std::endl;
                return -1;
              }
            }
            point->SetTargetLabEnergy(point->GetLabEnergy());
            point->SetTargetThickness(targetEffect->HasDensityFit()
                                          ? cmConversion * targetEffect->TargetThickness(point->GetLabEnergy(), configure)
                                          : targetThickness);
          }
          if (thinTarget) {
            point->SetTargetFallbackDepth(cmConversion * targetEffect->MeanEnergyLoss(point->GetLabEnergy(), configure,
                                                                                      targetEffect->GetDensity()));
            if (blendWeight <= 0.) {
              point->SetTargetFallbackOnly(true);
              AddTargetFallbackNodes(*this, &*point, targetEffect, cmConversion, configure, {&*segment, 0});
              continue;
            }
          }
          // How far the grid reaches above and below the point for this effect (beam spread, target, straggling).
          targetEffect->GridExtent(point->GetCMEnergy(), targetThickness, configure, cmConversion, bohrZ1, bohrZ2, forwardDepth, backwardDepth);
          // A THM (HOES) segment's Gaussian fold reaches thmConvolutionRange sigma (5; convolutionRange is 5 too
          // since dev 7d24c88).
          if (segment->IsTHM() && targetEffect->IsGaussian() && !targetEffect->IsTargetIntegration())
            forwardDepth = backwardDepth =
                TargetEffect::thmConvolutionRange * targetEffect->SigmaCM(point->GetCMEnergy(), configure, cmConversion);

          // Generate adaptive integration grid
          double startEnergy = point->GetCMEnergy() + forwardDepth;
          double endEnergy;

          /*
          THM (HOES) segments carry data below the entrance threshold, so their
          sub-point grids may extend to negative energies as long as the
          compound-system energy Ecm+SepE+ExE stays positive.
          */
          double minSubEnergy = segment->IsTHM() ? TargetEffect::minIntegrationEnergy - (entrancePair->GetSepE() + entrancePair->GetExE()) : TargetEffect::minIntegrationEnergy;

          if (targetEffect->IsTargetIntegration()) {
            // the target (or, for a thick target, everything down to 1 keV), plus the spread below it
            // (never a THM segment: refused above)
            endEnergy = point->GetCMEnergy() - backwardDepth;
            // Safety check: if backwardDepth > energy, set endEnergy to minimum
            if (endEnergy < minSubEnergy) {
              endEnergy = minSubEnergy;
            }
          } else {
            // Convolution or ConvCoefficients - use backward depth
            endEnergy = point->GetCMEnergy() - backwardDepth;
            // Safety check: if backwardDepth > energy, set endEnergy to minimum
            if (endEnergy < minSubEnergy) {
              endEnergy = minSubEnergy;
            }
          }

          if (targetEffect->IsBeamProfile()) {
            // The kernel is an absolute beam profile: sample where the beam
            // is, narrowed to the point's energy window (plus the resolution
            // tails) when it has one, not around the point's own energy.
            double low, high;
            targetEffect->BeamProfileSupport(low, high);  // lab
            low *= cmConversion; high *= cmConversion;
            double s = cmConversion * targetEffect->GetBeamTpcSigma();
            // The resolution tails: 4 s, or for a THM (HOES) point the 5 s of
            // its Gaussian fold -- a narrow peak 4-5 s from a point dominates
            // the tail there (-2 % at 81 keV from a 289 eV level, 18O(p,a) THM).
            double tails = segment->IsTHM() ? TargetEffect::thmConvolutionRange : 4.0;
            if (point->HasBinWindow()) {
              low = std::max(low, point->GetBinLowCM() - tails * s);
              high = std::min(high, point->GetBinHighCM() + tails * s);
            }
            // The floor of the Gaussian fold: a THM point may sit below the entrance
            // threshold (this was +1 keV, and the fold of a point at -45 keV came
            // out 66 times too large).
            if (low < minSubEnergy) low = minSubEnergy;
            if (high <= low) high = low + 0.001;
            startEnergy = high;
            endEnergy = low;
            // Q = E_gamma - E_cm: the entrance pair's threshold (S plus its total excitation, e1 + e2) minus the
            // energy of the final level the gamma pair feeds (its heavy excitation).  Before 2026-11-05 the
            // entrance pair's excitation was left out.
            point->SetPhotoKinematics(entrancePair->GetSepE() + entrancePair->GetExE() - exitPair->GetHeavyExE(),
                                      (entrancePair->GetM(1) + entrancePair->GetM(2)) * uconv);
          }
          if (targetEffect->IsUdr()) {
            // A delay tau maps the point's nominal flight time t0 to a true
            // flight time t0 - tau: the kernel's support in delay becomes a
            // window in true energy, extending above the point for positive
            // delays (moderation) and below it for negative ones.
            double mass = entrancePair->GetM(1) * uconv;
            std::shared_ptr<const UdrKernel> kernel = targetEffect->BuildUdrKernel(point->GetLabEnergy(), mass);
            DataSet::PointSetUp &setUp = Draft().SetUpFor(*point);
            setUp.udr = kernel;
            setUp.udrMass = mass;
            double t0 = UdrKinematics::TimeOfFlight(point->GetLabEnergy(), mass, targetEffect->GetUdrFlightPath());
            double eHighLab = UdrKinematics::EnergyFromTime(t0 - kernel->tauMax, mass, targetEffect->GetUdrFlightPath());
            double eLowLab = UdrKinematics::EnergyFromTime(t0 - kernel->tauMin, mass, targetEffect->GetUdrFlightPath());
            double cmPerLab = point->GetCMEnergy() / point->GetLabEnergy();
            if (eHighLab <= point->GetLabEnergy()) eHighLab = point->GetLabEnergy() * 1.0001;
            if (eLowLab <= 0.0 || eLowLab >= point->GetLabEnergy()) eLowLab = point->GetLabEnergy() * 0.9999;
            startEnergy = eHighLab * cmPerLab;
            endEnergy = std::max(eLowLab * cmPerLab, 0.001);
          }
          // Free-gas Doppler broadening comes first: every sub-point the kernel above reaches needs its own
          // Doppler window inside the grid, so the window is widened by the Doppler width at each end.
          if (targetEffect->IsDoppler()) {
            double lo, hi, unused;
            targetEffect->DopplerWindow(startEnergy, cmConversion, unused, hi);
            targetEffect->DopplerWindow(endEnergy, cmConversion, lo, unused);
            startEnergy = std::max(startEnergy, hi);
            endEnergy = std::max(std::min(endEnergy, lo), 0.001);
          }
          std::vector<double> energyGrid;
          int numPoints = targetEffect->NumSubPoints();
          if (configure.useAdaptiveGrid) {
            AdaptiveIntegrationGrid::GridConfig gridConfig;
            gridConfig.maxPoints = numPoints;
            gridConfig.entranceKey = segment->GetEntranceKey();
            gridConfig.inputWidthsArePhysical = (configure.paramMask & Config::TRANSFORM_PARAMETERS) && !compound->IsTransformedIn();
            gridConfig.baseEnergyStep = (startEnergy - endEnergy) / numPoints;
            gridConfig.resonanceWidthMultiplier = targetEffect->GetResonanceWidthMultiplier();
            gridConfig.pointsPerWidth = targetEffect->GetPointsPerWidth();
            AdaptiveIntegrationGrid gridGenerator(gridConfig);
            energyGrid = gridGenerator.GenerateGrid(startEnergy, endEnergy, compound);
            Draft().SetUpFor(*point).grid = MakeSubGridSpec(gridConfig, startEnergy, endEnergy, targetEffect, cmConversion);
          } else {
            double step = (startEnergy - endEnergy) / numPoints;
            for (int i = 0; i <= numPoints; i++)
              energyGrid.push_back(startEnergy - i * step);
          }

          for (size_t i = 0; i < energyGrid.size(); i++) {
            double subEnergy = energyGrid[i];
            Unit subPoint(point->GetCMAngle(), subEnergy, &*segment, 0);
            subPoint.BelongToRowOf(*point);
            if (targetEffect->IsTargetIntegration()) {
              double stoppingPower = cmConversion * targetEffect->GetStoppingPowerEq()->Evaluate(configure, subEnergy / cmConversion);
              subPoint.SetStoppingPower(stoppingPower);
            }
            this->AddSubPoint(*point, subPoint);
          }
        }
      }
    }
    {
      for (int t = 1; t < segment->NumTerms(); t++) {
        const TermRef component{&*segment, t};
        // The component's grid, kinematics and lab<->c.m. factor belong to the
        // pair the component calculates -- for a natural target another
        // compound nucleus with its own masses -- not to the base segment's.
        // The points were copied from the base segment with its c.m.
        // energies; convert them for this pair before the grid is laid out
        // (WorkingCopy::Initialize repeats this, idempotently, from the lab energy).
        PPair *compEntrance = compound->GetPair(compound->GetPairNumFromKey(component.GetEntranceKey()));
        PPair *compExit = compound->GetPairFromKey(component.GetExitKey());
        const double compConversion = Segment::CMEnergyFactor(compEntrance, compExit);
        for (int i = 0; i < component.NumPoints(); i++) segment->ConvertRowEnergyToCM(*this, t, i, compEntrance, compExit);
        if (segment->IsTargetEffect()) {
          TargetEffect *targetEffect = Draft().GetTargetEffect(component.GetTargetEffectNum());
          const int bohrZ1 = (compEntrance->GetPType() == 20) ? 0 : compEntrance->GetZ(1);
          const int bohrZ2 = compEntrance->GetZ(2);
          if (targetEffect->IsUdr()) {
            double eMin = 1.0e300, eMax = 0.0;
            for (UnitIterator point = component.GetPoints(*this).begin(); point < component.GetPoints(*this).end(); point++) {
              eMin = std::min(eMin, point->GetLabEnergy());
              eMax = std::max(eMax, point->GetLabEnergy());
            }
            if (eMax > 0.0 && !targetEffect->LoadUdrTable(eMin, eMax, configure.outStream, configure.cacheInputFiles)) return -1;
          }

          if (!CheckTargetWindows(targetEffect, segment->GetSegmentKey(), configure.outStream)) return -1;
          const bool thinTarget = targetEffect->IsTargetIntegration() && !targetEffect->IsThickTarget();
          // a component's windows come from its own entrance pair's level scheme
          const std::vector<std::pair<double, double> > windows =
              TargetWindowsFor(*this, targetEffect, compound, component, compConversion, configure, segment->GetSegmentKey(),
                               false);
          for (UnitIterator point = component.GetPoints(*this).begin(); point < component.GetPoints(*this).end(); point++) {
            double blendWeight = targetEffect->BlendWeight(point->GetLabEnergy(), windows);
            if (blendWeight <= 0. && !thinTarget) continue;  // see the base segment's loop above
            point->SetTargetEffectNum(component.GetTargetEffectNum());
            point->SetTargetBlendWeight(std::max(blendWeight, 0.0));
            point->SetTargetKinematics(compConversion, bohrZ1, bohrZ2);

            if (targetEffect->IsSubPointEffect()) {
              double forwardDepth = 0.0;
              double backwardDepth = 0.0;

              double targetThickness = 0.0;
              if (targetEffect->IsTargetIntegration()) {
                if (targetEffect->IsThickTarget()) {
                  targetThickness = std::max(point->GetCMEnergy() - 0.001, 0.0);
                } else {
                  targetThickness = compConversion * targetEffect->TargetThickness(point->GetLabEnergy(), configure,
                                                                                   targetEffect->GetLayoutDensity());
                  if (!(targetThickness > 0.0)) {
                    configure.outStream << "ERROR: target integration of a component of segment #" << segment->GetSegmentKey()
                                        << " at E_lab = " << point->GetLabEnergy() << " MeV: the stopping cross section times"
                                        << " the density is " << targetThickness / compConversion << " (not positive)." << std::endl;
                    return -1;
                  }
                }
                point->SetTargetLabEnergy(point->GetLabEnergy());
                point->SetTargetThickness(targetEffect->HasDensityFit()
                                              ? compConversion * targetEffect->TargetThickness(point->GetLabEnergy(), configure)
                                              : targetThickness);
              }
              if (thinTarget) {
                point->SetTargetFallbackDepth(compConversion * targetEffect->MeanEnergyLoss(point->GetLabEnergy(), configure,
                                                                                            targetEffect->GetDensity()));
                if (blendWeight <= 0.) {
                  point->SetTargetFallbackOnly(true);
                  AddTargetFallbackNodes(*this, &*point, targetEffect, compConversion, configure, component);
                  continue;
                }
              }
              // How far the grid reaches above and below the point for this effect (beam spread, target, straggling).
              targetEffect->GridExtent(point->GetCMEnergy(), targetThickness, configure, compConversion, bohrZ1, bohrZ2, forwardDepth, backwardDepth);
              // A THM (HOES) segment's Gaussian fold reaches thmConvolutionRange sigma (5), as above.
              if (segment->IsTHM() && targetEffect->IsGaussian() && !targetEffect->IsTargetIntegration())
                forwardDepth = backwardDepth =
                    TargetEffect::thmConvolutionRange * targetEffect->SigmaCM(point->GetCMEnergy(), configure, compConversion);

              // Generate adaptive integration grid for the leaf
              double startEnergy = point->GetCMEnergy() + forwardDepth;
              double endEnergy;

              /*
              THM (HOES) segments carry data below the entrance threshold, so their
              sub-point grids may extend to negative energies as long as the
              compound-system energy Ecm+SepE+ExE stays positive. Note this uses the
              component's own entrance pair, since an advanced (sum/ratio) component
              can have a different entrance channel than its master segment.
              */
              PPair *componentEntrancePair = compound->GetPair(compound->GetPairNumFromKey(component.GetEntranceKey()));
              double minSubEnergy = segment->IsTHM() ? TargetEffect::minIntegrationEnergy - (componentEntrancePair->GetSepE() + componentEntrancePair->GetExE()) : TargetEffect::minIntegrationEnergy;

              if (targetEffect->IsTargetIntegration()) {
                endEnergy = point->GetCMEnergy() - backwardDepth;
                // Safety check: if backwardDepth > energy, set endEnergy to minimum
                if (endEnergy < minSubEnergy) {
                  endEnergy = minSubEnergy;
                }
              } else {
                // Convolution or ConvCoefficients - use backward depth
                endEnergy = point->GetCMEnergy() - backwardDepth;
                // Safety check: if backwardDepth > energy, set endEnergy to minimum
                if (endEnergy < minSubEnergy) {
                  endEnergy = minSubEnergy;
                }
              }

              if (targetEffect->IsBeamProfile()) {
                double low, high;
                targetEffect->BeamProfileSupport(low, high);  // lab
                low *= compConversion; high *= compConversion;
                double s = compConversion * targetEffect->GetBeamTpcSigma();
                double tails = segment->IsTHM() ? TargetEffect::thmConvolutionRange : 4.0;
                if (point->HasBinWindow()) {
                  low = std::max(low, point->GetBinLowCM() - tails * s);
                  high = std::min(high, point->GetBinHighCM() + tails * s);
                }
                if (low < minSubEnergy) low = minSubEnergy;  // as above
                if (high <= low) high = low + 0.001;
                startEnergy = high;
                endEnergy = low;
                point->SetPhotoKinematics(compEntrance->GetSepE() + compEntrance->GetExE() - compExit->GetHeavyExE(),
                                          (compEntrance->GetM(1) + compEntrance->GetM(2)) * uconv);
              }
              if (targetEffect->IsUdr()) {
                double mass = compEntrance->GetM(1) * uconv;
                std::shared_ptr<const UdrKernel> kernel = targetEffect->BuildUdrKernel(point->GetLabEnergy(), mass);
                DataSet::PointSetUp &setUp = Draft().SetUpFor(*point);
                setUp.udr = kernel;
                setUp.udrMass = mass;
                double t0 = UdrKinematics::TimeOfFlight(point->GetLabEnergy(), mass, targetEffect->GetUdrFlightPath());
                double eHighLab = UdrKinematics::EnergyFromTime(t0 - kernel->tauMax, mass, targetEffect->GetUdrFlightPath());
                double eLowLab = UdrKinematics::EnergyFromTime(t0 - kernel->tauMin, mass, targetEffect->GetUdrFlightPath());
                double cmPerLab = point->GetCMEnergy() / point->GetLabEnergy();
                if (eHighLab <= point->GetLabEnergy()) eHighLab = point->GetLabEnergy() * 1.0001;
                if (eLowLab <= 0.0 || eLowLab >= point->GetLabEnergy()) eLowLab = point->GetLabEnergy() * 0.9999;
                startEnergy = eHighLab * cmPerLab;
                endEnergy = std::max(eLowLab * cmPerLab, 0.001);
              }
              // Free-gas Doppler broadening comes first: every sub-point the kernel above reaches needs its own
              // Doppler window inside the grid, so the window is widened by the Doppler width at each end.
              if (targetEffect->IsDoppler()) {
                double lo, hi, unused;
                targetEffect->DopplerWindow(startEnergy, compConversion, unused, hi);
                targetEffect->DopplerWindow(endEnergy, compConversion, lo, unused);
                startEnergy = std::max(startEnergy, hi);
                endEnergy = std::max(std::min(endEnergy, lo), 0.001);
              }
              std::vector<double> energyGrid;
              int numPoints = targetEffect->NumSubPoints();
              if (configure.useAdaptiveGrid) {
                AdaptiveIntegrationGrid::GridConfig gridConfig;
                gridConfig.maxPoints = numPoints;
                // The component's own entrance pair: its resonances, on its energy axis, are what the grid must resolve.
                gridConfig.entranceKey = component.GetEntranceKey();
                gridConfig.inputWidthsArePhysical = (configure.paramMask & Config::TRANSFORM_PARAMETERS) && !compound->IsTransformedIn();
                gridConfig.baseEnergyStep = (startEnergy - endEnergy) / numPoints;
                gridConfig.resonanceWidthMultiplier = targetEffect->GetResonanceWidthMultiplier();
                gridConfig.pointsPerWidth = targetEffect->GetPointsPerWidth();
                AdaptiveIntegrationGrid gridGenerator(gridConfig);
                energyGrid = gridGenerator.GenerateGrid(startEnergy, endEnergy, compound);
                Draft().SetUpFor(*point).grid =
                    MakeSubGridSpec(gridConfig, startEnergy, endEnergy, targetEffect, compConversion);
              } else {
                double step = (startEnergy - endEnergy) / numPoints;
                for (int i = 0; i <= numPoints; i++)
                  energyGrid.push_back(startEnergy - i * step);
              }

              for (size_t i = 0; i < energyGrid.size(); i++) {
                double subEnergy = energyGrid[i];
                Unit subPoint(point->GetCMAngle(), subEnergy, &*segment, t);
                subPoint.BelongToRowOf(*point);
                if (targetEffect->IsTargetIntegration()) {
                  double stoppingPower = compConversion * targetEffect->GetStoppingPowerEq()->Evaluate(configure, subEnergy / compConversion);
                  subPoint.SetStoppingPower(stoppingPower);
                }
                this->AddSubPoint(*point, subPoint);
              }
            }
          }
        }
      }
    }
  }
  this->CommitSubPoints();  // the sub-points given above take their places in the units
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
