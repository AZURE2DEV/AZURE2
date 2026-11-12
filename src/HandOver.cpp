// The hand-over of a parameter vector to the data: normalizations, energy shifts (with the re-energising of the
// shifted points), sample thicknesses, target densities, Doppler temperatures and screening potentials go into
// the segments, target effects and points, and the prior penalties of the target-effect parameters.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "CNuc.h"
#include "Config.h"
#include "WorkingCopy.h"

void WorkingCopy::FillThicknessFromParams(const vector_r &p) {
  if (set_->thicknessSegments_.empty() || set_->thicknessParamOffset_ < 0) return;
  for (int k = 1; k <= Data().NumThicknessParams(); k++) {
    const int idx = set_->thicknessParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the thickness block
    const Segment &seg = set_->segments_[set_->thicknessSegments_[k - 1] - 1];
    if (seg.IsThicknessGroupFollower()) continue;
    State(seg).thickness = std::max(0.0, p[idx]);
  }
  for (const Segment &seg : set_->segments_)
    if (seg.IsThicknessGroupFollower()) {
      const Segment &leader = set_->segments_[seg.GetThicknessLeaderPos()];
      State(seg).thickness = leader.GetSampleThickness(State(leader));
    }
}

double WorkingCopy::ThicknessPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= Data().NumThicknessParams(); k++) {
    const Segment &seg = set_->segments_[set_->thicknessSegments_[k - 1] - 1];
    chi2 += seg.ThicknessPenalty(State(seg));
  }
  return chi2;
}

void WorkingCopy::FillDopplerFromParams(const vector_r &p, const Config *configure) {
  if (set_->dopplerEffects_.empty() || set_->dopplerParamOffset_ < 0) return;
  for (int k = 1; k <= Data().NumDopplerParams(); k++) {
    const int idx = set_->dopplerParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the Doppler block
    const TargetEffect *te = Data().GetDopplerTargetEffect(k);
    double v = p[idx];
    // The grids and the integration windows span the layout temperature; below T0/4 the Gauss-Legendre pieces
    // (laid out for T0) would be too wide for the kernel.
    const double vMax = te->GetDopplerLayoutTemperature(), vMin = te->GetDopplerMinimum();
    if (v > vMax || !(v >= vMin)) {
      static std::atomic<bool> warned{false};
      if (!warned.exchange(true))
        (configure ? configure->outStream : std::cout)
            << "WARNING: doppler_" << Data().GetDopplerEffect(k) << " = " << v << " K is outside [" << vMin << ", "
            << vMax << "] K, the range its sub-point grids and integration pieces span; the calculation uses the"
            << " nearest end.  Give the line a wider range with 'doppler T vary sigma Tmax'. (Printed once.)"
            << std::endl;
      v = (v > vMax) ? vMax : vMin;
    }
    for (int n : set_->dopplerCopies_[k - 1]) EffectValues(n).dopplerTemperature = v;
  }
}

double WorkingCopy::DopplerPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= Data().NumDopplerParams(); k++)
    chi2 += Data().GetDopplerTargetEffect(k)->DopplerPenalty(EffectValues(Data().DopplerCopies(k).front()).dopplerTemperature);
  return chi2;
}

void WorkingCopy::FillDensitiesFromParams(const vector_r &p, const Config &configure) {
  if (set_->densityEffects_.empty() || set_->densityParamOffset_ < 0) return;
  std::vector<char> changed(Data().NumTargetEffects(), 0);
  bool any = false;
  for (int k = 1; k <= Data().NumDensityParams(); k++) {
    const int idx = set_->densityParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the density block
    const TargetEffect *lead = Data().GetDensityTargetEffect(k);
    double v = p[idx];
    // The sub-point grids span the layout density: beyond it the depth integral would run off the grid.
    const double vMax = lead->GetLayoutDensity(), vMin = 0.01 * lead->GetDensityNominal();
    if (v > vMax || !(v >= vMin)) {
      static std::atomic<bool> warned{false};
      if (!warned.exchange(true))
        configure.outStream << "WARNING: " << lead->GetTargetName() << ".areal_density = " << v << " is outside ["
                            << vMin << ", " << vMax << "], the range its sub-point grids span; the calculation uses the"
                            << " nearest end.  Give the target's areal_density a larger grid_max."
                            << " (Printed once.)" << std::endl;
      v = (v > vMax) ? vMax : vMin;
    }
    for (int n = 1; n <= Data().NumTargetEffects(); n++) {
      if (set_->densityParamOfEffect_[n - 1] != k) continue;
      TargetEffectValues &values = EffectValues(n);
      if (values.density == v) continue;
      values.density = v;
      changed[n - 1] = 1;
      any = true;
    }
  }
  if (!any) return;
  auto update = [&](const TermRef &seg) {
    if (!seg.segment->IsTargetEffect()) return;
    const int n = seg.GetTargetEffectNum();
    if (n < 1 || n > (int)changed.size() || !changed[n - 1]) return;
    const TargetEffect *te = Data().GetTargetEffect(n);
    const double density = EffectValues(n).density;
    for (UnitIterator point = seg.GetPoints(*this).begin(); point < seg.GetPoints(*this).end(); point++) {
      if (point->GetTargetEffectNum() != n || !(point->GetTargetThickness() > 0.0)) continue;
      const double scale = EffectValues(n).stoppingScale;
      point->SetTargetThickness(point->GetCMFactor() * te->TargetThickness(point->GetTargetLabEnergy(), configure, density, scale));
    }
  };
  for (const Segment &seg : set_->segments_) update({&seg, 0});
  for (const LeafRef &leaf : set_->leafRefs_) update({&set_->segments_[leaf.first], leaf.second});
}

double WorkingCopy::DensityPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= Data().NumDensityParams(); k++)
    chi2 += Data().GetDensityTargetEffect(k)->DensityPenalty(EffectValues(Data().GetDensityLeadCopy(k)).density);
  return chi2;
}

void WorkingCopy::ApplyScreening() {
  if (set_->screeningEffects_.empty()) return;
  std::vector<int> paramOf(Data().NumTargetEffects() + 1, 0);
  for (int k = 1; k <= Data().NumScreeningParams(); k++)
    for (int n : set_->screeningCopies_[k - 1]) paramOf[n] = k;
  auto apply = [&](const TermRef &seg) {
    const int te = seg.GetTargetEffectNum();
    if (te <= 0 || te >= (int)paramOf.size() || paramOf[te] == 0) return;
    const double ueMeV = EffectValues(te).screeningUe * 1.0e-6;
    for (int i = 1; i <= seg.NumPoints(); i++) {
      Unit *point = seg.GetPoint(*this, i);
      point->SetScreening(ueMeV, paramOf[te]);
      for (int j = 1; j <= point->NumSubPoints(); j++) point->GetSubPoint(*this, j)->SetScreening(ueMeV, paramOf[te]);
    }
  };
  for (const Segment &seg : set_->segments_) apply({&seg, 0});
  for (const LeafRef &leaf : set_->leafRefs_) apply({&set_->segments_[leaf.first], leaf.second});
}

void WorkingCopy::FillScreeningFromParams(const vector_r &p) {
  if (set_->screeningEffects_.empty() || set_->screeningParamOffset_ < 0) return;
  for (int k = 1; k <= Data().NumScreeningParams(); k++) {
    const int idx = set_->screeningParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the screening block
    for (int n : set_->screeningCopies_[k - 1]) EffectValues(n).screeningUe = p[idx];
  }
  this->ApplyScreening();
}

double WorkingCopy::ScreeningPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= Data().NumScreeningParams(); k++)
    chi2 += Data().GetScreeningTargetEffect(k)->ScreeningPenalty(EffectValues(Data().ScreeningCopies(k).front()).screeningUe);
  return chi2;
}

/*!
 * Fills the Normalizations from the Minuit parameter array.
 */

void WorkingCopy::FillNormsFromParams(const vector_r &p) {
  int i = Data().GetNormParamOffset();
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    if (segment->IsVaryNorm() && !segment->IsProfiledNorm()) {
      State(*segment).norm = p[i];
      i++;
    }
    // With component segments, no need to skip - total capture is handled within the segment
  }
  // Norm-group followers take their leader's norm (a follower's own parameter is fixed).
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++)
    if (segment->IsNormGroupFollower()) State(*segment).norm = State(Data().GetSegments()[segment->GetNormLeaderPos()]).norm;
}

/*!
 * Fills the Energy Shifts from the Minuit parameter array.
 */

void WorkingCopy::FillEnergyShiftsFromParams(const vector_r &p, azure2::Results &results, CNuc *theCNuc,
                                       const Config *configure) {
  this->FillEnergyShiftsFromParamsImpl(p, theCNuc, configure);
  this->BeginEvaluation(results);
}

void WorkingCopy::FillEnergyShiftsFromParamsImpl(const vector_r &p, CNuc *theCNuc, const Config *configure) {
  // The THM coherent backgrounds, the block after the shifts (none without
  // cbackground=): every caller fills the shifts right after the norms, with
  // the compound the model will run on.
  FillThmCoherentFromParams(p, theCNuc);
  // Screening first: sub-points rebuilt below (ReplaceSubPointGrids) inherit their parent's U_e.
  this->FillScreeningFromParams(p);
  // Fitted target densities: the thickness of each point of their lines (the grids stay put).
  if (configure) this->FillApparatusFromParams(p, *configure);
  if (configure) this->FillDensitiesFromParams(p, *configure);
  // Fitted Doppler temperatures: only the kernel values move (the grids and pieces stay put).
  this->FillDopplerFromParams(p, configure);
  // Fitted sample thicknesses: the n of the transmission / yield transforms.
  this->FillThicknessFromParams(p);
  int i = Data().GetEnergyShiftParamOffset();
  int j = Data().GetEnergyShiftSqrtParamOffset();

  // Re-energising a point (Segment::UpdatePointEnergyWithShift: Coulomb
  // functions of the point and its sub-points) touches only that point, so
  // the shift values are set here in the original order and the points of
  // the segments that moved are re-energised afterwards in parallel, each
  // segment once, in its final state.  A line search moves every shift at once; serially that cost
  // ~195 s per evaluation on the 11B+a mixed grid against 0.5 s for the
  // points.  Inside an enclosing parallel region (Minuit's numerical
  // gradient) the loop stays serial, as before.
  // The terms to re-energise: a segment's own and, with it, its leaves, which take the segment's shift.
  std::vector<TermRef> toUpdate;
  auto queueTerm = [&toUpdate](const TermRef &t) {
    if (std::find(toUpdate.begin(), toUpdate.end(), t) == toUpdate.end()) toUpdate.push_back(t);
  };
  auto queueUpdate = [&queueTerm](const Segment *s) {
    for (int t = 0; t < s->NumTerms(); t++) queueTerm({s, t});
  };
  // The terms of the segments whose points carry an energy shift, moved or not.
  std::vector<TermRef> shiftActiveSegments;
  for (SegmentIterator segment = Data().GetSegments().begin(); segment < Data().GetSegments().end(); segment++) {
    // Both terms are present for every segment: the constant shift at p[i]
    // and the sqrt(E) coefficient at p[j].  Either moving the points means
    // the segment must be re-shifted (the mapping E -> E + a + b*sqrt(E) is
    // applied from the original energies, so one update covers both).
    // A shift-group follower reads its group leader's constant shift.
    const int lead = segment->GetShiftLeaderPos();
    const double shiftConst = (lead >= 0) ? p[Data().GetEnergyShiftParamOffset() + lead] : p[i];
    const double shiftSqrt = (j >= 0 && j < (int)p.size()) ? p[j] : 0.0;
    // A shift function's parameters sit in the shift-function block, at the
    // owning segment's (this one or its group leader) start.
    SegmentState &state = State(*segment);
    std::vector<double> fnValues;
    if (segment->HasShiftFn()) {
      const int start = Data().GetSegment(segment->GetShiftFnSourcePos() + 1)->GetShiftFnParamStart();
      for (int q = 0; q < segment->NumShiftFnParams(); q++) {
        const int idx = Data().GetEnergyShiftFnParamOffset() + start + q;
        fnValues.push_back((start >= 0 && idx < (int)p.size()) ? p[idx] : state.shiftFnValues[q]);
      }
    }
    // A segment whose points carry a shift of an earlier evaluation (last shift not NaN) stays active, so a fixed
    // shift set back to 0 moves them back (upgrade_history: such a copy gives a fresh copy's bits).
    const bool shiftActive = segment->IsVaryEnergyShift() || shiftConst != 0.0 || lead >= 0 || segment->HasShiftFn() ||
                             segment->IsVaryEnergyShiftSqrt() || shiftSqrt != 0.0 || !std::isnan(state.lastEnergyShift);
    if (shiftActive) {
      for (int t = 0; t < segment->NumTerms(); t++) shiftActiveSegments.push_back({&*segment, t});
      if (state.lastEnergyShift == shiftConst && state.lastEnergyShiftSqrt == shiftSqrt &&
          (!segment->HasShiftFn() || state.lastShiftFnValues == fnValues)) {
        i++;
        j++;
        continue;
      }

      state.energyShift = shiftConst;
      state.lastEnergyShift = shiftConst;
      state.energyShiftSqrt = shiftSqrt;
      state.lastEnergyShiftSqrt = shiftSqrt;
      for (int q = 0; q < (int)fnValues.size(); q++) segment->SetShiftFnValue(state, q, fnValues[q]);
      if (segment->HasShiftFn()) state.lastShiftFnValues = fnValues;
      queueUpdate(&*segment);

    }

    i++;
    j++;
  }

  // The sub-point grids, for the levels and shifts just set.  Sub-points of
  // shifted segments are moved and energised by the step below.
  // The plan tasks whose units move (re-placed grids, shifted points): their slots are looked up again at the end.
  std::vector<int> movedTasks;
  if (theCNuc && configure) {
    std::vector<TermRef> reshift;
    this->ReplaceSubPointGrids(theCNuc, *configure, shiftActiveSegments, reshift, movedTasks);
    // Slots placed again: the row sets of their readers at the new energies (RowSet.h).
    if (!movedTasks.empty()) this->BuildRowSets(theCNuc, *configure);
    for (const TermRef &t : reshift) queueTerm(t);
  }

  const bool reenergise = !toUpdate.empty();
  if (reenergise) {
    this->ReleaseOverflowRows();
    for (const TermRef &s : toUpdate)
      for (int i = 0; i < s.NumPoints(); i++)
        if (const Unit *pt = s.GetPoint(*this, i + 1))
          if (pt->PlanTask() >= 0) movedTasks.push_back(pt->PlanTask());
    if (theCNuc && configure)
      for (const TermRef &t : toUpdate) t.segment->BeginShiftRefill(*this, t.term, theCNuc, *configure);
    // One work item per point (with its sub-points), largest first, so a
    // segment with thousands of sub-points per point does not serialize.
    struct Item { double cost; TermRef seg; int i; };
    std::vector<Item> work;
    for (const TermRef &s : toUpdate)
      for (int i = 0; i < s.NumPoints(); i++) {
        Unit *pt = s.GetPoint(*this, i + 1);
        work.push_back({1.0 + (pt ? pt->NumSubPoints() : 0), s, i});
      }
    std::stable_sort(work.begin(), work.end(), [](const Item &x, const Item &y) { return x.cost > y.cost; });
#ifdef _OPENMP
    const bool serial = omp_in_parallel();
#else
    const bool serial = true;
#endif
    // No parallel construct at all when already inside a parallel region: an
    // inactive nested team renumbers its thread 0, and per-thread state
    // indexed by thread number would then be shared between the outer threads.
    if (serial) {
      for (size_t q = 0; q < work.size(); q++)
        work[q].seg.segment->UpdatePointEnergyWithShift(*this, work[q].seg.term, work[q].i, theCNuc, configure);
    } else {
#pragma omp parallel for schedule(dynamic, 1)
      for (int q = 0; q < (int)work.size(); q++)
        work[q].seg.segment->UpdatePointEnergyWithShift(*this, work[q].seg.term, work[q].i, theCNuc, configure);
    }
    if (theCNuc && configure)
      for (const TermRef &t : toUpdate) t.segment->CommitShiftTable(*this, t.term, theCNuc, *configure);
  }
  this->PlaceOverflowNodes(theCNuc, movedTasks, reenergise);
  if (azure2::NodeCheckRequested() && theCNuc) {
    // Every unit of every plan task on a node of the tables, as the hand-over leaves them (the point pass checks its
    // own tasks again; an evaluation inside an enclosing parallel region has no point pass).
    std::vector<int> all(plan_.NumTasks());
    for (int t = 0; t < plan_.NumTasks(); t++) all[t] = t;
    long differ;
#pragma omp critical(azure_node_check)
    differ = this->CheckNodeTable(theCNuc, all);
    if (differ) throw std::logic_error("AZURE_NODE_CHECK: " + std::to_string(differ) + " differ (listed above)");
  }
}
