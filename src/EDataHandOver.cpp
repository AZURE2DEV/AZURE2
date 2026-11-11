// The hand-over of a parameter vector to the data: normalizations, energy shifts (with the re-energising of the
// shifted points), sample thicknesses, target densities, Doppler temperatures and screening potentials go into
// the segments, target effects and points, and the prior penalties of the target-effect parameters.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "CNuc.h"
#include "Config.h"
#include "EData.h"

void EData::FillThicknessFromParams(const vector_r &p) {
  if (thicknessSegments_.empty() || thicknessParamOffset_ < 0) return;
  for (int k = 1; k <= this->NumThicknessParams(); k++) {
    const int idx = thicknessParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the thickness block
    ESegment &seg = segments_[thicknessSegments_[k - 1] - 1];
    if (seg.IsThicknessGroupFollower()) continue;
    seg.SetSampleThickness(std::max(0.0, p[idx]));
  }
  for (ESegment &seg : segments_)
    if (seg.IsThicknessGroupFollower())
      seg.SetSampleThickness(segments_[seg.GetThicknessLeaderPos()].GetSampleThickness());
}

double EData::ThicknessPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= this->NumThicknessParams(); k++) chi2 += segments_[thicknessSegments_[k - 1] - 1].ThicknessPenalty();
  return chi2;
}

void EData::FillDopplerFromParams(const vector_r &p, const Config *configure) {
  if (dopplerEffects_.empty() || dopplerParamOffset_ < 0) return;
  for (int k = 1; k <= this->NumDopplerParams(); k++) {
    const int idx = dopplerParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the Doppler block
    TargetEffect *te = this->GetDopplerTargetEffect(k);
    double v = p[idx];
    // The grids and the integration windows span the layout temperature; below T0/4 the Gauss-Legendre pieces
    // (laid out for T0) would be too wide for the kernel.
    const double vMax = te->GetDopplerLayoutTemperature(), vMin = te->GetDopplerMinimum();
    if (v > vMax || !(v >= vMin)) {
      static std::atomic<bool> warned{false};
      if (!warned.exchange(true))
        (configure ? configure->outStream : std::cout)
            << "WARNING: doppler_" << this->GetDopplerEffect(k) << " = " << v << " K is outside [" << vMin << ", "
            << vMax << "] K, the range its sub-point grids and integration pieces span; the calculation uses the"
            << " nearest end.  Give the line a wider range with 'doppler T vary sigma Tmax'. (Printed once.)"
            << std::endl;
      v = (v > vMax) ? vMax : vMin;
    }
    for (int n : dopplerCopies_[k - 1]) this->GetTargetEffect(n)->SetDopplerTemperature(v);
  }
}

double EData::DopplerPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= this->NumDopplerParams(); k++)
    chi2 += this->GetDopplerTargetEffect(k)->DopplerPenalty();
  return chi2;
}

void EData::FillDensitiesFromParams(const vector_r &p, const Config &configure) {
  if (densityEffects_.empty() || densityParamOffset_ < 0) return;
  std::vector<char> changed(this->NumTargetEffects(), 0);
  bool any = false;
  for (int k = 1; k <= this->NumDensityParams(); k++) {
    const int idx = densityParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the density block
    const TargetEffect *lead = this->GetDensityTargetEffect(k);
    double v = p[idx];
    // The sub-point grids span the layout density: beyond it the depth integral would run off the grid.
    const double vMax = lead->GetLayoutDensity(), vMin = 0.01 * lead->GetDensityNominal();
    if (v > vMax || !(v >= vMin)) {
      static std::atomic<bool> warned{false};
      if (!warned.exchange(true))
        configure.outStream << "WARNING: density_" << this->GetDensityEffect(k) << " = " << v << " is outside ["
                            << vMin << ", " << vMax << "], the range its sub-point grids span; the calculation uses the"
                            << " nearest end.  Give the line a wider range with 'densityfit vary sigma max'."
                            << " (Printed once.)" << std::endl;
      v = (v > vMax) ? vMax : vMin;
    }
    for (int n = 1; n <= this->NumTargetEffects(); n++) {
      if (densityParamOfEffect_[n - 1] != k) continue;
      TargetEffect *te = this->GetTargetEffect(n);
      if (te->GetDensity() == v) continue;
      te->SetDensity(v);
      changed[n - 1] = 1;
      any = true;
    }
  }
  if (!any) return;
  auto update = [&](const TermRef &seg) {
    if (!seg.segment->IsTargetEffect()) return;
    const int n = seg.GetTargetEffectNum();
    if (n < 1 || n > (int)changed.size() || !changed[n - 1]) return;
    TargetEffect *te = this->GetTargetEffect(n);
    for (EPointIterator point = seg.GetPoints(*this).begin(); point < seg.GetPoints(*this).end(); point++) {
      if (point->GetTargetEffectNum() != n || !(point->GetTargetThickness() > 0.0)) continue;
      point->SetTargetThickness(point->GetCMFactor() * te->TargetThickness(point->GetTargetLabEnergy(), configure));
      // E_bar of the sigma(E_bar) fallback moves with the density too (EPoint::SetTargetFallbackOnly)
      if (!te->IsThickTarget())
        point->SetTargetFallbackDepth(point->GetCMFactor() *
                                      te->MeanEnergyLoss(point->GetTargetLabEnergy(), configure, te->GetDensity()));
    }
  };
  for (ESegment &seg : segments_) update({&seg, 0});
  for (const LeafRef &leaf : leafRefs_) update({&segments_[leaf.first], leaf.second});
}

double EData::DensityPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= this->NumDensityParams(); k++)
    chi2 += this->GetDensityTargetEffect(k)->DensityPenalty();
  return chi2;
}

void EData::ApplyScreening() {
  if (screeningEffects_.empty()) return;
  std::vector<int> paramOf(this->NumTargetEffects() + 1, 0);
  for (int k = 1; k <= this->NumScreeningParams(); k++)
    for (int n : screeningCopies_[k - 1]) paramOf[n] = k;
  auto apply = [&](const TermRef &seg) {
    const int te = seg.GetTargetEffectNum();
    if (te <= 0 || te >= (int)paramOf.size() || paramOf[te] == 0) return;
    const double ueMeV = this->GetTargetEffect(te)->GetScreeningUe() * 1.0e-6;
    for (int i = 1; i <= seg.NumPoints(); i++) {
      EPoint *point = seg.GetPoint(*this, i);
      point->SetScreening(ueMeV, paramOf[te]);
      for (int j = 1; j <= point->NumSubPoints(); j++) point->GetSubPoint(*this, j)->SetScreening(ueMeV, paramOf[te]);
    }
  };
  for (ESegment &seg : segments_) apply({&seg, 0});
  for (const LeafRef &leaf : leafRefs_) apply({&segments_[leaf.first], leaf.second});
}

void EData::FillScreeningFromParams(const vector_r &p) {
  if (screeningEffects_.empty() || screeningParamOffset_ < 0) return;
  for (int k = 1; k <= this->NumScreeningParams(); k++) {
    const int idx = screeningParamOffset_ + k - 1;
    if (idx < 0 || idx >= (int)p.size()) return;  // a vector without the screening block
    for (int n : screeningCopies_[k - 1]) this->GetTargetEffect(n)->SetScreeningUe(p[idx]);
  }
  this->ApplyScreening();
}

double EData::ScreeningPenalty() const {
  double chi2 = 0.0;
  for (int k = 1; k <= this->NumScreeningParams(); k++)
    chi2 += this->GetScreeningTargetEffect(k)->ScreeningPenalty();
  return chi2;
}

/*!
 * Fills the Normalizations from the Minuit parameter array.
 */

void EData::FillNormsFromParams(const vector_r &p) {
  int i = GetNormParamOffset();
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++) {
    if (segment->IsVaryNorm() && !segment->IsProfiledNorm()) {
      segment->SetNorm(p[i]);
      i++;
    }
    // With component segments, no need to skip - total capture is handled within the segment
  }
  // Norm-group followers take their leader's norm (a follower's own parameter is fixed).
  for (ESegmentIterator segment = GetSegments().begin(); segment < GetSegments().end(); segment++)
    if (segment->IsNormGroupFollower()) segment->SetNorm(GetSegments()[segment->GetNormLeaderPos()].GetNorm());
}

/*!
 * Fills the Energy Shifts from the Minuit parameter array.
 */

void EData::FillEnergyShiftsFromParams(const vector_r &p, azure2::Results &results, CNuc *theCNuc,
                                       const Config *configure) {
  this->FillEnergyShiftsFromParamsImpl(p, theCNuc, configure);
  this->BeginEvaluation(results);
}

void EData::FillEnergyShiftsFromParamsImpl(const vector_r &p, CNuc *theCNuc, const Config *configure) {
  // The THM coherent backgrounds, the block after the shifts (none without
  // cbackground=): every caller fills the shifts right after the norms, with
  // the compound the model will run on.
  FillThmCoherentFromParams(p, theCNuc);
  // Screening first: sub-points rebuilt below (ReplaceSubPointGrids) inherit their parent's U_e.
  this->FillScreeningFromParams(p);
  // Fitted target densities: the thickness of each point of their lines (the grids stay put).
  if (configure) this->FillDensitiesFromParams(p, *configure);
  // Fitted Doppler temperatures: only the kernel values move (the grids and pieces stay put).
  this->FillDopplerFromParams(p, configure);
  // Fitted sample thicknesses: the n of the transmission / yield transforms.
  this->FillThicknessFromParams(p);
  int i = GetEnergyShiftParamOffset();
  int j = GetEnergyShiftSqrtParamOffset();
  bool anyEnergyShifted = false;

  // Re-energising a point (ESegment::UpdatePointEnergyWithShift: Coulomb
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
  auto queueUpdate = [&queueTerm](ESegment *s) {
    for (int t = 0; t < s->NumTerms(); t++) queueTerm({s, t});
  };
  // The terms of the segments whose points carry an energy shift, moved or not.
  std::vector<TermRef> shiftActiveSegments;
  for (ESegmentIterator segment = this->GetSegments().begin(); segment < this->GetSegments().end(); segment++) {
    // Both terms are present for every segment: the constant shift at p[i]
    // and the sqrt(E) coefficient at p[j].  Either moving the points means
    // the segment must be re-shifted (the mapping E -> E + a + b*sqrt(E) is
    // applied from the original energies, so one update covers both).
    // A shift-group follower reads its group leader's constant shift.
    const int lead = segment->GetShiftLeaderPos();
    const double shiftConst = (lead >= 0) ? p[GetEnergyShiftParamOffset() + lead] : p[i];
    const double shiftSqrt = (j >= 0 && j < (int)p.size()) ? p[j] : 0.0;
    // A shift function's parameters sit in the shift-function block, at the
    // owning segment's (this one or its group leader) start.
    std::vector<double> fnValues;
    if (segment->HasShiftFn()) {
      const int start = this->GetSegment(segment->GetShiftFnSourcePos() + 1)->GetShiftFnParamStart();
      for (int q = 0; q < segment->NumShiftFnParams(); q++) {
        const int idx = GetEnergyShiftFnParamOffset() + start + q;
        fnValues.push_back((start >= 0 && idx < (int)p.size()) ? p[idx] : segment->GetShiftFnValues()[q]);
      }
    }
    // A segment whose points carry a shift of an earlier evaluation (last shift not NaN) stays active, so a fixed
    // shift set back to 0 moves them back (upgrade_history: such a copy gives a fresh copy's bits).
    const bool shiftActive = segment->IsVaryEnergyShift() || shiftConst != 0.0 || lead >= 0 || segment->HasShiftFn() ||
                             segment->IsVaryEnergyShiftSqrt() || shiftSqrt != 0.0 || !std::isnan(segment->GetLastEnergyShift());
    if (shiftActive) {
      for (int t = 0; t < segment->NumTerms(); t++) shiftActiveSegments.push_back({&*segment, t});
      if (segment->GetLastEnergyShift() == shiftConst && segment->GetLastEnergyShiftSqrt() == shiftSqrt &&
          (!segment->HasShiftFn() || segment->GetLastShiftFnValues() == fnValues)) {
        i++;
        j++;
        continue;
      }

      segment->SetEnergyShift(shiftConst);
      segment->SetLastEnergyShift(shiftConst);
      segment->SetEnergyShiftSqrt(shiftSqrt);
      segment->SetLastEnergyShiftSqrt(shiftSqrt);
      for (int q = 0; q < (int)fnValues.size(); q++) segment->SetShiftFnValue(q, fnValues[q]);
      if (segment->HasShiftFn()) segment->SetLastShiftFnValues(fnValues);
      queueUpdate(&*segment);
      anyEnergyShifted = true;

    }

    i++;
    j++;
  }

  // The sub-point grids, for the levels and shifts just set.  Sub-points of
  // shifted segments are moved and energised by the step below.
  if (theCNuc && configure) {
    std::vector<TermRef> reshift;
    if (this->ReplaceSubPointGrids(theCNuc, *configure, shiftActiveSegments, reshift) > 0) anyEnergyShifted = true;
    for (const TermRef &t : reshift) queueTerm(t);
  }

  if (!toUpdate.empty()) {
    if (theCNuc && configure)
      for (const TermRef &t : toUpdate) t.segment->BeginShiftRefill(t.term, theCNuc, *configure);
    // One work item per point (with its sub-points), largest first, so a
    // segment with thousands of sub-points per point does not serialize.
    // A point mapped onto another (EData::MapData) copies its values from
    // that host, so it is re-energised with it, on the same thread, after
    // it (ESegment::UpdatePointGroupEnergyWithShift); it gets no item of its
    // own.  A mapped point costs little (a copy, its Coulomb amplitudes).
    struct Item { double cost; TermRef seg; int i; };
    std::vector<Item> work;
    for (const TermRef &s : toUpdate)
      for (int i = 0; i < s.NumPoints(); i++) {
        if (s.segment->ShiftHost(*this, s.term, i)) continue;
        EPoint *pt = s.GetPoint(*this, i + 1);
        const double nsp = pt ? pt->NumSubPoints() : 0;
        const double nMapped = pt ? pt->PlanMapped(this).size() : 0;
        work.push_back({(1.0 + nsp) * (1.0 + 0.1 * nMapped), s, i});
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
        work[q].seg.segment->UpdatePointGroupEnergyWithShift(*this, work[q].seg.term, work[q].i, theCNuc, configure);
    } else {
#pragma omp parallel for schedule(dynamic, 1)
      for (int q = 0; q < (int)work.size(); q++)
        work[q].seg.segment->UpdatePointGroupEnergyWithShift(*this, work[q].seg.term, work[q].i, theCNuc, configure);
    }
    if (theCNuc && configure)
      for (const TermRef &t : toUpdate) t.segment->CommitShiftTable(t.term, theCNuc, *configure);
  }

  // The energy-based point mapping (MapData) follows the shifted energies.
  if (anyEnergyShifted) {
    this->MapData(theCNuc);
  }
}
