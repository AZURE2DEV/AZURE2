// AZURECalc, the residual rows: the residual vector of the objective and its Jacobian (ResidualJacobian: the
// adjoint rows of AZUREGrad, the energy-shift columns, the whitening of a data covariance, the model data, the
// prior covariance, Park's wall and the abundance sum as rows), the rows alone (ResidualsOnly), and the
// energy-shift columns by central differences (FillShiftColumns, which a pyazr session's Jacobian calls too).  One
// unit because the three build the same rows in the same order: the least-squares fitters, G2 and the
// Jacobian-path gradient read them.
#include "AZURECalc.h"
#include "RunDiagnostics.h"
#include "azure2/support/Environment.h"
#include "ModelPass.h"
#include "Objective.h"
#include "ShiftedFold.h"
#include "ShiftParamTarget.h"
#include "InferenceTypes.h"
#include "Profile.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "Config.h"
#include "CompoundSystem.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "Unit.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "AZUREGrad.h"

ShiftParamTarget ShiftTargets(WorkingCopy *ld, const ParamDesc &d) {
  ShiftParamTarget t{d.kind, d.fnTerm, {d.segment}};
  const int pos = d.segment - 1;
  for (int s = 1; s <= ld->Data().NumSegments(); s++) {
    if (s == d.segment) continue;
    const Segment *seg = ld->Data().GetSegment(s);
    if ((d.kind == ParamKind::EnergyShift && seg->GetShiftLeaderPos() == pos) ||
        (d.kind == ParamKind::EnergyShiftFn && seg->GetShiftFnSourcePos() == pos))
      t.segs.push_back(s);
  }
  return t;
}

ShiftColumnCount FillShiftColumns(const CompoundSystem *lc, LevelState &state, const LevelInputs &inputs, WorkingCopy *ld, azure2::Results &results, const Config &config,
                                  const ParamIndexMap &pmap,
                                  const std::vector<bool> &fixed, const vector_r &full,
                                  const std::vector<char> *skipPacked, size_t nRes, int nCols, vector_r &jac) {
  ShiftColumnCount count;

  // A shift of a self-contained segment (no components, not a total-capture
  // sum) moves only that segment's rows, so
  // its column is differenced on those rows alone; every other row of a full
  // re-evaluation would difference to exactly zero.  With 55 free shifts on
  // 11B+a the full path took over an hour per gradient.  Same criterion as the
  // segment-local shift gradient in Gradient().
  bool rowLayoutUnknown = false;
  std::vector<int> rowStart(ld->Data().NumSegments() + 2, 0);
  std::vector<const Segment *> segments;
  {
    int row = 0;
    for (int s = 1; s <= ld->Data().NumSegments(); s++) {
      rowStart[s] = row;
      const Segment *seg = ld->Data().GetSegment(s);
      if (!seg) continue;
      segments.push_back(seg);
      for (int pid = 0; pid < seg->NumPoints(); pid++) {
        Unit *pt = seg->GetPoint(*ld, pid + 1);
        if (!pt) continue;
        row++;
      }
    }
    rowStart[ld->Data().NumSegments() + 1] = row;
    if ((size_t)row != nRes) rowLayoutUnknown = true;  // row layout not understood: use the full path
  }
  // Residuals of one segment's rows, in row order, exactly as
  // ComputeResidualJacobian forms them.
  auto segmentResiduals = [&](const Segment *seg, vector_r &out, azure2::LeafValues leaves) {
    const double norm = ld->State(*seg).norm;
    out.clear();
    for (int pid = 0; pid < seg->NumPoints(); pid++) {
      if (!seg->GetPoint(*ld, pid + 1)) continue;
      azure2::EvaluateModelPoint(*ld, *seg, pid, lc, state, config, results, leaves);
      out.push_back(azure2::PointResidual(*ld, *seg, pid, results, norm).Pull());
    }
  };
  // The residuals of the segments `segs` (each into its own vector): folded again on the fixed nodes (ShiftedFold.h), or
  // one point pass over them (their energy tasks, then their folds), and each point from the records it wrote; after
  // a failed point segment by segment, point by point.
  auto groupResiduals = [&](const std::vector<int> &segs, std::vector<vector_r> &out) {
    std::vector<const Segment *> group;
    for (int s : segs) group.push_back(ld->Data().GetSegment(s));
    const bool pre = azure2::RefoldShiftedSegments(*ld, group, config, results) ||
                     ld->CalculatePointsInParallel(group, lc, state, config, results) == 0;
    for (size_t m = 0; m < group.size(); m++)
      segmentResiduals(group[m], out[m], pre ? azure2::LeafValues::Ready : azure2::LeafValues::Calculate);
  };
  // Every row, forward only: one parallel pass over the points (as a
  // chi-squared evaluation), then the residuals in row order.  If a point
  // failed in the pass, the combination recalculates point by point.
  auto allResiduals = [&](vector_r &out) {
    const bool pre = ld->CalculatePointsInParallel(segments, lc, state, config, results) == 0;
    out.clear();
    for (const Segment *seg : segments) {
      const double norm = ld->State(*seg).norm;
      for (int pid = 0; pid < seg->NumPoints(); pid++) {
        Unit *pt = seg->GetPoint(*ld, pid + 1);
        if (!pt) continue;
        azure2::EvaluateModelPoint(*ld, *seg, pid, lc, state, config, results, pre ? azure2::LeafValues::Ready : azure2::LeafValues::Calculate);
        out.push_back(azure2::PointResidual(*ld, *seg, pid, results, norm).Pull());
      }
    }
  };

  // Local shift columns: a parameter whose segments (its own, plus its shift
  // group's followers for a group's constant shift or a shift-function
  // parameter) are all self-contained is differenced on their rows alone.
  struct LocalShift { int f, packed; ShiftParamTarget target; };
  std::vector<LocalShift> localCols;
  std::vector<int> fullShifts;
  for (int f = 0; f < pmap.NumFull() && f < (int)full.size() && f < (int)fixed.size(); f++) {
    if (fixed[f]) continue;
    const ParamKind kind = pmap.Desc(f).kind;
    if (kind != ParamKind::EnergyShift && kind != ParamKind::EnergyShiftSqrt && kind != ParamKind::EnergyShiftFn &&
        kind != ParamKind::Density && kind != ParamKind::Apparatus) continue;
    const int packed = pmap.FullToPacked(f);
    if (packed < 0 || packed >= nCols) continue;
    // Columns the caller holds fixed: not differenced.
    if (skipPacked && packed < (int)skipPacked->size() && (*skipPacked)[packed]) continue;
    // A fitted target density, or an apparatus parameter (a stopping scale moves the depth window of every point of
    // its lines, a kernel's width every fold of its effects): full re-evaluation.
    if (kind == ParamKind::Density || kind == ParamKind::Apparatus) {
      fullShifts.push_back(f);
      continue;
    }
    const int sNum = pmap.Desc(f).segment;
    if (!rowLayoutUnknown && sNum >= 1 && sNum <= ld->Data().NumSegments()) {
      const ShiftParamTarget target = ShiftTargets(ld, pmap.Desc(f));
      if (target.Local(ld)) {
        localCols.push_back({f, packed, target});
        continue;
      }
    }
    fullShifts.push_back(f);
  }
  // The points of `segs` at the shift set now, with their slots on the nodes they are at, as the hand-over leaves them
  // (WorkingCopy::MoveSegmentsToShift).
  auto moveSegments = [&](const std::vector<int> &segs) {
    std::vector<const Segment *> group;
    for (int s : segs) group.push_back(ld->Data().GetSegment(s));
    ld->MoveSegmentsToShift(group, lc, config);
  };
  for (const LocalShift &ls : localCols) {
    AZUREProfile::Timer profileColumn(config.profile, AZUREProfile::kShiftLocal);
    const std::vector<int> &segs = ls.target.segs;
    std::vector<vector_r> rP(segs.size()), rM(segs.size());
    const double x0 = full[ls.f];
    // A step that keeps the points inside their planned nodes (ShiftedFold.h, question 5); none: the full path below
    const double h = azure2::InsideShiftStep(*ld, segs, x0, 1.0e-6 * (std::fabs(x0) + 1.0),
                                             [&](double v) { ls.target.Set(ld, v); moveSegments(segs); });
    if (h == 0.0) { fullShifts.push_back(ls.f); continue; }
    const double d0 = ls.target.Get(ld);
    ls.target.Set(ld, x0 + h);
    moveSegments(segs);
    groupResiduals(segs, rP);
    ls.target.Set(ld, x0 - h);
    moveSegments(segs);
    groupResiduals(segs, rM);
    ls.target.Set(ld, d0);
    moveSegments(segs);
    for (size_t m = 0; m < segs.size(); m++) {
      const int sNum = segs[m];
      const size_t r0 = (size_t)rowStart[sNum], nSegRows = (size_t)(rowStart[sNum + 1] - rowStart[sNum]);
      if (rP[m].size() == nSegRows && rM[m].size() == nSegRows)
        for (size_t r = 0; r < nSegRows; r++)
          jac[(r0 + r) * (size_t)nCols + (size_t)ls.packed] = (rP[m][r] - rM[m][r]) / (2.0 * h);
    }
  }
  count.local = (int)localCols.size();

  vector_r rPlus, rMinus;
  for (int f : fullShifts) {
    AZUREProfile::Timer profileColumn(config.profile, AZUREProfile::kShiftFull);
    const int packed = pmap.FullToPacked(f);
    const double x0 = full[f];
    // A few eV: far inside any quoted beam-energy calibration uncertainty,
    // far outside the noise of the forward model.
    const double h = 1.0e-6 * (std::fabs(x0) + 1.0);
    count.full++;
    auto residualsAt = [&](double value, vector_r &out) {
      vector_r fk = full;
      fk[f] = value;
      lc->FillCompoundFromParams(fk, state);
      ld->FillNormsFromParams(fk);
      ld->FillEnergyShiftsFromParams(fk, results, lc, &state, &inputs, &config);
      allResiduals(out);
    };
    residualsAt(x0 + h, rPlus);
    residualsAt(x0 - h, rMinus);
    if (rPlus.size() == nRes && rMinus.size() == nRes)
      for (size_t r = 0; r < nRes; r++)
        jac[r * (size_t)nCols + (size_t)packed] = (rPlus[r] - rMinus[r]) / (2.0 * h);
  }
  // Leave the model at the parameters that were asked for.
  lc->FillCompoundFromParams(full, state);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, results, lc, &state, &inputs, &config);
  return count;
}

bool AZURECalc::ResidualJacobian(const vector_r &full, vector_r &residuals,
                                 vector_r &jac, std::vector<int> &packedToFull) const {
  if (!FiniteVector(full)) return false;  // a non-finite parameter: a failed evaluation (FiniteVector)
  // Pooled working copies keep their applied energy shifts and mapping between calls; a fresh
  // Clone() re-shifted every shifted segment and rebuilt the mapping each time (~1 s here).
  const CompoundSystem *lc = compound();
  std::unique_ptr<LevelState> pooled = GetPooledState();
  LevelState &ls = *pooled;
  WorkingData *working = GetPooledData();
  WorkingCopy *ld = working->data.get();
  azure2::Results &results = working->results;
  lc->FillCompoundFromParams(full, ls);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, results, lc, &ls, &inputs_, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) {
    lc->CalcShiftFunctions(configure(), ls);
    shiftDeriv = BuildShiftDerivTable(lc, ls, configure());
    sdp = &shiftDeriv;
  }

  // Fixed-parameter mask, in the Minuit parameter order: the minimizer's.
  const std::vector<bool> fixed = FixedMask();

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed);
  int nCols = 0;
  const auto tj0 = std::chrono::steady_clock::now();
  bool ok = ComputeResidualJacobian(lc, ls, ld, results, configure(), pmap, sdp, residuals, jac, nCols);
  const auto tj1 = std::chrono::steady_clock::now();
  int nShiftLocal = 0, nShiftFull = 0;

  // Energy-shift columns by central differences (FillShiftColumns, shared with
  // AZUREAPI::CalculateResidualJacobianRWA).  The adjoint does not cover a
  // shift (it moves the energy of every quantity of the segment), and left as
  // zero the Levenberg-Marquardt and GSL drivers flag the shifts as
  // data-insensitive and never move them.
  if (ok) {
    const ShiftColumnCount cols = FillShiftColumns(lc, ls, inputs_, ld, results, configure(), pmap, fixed, full, nullptr,
                                                   residuals.size(), nCols, jac);
    nShiftLocal = cols.local;
    nShiftFull = cols.full;
  }
  if (azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kGradTiming)) {
    const int c = ++configure().diagnostics->residualJacobianTimings;
    if (c <= 5 || c % 50 == 0) {
      auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
      std::cerr << "[grad-timing] residual Jacobian " << c << ": adjoint rows " << ms(tj0, tj1) << " ms, energy-shift columns "
                << ms(tj1, std::chrono::steady_clock::now()) << " ms (" << nShiftLocal << " segment-local, " << nShiftFull
                << " full)" << std::endl;
    }
  }

  // Data covariance: whiten the data rows (and their Jacobian rows, energy-shift columns
  // included), w = W r per group (W^T W = M^-1), so sum w^2 = r^T M^-1 r.
  if (ok && ld->Data().HasCovariance()) ld->WhitenResidualRows(residuals, &jac, nCols);

  // The model data (T9): their rows after the data rows, with their Jacobian rows (no energy-shift column).
  if (ok) ok = Objective::Engine().AppendModelDataRows(ObjectiveContext{lc, ls, inputs_, ld, results, full, configure()}, pmap, sdp,
                                                       nCols, residuals, jac);

  // Prior parameter covariance: its n whitened rows L^-1 (p - mu) last, Jacobian L^-1 on its parameters' columns.
  // The least-squares drivers (LM, GSL), the Jacobian-path gradient and G2 take the prior from these rows.
  if (ok)
    if (const PriorCovariance *pc = ld->Data().GetPriorCovariance())
      pc->AppendRows(full, [&](int f) { return f < pmap.NumFull() ? pmap.FullToPacked(f) : -1; }, nCols, residuals,
                     &jac);
  // The terms' own rows (ObjectiveTerm::AppendRows: Park's J > 0 wall, the gamma_decay records, the abundance sums;
  // objective.md D2, D7), last.
  if (ok)
    Objective::Engine().AppendModelRows(ObjectiveContext{lc, ls, inputs_, ld, results, full, configure()}, &pmap, nCols,
                                        residuals, &jac);

  packedToFull.resize(pmap.NumPacked());
  for (int a = 0; a < pmap.NumPacked(); a++) packedToFull[a] = pmap.PackedToFull(a);

  ReturnPooledState(std::move(pooled));
  ReturnPooledData(working);
  return ok;
}

bool AZURECalc::ResidualsOnly(const vector_r &full, vector_r &residuals) const {
  if (!FiniteVector(full)) return false;  // a non-finite parameter: a failed evaluation (FiniteVector)
  // Forward-only counterpart of ResidualJacobian (no adjoint): same row order
  // and residual r = (model/J - data*n)/(err*n) as ComputeResidualJacobian.
  // Pooled working copies keep their applied energy shifts and mapping between calls; a fresh
  // Clone() re-shifted every shifted segment and rebuilt the mapping each time (~1 s here).
  const CompoundSystem *lc = compound();
  std::unique_ptr<LevelState> pooled = GetPooledState();
  LevelState &ls = *pooled;
  WorkingData *working = GetPooledData();
  WorkingCopy *ld = working->data.get();
  azure2::Results &results = working->results;
  lc->FillCompoundFromParams(full, ls);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, results, lc, &ls, &inputs_, &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure(), ls, CompoundSystem::ShiftOrders::kFirstDerivativesOnly);

  residuals.clear();
  for (int i = 1; i <= ld->Data().NumSegments(); i++) {
    const Segment *segment = ld->Data().GetSegment(i);
    if (!segment) continue;
    const double norm = ld->State(*segment).norm;
    for (int pid = 0; pid < segment->NumPoints(); pid++) {
      azure2::EvaluateModelPoint(*ld, *segment, pid, lc, ls, configure(), results, azure2::LeafValues::Calculate);
      residuals.push_back(azure2::PointResidual(*ld, *segment, pid, results, norm).Pull());
    }
  }
  if (ld->Data().HasCovariance()) ld->WhitenResidualRows(residuals, nullptr, 0);
  // The model data (T9), as in ResidualJacobian.
  Objective::Engine().ModelDataResiduals(ObjectiveContext{lc, ls, inputs_, ld, results, full, configure()}, residuals);
  // Prior parameter covariance, as the last rows of ResidualJacobian.
  if (const PriorCovariance *pc = ld->Data().GetPriorCovariance())
    pc->AppendRows(full, [](int) { return -1; }, 0, residuals, nullptr);
  // The terms' own rows, as the last rows of ResidualJacobian.
  Objective::Engine().AppendModelRows(ObjectiveContext{lc, ls, inputs_, ld, results, full, configure()}, nullptr, 0,
                                      residuals, nullptr);
  ReturnPooledState(std::move(pooled));
  ReturnPooledData(working);
  return true;
}
