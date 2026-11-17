// AZUREAPI, the derivatives of the objective: the analytic gradient of the chi-squared (Chi2GradEGammaNorm,
// CalculateChi2GradRWA) and of the whole objective (CalculateObjectiveGradRWA), the residual Jacobian with its
// energy-shift and finite-difference columns (CalculateResidualJacobianRWA), the gradient of each point's model
// value (CalculateModelGradientsRWA), and the rows of the THM experiments they share (ComputeTHMRows).  One unit
// because all of them differentiate the value of AZUREAPIChi2.cpp through the same adjoint (AZUREGrad.h) and
// the same index map.
#include <algorithm>
#include <stdexcept>
#include "Objective.h"
#include "ParameterLimitsManager.h"
#include "AZUREAPI.h"
#include "CompoundSystemSetUp.h"
#include "AdaptiveExtrapolation.h"
#include "ChannelFunc.h"
#include "ModelPass.h"
#include "ScatteringLength.h"
#include "AZUREParams.h"

#include "GSLException.h"
#include "ECAmplitudeCache.h"

#include "Config.h"
#include "CompoundSystem.h"
#include "LevelReport.h"
#include "CaptureAyRequests.h"
#include "ProjectSetup.h"
#include "PPair.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "Unit.h"
#include "TargetEffect.h"
#include "AMatrixFunc.h"
#include "AZUREGrad.h"
#include "CovarianceBand.h"
#include "PriorCovarianceFile.h"
#include "AZURECalc.h"
#include "CoulFunc.h"
#include "ECIntegral.h"
#include "Profile.h"
#include "InputFiles.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <new>
#include <cmath>
#include "AZUREAPIInternal.h"

std::vector<AZUREAPI::THMRows> AZUREAPI::ComputeTHMRows(const vector_r &full,
                                                        const ParamIndexMap &pmap) {
  return ::ComputeTHMRows(compound(), state(), inputs(), data(), results(), configure(), full, pmap);
}

bool AZUREAPI::Chi2GradEGammaNorm(const vector_r &full, vector_r &gradFull,
                                  double &chi2Out) {
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);
  // Only one API request runs at a time, so operate on the canonical
  // compound/data in place (like UpdateSegments) rather than cloning; every
  // entry point re-fills from its own parameters before use.
  const CompoundSystem *lc = compound();
  WorkingCopy *ld = data();
  lc->FillCompoundFromParams(full, state());
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, results(), lc, &state(), &inputs(), &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (brune) {
    lc->CalcShiftFunctions(configure(), state());
    shiftDeriv = BuildShiftDerivTable(lc, state(), configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed_);
  GradAccum accum;
  accum.Init(lc, state());

  // chi2 = sum (fit/J - data*n)^2/(err*n)^2 (ModelPass.h DataResidual).  d(chi2)/d(model) = 2 r / err^2 / J;
  // the same model gives the data-term norm gradient, accumulated per segment.
  std::vector<double> normData(ld->Data().NumSegments() + 1, 0.0);
  double chi2 = 0.0;
  // A data covariance couples the points: their cotangent 2 (M^-1 r)_i / (err n J) needs
  // every residual of the group first (one forward pass of those points), and the
  // group's chi2 is r^T M^-1 r in place of the sum of r_i^2 accumulated below.
  if (ld->Data().HasCovariance()) {
    ld->PrepareCovarianceSeeds(lc, state(), configure(), results());
    chi2 += ld->CovarianceChi2Correction(results());
  }
  FitBarFn fb = [&](const Segment *seg, int i, int pid, double model) -> double {
    Unit *pt = seg->GetPoint(*ld, pid + 1);
    if (!pt) return 0.0;
    double norm = ld->State(*seg).norm;
    const double dataval = seg->LabValue(pid), labErr = seg->LabError(pid), jac = Segment::Jacobian(ld->State(*seg), 0, pid);
    const azure2::Residual res = azure2::DataResidual(model, jac, dataval, labErr, norm);
    const double r = res.difference, err = res.error;
    if (err == 0.0) return 0.0;
    // chi2 falls out of the same residual the gradient uses (so value and
    // gradient stay consistent). Runs in the parallel point loop, so guard it.
#pragma omp atomic
    chi2 += (r * r) / (err * err);
    if (seg->GetCovarianceGroup() >= 0) {
      // d chi2/d r_i = 2 z_i, z = M^-1 r; d r_i/d model = 1/(err n J), d r_i/d n = -(model/J)/(err n^2).
      const double z = pt->GetCovarianceSeed(results());
      const int o = ld->Data().NormOwner(i);
      if (norm != 0.0 && o >= 1 && o < (int)normData.size() && ld->Data().GetSegment(o)->IsVaryNorm()) {
        const double dNorm = -2.0 * z * res.model / (labErr * norm * norm);
#pragma omp atomic
        normData[o] += dNorm;
      }
      return 2.0 * z / err / res.jacobian;
    }
    // A norm-group follower's points move its leader's norm (o), the one parameter.
    const int o = ld->Data().NormOwner(i);
    if (norm != 0.0 && o >= 1 && o < (int)normData.size() && ld->Data().GetSegment(o)->IsVaryNorm()) {
      double e2 = labErr * labErr;
      // fitBarFn runs inside the parallel point loop of AccumulateEGammaGradient,
      // and all points of a segment share the same normData[o], so guard the
      // accumulation.
      double dNorm = -2.0 * r * dataval / (e2 * norm * norm) - 2.0 * r * r / (e2 * norm * norm * norm);
#pragma omp atomic
      normData[o] += dNorm;
    }
    return 2.0 * r / (err * err) / res.jacobian;
  };

  bool ok = AccumulateEGammaGradient(lc, state(), ld, results(), configure(), pmap, sdp, fb, accum);
  // The gradients of the terms of the model alone (ObjectiveTerm::AddGradient; false -> finite differences).
  if (ok)
    ok = Objective::Engine().AddGradient(ObjectiveContext{lc, state(), inputs(), ld, results(), full, configure()}, pmap, sdp,
                                         accum, TermSet::kModelOnly);
  if (ok) {
    // Their values: the model terms of the objective (T6, T7, T9, T11; before 2026-11-09 added here in the order 9, 6, 7).
    chi2 += Objective::Engine().ModelChi2(ObjectiveContext{lc, state(), inputs(), ld, results(), full, configure()});
    // The normalizations first, assigned: the scatter after them adds a term's part on a norm (GradAccum::direct).
    for (int s = 1; s <= ld->Data().NumSegments(); s++) {
      const Segment *seg = ld->Data().GetSegment(s);
      if (!seg || !seg->IsVaryNorm() || seg->IsNormGroupFollower()) continue;  // a follower's norm is fixed
      int idx = pmap.NormIndex(s);
      if (idx >= 0 && idx < (int)gradFull.size()) gradFull[idx] = normData[s];
    }
    accum.Scatter(pmap, gradFull);
    chi2Out = chi2;
  }
  return ok;
}

vector_r AZUREAPI::CalculateChi2GradRWA(const vector_r &params) {
  vector_r full = MapPackedToFull(params, all_rwa_, fixed_);

  // The analytic gradient pass runs the full forward model at every point, so
  // chi2 comes out as a byproduct -- no separate forward pass. Only fall back to
  // a standalone chi2 evaluation if the analytic path bails.
  vector_r gradFull(all_rwa_.size(), 0.0);
  double chi2 = 0.0;
  bool eg = Chi2GradEGammaNorm(full, gradFull, chi2);
  if (!eg) chi2 = CalculateChi2RWA(params);

  ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), fixed_);

  // THM (HOES) segments are skipped by the adjoint, so neither their chi2 nor
  // their gradient is in what it returned.  Add both: the chi2 with the norm
  // profiled as the CLI does, the gradient as 2 J^T r over the THM rows (whose
  // J carries the dependence of the profiled scale; at the profile optimum
  // that term drops out of the gradient anyway, sum_i r_i m_i/e_i = 0).
  if (eg) {
    std::vector<THMRows> thm = ComputeTHMRows(full, pmap);
    const int nCols = pmap.NumPacked();
    for (const THMRows &tr : thm) {
      const Segment *seg = data()->Data().GetSegment(tr.segment);
      double c = 0.0;
      if (data()->Data().ThmGroupOf(tr.segment) >= 0) {
        // A THM experiment, profiled by ComputeTHMRows: its residuals.
        for (double r : tr.r) c += r * r;
      } else if (seg->IsProfiledNorm()) {
        c = seg->ProfileNormChiSquared(*data(), results(), results().objective.norm[tr.segment]);
      } else {
        for (double r : tr.r) c += r * r;
      }
      results().objective.segmentChi2[tr.segment] = c;
      chi2 += c;
      for (int col = 0; col < nCols; col++) {
        double g = 0.0;
        for (size_t i = 0; i < tr.r.size(); i++) g += tr.r[i] * tr.J[i * nCols + col];
        gradFull[pmap.PackedToFull(col)] += 2.0 * g;
      }
    }
  }

  // Finite differences for energy shifts (and the whole block if the analytic
  // path bailed), using the scalar chi-squared.
  for (int f = 0; f < (int)full.size() && f < pmap.NumFull(); f++) {
    if (fixed_[f]) continue;
    ParamKind kind = pmap.Desc(f).kind;
    // THM coherent backgrounds (cbackground=) are columns of the THM rows above.
    if (eg && (IsLevelSlot(kind) || kind == ParamKind::Norm ||
               kind == ParamKind::Screening || kind == ParamKind::Doppler || kind == ParamKind::Thickness ||
               kind == ParamKind::ThmCoherent))
      continue;
    int packed = pmap.FullToPacked(f);
    if (packed < 0 || packed >= (int)params.size()) continue;
    double x0 = params[packed];
    double h = 1.0e-6 * (std::fabs(x0) + 1.0);
    vector_r pp = params;
    pp[packed] = x0 + h;
    vector_r pm = params;
    pm[packed] = x0 - h;
    gradFull[f] = (CalculateChi2RWA(pp) - CalculateChi2RWA(pm)) / (2.0 * h);
  }

  vector_r out;
  out.reserve(1 + params.size());
  out.push_back(chi2);
  for (int i = 0; i < (int)all_rwa_.size(); i++)
    if (!fixed_[i]) out.push_back(gradFull[i]);
  return out;
}

vector_r AZUREAPI::CalculateObjectiveGradRWA(const vector_r &params) {
  // The gradient of CalculateChi2RWA's terms (data, Park, abundances, scattering lengths) ...
  vector_r out = CalculateChi2GradRWA(params);
  const vector_r full = MapPackedToFull(params, all_rwa_, fixed_);
  vector_r gradFull(all_rwa_.size(), 0.0);
  // ... plus the diagonal priors (abundances are already in), the nuisance priors and the prior covariance.
  const ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), fixed_);
  ForEachDiagonalPrior(data(), pmap, [&](int idx, double nominal, double sigma) {
    if (idx < 0 || idx >= (int)gradFull.size() || pmap.Desc(idx).kind == ParamKind::Abundance) return;
    gradFull[idx] += GaussianPriorGradient(full[idx], nominal, sigma);
  });
  AddNuisanceGradient(full, gradFull);
  if (const PriorCovariance *pc = data()->Data().GetPriorCovariance()) pc->AddGradient(full, gradFull);
  for (int i = 0, k = 1; i < (int)all_rwa_.size(); i++)
    if (!fixed_[i]) out[k++] += gradFull[i];
  out[0] = CalculateObjectiveRWA(params);
  // ... plus the gradients the other terms state (ObjectiveTerm::AddGradient of a term not of the model alone; none of
  // the engine's own terms has one), at x, where CalculateObjectiveRWA has just left the session.
  GradAccum accum;
  accum.Init(compound(), state());
  bool added = false;
  if (!Objective::Engine().AddGradient(ObjectiveContext{compound(), state(), inputs(), data(), results(), full, configure()},
                                       pmap, nullptr, accum, TermSet::kOthers, &added))
    throw std::runtime_error("CalculateObjectiveGradRWA: a term of the objective has no gradient here");
  if (added) {
    vector_r termsFull(all_rwa_.size(), 0.0);
    accum.Scatter(pmap, termsFull);
    for (int i = 0, k = 1; i < (int)all_rwa_.size(); i++)
      if (!fixed_[i]) out[k++] += termsFull[i];
  }
  return out;
}

vector_r AZUREAPI::CalculateResidualJacobianRWA(const vector_r &params, const std::vector<int> &skipColumns) {
  AZUREProfile::ScopedReport profileReport(configure().profile, "api jacobian");
  vector_r full = MapPackedToFull(params, all_rwa_, fixed_);

  // One request at a time: operate on the canonical compound/data in place
  // (re-filled here from these parameters) instead of cloning.
  const CompoundSystem *lc = compound();
  WorkingCopy *ld = data();
  lc->FillCompoundFromParams(full, state());
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, results(), lc, &state(), &inputs(), &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) {
    lc->CalcShiftFunctions(configure(), state());
    shiftDeriv = BuildShiftDerivTable(lc, state(), configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed_);
  vector_r residuals, jacobian;
  int nCols = 0;
  // THM rows are left zero by the adjoint and filled here: residuals with the
  // profiled norm, Jacobian from the HOES adjoint plus the
  // analytic dependence of the profiled scale (see ComputeTHMRows).
  bool ok = ComputeResidualJacobian(lc, state(), ld, results(), configure(), pmap, sdp, residuals, jacobian, nCols,
                                    /*skipTHM=*/true);

  if (!ok) return vector_r{-1.0};

  {
    std::vector<THMRows> thm = ComputeTHMRows(full, pmap);
    for (const THMRows &tr : thm) {
      for (size_t i = 0; i < tr.r.size(); i++) {
        const size_t row = (size_t)tr.firstRow + i;
        if (row >= residuals.size()) continue;
        residuals[row] = tr.r[i];
        for (int c = 0; c < nCols; c++) jacobian[row * (size_t)nCols + c] = tr.J[i * nCols + c];
      }
    }
  }

  // ---- Energy-shift columns, by central differences -----------------------
  //
  // The adjoint covers level energies, reduced widths and normalizations. An
  // energy shift is a different derivative: it translates the energy axis of a
  // whole segment, so what is wanted is d(model)/dE rather than d(model)/d(a
  // parameter). Every energy-dependent quantity moves with it --- the level
  // matrix, the penetrabilities and shift functions, the Coulomb and
  // hard-sphere phases, the external-capture amplitudes, and, where a segment
  // carries target integration, the sub-point grid itself (which brings
  // boundary terms with it). AZURE2 applies a shift by rebuilding all of that
  // in UpdatePointEnergiesWithShift, so there is no cheap analytic route
  // through the existing forward code.
  //
  // Returning these columns as zero, which is what this function used to do,
  // is the worst of the options: a least-squares driver handed a zero column
  // simply never moves that parameter, converges, and reports success.
  //
  // FillShiftColumns (dev, shared with the CLI's AZURECalc::ResidualJacobian)
  // differences them forward only: a self-contained segment's shift on that
  // segment's rows alone, anything else with one parallel point pass per side.
  // It forms the plain residuals, so a project with THM segments (rows with
  // their profiled norm, a THM experiment's shared profile) takes the forward
  // residuals of EvaluateFilledChi2 instead, two evaluations per column.
  {
    // Columns the caller holds fixed (skipColumns): not differenced, left zero.
    std::vector<char> skip(params.size(), 0);
    for (int c : skipColumns)
      if (c >= 0 && c < (int)params.size()) skip[c] = 1;
    bool anyTHM = false;
    for (int i = 1; i <= ld->Data().NumSegments(); i++) anyTHM = anyTHM || (ld->Data().GetSegment(i) && ld->Data().GetSegment(i)->IsTHM());
    if (!anyTHM) {
      FillShiftColumns(lc, state(), inputs(), ld, results(), configure(), pmap, fixed_, full, &skip, residuals.size(), nCols, jacobian);
    } else {
      const size_t nRes = residuals.size();
      vector_r rPlus, rMinus;
      for (int f = 0; f < pmap.NumFull(); f++) {
        if (fixed_[f]) continue;
        const ParamKind kind = pmap.Desc(f).kind;
        if (kind != ParamKind::EnergyShift && kind != ParamKind::EnergyShiftSqrt && kind != ParamKind::EnergyShiftFn &&
            kind != ParamKind::Density && kind != ParamKind::Apparatus)
          continue;
        const int packed = pmap.FullToPacked(f);
        if (packed < 0 || packed >= (int)params.size() || skip[packed]) continue;

        const double x0 = params[packed];
        // A few eV for a shift: far inside any quoted beam-energy calibration
        // uncertainty, far outside the noise of the forward model; a density
        // (atoms/cm^2) relative to its value.
        const double h = (kind == ParamKind::Density) ? 1.0e-6 * std::fabs(x0) : 1.0e-6 * (std::fabs(x0) + 1.0);
        if (!(h > 0.0)) continue;

        // The forward residuals: the same rows as ComputeResidualJacobian's,
        // THM rows included with their profiled norm, and no adjoint to pay for.
        auto residualsAt = [&](double value, vector_r &out) -> bool {
          vector_r pk = params;
          pk[packed] = value;
          FillFromFullRWA(MapPackedToFull(pk, all_rwa_, fixed_));
          EvaluateFilledChi2(&out);
          return true;
        };

        const bool okp = residualsAt(x0 + h, rPlus);
        const bool okm = residualsAt(x0 - h, rMinus);
        if (okp && okm && rPlus.size() == nRes && rMinus.size() == nRes) {
          for (size_t r = 0; r < nRes; r++)
            jacobian[r * (size_t)nCols + (size_t)packed] = (rPlus[r] - rMinus[r]) / (2.0 * h);
        }
      }
      // Leave the model at the parameters that were asked for.
      FillFromFullRWA(full);
      EvaluateFilledChi2(nullptr);
    }
  }

  // Data covariance: whiten the data rows, w = W r per group (W^T W = M^-1, sum w^2 = r^T M^-1 r).
  if (ld->Data().HasCovariance()) ld->WhitenResidualRows(residuals, &jacobian, nCols);

  // The model data (T9): rows after the data rows, as AZURECalc::ResidualJacobian.
  if (!Objective::Engine().AppendModelDataRows(ObjectiveContext{lc, state(), inputs(), ld, results(), full, configure()}, pmap, sdp, nCols,
                                               residuals, jacobian))
    return vector_r{-1.0};

  // Prior parameter covariance: its whitened rows L^-1 (p - mu) last, Jacobian L^-1 on its parameters' columns.
  if (const PriorCovariance *pc = ld->Data().GetPriorCovariance())
    pc->AppendRows(full, [&](int f) { return f < pmap.NumFull() ? pmap.FullToPacked(f) : -1; }, nCols, residuals,
                   &jacobian);

  vector_r out;
  out.reserve(2 + residuals.size() + jacobian.size());
  out.push_back((double)residuals.size());
  out.push_back((double)nCols);
  out.insert(out.end(), residuals.begin(), residuals.end());
  out.insert(out.end(), jacobian.begin(), jacobian.end());
  return out;
}

vector_r AZUREAPI::CalculateModelGradientsRWA(const vector_r &params) {
  vector_r full = MapPackedToFull(params, all_rwa_, fixed_);

  // Same preparation as CalculateResidualJacobianRWA: one request at a time, so
  // the canonical compound/data are re-filled in place from these parameters.
  const CompoundSystem *lc = compound();
  WorkingCopy *ld = data();
  lc->FillCompoundFromParams(full, state());
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, results(), lc, &state(), &inputs(), &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) {
    lc->CalcShiftFunctions(configure(), state());
    shiftDeriv = BuildShiftDerivTable(lc, state(), configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed_);

  // A band is sensitive only to the R-matrix (and THM coherent-background)
  // parameters, and covariance.dat spans exactly those columns
  // (BandPackedColumns) -- so reduce each full packed row to them here
  // rather than shipping zero columns for every normalization.
  const std::vector<int> rc = BandPackedColumns(pmap);
  const int nCols = (int)rc.size();

  // THM (HOES) rows come from the HOES adjoint (the T-matrix adjoint does not
  // cover them); see ComputeTHMRows.
  std::map<Unit *, vector_r> grad;
  if (!ComputeModelGradients(lc, state(), ld, results(), configure(), pmap, sdp, grad, /*skipTHM=*/true))
    return vector_r{-1.0};
  {
    std::vector<THMRows> thm = ComputeTHMRows(full, pmap);
    const int nPacked = pmap.NumPacked();
    for (const THMRows &tr : thm) {
      const Segment *seg = ld->Data().GetSegment(tr.segment);
      size_t i = 0;
      for (int pid = 0; pid < seg->NumPoints(); pid++) {
        Unit *pt = seg->GetPoint(*ld, pid + 1);
        if (!pt) continue;
        std::map<Unit *, vector_r>::iterator it = grad.find(pt);
        if (it != grad.end() && i < tr.m.size())
          for (int c = 0; c < nPacked && c < (int)it->second.size(); c++)
            it->second[c] = tr.Jm[i * nPacked + c];
        i++;
      }
    }
  }

  // Walk the segments exactly as UpdateSegments does, so row k of segment s
  // lines up with GET_CALCULATED_SEGMENT's point k: segments sharing a segment
  // key collapse to one calculated segment, and only the first of them is used.
  std::vector<int> counts;
  vector_r rows;
  int prevKey = -1;
  const std::vector<Segment> &segments = ld->Data().GetSegments();
  for (int i = 0; i < (int)segments.size(); ++i) {
    const int newKey = segments[i].GetSegmentKey();
    if (prevKey == newKey) continue;
    prevKey = newKey;

    const UnitRange<Unit> points = segments[i].GetPoints(*ld);
    int n = 0;
    for (int k = 0; k < (int)points.size(); ++k) {
      std::map<Unit *, vector_r>::const_iterator it = grad.find(&points[k]);
      if (it == grad.end()) continue;
      const vector_r &f = it->second;
      for (int c = 0; c < nCols; ++c)
        rows.push_back(rc[c] < (int)f.size() ? f[rc[c]] : 0.0);
      ++n;
    }
    counts.push_back(n);
  }

  vector_r out;
  out.reserve(2 + counts.size() + rows.size());
  out.push_back((double)counts.size());
  out.push_back((double)nCols);
  for (int i = 0; i < (int)counts.size(); ++i) out.push_back((double)counts[i]);
  out.insert(out.end(), rows.begin(), rows.end());
  return out;
}
