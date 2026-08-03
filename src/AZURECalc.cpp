#include "AZURECalc.h"
#include "ParameterLabel.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "ESegment.h"
#include "EPoint.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "AZUREGrad.h"
#include "CovarianceBand.h"
#include "GSLException.h"
#include <iostream>
#include <iomanip>
#include <thread>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// chi^2 contribution of a single segment: its data term plus its own
// energy-shift nuisance penalty, evaluated at the segment's current shift.
// Mirrors exactly the per-segment math in AZURECalc::Chi2Value (data term +
// lines guarded by IsVaryEnergyShift), so a finite difference of this over the
// segment's shift equals a finite difference of the whole chi^2 -- *provided*
// no other segment's contribution depends on this shift.  That holds when the
// segment is self-contained: no components and no cross-segment energy mapping
// (see the eligibility guard in Gradient()).  All other segments then produce
// identical +h/-h contributions that cancel in the central difference.
double SegmentLocalChi2(ESegment* segment, CNuc* lc, const Config& config,
                        EData* ld) {
  double chi = 0.0;
  const double norm = segment->GetNorm();
  for(int pid = 0; pid < segment->NumPoints(); pid++) {
    double th = segment->CalculateTheoreticalCrossSection(pid, lc, config, ld);
    EPoint* pt = segment->GetPoint(pid + 1);
    if(!pt) continue;
    pt->SetFitCrossSection(th);
    double r = th - pt->GetCMCrossSection() * norm;
    double err = pt->GetCMCrossSectionError() * norm;
    if(err != 0.0) chi += (r * r) / (err * err);
  }
  if(segment->IsVaryEnergyShift()) {
    double sh = segment->GetEnergyShift();
    double shn = segment->GetNominalEnergyShift();
    double she = segment->GetEnergyShiftError();
    if(she != 0.0) chi += pow((sh - shn) / she, 2.0);
  }
  return chi;
}

}  // namespace

double AZURECalc::operator()(const vector_r& p) const {

  int thisIteration=data()->Iterations();
  data()->Iterate();
  bool isFit=data()->IsFit();

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  if(isFit) {

    // New multithreading with object pools
    if (!pools_initialized_) {
      InitializePools();
    }
    
    // Get objects from pool
    localCompound = GetPooledCNuc();
    localData = GetPooledEData();

    // Old multithreading
    //localCompound = compound()->Clone();
    //localData = data()->Clone();
  } else {
    localCompound = compound();
    localData = data();
  }

  //Fill Compound Nucleus From Minuit Parameters
  localCompound->FillCompoundFromParams(p);
  localData->FillNormsFromParams(p);
  localData->FillEnergyShiftsFromParams(p,localData,localCompound,&configure());
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  // Process segments with components - use new integrated calculation method
  double chiSquared=0.0;
  for(int i = 1; i <= localData->NumSegments(); i++) {
    ESegment* segment = localData->GetSegment(i);
    if(segment) {
      // Recalculate points using the new combined calculation method
      for(int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        double theoreticalValue = segment->CalculateTheoreticalCrossSection(pointIdx, localCompound, configure(), localData);
        EPoint* point = segment->GetPoint(pointIdx + 1);
        if(point) {
          point->SetFitCrossSection(theoreticalValue);
        }
      }
      
      // Recalculate chi-squared for this segment with components
      double segmentChiSquared = 0.0;
      for(int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        EPoint* point = segment->GetPoint(pointIdx + 1);
        if(point) {
          double residual = point->GetFitCrossSection() - point->GetCMCrossSection() * segment->GetNorm();
          double error = point->GetCMCrossSectionError() * segment->GetNorm();
          if(error != 0.0) {
            segmentChiSquared += (residual * residual) / (error * error);
          }
        }
      }

      // Add normalization chi-squared contribution
      double dataNorm=segment->GetNorm();
      double dataNormNominal=segment->GetNominalNorm();
      double dataNormError=dataNormNominal/100.*segment->GetNormError();
      if(dataNormError!=0.){
	      chiSquared += pow((dataNorm-dataNormNominal)/dataNormError,2.0);
      }

      // Add energy shift chi-squared contribution
      if(segment->IsVaryEnergyShift()) {
        double energyShift=segment->GetEnergyShift();
        double energyShiftNominal=segment->GetNominalEnergyShift();
        double energyShiftError=segment->GetEnergyShiftError();
        if(energyShiftError!=0.){
          chiSquared += pow((energyShift-energyShiftNominal)/energyShiftError,2.0);
        }
      }

      segment->SetSegmentChiSquared(segmentChiSquared);
      chiSquared += segmentChiSquared;
    }
  }

  // Add nuisance parameter chi-squared contributions
  if(limitsManager_) {
    chiSquared += CalculateNuisanceChiSquared(p);
  }

  if(!localData->IsErrorAnalysis()&&thisIteration!=0) {
    if(thisIteration%10==0) configure().outStream
			       << "\r\tIteration: " << std::setw(6) << thisIteration
			       << " Chi-Squared: " << chiSquared;  configure().outStream.flush();

    if(thisIteration%kOutputInterval==0) WriteIterationOutput(p);
  }
  if(isFit) {

    // New multithreading with object pools
    ReturnPooledCNuc(localCompound);
    ReturnPooledEData(localData);
    
    // Old multithreading
    //delete localCompound;
    //delete localData;
  }

  // Make a check if chiSquared is NaN
  if(std::isnan(chiSquared)) {
      // In that case return infinite since MINUIT2 can have issues with NaN values
      return std::numeric_limits<double>::infinity();
  }

  if(configure().stopFlag&&isFit) return 0.;
  else return chiSquared;
}

void AZURECalc::WriteIterationOutput(const vector_r& p) const {
  // Work on clones: this runs in the middle of a fit and must not disturb the
  // objects the minimizer is stepping.
  CNuc* lc = compound()->Clone();
  EData* ld = data()->Clone();

  try {
    lc->FillCompoundFromParams(p);
    ld->FillNormsFromParams(p);
    ld->FillEnergyShiftsFromParams(p, ld, lc, &configure());
    if(configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure());

    for(int i = 1; i <= ld->NumSegments(); i++) {
      ESegment* segment = ld->GetSegment(i);
      if(!segment) continue;
      for(int pid = 0; pid < segment->NumPoints(); pid++) {
        EPoint* pt = segment->GetPoint(pid + 1);
        if(pt) pt->SetFitCrossSection(
            segment->CalculateTheoreticalCrossSection(pid, lc, configure(), ld));
      }
    }

    AZUREParams params;
    lc->FillMnParams(params.GetMinuitParams(), &configure());
    ld->FillMnParams(params.GetMinuitParams());
    WriteParameters(params, configure());
    ld->WriteOutputFiles(configure(), true);
    lc->TransformOut(configure());
    lc->PrintTransformParams(configure());
  } catch(...) {
    // An intermediate snapshot is a convenience, never a reason to abort a fit.
  }

  delete lc;
  delete ld;
}

double AZURECalc::Chi2Value(const vector_r& p) const {
  // Side-effect-free chi-squared (mirrors the chi-squared math of operator()).
  CNuc* lc = compound()->Clone();
  EData* ld = data()->Clone();
  lc->FillCompoundFromParams(p);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, ld, lc, &configure());
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure());

  double chiSquared = 0.0;
  for(int i = 1; i <= ld->NumSegments(); i++) {
    ESegment* segment = ld->GetSegment(i);
    if(!segment) continue;
    double segChi = 0.0;
    for(int pid = 0; pid < segment->NumPoints(); pid++) {
      double th = segment->CalculateTheoreticalCrossSection(pid, lc, configure(), ld);
      EPoint* pt = segment->GetPoint(pid + 1);
      if(pt) {
        pt->SetFitCrossSection(th);
        double r = th - pt->GetCMCrossSection() * segment->GetNorm();
        double err = pt->GetCMCrossSectionError() * segment->GetNorm();
        if(err != 0.0) segChi += (r * r) / (err * err);
      }
    }
    double dataNorm = segment->GetNorm();
    double nom = segment->GetNominalNorm();
    double nerr = nom / 100.0 * segment->GetNormError();
    if(nerr != 0.0) chiSquared += pow((dataNorm - nom) / nerr, 2.0);
    if(segment->IsVaryEnergyShift()) {
      double sh = segment->GetEnergyShift();
      double shn = segment->GetNominalEnergyShift();
      double she = segment->GetEnergyShiftError();
      if(she != 0.0) chiSquared += pow((sh - shn) / she, 2.0);
    }
    chiSquared += segChi;
  }
  if(limitsManager_) chiSquared += CalculateNuisanceChiSquared(p);

  delete lc;
  delete ld;
  return chiSquared;
}

std::vector<double> AZURECalc::Gradient(const std::vector<double>& p) const {
  std::vector<double> grad(p.size(), 0.0);
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);

  // --- Analytic energy / reduced-width block via the shared adjoint engine. ---
  CNuc* lc = compound()->Clone();
  EData* ld = data()->Clone();
  lc->FillCompoundFromParams(p);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r* sdp = nullptr;
  if(brune) {
    lc->CalcShiftFunctions(configure());
    shiftDeriv = BuildShiftDerivTable(lc, configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, std::vector<bool>());
  GradAccum accum;
  accum.Init(lc);

  // chi^2 = sum (fit - data*n)^2 / (err*n)^2 + penalties.  The data-term model
  // cotangent is d(chi2)/d(model) = 2 r / err^2.  The same per-point model also
  // gives the analytic normalization gradient, accumulated per segment here so
  // it needs no extra forward pass:
  //   d(chi2)/d(n_s) = sum_pts [ -2 r data/(cmErr^2 n^2) - 2 r^2/(cmErr^2 n^3) ].
  std::vector<double> normData(ld->NumSegments() + 1, 0.0);
  FitBarFn fitBarFn = [&](ESegment* seg, int i, int pid, double model) -> double {
    EPoint* pt = seg->GetPoint(pid + 1);
    if(!pt) return 0.0;
    double norm = seg->GetNorm();
    double dataval = pt->GetCMCrossSection();
    double cmErr = pt->GetCMCrossSectionError();
    double r = model - dataval * norm;
    double err = cmErr * norm;
    if(err == 0.0) return 0.0;
    if(seg->IsVaryNorm() && norm != 0.0 && i >= 1 && i < (int)normData.size()) {
      double e2 = cmErr * cmErr;
      // fitBarFn runs inside the parallel point loop of AccumulateEGammaGradient,
      // and all points of a segment share the same normData[i], so guard the
      // accumulation.
      double dNorm = -2.0 * r * dataval / (e2 * norm * norm)
                     - 2.0 * r * r / (e2 * norm * norm * norm);
#pragma omp atomic
      normData[i] += dNorm;
    }
    return 2.0 * r / (err * err);
  };

  bool eg = AccumulateEGammaGradient(lc, ld, configure(), pmap, sdp, fitBarFn, accum);
  if(!eg && std::getenv("AZURE_GRAD_DEBUG")) {
    std::cerr << "[grad] analytic energy/gamma adjoint bailed -> full finite "
                 "differences (no speed-up). An unsupported segment/config is "
                 "present (RMC, or a not-yet-handled case)." << std::endl;
  }
  if(eg) {
    accum.Scatter(pmap, grad);                 // energies + reduced widths
    AddNuisanceGradient(p, grad);              // nuisance penalty
    for(int s = 1; s <= ld->NumSegments(); s++) {   // normalizations
      ESegment* seg = ld->GetSegment(s);
      if(!seg || !seg->IsVaryNorm()) continue;
      int idx = pmap.NormIndex(s);
      if(idx < 0 || idx >= (int)grad.size()) continue;
      double g = normData[s];
      double n0 = seg->GetNominalNorm();
      double nerr = n0 / 100.0 * seg->GetNormError();
      if(nerr != 0.0) g += 2.0 * (seg->GetNorm() - n0) / (nerr * nerr);
      grad[idx] = g;
    }
  }

  // Fixed-parameter mask: Minuit ignores the gradient of fixed parameters, so
  // do not waste finite differences on them (large fits fix most energy shifts).
  AZUREParams fp;
  compound()->FillMnParams(fp.GetMinuitParams(), &configure());
  data()->FillMnParams(fp.GetMinuitParams());
  const int nMn = fp.GetMinuitParams().Params().size();

  // Segment-local energy-shift gradient is exact only when a shift's effect is
  // confined to its own segment's data term.  A dataset that maps points across
  // segments breaks that confinement, so detect any mapping once and, if
  // present, keep every shift on the (correct) full-dataset finite difference.
  // Likewise, nuisance penalties handled by the limits manager can couple a
  // shift to the global chi^2, so disable the local path when one is active.
  bool datasetHasMapping = false;
  for(int s = 1; s <= ld->NumSegments() && !datasetHasMapping; s++) {
    ESegment* seg = ld->GetSegment(s);
    if(!seg) continue;
    for(int pid = 0; pid < seg->NumPoints(); pid++) {
      EPoint* pt = seg->GetPoint(pid + 1);
      if(pt && pt->IsMapped()) { datasetHasMapping = true; break; }
    }
  }
  const bool localShiftOk = !datasetHasMapping && (limitsManager_ == nullptr);

  // --- Finite differences only for what is left: non-fixed energy shifts, and
  //     (if the analytic path bailed) the energy/gamma and norm blocks too. ---
  for(int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++) {
    if(idx < nMn && fp.GetMinuitParams().Parameter(idx).IsFixed()) continue;
    ParamKind kind = pmap.Desc(idx).kind;
    if(eg && (kind == ParamKind::LevelEnergy || kind == ParamKind::Gamma ||
              kind == ParamKind::Norm)) continue;  // analytic
    double x0 = p[idx];
    double h = 1.0e-6 * (std::fabs(x0) + 1.0);

    // Fast path: a free energy shift on a self-contained segment only moves that
    // segment's points, so finite-difference just that segment's chi^2 on the
    // already-filled clone instead of re-solving the entire dataset twice.
    if(kind == ParamKind::EnergyShift && localShiftOk) {
      int s = pmap.Desc(idx).segment;
      ESegment* seg = (s >= 1) ? ld->GetSegment(s) : nullptr;
      if(seg && !seg->HasComponents() && !seg->IsTotalCapture()) {
        double d0 = seg->GetEnergyShift();
        seg->SetEnergyShift(x0 + h);
        seg->UpdatePointEnergiesWithShift(lc, &configure());
        double chiP = SegmentLocalChi2(seg, lc, configure(), ld);
        seg->SetEnergyShift(x0 - h);
        seg->UpdatePointEnergiesWithShift(lc, &configure());
        double chiM = SegmentLocalChi2(seg, lc, configure(), ld);
        seg->SetEnergyShift(d0);                 // restore base state on the clone
        seg->UpdatePointEnergiesWithShift(lc, &configure());
        grad[idx] = (chiP - chiM) / (2.0 * h);

        // Optional self-check: compare against the exact full-dataset finite
        // difference for this shift.  Enable with AZURE_SHIFT_GRAD_CHECK=1 to
        // validate the fast path on real data (one gradient evaluation), then
        // disable it for the actual fit.
        if(std::getenv("AZURE_SHIFT_GRAD_CHECK")) {
          vector_r pp = p; pp[idx] = x0 + h;
          vector_r pm = p; pm[idx] = x0 - h;
          double gFull = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * h);
          double denom = std::max(1.0, std::fabs(gFull));
          std::cerr << "[shift-grad-check] seg " << s << " param " << idx
                    << ": local=" << grad[idx] << " full=" << gFull
                    << " reldiff=" << std::fabs(grad[idx] - gFull) / denom
                    << std::endl;
        }
        continue;
      }
    }

    vector_r pp = p; pp[idx] = x0 + h;
    vector_r pm = p; pm[idx] = x0 - h;
    grad[idx] = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * h);
  }

  delete lc;
  delete ld;

  return grad;
}

bool AZURECalc::ResidualJacobian(const vector_r& full, vector_r& residuals,
                                 vector_r& jac, std::vector<int>& packedToFull) const {
  CNuc* lc = compound()->Clone();
  EData* ld = data()->Clone();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r* sdp = nullptr;
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) {
    lc->CalcShiftFunctions(configure());
    shiftDeriv = BuildShiftDerivTable(lc, configure());
    sdp = &shiftDeriv;
  }

  // Fixed-parameter mask, in the Minuit parameter order.
  AZUREParams tp;
  compound()->FillMnParams(tp.GetMinuitParams(), &configure());
  data()->FillMnParams(tp.GetMinuitParams());
  int nMn = tp.GetMinuitParams().Params().size();
  std::vector<bool> fixed(nMn);
  for(int i = 0; i < nMn; i++) fixed[i] = tp.GetMinuitParams().Parameter(i).IsFixed();

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed);
  int nCols = 0;
  bool ok = ComputeResidualJacobian(lc, ld, configure(), pmap, sdp, residuals, jac, nCols);

  packedToFull.resize(pmap.NumPacked());
  for(int a = 0; a < pmap.NumPacked(); a++) packedToFull[a] = pmap.PackedToFull(a);

  delete lc;
  delete ld;
  return ok;
}

double AZURECalc::RunLevenbergMarquardt(AZUREParams& params, int maxIter,
                                        BandCovariance* bandCovOut) const {
  ROOT::Minuit2::MnUserParameters& mp = params.GetMinuitParams();
  const int nFull = mp.Params().size();
  vector_r full(nFull);
  for(int i = 0; i < nFull; i++) full[i] = mp.Value(i);

  // First Jacobian fixes the column <-> parameter mapping; also collects the
  // free-parameter start values, projection limits, and Gaussian priors.
  vector_r r, J, x;
  std::vector<int> p2f;
  std::vector<double> lo, hi, pen_nom, pen_inv2;
  const int nFree = PrepareFreeParams(full, mp, r, J, p2f, x, lo, hi, pen_nom, pen_inv2);
  if(nFree < 0) return -1.0;          // analytic Jacobian unsupported -> caller falls back
  if(nFree == 0) return Chi2Value(full);

  double cost = Chi2Value(full);
  double lambda = 1.0e-3;

  std::vector<double> H(nFree * nFree), g(nFree), Hl(nFree * nFree), dx(nFree);

  for(int iter = 0; iter < maxIter; iter++) {
    if(iter > 0 && !ResidualJacobian(full, r, J, p2f)) break;
    // MIGRAD reaches WriteIterationOutput through operator(); this path
    // evaluates through the side-effect-free Chi2Value and so must drive the
    // snapshots itself.
    if(iter > 0 && iter % kIterationOutputInterval == 0) WriteIterationOutput(full);
    const int nRes = (int)r.size();

    // Gauss-Newton normal equations: H = J^T J (+ penalty diag), g = J^T r (+ penalty grad).
    std::fill(H.begin(), H.end(), 0.0);
    std::fill(g.begin(), g.end(), 0.0);
    for(int i = 0; i < nRes; i++) {
      const double* Ji = &J[(size_t)i * nFree];
      double ri = r[i];
      for(int a = 0; a < nFree; a++) {
        g[a] += Ji[a] * ri;
        double Jia = Ji[a];
        double* Ha = &H[(size_t)a * nFree];
        for(int b = a; b < nFree; b++) Ha[b] += Jia * Ji[b];
      }
    }
    for(int a = 0; a < nFree; a++) {
      for(int b = 0; b < a; b++) H[(size_t)a*nFree + b] = H[(size_t)b*nFree + a];  // symmetrize
      if(pen_inv2[a] != 0.0) {
        H[(size_t)a*nFree + a] += pen_inv2[a];
        g[a] += (x[a] - pen_nom[a]) * pen_inv2[a];
      }
    }

    // Convergence on the (projected) gradient.
    double gmax = 0.0;
    for(int a = 0; a < nFree; a++) gmax = std::max(gmax, std::fabs(g[a]));
    if(gmax < 1.0e-8) break;

    // Inner loop: damp until a step decreases the cost.
    bool accepted = false;
    for(int tries = 0; tries < 40 && !accepted; tries++) {
      // (H + lambda*diag(H)) dx = -g
      Hl = H;
      for(int a = 0; a < nFree; a++) Hl[(size_t)a*nFree + a] += lambda * H[(size_t)a*nFree + a];

      gsl_matrix* A = gsl_matrix_alloc(nFree, nFree);
      for(int a = 0; a < nFree; a++)
        for(int b = 0; b < nFree; b++) gsl_matrix_set(A, a, b, Hl[(size_t)a*nFree + b]);
      gsl_vector* bvec = gsl_vector_alloc(nFree);
      for(int a = 0; a < nFree; a++) gsl_vector_set(bvec, a, -g[a]);
      gsl_vector* sol = gsl_vector_alloc(nFree);

      gsl_set_error_handler_off();
      int status = gsl_linalg_cholesky_decomp1(A);
      if(status == 0) status = gsl_linalg_cholesky_solve(A, bvec, sol);

      bool solved = (status == 0);
      if(solved) for(int a = 0; a < nFree; a++) dx[a] = gsl_vector_get(sol, a);
      gsl_matrix_free(A); gsl_vector_free(bvec); gsl_vector_free(sol);

      if(!solved) { lambda *= 4.0; if(lambda > 1.0e14) { accepted = false; break; } continue; }

      vector_r xnew(nFree), fullNew = full;
      for(int a = 0; a < nFree; a++) {
        xnew[a] = std::min(std::max(x[a] + dx[a], lo[a]), hi[a]);  // clamp to limits
        fullNew[p2f[a]] = xnew[a];
      }
      double costNew = Chi2Value(fullNew);

      if(costNew < cost && std::isfinite(costNew)) {
        double rel = (cost - costNew) / std::max(std::fabs(cost), 1.0);
        full = fullNew; x = xnew; cost = costNew;
        lambda = std::max(lambda * 0.3, 1.0e-12);
        accepted = true;
        configure().outStream << "\r\tLM iteration: " << std::setw(4) << iter + 1
                              << "  Chi-Squared: " << cost << "        ";
        configure().outStream.flush();
        if(rel < 1.0e-9) iter = maxIter;   // converged
      } else {
        lambda *= 4.0;
        if(lambda > 1.0e14) break;
      }
    }
    if(!accepted) break;
  }
  configure().outStream << std::endl;

  // Write best-fit values back, and parameter errors from (J^T J + penalties)^{-1}.
  for(int i = 0; i < nFull; i++) mp.SetValue(i, full[i]);
  FinalizeLeastSquaresCovariance(full, p2f, pen_inv2, mp, bandCovOut);

  return cost;
}

int AZURECalc::PrepareFreeParams(const vector_r& full,
                                 const ROOT::Minuit2::MnUserParameters& mp,
                                 vector_r& r, vector_r& J, std::vector<int>& p2f,
                                 vector_r& x, std::vector<double>& lo,
                                 std::vector<double>& hi, std::vector<double>& pen_nom,
                                 std::vector<double>& pen_inv2) const {
  // First Jacobian fixes the column <-> parameter mapping.
  if(!ResidualJacobian(full, r, J, p2f)) return -1;   // unsupported -> caller falls back
  const int nFree = (int)p2f.size();

  // Free-parameter values, limits, and quadratic-penalty (prior) terms.
  x.assign(nFree, 0.0);
  lo.assign(nFree, -std::numeric_limits<double>::infinity());
  hi.assign(nFree,  std::numeric_limits<double>::infinity());
  pen_nom.assign(nFree, 0.0);
  pen_inv2.assign(nFree, 0.0);
  for(int a = 0; a < nFree; a++) {
    int f = p2f[a];
    x[a] = full[f];
    const auto& par = mp.Parameter(f);
    if(par.HasLowerLimit()) lo[a] = par.LowerLimit();
    if(par.HasUpperLimit()) hi[a] = par.UpperLimit();
    // Gaussian penalties contribute (1/sigma^2) to the diagonal of J^T J and
    // (x-nominal)/sigma^2 to the gradient.  Norm / shift / nuisance penalties
    // are identified the same way operator()/CalculateNuisanceChiSquared do.
    std::string name = par.GetName();
    if(name.find("norm") != std::string::npos) {
      for(int s = 1; s <= data()->NumSegments(); s++) {
        ESegment* seg = data()->GetSegment(s);
        if(!seg || !seg->IsVaryNorm()) continue;
        char vn[64]; snprintf(vn, sizeof(vn), "segment_%d_norm", seg->GetSegmentKey());
        if(name == vn) {
          double n0 = seg->GetNominalNorm();
          double sig = n0 / 100.0 * seg->GetNormError();
          if(sig != 0.0) { pen_nom[a] = n0; pen_inv2[a] = 1.0 / (sig * sig); }
          break;
        }
      }
    } else if(name.find("shift") != std::string::npos) {
      for(int s = 1; s <= data()->NumSegments(); s++) {
        ESegment* seg = data()->GetSegment(s);
        if(!seg || !seg->IsVaryEnergyShift()) continue;
        char vn[64]; snprintf(vn, sizeof(vn), "segment_%d_energy_shift", seg->GetSegmentKey());
        if(name == vn) {
          double sig = seg->GetEnergyShiftError();
          if(sig != 0.0) { pen_nom[a] = seg->GetNominalEnergyShift(); pen_inv2[a] = 1.0/(sig*sig); }
          break;
        }
      }
    } else if(limitsManager_ && limitsManager_->IsNuisanceParameterByIndex(a)) {
      double sig = limitsManager_->GetConvertedErrorByIndex(a);
      if(sig > 0.0) { pen_nom[a] = limitsManager_->GetConvertedNominalValueByIndex(a);
                      pen_inv2[a] = 1.0 / (sig * sig); }
    }
  }
  return nFree;
}

void AZURECalc::FinalizeLeastSquaresCovariance(const vector_r& full,
                                               const std::vector<int>& p2f,
                                               const std::vector<double>& pen_inv2,
                                               ROOT::Minuit2::MnUserParameters& mp,
                                               BandCovariance* bandCovOut) const {
  const int nFree = (int)p2f.size();
  vector_r r, J;
  std::vector<int> p2fLocal;
  if(!ResidualJacobian(full, r, J, p2fLocal)) return;

  std::vector<double> H(nFree * nFree, 0.0);
  const int nRes = (int)r.size();
  for(int i = 0; i < nRes; i++) {
    const double* Ji = &J[(size_t)i * nFree];
    for(int a = 0; a < nFree; a++) {
      double* Ha = &H[(size_t)a * nFree];
      for(int b = a; b < nFree; b++) Ha[b] += Ji[a] * Ji[b];
    }
  }
  // Pure J^T J diagonal = each parameter's data sensitivity (column-norm^2 in
  // J), captured before priors are folded into H.  Near zero => unconstrained.
  std::vector<double> jtjDiag(nFree);
  for(int a = 0; a < nFree; a++) {
    for(int b = 0; b < a; b++) H[(size_t)a*nFree + b] = H[(size_t)b*nFree + a];
    jtjDiag[a] = H[(size_t)a*nFree + a];
    if(pen_inv2[a] != 0.0) H[(size_t)a*nFree + a] += pen_inv2[a];
  }

  // The Minuit name alone ("width_1_2") does not say which level or channel is
  // meant, so pair it with the level's J^pi and energy and the channel's pair,
  // L and S.
  auto describe = [&](int a) {
    const int minuitIndex = p2f[a];
    std::string text = mp.GetName(minuitIndex);
    const std::string label = AZURELabel::Parameter(compound(), data(), minuitIndex);
    if(!label.empty()) text += "  =  " + label;
    return text;
  };

  // Name the free parameters the data barely constrain (near-zero J column).
  auto reportInsensitive = [&]() {
    double mx = 0.0;
    for(int a = 0; a < nFree; a++) mx = std::max(mx, jtjDiag[a]);
    int nInsensitive = 0;
    for(int a = 0; a < nFree; a++) {
      if(mx <= 0.0 || jtjDiag[a] <= 1.e-10 * mx) {
        if(!nInsensitive)
          configure().outStream << "  Data-insensitive parameter(s) -- fix, constrain, or add "
                                   "a prior to one of these:" << std::endl;
        configure().outStream << "    " << describe(a) << std::endl;
        nInsensitive++;
      }
    }
    if(!nInsensitive)
      configure().outStream << "No single parameter is insensitive; the flat direction is a "
                               "degeneracy between two or more correlated parameters." << std::endl;
  };

  // Covariance = (J^T J + priors)^{-1}.  If it is singular (unconstrained or
  // degenerate parameter), retry with a small ridge so the flat direction gets
  // a large finite variance instead of losing the errors and band; warn.
  double maxDiag = 0.0;
  for(int a = 0; a < nFree; a++) maxDiag = std::max(maxDiag, H[(size_t)a*nFree + a]);
  const double ridges[4] = {0.0, maxDiag*1.e-9, maxDiag*1.e-6, maxDiag*1.e-3};
  gsl_matrix* A = gsl_matrix_alloc(nFree, nFree);
  gsl_set_error_handler_off();
  int used = -1;
  for(int t = 0; t < 4; t++) {
    for(int a = 0; a < nFree; a++)
      for(int b = 0; b < nFree; b++)
        gsl_matrix_set(A, a, b, H[(size_t)a*nFree + b] + (a==b ? ridges[t] : 0.0));
    if(gsl_linalg_cholesky_decomp1(A) == 0 && gsl_linalg_cholesky_invert(A) == 0) { used = t; break; }
  }
  if(used < 0) {
    configure().outStream << "Warning: the fit curvature matrix (J^T J) is singular; parameter errors "
                             "and the uncertainty band are unavailable." << std::endl;
    reportInsensitive();
  } else {
    if(used > 0) {
      configure().outStream << "Warning: the fit curvature matrix was near-singular; regularized it "
                               "with a small ridge to obtain the covariance. Uncertainties along the "
                               "poorly-constrained direction(s) are large and approximate." << std::endl;
      reportInsensitive();
      // Strongly correlated pairs (|rho| near 1) are the degenerate combinations.
      struct CorrPair { double c; int a, b; };
      std::vector<CorrPair> pairs;
      for(int a = 0; a < nFree; a++) {
        double va = gsl_matrix_get(A, a, a);
        if(va <= 0.0) continue;
        for(int b = a+1; b < nFree; b++) {
          double vb = gsl_matrix_get(A, b, b);
          if(vb <= 0.0) continue;
          double c = gsl_matrix_get(A, a, b) / std::sqrt(va*vb);
          if(std::fabs(c) > 0.95) pairs.push_back({c, a, b});
        }
      }
      std::sort(pairs.begin(), pairs.end(),
                [](const CorrPair& x, const CorrPair& y){ return std::fabs(x.c) > std::fabs(y.c); });
      for(size_t k = 0; k < pairs.size() && k < 5; k++)
        configure().outStream << "  Correlated (rho=" << pairs[k].c << "):" << std::endl
                              << "    " << describe(pairs[k].a) << std::endl
                              << "    " << describe(pairs[k].b) << std::endl;
    }
    for(int a = 0; a < nFree; a++) {
      double v = gsl_matrix_get(A, a, a);
      if(v > 0.0) mp.SetError(p2f[a], std::sqrt(v));
    }
    // Name free parameters whose error dwarfs their value (they dominate the band).
    {
      int n = 0;
      for(int a = 0; a < nFree; a++) {
        double v = gsl_matrix_get(A, a, a);
        if(v <= 0.0) continue;
        double val = mp.Value(p2f[a]);
        if(std::sqrt(v) > 10.0 * std::max(std::fabs(val), 1.e-30)) {
          if(!n)
            configure().outStream << "Note: weakly-determined parameter(s) (error > 1000% of value) "
                                     "dominate the uncertainty band. Constrain or fix these to "
                                     "tighten it:" << std::endl;
          configure().outStream << "    " << describe(a) << std::endl;
          n++;
        }
      }
    }
    // Export the full covariance for cross-section bands.  Columns are the
    // packed free parameters, in the same order as p2f, so we tag them with
    // their parameter identities via a matching index map.
    if(bandCovOut) {
      AZUREParams tp;
      compound()->FillMnParams(tp.GetMinuitParams(), &configure());
      data()->FillMnParams(tp.GetMinuitParams());
      int nMn = tp.GetMinuitParams().Params().size();
      std::vector<bool> fixed(nMn);
      for(int i = 0; i < nMn; i++) fixed[i] = tp.GetMinuitParams().Parameter(i).IsFixed();
      ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), fixed);
      if(pmap.NumPacked() == nFree) {
        // Keep only the R-matrix sub-block: the band is insensitive to norms
        // and energy shifts, so they are dropped from the saved covariance.
        const std::vector<int> rc = RMatrixPackedColumns(pmap);
        const int m = (int)rc.size();
        bandCovOut->cols.resize(m);
        for(int a = 0; a < m; a++) bandCovOut->cols[a] = pmap.Desc(pmap.PackedToFull(rc[a]));
        bandCovOut->M.assign(m, std::vector<double>(m, 0.0));
        for(int a = 0; a < m; a++)
          for(int b = 0; b < m; b++) bandCovOut->M[a][b] = gsl_matrix_get(A, rc[a], rc[b]);
      }
    }
  }
  gsl_matrix_free(A);
}

bool AZURECalc::ResidualsOnly(const vector_r& full, vector_r& residuals) const {
  // Forward-only counterpart of ResidualJacobian (no adjoint): same row order
  // and residual r = (model - data*n)/(cmErr*n) as ComputeResidualJacobian.
  CNuc* lc = compound()->Clone();
  EData* ld = data()->Clone();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure());

  residuals.clear();
  for(int i = 1; i <= ld->NumSegments(); i++) {
    ESegment* segment = ld->GetSegment(i);
    if(!segment) continue;
    const double norm = segment->GetNorm();
    for(int pid = 0; pid < segment->NumPoints(); pid++) {
      EPoint* pt = segment->GetPoint(pid + 1);
      if(!pt) continue;
      double th = segment->CalculateTheoreticalCrossSection(pid, lc, configure(), ld);
      double denom = pt->GetCMCrossSectionError() * norm;
      residuals.push_back(denom != 0.0 ? (th - pt->GetCMCrossSection() * norm) / denom : 0.0);
    }
  }
  delete lc;
  delete ld;
  return true;
}

namespace {

// State shared with the GSL least-squares callbacks; owned by RunGSLNonlinear.
struct GSLNlsData {
  const AZURECalc* self;
  int nFree;
  int nData;                          // number of data residual rows
  const std::vector<int>* p2f;        // packed free param -> full-vector index
  const std::vector<double>* lo;
  const std::vector<double>* hi;
  const std::vector<double>* pen_nom;
  const std::vector<double>* pen_inv2;
  const std::vector<int>* penRows;    // free-param indices that carry a prior
  vector_r baseFull;                  // full vector holding the fixed params
};

// Project x into [lo,hi] and scatter it into the full parameter vector.
void gslExpandFull(const GSLNlsData* d, const gsl_vector* x, vector_r& full) {
  full = d->baseFull;
  for(int a = 0; a < d->nFree; a++) {
    double v = gsl_vector_get(x, a);
    v = std::min(std::max(v, (*d->lo)[a]), (*d->hi)[a]);
    full[(*d->p2f)[a]] = v;
  }
}

// Residual vector: data rows plus Gaussian-prior rows sqrt(1/sigma^2)*(x-nom).
int gslResidualF(const gsl_vector* x, void* params, gsl_vector* f) {
  const GSLNlsData* d = static_cast<const GSLNlsData*>(params);
  vector_r full;
  gslExpandFull(d, x, full);
  vector_r res;
  if(!d->self->ResidualsOnly(full, res) || (int)res.size() != d->nData) return GSL_EBADFUNC;
  for(int i = 0; i < d->nData; i++) gsl_vector_set(f, i, res[i]);
  for(size_t k = 0; k < d->penRows->size(); k++) {
    int a = (*d->penRows)[k];
    double s = std::sqrt((*d->pen_inv2)[a]);
    gsl_vector_set(f, d->nData + (int)k, s * (full[(*d->p2f)[a]] - (*d->pen_nom)[a]));
  }
  return GSL_SUCCESS;
}

// Residual Jacobian: analytic data-row block plus the constant prior rows.
int gslResidualDf(const gsl_vector* x, void* params, gsl_matrix* Jm) {
  const GSLNlsData* d = static_cast<const GSLNlsData*>(params);
  vector_r full;
  gslExpandFull(d, x, full);
  vector_r r, Jflat;
  std::vector<int> p2f;
  if(!d->self->ResidualJacobian(full, r, Jflat, p2f)) return GSL_EBADFUNC;
  if((int)p2f.size() != d->nFree || (int)r.size() != d->nData) return GSL_EBADFUNC;
  const int nFree = d->nFree;
  for(int i = 0; i < d->nData; i++)
    for(int a = 0; a < nFree; a++)
      gsl_matrix_set(Jm, i, a, Jflat[(size_t)i * nFree + a]);
  for(size_t k = 0; k < d->penRows->size(); k++) {
    int row = d->nData + (int)k;
    for(int a = 0; a < nFree; a++) gsl_matrix_set(Jm, row, a, 0.0);
    int a = (*d->penRows)[k];
    gsl_matrix_set(Jm, row, a, std::sqrt((*d->pen_inv2)[a]));
  }
  return GSL_SUCCESS;
}

// Per-iteration progress line (chi^2 = ||residual||^2, priors included).
void gslProgress(const size_t iter, void* params,
                 const gsl_multifit_nlinear_workspace* w) {
  const GSLNlsData* d = static_cast<const GSLNlsData*>(params);
  gsl_vector* f = gsl_multifit_nlinear_residual(w);
  double chi2 = 0.0;
  gsl_blas_ddot(f, f, &chi2);
  d->self->configure().outStream << "\r\tGSL-LM iteration: " << std::setw(4) << iter
                                 << "  Chi-Squared: " << chi2 << "        ";
  d->self->configure().outStream.flush();

  // Intermediate output, so a long GSL fit also leaves usable files behind
  // before it converges rather than only at the end.
  if(iter > 0 && iter % AZURECalc::kIterationOutputInterval == 0) {
    vector_r full = d->baseFull;
    const gsl_vector* x = gsl_multifit_nlinear_position(w);
    for(int a = 0; a < d->nFree; a++) full[(*d->p2f)[a]] = gsl_vector_get(x, a);
    d->self->WriteIterationOutput(full);
  }
}

}  // namespace

double AZURECalc::RunGSLNonlinear(AZUREParams& params, int maxIter,
                                  BandCovariance* bandCovOut) const {
  ROOT::Minuit2::MnUserParameters& mp = params.GetMinuitParams();
  const int nFull = mp.Params().size();
  vector_r full(nFull);
  for(int i = 0; i < nFull; i++) full[i] = mp.Value(i);

  vector_r r0, J0, x;
  std::vector<int> p2f;
  std::vector<double> lo, hi, pen_nom, pen_inv2;
  const int nFree = PrepareFreeParams(full, mp, r0, J0, p2f, x, lo, hi, pen_nom, pen_inv2);
  if(nFree < 0) return -1.0;              // analytic Jacobian unsupported -> caller falls back
  if(nFree == 0) return Chi2Value(full);
  const int nData = (int)r0.size();

  // Gaussian priors on norms / shifts / nuisance params become extra residual rows.
  std::vector<int> penRows;
  for(int a = 0; a < nFree; a++) if(pen_inv2[a] > 0.0) penRows.push_back(a);
  const int nRes = nData + (int)penRows.size();

  GSLNlsData d;
  d.self = this;  d.nFree = nFree;  d.nData = nData;
  d.p2f = &p2f;  d.lo = &lo;  d.hi = &hi;
  d.pen_nom = &pen_nom;  d.pen_inv2 = &pen_inv2;  d.penRows = &penRows;
  d.baseFull = full;

  gsl_multifit_nlinear_fdf fdf;
  fdf.f = gslResidualF;
  fdf.df = gslResidualDf;
  fdf.fvv = nullptr;            // second directional derivative via finite differences of f
  fdf.n = nRes;
  fdf.p = nFree;
  fdf.params = &d;

  // Geodesic-accelerated LM in a trust region; the default solver factorizes J
  // directly (QR) rather than forming J^T J.
  gsl_multifit_nlinear_parameters gp = gsl_multifit_nlinear_default_parameters();
  gp.trs = gsl_multifit_nlinear_trs_lmaccel;

  const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace* w = gsl_multifit_nlinear_alloc(T, &gp, nRes, nFree);
  if(!w) return -1.0;

  gsl_vector* x0 = gsl_vector_alloc(nFree);
  for(int a = 0; a < nFree; a++) gsl_vector_set(x0, a, x[a]);

  gsl_set_error_handler_off();
  gsl_multifit_nlinear_init(x0, &fdf, w);

  const double xtol = 1.0e-9, gtol = 1.0e-9, ftol = 1.0e-9;
  int info = 0;
  int status = gsl_multifit_nlinear_driver(maxIter, xtol, gtol, ftol,
                                           gslProgress, &d, &info, w);
  configure().outStream << std::endl;

  std::string reason = "reached the maximum number of iterations";
  if(status == GSL_SUCCESS) reason = (info == 1) ? "converged (small step)" : "converged (small gradient)";
  else if(status != GSL_EMAXITER) reason = std::string("stopped: ") + gsl_strerror(status);
  configure().outStream << "GSL trust-region solver " << reason << " after "
                        << gsl_multifit_nlinear_niter(w) << " iterations." << std::endl;

  // Read back the (projected) best-fit parameters.
  gsl_vector* xf = gsl_multifit_nlinear_position(w);
  for(int a = 0; a < nFree; a++) {
    double v = gsl_vector_get(xf, a);
    v = std::min(std::max(v, lo[a]), hi[a]);
    full[p2f[a]] = v;
  }
  for(int i = 0; i < nFull; i++) mp.SetValue(i, full[i]);

  gsl_vector_free(x0);
  gsl_multifit_nlinear_free(w);

  double cost = Chi2Value(full);

  // Errors and band covariance from (J^T J + priors)^{-1}, as the LM path does.
  FinalizeLeastSquaresCovariance(full, p2f, pen_inv2, mp, bandCovOut);

  return cost;
}

void AZURECalc::AddNuisanceGradient(const vector_r& p, std::vector<double>& grad) const {
  if(!limitsManager_) return;

  AZUREParams tempParams;
  compound()->FillMnParams(tempParams.GetMinuitParams(), &configure());
  data()->FillMnParams(tempParams.GetMinuitParams());

  std::vector<int> nonFixedToActualIndex;
  for(int i = 0; i < tempParams.GetMinuitParams().Params().size(); i++) {
    if(!tempParams.GetMinuitParams().Parameter(i).IsFixed() ||
       tempParams.GetMinuitParams().Parameter(i).GetName().find("segment") != std::string::npos) {
      nonFixedToActualIndex.push_back(i);
    }
  }

  for(int nf = 0; nf < (int)nonFixedToActualIndex.size() && nf < (int)p.size(); nf++) {
    int actualIndex = nonFixedToActualIndex[nf];
    std::string paramName = tempParams.GetMinuitParams().Parameter(actualIndex).GetName();
    if(paramName.find("norm") != std::string::npos || paramName.find("shift") != std::string::npos) continue;
    if(!limitsManager_->IsNuisanceParameterByIndex(nf)) continue;
    double nominalValue = limitsManager_->GetConvertedNominalValueByIndex(nf);
    double paramError = limitsManager_->GetConvertedErrorByIndex(nf);
    if(paramError > 0.0) {
      // d/d p[nf] of ((p[nf]-nominal)/error)^2.
      grad[nf] += 2.0 * (p[nf] - nominalValue) / (paramError * paramError);
    }
  }
}

double AZURECalc::CalculateNuisanceChiSquared(const vector_r& p) const {
  double nuisanceChiSquared = 0.0;
  
  // Create temporary AZUREParams to get parameter names
  AZUREParams tempParams;
  compound()->FillMnParams(tempParams.GetMinuitParams(), &configure());
  data()->FillMnParams(tempParams.GetMinuitParams());
  
  // Build mapping from non-fixed parameter index to actual parameter index
  std::vector<int> nonFixedToActualIndex;
  for(int i = 0; i < tempParams.GetMinuitParams().Params().size(); i++) {
    if(!tempParams.GetMinuitParams().Parameter(i).IsFixed() || tempParams.GetMinuitParams().Parameter(i).GetName().find("segment") != std::string::npos) {
      nonFixedToActualIndex.push_back(i);
    }
  }
  
  // Check each non-fixed parameter to see if it's marked as nuisance
  for(int nonFixedIndex = 0; nonFixedIndex < nonFixedToActualIndex.size() && nonFixedIndex < p.size(); nonFixedIndex++) {
    int actualIndex = nonFixedToActualIndex[nonFixedIndex];
    std::string paramName = tempParams.GetMinuitParams().Parameter(actualIndex).GetName();

    // If norm or shift in param name, skip
    if(paramName.find("norm") != std::string::npos || paramName.find("shift") != std::string::npos) {
      continue;
    }
    
    // First check if this parameter is marked as nuisance (fast check)
    if(!limitsManager_->IsNuisanceParameterByIndex(nonFixedIndex)) {
      continue; // Skip if not a nuisance parameter
    }
    
    // Only do expensive conversions if parameter is marked as nuisance
    double nominalValue = limitsManager_->GetConvertedNominalValueByIndex(nonFixedIndex);
    double paramError = limitsManager_->GetConvertedErrorByIndex(nonFixedIndex);
    
    // If we got valid values (non-zero error means this parameter has valid nuisance settings)
    if(paramError > 0.0) {
      double paramValue = p[nonFixedIndex];      
      double deviation = (paramValue - nominalValue) / paramError;
      nuisanceChiSquared += deviation * deviation;
    }
  }
  
  return nuisanceChiSquared;
}

/*!
 * Initialize object pools with pre-allocated CNuc and EData objects
 */
void AZURECalc::InitializePools() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);
  if (pools_initialized_) return;
  
  // Calculate pool size based on OpenMP threads (fixes interaction with OpenMP)
  int pool_size = 4; // default minimum
#ifdef _OPENMP
  pool_size = std::max(4, omp_get_max_threads());
#else
  pool_size = std::max(4, static_cast<int>(std::thread::hardware_concurrency()));
#endif
  
  // Pre-allocate CNuc objects by cloning once
  for (int i = 0; i < pool_size; ++i) {
    cnuc_pool_.push(std::unique_ptr<CNuc>(compound()->Clone()));
  }
  
  // Pre-allocate EData objects by cloning once
  for (int i = 0; i < pool_size; ++i) {
    edata_pool_.push(std::unique_ptr<EData>(data()->Clone()));
  }
  
  pools_initialized_ = true;
}

/*!
 * Get a CNuc object from the pool, creating new if pool is empty
 */
CNuc* AZURECalc::GetPooledCNuc() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);
  
  if (!cnuc_pool_.empty()) {
    auto obj = std::move(cnuc_pool_.top());
    cnuc_pool_.pop();
    return obj.release();
  }
  
  // Fallback: create new if pool is empty (shouldn't happen often)
  return compound_->Clone();
}

/*!
 * Get an EData object from the pool, creating new if pool is empty
 */
EData* AZURECalc::GetPooledEData() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);
  
  if (!edata_pool_.empty()) {
    auto obj = std::move(edata_pool_.top());
    edata_pool_.pop();
    return obj.release();
  }
  
  // Fallback: create new if pool is empty (shouldn't happen often)
  return data_->Clone();
}

/*!
 * Return a CNuc object to the pool for reuse
 */
void AZURECalc::ReturnPooledCNuc(CNuc* obj) const {
  if (!obj) return;
  
  std::lock_guard<std::mutex> lock(pool_mutex_);
  cnuc_pool_.push(std::unique_ptr<CNuc>(obj));
}

/*!
 * Return an EData object to the pool for reuse  
 */
void AZURECalc::ReturnPooledEData(EData* obj) const {
  if (!obj) return;
  
  std::lock_guard<std::mutex> lock(pool_mutex_);
  edata_pool_.push(std::unique_ptr<EData>(obj));
}

/*!
 * Write parameters to file
 */
void AZURECalc::WriteParameters(AZUREParams& params, const Config& configure) const {
  char filename[256];
  snprintf(filename, sizeof(filename),"%sparam.fit",configure.outputdir.c_str());
  std::ofstream out;
  out.open(filename);
  if(out) {
    out.precision(7);
    for(int i=0;i<params.GetMinuitParams().Params().size();i++) {
      out << std::setw(20) << params.GetMinuitParams().GetName(i)
	  << std::scientific << std::setw(20) <<  params.GetMinuitParams().Value(i)
	  << std::scientific << std::setw(20) <<  params.GetMinuitParams().Error(i) << std::endl;
    }
    out.flush();
    out.close();
  } else configure.outStream << "Could not save param.fit file." << std::endl;
}