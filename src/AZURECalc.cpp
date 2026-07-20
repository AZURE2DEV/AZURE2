#include "AZURECalc.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "ESegment.h"
#include "EPoint.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "AZUREGrad.h"
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

#ifdef _OPENMP
#include <omp.h>
#endif

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

    if(thisIteration%100==0) {
      AZUREParams params;
      localCompound->FillMnParams(params.GetMinuitParams(), &configure());
      localData->FillMnParams(params.GetMinuitParams());
      WriteParameters(params,configure());
      localData->WriteOutputFiles(configure(),isFit);
      localCompound->TransformOut(configure());
      localCompound->PrintTransformParams(configure());
    }
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

double AZURECalc::Chi2Value(const vector_r& p, bool thmOnly) const {
  // Side-effect-free chi-squared (mirrors the chi-squared math of operator()).
  // With thmOnly=true, only THM (HOES) segments are computed and summed; this is
  // used when invoking Gradient() to finite-difference the THM contribution (in this
  // way the computing cost scales with the number of THM points alone, rather than the whole dataset).
  // The nuisance penalty is excluded there because it is per-parameter and it is analitically treated by Gradient().
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
    if(thmOnly && !segment->IsTHM()) continue;
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
  if(!thmOnly && limitsManager_) chiSquared += CalculateNuisanceChiSquared(p);

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
      // THM segments never enter the analytic accumulation (their normData
      // stays 0); their norm gradient -- data term AND penalty -- comes from
      // the THM-only finite differences below.
      if(seg->IsTHM()) continue;
      int idx = pmap.NormIndex(s);
      if(idx < 0 || idx >= (int)grad.size()) continue;
      double g = normData[s];
      double n0 = seg->GetNominalNorm();
      double nerr = n0 / 100.0 * seg->GetNormError();
      if(nerr != 0.0) g += 2.0 * (seg->GetNorm() - n0) / (nerr * nerr);
      grad[idx] = g;
    }
  }

  delete lc;
  delete ld;

  // Fixed-parameter mask: Minuit ignores the gradient of fixed parameters, so
  // do not waste finite differences on them (large fits fix most energy shifts).
  AZUREParams fp;
  compound()->FillMnParams(fp.GetMinuitParams(), &configure());
  data()->FillMnParams(fp.GetMinuitParams());
  const int nMn = fp.GetMinuitParams().Params().size();

  // THM (HOES) segments are excluded from the analytic adjoint (it
  // differentiates the T-matrix observable, not HOES); their energy/gamma/norm
  // contribution is added below by central differences of the THM-only chi^2.
  bool haveTHM = false;
  for(int s = 1; s <= data()->NumSegments() && !haveTHM; s++) {
    ESegment* seg = data()->GetSegment(s);
    if(seg && seg->IsTHM()) haveTHM = true;
  }

  // --- Finite differences only for what is left: non-fixed energy shifts, the
  //     THM-segment contribution to energies/gammas/norms, and (if the
  //     analytic path bailed) the full energy/gamma and norm blocks too. ---
  for(int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++) {
    if(idx < nMn && fp.GetMinuitParams().Parameter(idx).IsFixed()) continue;
    ParamKind kind = pmap.Desc(idx).kind;
    if(eg && (kind == ParamKind::LevelEnergy || kind == ParamKind::Gamma ||
              kind == ParamKind::Norm)) {
      if(!haveTHM) continue;                   // fully analytic
      if(kind == ParamKind::Norm) {
        // A norm belongs to exactly one segment: non-THM norms are fully
        // analytic (handled above); THM norms are fully finite-difference.
        ESegment* seg = data()->GetSegment(pmap.Desc(idx).segment);
        if(!seg || !seg->IsTHM()) continue;
      }
      double x0 = p[idx];
      double h = 1.0e-6 * (std::fabs(x0) + 1.0);
      vector_r pp = p; pp[idx] = x0 + h;
      vector_r pm = p; pm[idx] = x0 - h;
      // += : energies/gammas already carry the analytic non-THM part (a THM
      // norm's grad entry is still 0). Penalty terms inside the THM-only chi^2
      // that do not depend on p[idx] cancel in the central difference.
      grad[idx] += (Chi2Value(pp, true) - Chi2Value(pm, true)) / (2.0 * h);
      continue;
    }
    double x0 = p[idx];
    double h = 1.0e-6 * (std::fabs(x0) + 1.0);
    vector_r pp = p; pp[idx] = x0 + h;
    vector_r pm = p; pm[idx] = x0 - h;
    grad[idx] = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * h);
  }

  // Env-gated self-check (first call only): compare the assembled gradient
  // against full central finite differences of the chi-squared. Slow;
  // diagnostics for the analytic/THM-hybrid path only.
  if(std::getenv("AZURE_GRAD_CHECK")) {
    static bool checked = false;
    if(!checked) {
      checked = true;
      double worst = 0.0; int worstIdx = -1; double worstGrad = 0.0, worstFd = 0.0;
      for(int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++) {
        if(idx < nMn && fp.GetMinuitParams().Parameter(idx).IsFixed()) continue;
        double x0 = p[idx];
        double h = 1.0e-6 * (std::fabs(x0) + 1.0);
        vector_r pp = p; pp[idx] = x0 + h;
        vector_r pm = p; pm[idx] = x0 - h;
        double fd = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * h);
        double d = std::fabs(grad[idx] - fd) / (std::fabs(fd) + 1.0);
        if(d > worst) { worst = d; worstIdx = idx; worstGrad = grad[idx]; worstFd = fd; }
      }
      std::cerr << "[grad-check] worst |grad-fd|/(|fd|+1) = " << worst
                << " at param " << worstIdx << " (grad " << worstGrad
                << ", fd " << worstFd << "), analytic path "
                << (eg ? "ACTIVE" : "BAILED") << std::endl;
    }
  }

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

double AZURECalc::RunLevenbergMarquardt(AZUREParams& params, int maxIter) const {
  ROOT::Minuit2::MnUserParameters& mp = params.GetMinuitParams();
  const int nFull = mp.Params().size();
  vector_r full(nFull);
  for(int i = 0; i < nFull; i++) full[i] = mp.Value(i);

  // First Jacobian fixes the column <-> parameter mapping.
  vector_r r, J;
  std::vector<int> p2f;
  if(!ResidualJacobian(full, r, J, p2f)) return -1.0;   // unsupported -> caller falls back
  const int nFree = (int)p2f.size();
  if(nFree == 0) return Chi2Value(full);

  // Free-parameter values, limits, and quadratic-penalty (prior) terms.
  vector_r x(nFree);
  std::vector<double> lo(nFree, -std::numeric_limits<double>::infinity());
  std::vector<double> hi(nFree,  std::numeric_limits<double>::infinity());
  std::vector<double> pen_nom(nFree, 0.0), pen_inv2(nFree, 0.0);
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
        char vn[64]; sprintf(vn, "segment_%d_norm", seg->GetSegmentKey());
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
        char vn[64]; sprintf(vn, "segment_%d_energy_shift", seg->GetSegmentKey());
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

  double cost = Chi2Value(full);
  double lambda = 1.0e-3;

  std::vector<double> H(nFree * nFree), g(nFree), Hl(nFree * nFree), dx(nFree);

  for(int iter = 0; iter < maxIter; iter++) {
    if(iter > 0 && !ResidualJacobian(full, r, J, p2f)) break;
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

  if(ResidualJacobian(full, r, J, p2f)) {
    const int nRes = (int)r.size();
    std::fill(H.begin(), H.end(), 0.0);
    for(int i = 0; i < nRes; i++) {
      const double* Ji = &J[(size_t)i * nFree];
      for(int a = 0; a < nFree; a++) {
        double* Ha = &H[(size_t)a * nFree];
        for(int b = a; b < nFree; b++) Ha[b] += Ji[a] * Ji[b];
      }
    }
    for(int a = 0; a < nFree; a++) {
      for(int b = 0; b < a; b++) H[(size_t)a*nFree + b] = H[(size_t)b*nFree + a];
      if(pen_inv2[a] != 0.0) H[(size_t)a*nFree + a] += pen_inv2[a];
    }
    gsl_matrix* A = gsl_matrix_alloc(nFree, nFree);
    for(int a = 0; a < nFree; a++)
      for(int b = 0; b < nFree; b++) gsl_matrix_set(A, a, b, H[(size_t)a*nFree + b]);
    gsl_set_error_handler_off();
    if(gsl_linalg_cholesky_decomp1(A) == 0 && gsl_linalg_cholesky_invert(A) == 0) {
      for(int a = 0; a < nFree; a++) {
        double v = gsl_matrix_get(A, a, a);
        if(v > 0.0) mp.SetError(p2f[a], std::sqrt(v));
      }
    }
    gsl_matrix_free(A);
  }

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
  sprintf(filename,"%sparam.fit",configure.outputdir.c_str());
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