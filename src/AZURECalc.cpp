#include "AZURECalc.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <string>
#include <set>
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

// setenv/unsetenv are POSIX and missing from the Windows (MinGW) runtime.
// There, _putenv("NAME=value") sets a variable for getenv and "NAME=" removes it.
namespace {
void setEnvVar(const char *name, const char *value) {
#ifdef _WIN32
  static thread_local std::string entry;
  entry = std::string(name) + "=" + value;
  _putenv(entry.c_str());
#else
  setenv(name, value, 1);
#endif
}
void unsetEnvVar(const char *name) {
#ifdef _WIN32
  static thread_local std::string entry;
  entry = std::string(name) + "=";
  _putenv(entry.c_str());
#else
  unsetenv(name);
#endif
}
}  // namespace
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
double SegmentLocalChi2(ESegment *segment, CNuc *lc, const Config &config,
                        EData *ld) {
  double chi = 0.0;
  const double norm = segment->GetNorm();
  for (int pid = 0; pid < segment->NumPoints(); pid++) {
    double th = segment->CalculateTheoreticalCrossSection(pid, lc, config, ld);
    EPoint *pt = segment->GetPoint(pid + 1);
    if (!pt) continue;
    pt->SetFitCrossSection(th);
    double r = th - pt->GetCMCrossSection() * norm;
    double err = pt->GetCMCrossSectionError() * norm;
    if (err != 0.0) chi += (r * r) / (err * err);
  }
  if (segment->IsVaryEnergyShift()) {
    double sh = segment->GetEnergyShift();
    double shn = segment->GetNominalEnergyShift();
    double she = segment->GetEnergyShiftError();
    if (she != 0.0) chi += pow((sh - shn) / she, 2.0);
  }
  if (segment->IsVaryEnergyShiftSqrt()) {
    double sq = segment->GetEnergyShiftSqrt();
    double sqn = segment->GetNominalEnergyShiftSqrt();
    double sqe = segment->GetEnergyShiftSqrtError();
    if (sqe != 0.0) chi += pow((sq - sqn) / sqe, 2.0);
  }
  return chi;
}

}  // namespace

double AZURECalc::operator()(const vector_r &p) const {
  int thisIteration = data()->Iterations();
  data()->Iterate();
  bool isFit = data()->IsFit();

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  if (isFit) {
    // New multithreading with object pools
    if (!pools_initialized_) {
      InitializePools();
    }

    // Get objects from pool
    localCompound = GetPooledCNuc();
    localData = GetPooledEData();

    // Old multithreading
    // localCompound = compound()->Clone();
    // localData = data()->Clone();
  } else {
    localCompound = compound();
    localData = data();
  }

  // Fill Compound Nucleus From Minuit Parameters
  localCompound->FillCompoundFromParams(p);
  localData->FillNormsFromParams(p);
  localData->FillEnergyShiftsFromParams(p, localData, localCompound, &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  // Calculate all points in one parallel pass when not already inside a parallel
  // region (MIGRAD with the analytic gradient, line searches, plain calculations).
  const bool pre = PrecalculatePoints(localCompound, localData);

  // Process segments with components - use new integrated calculation method
  double chiSquared = 0.0;
  for (int i = 1; i <= localData->NumSegments(); i++) {
    ESegment *segment = localData->GetSegment(i);
    if (segment) {
      // Recalculate points using the new combined calculation method
      for (int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        double theoreticalValue = segment->CalculateTheoreticalCrossSection(pointIdx, localCompound, configure(), localData, pre);
        EPoint *point = segment->GetPoint(pointIdx + 1);
        if (point) {
          point->SetFitCrossSection(theoreticalValue);
        }
      }

      // Recalculate chi-squared for this segment with components
      double segmentChiSquared = 0.0;
      for (int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        EPoint *point = segment->GetPoint(pointIdx + 1);
        if (point) {
          double residual = point->GetFitCrossSection() - point->GetCMCrossSection() * segment->GetNorm();
          double error = point->GetCMCrossSectionError() * segment->GetNorm();
          if (error != 0.0) {
            segmentChiSquared += (residual * residual) / (error * error);
          }
        }
      }

      // Add normalization chi-squared contribution
      double dataNorm = segment->GetNorm();
      double dataNormNominal = segment->GetNominalNorm();
      double dataNormError = dataNormNominal / 100. * segment->GetNormError();
      if (dataNormError != 0.) {
        chiSquared += pow((dataNorm - dataNormNominal) / dataNormError, 2.0);
      }

      // Add energy shift chi-squared contribution
      if (segment->IsVaryEnergyShift()) {
        double energyShift = segment->GetEnergyShift();
        double energyShiftNominal = segment->GetNominalEnergyShift();
        double energyShiftError = segment->GetEnergyShiftError();
        if (energyShiftError != 0.) {
          chiSquared += pow((energyShift - energyShiftNominal) / energyShiftError, 2.0);
        }
      }
      // ... and the same for the sqrt(E) energy-shift coefficient
      if (segment->IsVaryEnergyShiftSqrt()) {
        double sq = segment->GetEnergyShiftSqrt();
        double sqn = segment->GetNominalEnergyShiftSqrt();
        double sqe = segment->GetEnergyShiftSqrtError();
        if (sqe != 0.) {
          chiSquared += pow((sq - sqn) / sqe, 2.0);
        }
      }

      segment->SetSegmentChiSquared(segmentChiSquared);
      chiSquared += segmentChiSquared;
    }
  }

  // Add nuisance parameter chi-squared contributions
  if (limitsManager_) {
    chiSquared += CalculateNuisanceChiSquared(p);
  }
  // Park formalism: keep every level's overlap J positive.
  if (configure().paramMask & Config::USE_PARK_FORMALISM) {
    double parkPenalty = localCompound->ParkNormPenalty();
    chiSquared += parkPenalty;
    localData->SetParkPenalty(parkPenalty);
  }

  if (!localData->IsErrorAnalysis() && thisIteration != 0) {
    if (thisIteration % 10 == 0) configure().outStream
        << "\r\tIteration: " << std::setw(6) << thisIteration
        << " Chi-Squared: " << chiSquared;
    configure().outStream.flush();

    if (thisIteration % kOutputInterval == 0) WriteIterationOutput(p);
  }
  if (isFit) {
    // New multithreading with object pools
    ReturnPooledCNuc(localCompound);
    ReturnPooledEData(localData);

    // Old multithreading
    // delete localCompound;
    // delete localData;
  }

  // Make a check if chiSquared is NaN
  if (std::isnan(chiSquared)) {
    // In that case return infinite since MINUIT2 can have issues with NaN values
    return std::numeric_limits<double>::infinity();
  }

  if (configure().stopFlag && isFit)
    return 0.;
  else
    return chiSquared;
}

void AZURECalc::WriteIterationOutput(const vector_r &p) const {
  // Work on private copies: this runs in the middle of a fit and must not
  // disturb the objects the minimizer is stepping.  The copies are made once
  // and reused (every snapshot refills them from p below, exactly as the pooled
  // objects are refilled in operator()), instead of being cloned and deleted at
  // every snapshot -- see output_data_ in AZURECalc.h.  The lock also keeps two
  // threads from writing the snapshot files at the same time.
  std::lock_guard<std::mutex> outputLock(output_mutex_);
  if (!output_compound_) output_compound_.reset(compound()->Clone());
  if (!output_data_) output_data_.reset(data()->Clone());
  CNuc *lc = output_compound_.get();
  EData *ld = output_data_.get();

  try {
    lc->FillCompoundFromParams(p);
    ld->FillNormsFromParams(p);
    ld->FillEnergyShiftsFromParams(p, ld, lc, &configure());
    if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure());

    for (int i = 1; i <= ld->NumSegments(); i++) {
      ESegment *segment = ld->GetSegment(i);
      if (!segment) continue;
      for (int pid = 0; pid < segment->NumPoints(); pid++) {
        EPoint *pt = segment->GetPoint(pid + 1);
        if (pt) pt->SetFitCrossSection(
            segment->CalculateTheoreticalCrossSection(pid, lc, configure(), ld));
      }
    }

    AZUREParams params;
    lc->FillMnParams(params.GetMinuitParams(), &configure());
    ld->FillMnParams(params.GetMinuitParams());
    // FillMnParams reads the levels' input values (GetE/GetGamma), not the
    // point being evaluated; copy p in so param.fit matches this snapshot.
    ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
    if (mp.Params().size() == p.size())
      for (unsigned int i = 0; i < p.size(); i++) mp.SetValue(i, p[i]);
    WriteParameters(params, configure());
    ld->WriteOutputFiles(configure(), true);
    lc->TransformOut(configure());
    lc->PrintTransformParams(configure());
  } catch (...) {
    // An intermediate snapshot is a convenience, never a reason to abort a fit.
    // A failure may leave the copies half-filled: drop them so the next
    // snapshot starts from fresh clones.
    output_compound_.reset();
    output_data_.reset();
  }
}

bool AZURECalc::PrecalculatePoints(CNuc *lc, EData *ld) const {
#ifdef _OPENMP
  if (omp_in_parallel()) return false;
#else
  return false;
#endif
  if (std::getenv("AZURE_SERIAL_CHI2")) return false;
  std::vector<EPoint *> work;
  std::set<ESegment *> seen;
  auto add = [&](ESegment *seg) {
    if (!seg || !seen.insert(seg).second) return;
    for (int p = 1; p <= seg->NumPoints(); p++) {
      EPoint *pt = seg->GetPoint(p);
      if (pt && !pt->IsMapped()) work.push_back(pt);
    }
  };
  for (int i = 1; i <= ld->NumSegments(); i++) {
    ESegment *seg = ld->GetSegment(i);
    if (!seg) continue;
    add(seg);
    for (ESegment *comp : seg->GetComponentSegments()) add(comp);
  }
  const Config &cfg = configure();
#pragma omp parallel for schedule(dynamic, 4)
  for (int k = 0; k < (int)work.size(); k++) {
    try {
      work[k]->Calculate(lc, cfg);
    } catch (...) {
      work[k]->SetFitCrossSection(0.0);  // as the per-point path, which returns 0 on failure
    }
  }
  return true;
}

double AZURECalc::Chi2On(CNuc *lc, EData *ld, const vector_r &p, std::vector<double> *segChis) const {
  lc->FillCompoundFromParams(p);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, ld, lc, &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure());
  const bool pre = PrecalculatePoints(lc, ld);

  double chiSquared = 0.0;
  if (segChis) segChis->assign(ld->NumSegments() + 1, 0.0);
  for (int i = 1; i <= ld->NumSegments(); i++) {
    ESegment *segment = ld->GetSegment(i);
    if (!segment) continue;
    double segChi = 0.0;
    for (int pid = 0; pid < segment->NumPoints(); pid++) {
      double th = segment->CalculateTheoreticalCrossSection(pid, lc, configure(), ld, pre);
      EPoint *pt = segment->GetPoint(pid + 1);
      if (pt) {
        pt->SetFitCrossSection(th);
        double r = th - pt->GetCMCrossSection() * segment->GetNorm();
        double err = pt->GetCMCrossSectionError() * segment->GetNorm();
        if (err != 0.0) segChi += (r * r) / (err * err);
      }
    }
    if (segChis) (*segChis)[i] = segChi;
    double dataNorm = segment->GetNorm();
    double nom = segment->GetNominalNorm();
    double nerr = nom / 100.0 * segment->GetNormError();
    if (nerr != 0.0) chiSquared += pow((dataNorm - nom) / nerr, 2.0);
    if (segment->IsVaryEnergyShift()) {
      double sh = segment->GetEnergyShift();
      double shn = segment->GetNominalEnergyShift();
      double she = segment->GetEnergyShiftError();
      if (she != 0.0) chiSquared += pow((sh - shn) / she, 2.0);
    }
    if (segment->IsVaryEnergyShiftSqrt()) {
      double sq = segment->GetEnergyShiftSqrt();
      double sqn = segment->GetNominalEnergyShiftSqrt();
      double sqe = segment->GetEnergyShiftSqrtError();
      if (sqe != 0.0) chiSquared += pow((sq - sqn) / sqe, 2.0);
    }
    chiSquared += segChi;
  }
  if (limitsManager_) chiSquared += CalculateNuisanceChiSquared(p);
  return chiSquared;
}

double AZURECalc::Chi2Value(const vector_r &p) const {
  // Side-effect-free chi-squared (mirrors the chi-squared math of operator()).
  // Pooled working copies keep their applied energy shifts and mapping between
  // calls; a fresh Clone() re-shifts every shifted segment and rebuilds the
  // mapping each time.
  CNuc *lc = GetPooledCNuc();
  EData *ld = GetPooledEData();
  std::vector<double> segPooled;
  const bool poolCheck = std::getenv("AZURE_POOL_CHECK") != nullptr;
  const double chi = Chi2On(lc, ld, p, poolCheck ? &segPooled : nullptr);
  ReturnPooledCNuc(lc);
  ReturnPooledEData(ld);
  if (poolCheck) {
    static std::atomic<int> pc{0};
    if (++pc <= 5) {
      CNuc *fc = compound()->Clone();
      EData *fd = data()->Clone();
      std::vector<double> segFresh;
      setEnvVar("AZURE_SERIAL_CHI2", "1");  // reference: per-point (serial) path on a fresh copy
      const double chiFresh = Chi2On(fc, fd, p, &segFresh);
      unsetEnvVar("AZURE_SERIAL_CHI2");
      delete fc;
      delete fd;
      std::cerr << "[pool-check] pooled+parallel " << chi << " fresh+serial " << chiFresh << std::endl;
      std::vector<std::pair<double, int>> d;
      for (size_t i = 1; i < segFresh.size() && i < segPooled.size(); i++)
        d.push_back({std::fabs(segPooled[i] - segFresh[i]), (int)i});
      std::sort(d.rbegin(), d.rend());
      for (int k = 0; k < 6 && k < (int)d.size() && d[k].first > 1e-6; k++)
        std::cerr << "[pool-check]   segment " << d[k].second << ": pooled " << segPooled[d[k].second] << " fresh "
                  << segFresh[d[k].second] << std::endl;
    }
  }
  return chi;
}

std::vector<double> AZURECalc::Gradient(const std::vector<double> &p) const {
  std::vector<double> grad(p.size(), 0.0);
  const auto tg0 = std::chrono::steady_clock::now();
  // With G2 enabled, one residual Jacobian gives both the gradient (2 J^T r +
  // penalties) and Minuit's G2, which it requests at the same point next.
  if (g2Enabled_ && !std::getenv("AZURE_GRAD_ADJOINT") && GradientFromJacobian(p, grad)) {
    if (std::getenv("AZURE_GRAD_TIMING")) {
      static std::atomic<int> jc{0};
      int c = ++jc;
      if (c <= 5 || c % 50 == 0)
        std::cerr << "[grad-timing] call " << c << " (Jacobian path, gradient + G2): "
                  << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tg0).count() << " ms" << std::endl;
    }
    // AZURE_GRAD_CHECK: on the first call also take the adjoint gradient (which
    // itself runs the sampled finite-difference check) and compare all components.
    static std::atomic<bool> jChecked{false};
    if (std::getenv("AZURE_GRAD_CHECK") && !jChecked.exchange(true)) {
      setEnvVar("AZURE_GRAD_ADJOINT", "1");
      std::vector<double> ga = Gradient(p);
      unsetEnvVar("AZURE_GRAD_ADJOINT");
      AZUREParams fpc;
      compound()->FillMnParams(fpc.GetMinuitParams(), &configure());
      data()->FillMnParams(fpc.GetMinuitParams());
      const int nMnc = fpc.GetMinuitParams().Params().size();
      double worst = 0.0; int worstIdx = -1, nFree = 0;
      for (size_t i = 0; i < grad.size() && i < ga.size(); i++) {
        if ((int)i < nMnc && fpc.GetMinuitParams().Parameter(i).IsFixed()) continue;  // Minuit ignores these
        nFree++;
        if (grad[i] == 0.0 && ga[i] == 0.0) continue;
        const double rel = std::fabs(grad[i] - ga[i]) / (std::fabs(grad[i]) + std::fabs(ga[i]) + 1.0);
        if (rel > worst) { worst = rel; worstIdx = (int)i; }
      }
      std::cerr << "[grad-check] Jacobian-path vs adjoint gradient over " << nFree << " free parameters: worst relative difference " << worst << " at idx "
                << worstIdx << (worstIdx >= 0 ? " (" + std::to_string(grad[worstIdx]) + " vs " + std::to_string(ga[worstIdx]) + ")" : "")
                << std::endl;
    }
    return grad;
  }
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);

  // --- Analytic energy / reduced-width block via the shared adjoint engine. ---
  // Pooled working copies keep their applied energy shifts and mapping between calls; a fresh
  // Clone() re-shifted every shifted segment and rebuilt the mapping each time (~1 s here).
  CNuc *lc = GetPooledCNuc();
  EData *ld = GetPooledEData();
  lc->FillCompoundFromParams(p);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (brune) {
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
  FitBarFn fitBarFn = [&](ESegment *seg, int i, int pid, double model) -> double {
    EPoint *pt = seg->GetPoint(pid + 1);
    if (!pt) return 0.0;
    double norm = seg->GetNorm();
    double dataval = pt->GetCMCrossSection();
    double cmErr = pt->GetCMCrossSectionError();
    double r = model - dataval * norm;
    double err = cmErr * norm;
    if (err == 0.0) return 0.0;
    if (seg->IsVaryNorm() && norm != 0.0 && i >= 1 && i < (int)normData.size()) {
      double e2 = cmErr * cmErr;
      // fitBarFn runs inside the parallel point loop of AccumulateEGammaGradient,
      // and all points of a segment share the same normData[i], so guard the
      // accumulation.
      double dNorm = -2.0 * r * dataval / (e2 * norm * norm) - 2.0 * r * r / (e2 * norm * norm * norm);
#pragma omp atomic
      normData[i] += dNorm;
    }
    return 2.0 * r / (err * err);
  };

  bool eg = AccumulateEGammaGradient(lc, ld, configure(), pmap, sdp, fitBarFn, accum);
  const auto tg1 = std::chrono::steady_clock::now();
  if (!eg && std::getenv("AZURE_GRAD_DEBUG")) {
    std::cerr << "[grad] analytic energy/gamma adjoint bailed -> full finite "
                 "differences (no speed-up). An unsupported segment/config is "
                 "present (RMC, or a not-yet-handled case)."
              << std::endl;
  }
  if (eg) {
    AddParkPenaltyGradient(lc, configure(), accum); // Park J > 0 wall
    accum.Scatter(pmap, grad);                      // energies + reduced widths
    AddNuisanceGradient(p, grad);                   // nuisance penalty
    for (int s = 1; s <= ld->NumSegments(); s++) {  // normalizations
      ESegment *seg = ld->GetSegment(s);
      if (!seg || !seg->IsVaryNorm()) continue;
      int idx = pmap.NormIndex(s);
      if (idx < 0 || idx >= (int)grad.size()) continue;
      double g = normData[s];
      double n0 = seg->GetNominalNorm();
      double nerr = n0 / 100.0 * seg->GetNormError();
      if (nerr != 0.0) g += 2.0 * (seg->GetNorm() - n0) / (nerr * nerr);
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
  for (int s = 1; s <= ld->NumSegments() && !datasetHasMapping; s++) {
    ESegment *seg = ld->GetSegment(s);
    if (!seg) continue;
    for (int pid = 0; pid < seg->NumPoints(); pid++) {
      EPoint *pt = seg->GetPoint(pid + 1);
      if (pt && pt->IsMapped()) {
        datasetHasMapping = true;
        break;
      }
    }
  }
  const bool localShiftOk = !datasetHasMapping && !HasNuisanceParameters();

  // --- Finite differences only for what is left: non-fixed energy shifts, and
  //     (if the analytic path bailed) the energy/gamma and norm blocks too. ---
  for (int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++) {
    if (idx < nMn && fp.GetMinuitParams().Parameter(idx).IsFixed()) continue;
    ParamKind kind = pmap.Desc(idx).kind;
    if (eg && (kind == ParamKind::LevelEnergy || kind == ParamKind::Gamma || kind == ParamKind::Norm)) continue;  // analytic
    double x0 = p[idx];
    double h = 1.0e-6 * (std::fabs(x0) + 1.0);

    // Fast path: a free energy shift on a self-contained segment only moves that
    // segment's points, so finite-difference just that segment's chi^2 on the
    // already-filled clone instead of re-solving the entire dataset twice.
    if ((kind == ParamKind::EnergyShift || kind == ParamKind::EnergyShiftSqrt) && localShiftOk) {
      int s = pmap.Desc(idx).segment;
      ESegment *seg = (s >= 1) ? ld->GetSegment(s) : nullptr;
      if (seg && !seg->HasComponents() && !seg->IsTotalCapture()) {
        // The constant shift and the sqrt(E) coefficient enter the same
        // E -> E + a + b*sqrt(E) mapping; difference whichever one this is.
        const bool isSqrt = (kind == ParamKind::EnergyShiftSqrt);
        auto setShift = [&](double v) {
          if (isSqrt) seg->SetEnergyShiftSqrt(v); else seg->SetEnergyShift(v);
        };
        double d0 = isSqrt ? seg->GetEnergyShiftSqrt() : seg->GetEnergyShift();
        setShift(x0 + h);
        seg->UpdatePointEnergiesWithShift(lc, &configure());
        double chiP = SegmentLocalChi2(seg, lc, configure(), ld);
        setShift(x0 - h);
        seg->UpdatePointEnergiesWithShift(lc, &configure());
        double chiM = SegmentLocalChi2(seg, lc, configure(), ld);
        setShift(d0);  // restore base state on the clone
        seg->UpdatePointEnergiesWithShift(lc, &configure());
        grad[idx] = (chiP - chiM) / (2.0 * h);

        // Optional self-check: compare against the exact full-dataset finite
        // difference for this shift.  Enable with AZURE_SHIFT_GRAD_CHECK=1 to
        // validate the fast path on real data (one gradient evaluation), then
        // disable it for the actual fit.
        if (std::getenv("AZURE_SHIFT_GRAD_CHECK")) {
          vector_r pp = p;
          pp[idx] = x0 + h;
          vector_r pm = p;
          pm[idx] = x0 - h;
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

    vector_r pp = p;
    pp[idx] = x0 + h;
    vector_r pm = p;
    pm[idx] = x0 - h;
    grad[idx] = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * h);
  }

  const auto tg2 = std::chrono::steady_clock::now();
  static std::atomic<int> gradCalls{0};
  if (std::getenv("AZURE_GRAD_CHECK") || std::getenv("AZURE_GRAD_TIMING")) {
    int c = ++gradCalls;
    if (c <= 5 || c % 50 == 0) {
      auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
      int nfd = 0;
      for (int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++) {
        if (idx < nMn && fp.GetMinuitParams().Parameter(idx).IsFixed()) continue;
        ParamKind kind = pmap.Desc(idx).kind;
        if (!(eg && (kind == ParamKind::LevelEnergy || kind == ParamKind::Gamma || kind == ParamKind::Norm))) nfd++;
      }
      std::cerr << "[grad-timing] call " << c << ": setup+adjoint " << ms(tg0, tg1) << " ms (analytic " << (eg ? "ok" : "BAILED")
                << "), finite-difference part " << ms(tg1, tg2) << " ms for " << nfd << " parameters" << std::endl;
    }
  }

  // Optional full self-check, AZURE_GRAD_CHECK=1 (or =exit to stop the run after
  // it): on the first call, compare every non-fixed gradient component with a
  // central finite difference of the side-effect-free chi-squared and list the
  // worst disagreements.  Minuit itself never checks (CheckGradient() is false).
  static std::atomic<bool> gradChecked{false};
  if (const char *gc = std::getenv("AZURE_GRAD_CHECK")) {
    if (!gradChecked.exchange(true)) {
      static const char *kindName[] = {"energy", "gamma", "norm", "shift", "shift_sqrt"};
      struct Row { double score, g, fd; int idx; ParamDesc d; };
      std::vector<Row> rows;
      int nBad[5] = {0, 0, 0, 0, 0}, nAll[5] = {0, 0, 0, 0, 0};
      const auto tc0 = std::chrono::steady_clock::now();
      const double f0 = Chi2Value(p);
      std::cerr << "[grad-check] one Chi2Value: " << std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count() << " s" << std::endl;
      std::vector<int> freeIdx;
      for (int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++)
        if (!(idx < nMn && fp.GetMinuitParams().Parameter(idx).IsFixed())) freeIdx.push_back(idx);
      const char *cn = std::getenv("AZURE_GRAD_CHECK_N");
      const int nCheck = cn ? std::atoi(cn) : 16;
      const int stride = std::max(1, (int)freeIdx.size() / std::max(1, nCheck));
      for (int q = 0; q < (int)freeIdx.size(); q += stride) {
        const int idx = freeIdx[q];
        double x0 = p[idx], h = 1.0e-6 * (std::fabs(x0) + 1.0);
        vector_r pp = p, pm = p;
        pp[idx] = x0 + h;
        pm[idx] = x0 - h;
        double fd = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * h);
        double score = std::fabs(grad[idx] - fd) / (std::fabs(fd) + std::fabs(grad[idx]) + 1.0e-6 * std::fabs(f0) + 1.0);
        ParamDesc d = pmap.Desc(idx);
        int k = (int)d.kind;
        nAll[k]++;
        if (score > 1.0e-2) nBad[k]++;
        rows.push_back({score, grad[idx], fd, idx, d});
      }
      std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b) { return a.score > b.score; });
      std::cerr << "[grad-check] chi2 " << f0 << ", analytic block " << (eg ? "used" : "BAILED (all finite differences)")
                << "; components with relative disagreement > 1e-2:";
      for (int k = 0; k < 5; k++)
        if (nAll[k]) std::cerr << " " << kindName[k] << " " << nBad[k] << "/" << nAll[k];
      std::cerr << std::endl;
      for (int r = 0; r < (int)rows.size() && r < 40; r++) {
        const Row &w = rows[r];
        std::cerr << "[grad-check] idx " << w.idx << " " << kindName[(int)w.d.kind] << " J" << w.d.jGroup << " L" << w.d.level
                  << " ch" << w.d.channel << " seg" << w.d.segment << "  analytic " << w.g << "  fd " << w.fd
                  << "  rel " << w.score << std::endl;
      }
      if (std::string(gc) == "exit") std::exit(0);
    }
  }

  ReturnPooledCNuc(lc);
  ReturnPooledEData(ld);

  return grad;
}

bool AZURECalc::GradientFromJacobian(const std::vector<double> &p, std::vector<double> &grad) const {
  if (configure().paramMask & Config::USE_PARK_FORMALISM) return false;  // Park wall penalty not in J
  vector_r r, jac;
  std::vector<int> packedToFull;
  if (!ResidualJacobian(p, r, jac, packedToFull) || r.empty()) return false;
  const size_t nCols = packedToFull.size(), nRes = r.size();
  std::vector<double> g2(p.size(), 0.0);
  grad.assign(p.size(), 0.0);
  for (size_t i = 0; i < nRes; i++)
    for (size_t a = 0; a < nCols; a++) {
      const double d = jac[i * nCols + a];
      const int f = packedToFull[a];
      if (f < 0 || f >= (int)p.size()) continue;
      grad[f] += 2.0 * r[i] * d;
      g2[f] += 2.0 * d * d;
    }
  // Penalties (x - x0)^2 / sigma^2: gradient 2 (x - x0) / sigma^2, curvature 2 / sigma^2.
  ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), std::vector<bool>());
  for (int s = 1; s <= data()->NumSegments(); s++) {
    ESegment *seg = data()->GetSegment(s);
    if (!seg) continue;
    int idx;
    if (seg->IsVaryNorm() && (idx = pmap.NormIndex(s)) >= 0 && idx < (int)p.size()) {
      const double n0 = seg->GetNominalNorm(), nerr = n0 / 100.0 * seg->GetNormError();
      if (nerr != 0.0) { grad[idx] += 2.0 * (p[idx] - n0) / (nerr * nerr); g2[idx] += 2.0 / (nerr * nerr); }
    }
    if (seg->IsVaryEnergyShift() && (idx = pmap.EnergyShiftIndex(s)) >= 0 && idx < (int)p.size()) {
      const double e = seg->GetEnergyShiftError();
      if (e != 0.0) { grad[idx] += 2.0 * (p[idx] - seg->GetNominalEnergyShift()) / (e * e); g2[idx] += 2.0 / (e * e); }
    }
    if (seg->IsVaryEnergyShiftSqrt() && (idx = pmap.EnergyShiftSqrtIndex(s)) >= 0 && idx < (int)p.size()) {
      const double e = seg->GetEnergyShiftSqrtError();
      if (e != 0.0) { grad[idx] += 2.0 * (p[idx] - seg->GetNominalEnergyShiftSqrt()) / (e * e); g2[idx] += 2.0 / (e * e); }
    }
  }
  g2CacheP_ = p;
  g2Cache_ = g2;
  return true;
}

std::vector<double> AZURECalc::G2(const std::vector<double> &p) const {
  if (!g2Cache_.empty() && g2CacheP_ == p) return g2Cache_;
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<double> g2(p.size(), 0.0);
  vector_r r, jac;
  std::vector<int> packedToFull;
  if (!ResidualJacobian(p, r, jac, packedToFull) || r.empty()) return std::vector<double>();
  const size_t nCols = packedToFull.size(), nRes = r.size();
  for (size_t i = 0; i < nRes; i++)
    for (size_t a = 0; a < nCols; a++) {
      const double d = jac[i * nCols + a];
      const int f = packedToFull[a];
      if (f >= 0 && f < (int)g2.size()) g2[f] += 2.0 * d * d;
    }
  // Penalty curvatures, (x - x0)^2 / sigma^2 -> 2 / sigma^2.
  ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), std::vector<bool>());
  for (int s = 1; s <= data()->NumSegments(); s++) {
    ESegment *seg = data()->GetSegment(s);
    if (!seg) continue;
    int idx;
    if (seg->IsVaryNorm() && (idx = pmap.NormIndex(s)) >= 0 && idx < (int)g2.size()) {
      const double nerr = seg->GetNominalNorm() / 100.0 * seg->GetNormError();
      if (nerr != 0.0) g2[idx] += 2.0 / (nerr * nerr);
    }
    if (seg->IsVaryEnergyShift() && (idx = pmap.EnergyShiftIndex(s)) >= 0 && idx < (int)g2.size() &&
        seg->GetEnergyShiftError() != 0.0)
      g2[idx] += 2.0 / (seg->GetEnergyShiftError() * seg->GetEnergyShiftError());
    if (seg->IsVaryEnergyShiftSqrt() && (idx = pmap.EnergyShiftSqrtIndex(s)) >= 0 && idx < (int)g2.size() &&
        seg->GetEnergyShiftSqrtError() != 0.0)
      g2[idx] += 2.0 / (seg->GetEnergyShiftSqrtError() * seg->GetEnergyShiftSqrtError());
  }
  static std::atomic<int> g2Calls{0};
  if (std::getenv("AZURE_GRAD_TIMING")) {
    int c = ++g2Calls;
    if (c <= 5 || c % 50 == 0)
      std::cerr << "[g2-timing] call " << c << ": " << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count()
                << " ms" << std::endl;
  }
  return g2;
}

bool AZURECalc::HasNuisanceParameters() const {
  // The CLI always hands AZURECalc a ParameterLimitsManager, so "a manager
  // exists" is not the question; whether it declares any nuisance parameter is.
  if (!limitsManager_) return false;
  AZUREParams fp;
  compound()->FillMnParams(fp.GetMinuitParams(), &configure());
  data()->FillMnParams(fp.GetMinuitParams());
  int nonFixed = 0;
  for (const auto &mp : fp.GetMinuitParams().Parameters())
    if (!mp.IsFixed()) nonFixed++;
  for (int i = 0; i < nonFixed; i++)
    if (limitsManager_->IsNuisanceParameterByIndex(i)) return true;
  return false;
}

bool AZURECalc::EnableG2(const std::vector<double> &p) {
  g2Enabled_ = false;
  if (HasNuisanceParameters()) return false;  // nuisance penalties are not in the Jacobian
  std::vector<double> g2 = G2(p);
  g2Enabled_ = !g2.empty();
  return g2Enabled_;
}

bool AZURECalc::ResidualJacobian(const vector_r &full, vector_r &residuals,
                                 vector_r &jac, std::vector<int> &packedToFull) const {
  // Pooled working copies keep their applied energy shifts and mapping between calls; a fresh
  // Clone() re-shifted every shifted segment and rebuilt the mapping each time (~1 s here).
  CNuc *lc = GetPooledCNuc();
  EData *ld = GetPooledEData();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) {
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
  for (int i = 0; i < nMn; i++) fixed[i] = tp.GetMinuitParams().Parameter(i).IsFixed();

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed);
  int nCols = 0;
  bool ok = ComputeResidualJacobian(lc, ld, configure(), pmap, sdp, residuals, jac, nCols);

  // Energy-shift columns by central differences, as AZUREAPI::
  // CalculateResidualJacobianRWA does for pyazr.  The adjoint does not cover a
  // shift (it moves the energy of every quantity of the segment), and left as
  // zero the Levenberg-Marquardt and GSL drivers flag the shifts as
  // data-insensitive and never move them.  Differencing the residual vector
  // of the same function keeps the rows aligned.
  if (ok) {
    const size_t nRes = residuals.size();
    vector_r rPlus, rMinus, jTmp;
    int nc = 0;
    for (int f = 0; f < pmap.NumFull() && f < (int)full.size(); f++) {
      if (fixed[f]) continue;
      if (pmap.Desc(f).kind != ParamKind::EnergyShift &&
          pmap.Desc(f).kind != ParamKind::EnergyShiftSqrt) continue;
      const int packed = pmap.FullToPacked(f);
      if (packed < 0 || packed >= nCols) continue;
      const double x0 = full[f];
      const double h = 1.0e-6 * (std::fabs(x0) + 1.0);
      auto residualsAt = [&](double value, vector_r &out) -> bool {
        vector_r fk = full;
        fk[f] = value;
        lc->FillCompoundFromParams(fk);
        ld->FillNormsFromParams(fk);
        ld->FillEnergyShiftsFromParams(fk, ld, lc, &configure());
        out.clear();
        jTmp.clear();
        return ComputeResidualJacobian(lc, ld, configure(), pmap, sdp, out, jTmp, nc);
      };
      const bool okp = residualsAt(x0 + h, rPlus);
      const bool okm = residualsAt(x0 - h, rMinus);
      if (okp && okm && rPlus.size() == nRes && rMinus.size() == nRes)
        for (size_t r = 0; r < nRes; r++)
          jac[r * (size_t)nCols + (size_t)packed] = (rPlus[r] - rMinus[r]) / (2.0 * h);
    }
    lc->FillCompoundFromParams(full);
    ld->FillNormsFromParams(full);
    ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());
  }

  packedToFull.resize(pmap.NumPacked());
  for (int a = 0; a < pmap.NumPacked(); a++) packedToFull[a] = pmap.PackedToFull(a);

  ReturnPooledCNuc(lc);
  ReturnPooledEData(ld);
  return ok;
}

double AZURECalc::RunLevenbergMarquardt(AZUREParams &params, int maxIter,
                                        BandCovariance *bandCovOut) const {
  ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  const int nFull = mp.Params().size();
  vector_r full(nFull);
  for (int i = 0; i < nFull; i++) full[i] = mp.Value(i);

  // First Jacobian fixes the column <-> parameter mapping; also collects the
  // free-parameter start values, projection limits, and Gaussian priors.
  vector_r r, J, x;
  std::vector<int> p2f;
  std::vector<double> lo, hi, pen_nom, pen_inv2;
  const int nFree = PrepareFreeParams(full, mp, r, J, p2f, x, lo, hi, pen_nom, pen_inv2);
  if (nFree < 0) return -1.0;  // analytic Jacobian unsupported -> caller falls back
  if (nFree == 0) return Chi2Value(full);

  double cost = Chi2Value(full);
  double lambda = 1.0e-3;

  std::vector<double> H(nFree * nFree), g(nFree), Hl(nFree * nFree), dx(nFree);

  for (int iter = 0; iter < maxIter; iter++) {
    if (iter > 0 && !ResidualJacobian(full, r, J, p2f)) break;
    // MIGRAD reaches WriteIterationOutput through operator(); this path
    // evaluates through the side-effect-free Chi2Value and so must drive the
    // snapshots itself.
    if (iter > 0 && iter % kIterationOutputInterval == 0) WriteIterationOutput(full);
    const int nRes = (int)r.size();

    // Gauss-Newton normal equations: H = J^T J (+ penalty diag), g = J^T r (+ penalty grad).
    std::fill(H.begin(), H.end(), 0.0);
    std::fill(g.begin(), g.end(), 0.0);
    for (int i = 0; i < nRes; i++) {
      const double *Ji = &J[(size_t)i * nFree];
      double ri = r[i];
      for (int a = 0; a < nFree; a++) {
        g[a] += Ji[a] * ri;
        double Jia = Ji[a];
        double *Ha = &H[(size_t)a * nFree];
        for (int b = a; b < nFree; b++) Ha[b] += Jia * Ji[b];
      }
    }
    for (int a = 0; a < nFree; a++) {
      for (int b = 0; b < a; b++) H[(size_t)a * nFree + b] = H[(size_t)b * nFree + a];  // symmetrize
      if (pen_inv2[a] != 0.0) {
        H[(size_t)a * nFree + a] += pen_inv2[a];
        g[a] += (x[a] - pen_nom[a]) * pen_inv2[a];
      }
    }

    // Convergence on the (projected) gradient.
    double gmax = 0.0;
    for (int a = 0; a < nFree; a++) gmax = std::max(gmax, std::fabs(g[a]));
    if (gmax < 1.0e-8) break;

    // Inner loop: damp until a step decreases the cost.
    bool accepted = false;
    for (int tries = 0; tries < 40 && !accepted; tries++) {
      // (H + lambda*diag(H)) dx = -g
      Hl = H;
      for (int a = 0; a < nFree; a++) Hl[(size_t)a * nFree + a] += lambda * H[(size_t)a * nFree + a];

      gsl_matrix *A = gsl_matrix_alloc(nFree, nFree);
      for (int a = 0; a < nFree; a++)
        for (int b = 0; b < nFree; b++) gsl_matrix_set(A, a, b, Hl[(size_t)a * nFree + b]);
      gsl_vector *bvec = gsl_vector_alloc(nFree);
      for (int a = 0; a < nFree; a++) gsl_vector_set(bvec, a, -g[a]);
      gsl_vector *sol = gsl_vector_alloc(nFree);

      gsl_set_error_handler_off();
      int status = gsl_linalg_cholesky_decomp1(A);
      if (status == 0) status = gsl_linalg_cholesky_solve(A, bvec, sol);

      bool solved = (status == 0);
      if (solved)
        for (int a = 0; a < nFree; a++) dx[a] = gsl_vector_get(sol, a);
      gsl_matrix_free(A);
      gsl_vector_free(bvec);
      gsl_vector_free(sol);

      if (!solved) {
        lambda *= 4.0;
        if (lambda > 1.0e14) {
          accepted = false;
          break;
        }
        continue;
      }

      vector_r xnew(nFree), fullNew = full;
      for (int a = 0; a < nFree; a++) {
        xnew[a] = std::min(std::max(x[a] + dx[a], lo[a]), hi[a]);  // clamp to limits
        fullNew[p2f[a]] = xnew[a];
      }
      double costNew = Chi2Value(fullNew);

      if (costNew < cost && std::isfinite(costNew)) {
        double rel = (cost - costNew) / std::max(std::fabs(cost), 1.0);
        full = fullNew;
        x = xnew;
        cost = costNew;
        lambda = std::max(lambda * 0.3, 1.0e-12);
        accepted = true;
        configure().outStream << "\r\tLM iteration: " << std::setw(4) << iter + 1
                              << "  Chi-Squared: " << cost << "        ";
        configure().outStream.flush();
        if (rel < 1.0e-9) iter = maxIter;  // converged
      } else {
        lambda *= 4.0;
        if (lambda > 1.0e14) break;
      }
    }
    if (!accepted) break;
  }
  configure().outStream << std::endl;

  // Write best-fit values back, and parameter errors from (J^T J + penalties)^{-1}.
  for (int i = 0; i < nFull; i++) mp.SetValue(i, full[i]);
  FinalizeLeastSquaresCovariance(full, p2f, pen_inv2, mp, bandCovOut);

  return cost;
}

int AZURECalc::PrepareFreeParams(const vector_r &full,
                                 const ROOT::Minuit2::MnUserParameters &mp,
                                 vector_r &r, vector_r &J, std::vector<int> &p2f,
                                 vector_r &x, std::vector<double> &lo,
                                 std::vector<double> &hi, std::vector<double> &pen_nom,
                                 std::vector<double> &pen_inv2) const {
  // First Jacobian fixes the column <-> parameter mapping.
  if (!ResidualJacobian(full, r, J, p2f)) return -1;  // unsupported -> caller falls back
  const int nFree = (int)p2f.size();

  // Free-parameter values, limits, and quadratic-penalty (prior) terms.
  x.assign(nFree, 0.0);
  lo.assign(nFree, -std::numeric_limits<double>::infinity());
  hi.assign(nFree, std::numeric_limits<double>::infinity());
  pen_nom.assign(nFree, 0.0);
  pen_inv2.assign(nFree, 0.0);
  for (int a = 0; a < nFree; a++) {
    int f = p2f[a];
    x[a] = full[f];
    const auto &par = mp.Parameter(f);
    if (par.HasLowerLimit()) lo[a] = par.LowerLimit();
    if (par.HasUpperLimit()) hi[a] = par.UpperLimit();
    // Gaussian penalties contribute (1/sigma^2) to the diagonal of J^T J and
    // (x-nominal)/sigma^2 to the gradient.  Norm / shift / nuisance penalties
    // are identified the same way operator()/CalculateNuisanceChiSquared do.
    std::string name = par.GetName();
    if (name.find("norm") != std::string::npos) {
      for (int s = 1; s <= data()->NumSegments(); s++) {
        ESegment *seg = data()->GetSegment(s);
        if (!seg || !seg->IsVaryNorm()) continue;
        char vn[64];
        snprintf(vn, sizeof(vn), "segment_%d_norm", seg->GetSegmentKey());
        if (name == vn) {
          double n0 = seg->GetNominalNorm();
          double sig = n0 / 100.0 * seg->GetNormError();
          if (sig != 0.0) {
            pen_nom[a] = n0;
            pen_inv2[a] = 1.0 / (sig * sig);
          }
          break;
        }
      }
    } else if (name.find("shift") != std::string::npos) {
      for (int s = 1; s <= data()->NumSegments(); s++) {
        ESegment *seg = data()->GetSegment(s);
        if (!seg) continue;
        char vn[64], vs[64];
        snprintf(vn, sizeof(vn), "segment_%d_energy_shift", seg->GetSegmentKey());
        snprintf(vs, sizeof(vs), "segment_%d_energy_shift_sqrt", seg->GetSegmentKey());
        if (name == vn && seg->IsVaryEnergyShift()) {
          double sig = seg->GetEnergyShiftError();
          if (sig != 0.0) {
            pen_nom[a] = seg->GetNominalEnergyShift();
            pen_inv2[a] = 1.0 / (sig * sig);
          }
          break;
        }
        if (name == vs && seg->IsVaryEnergyShiftSqrt()) {
          double sig = seg->GetEnergyShiftSqrtError();
          if (sig != 0.0) {
            pen_nom[a] = seg->GetNominalEnergyShiftSqrt();
            pen_inv2[a] = 1.0 / (sig * sig);
          }
          break;
        }
      }
    } else if (limitsManager_ && limitsManager_->IsNuisanceParameterByIndex(a)) {
      double sig = limitsManager_->GetConvertedErrorByIndex(a);
      if (sig > 0.0) {
        pen_nom[a] = limitsManager_->GetConvertedNominalValueByIndex(a);
        pen_inv2[a] = 1.0 / (sig * sig);
      }
    }
  }
  return nFree;
}

void AZURECalc::FinalizeLeastSquaresCovariance(const vector_r &full,
                                               const std::vector<int> &p2f,
                                               const std::vector<double> &pen_inv2,
                                               ROOT::Minuit2::MnUserParameters &mp,
                                               BandCovariance *bandCovOut) const {
  const int nFree = (int)p2f.size();
  vector_r r, J;
  std::vector<int> p2fLocal;
  if (!ResidualJacobian(full, r, J, p2fLocal)) return;

  std::vector<double> H(nFree * nFree, 0.0);
  const int nRes = (int)r.size();
  for (int i = 0; i < nRes; i++) {
    const double *Ji = &J[(size_t)i * nFree];
    for (int a = 0; a < nFree; a++) {
      double *Ha = &H[(size_t)a * nFree];
      for (int b = a; b < nFree; b++) Ha[b] += Ji[a] * Ji[b];
    }
  }
  // Pure J^T J diagonal = each parameter's data sensitivity (column-norm^2 in
  // J), captured before priors are folded into H.  Near zero => unconstrained.
  std::vector<double> jtjDiag(nFree);
  for (int a = 0; a < nFree; a++) {
    for (int b = 0; b < a; b++) H[(size_t)a * nFree + b] = H[(size_t)b * nFree + a];
    jtjDiag[a] = H[(size_t)a * nFree + a];
    if (pen_inv2[a] != 0.0) H[(size_t)a * nFree + a] += pen_inv2[a];
  }

  // The Minuit name alone ("width_1_2") does not say which level or channel is
  // meant, so pair it with the level's J^pi and energy and the channel's pair,
  // L and S.
  auto describe = [&](int a) {
    const int minuitIndex = p2f[a];
    std::string text = mp.GetName(minuitIndex);
    const std::string label = AZURELabel::Parameter(compound(), data(), minuitIndex);
    if (!label.empty()) text += "  =  " + label;
    return text;
  };

  // Name the free parameters the data barely constrain (near-zero J column).
  auto reportInsensitive = [&]() {
    double mx = 0.0;
    for (int a = 0; a < nFree; a++) mx = std::max(mx, jtjDiag[a]);
    int nInsensitive = 0;
    for (int a = 0; a < nFree; a++) {
      if (mx <= 0.0 || jtjDiag[a] <= 1.e-10 * mx) {
        if (!nInsensitive)
          configure().outStream << "  Data-insensitive parameter(s) -- fix, constrain, or add "
                                   "a prior to one of these:"
                                << std::endl;
        configure().outStream << "    " << describe(a) << std::endl;
        nInsensitive++;
      }
    }
    if (!nInsensitive)
      configure().outStream << "No single parameter is insensitive; the flat direction is a "
                               "degeneracy between two or more correlated parameters."
                            << std::endl;
  };

  // Covariance = (J^T J + priors)^{-1}.  If it is singular (unconstrained or
  // degenerate parameter), retry with a small ridge so the flat direction gets
  // a large finite variance instead of losing the errors and band; warn.
  double maxDiag = 0.0;
  for (int a = 0; a < nFree; a++) maxDiag = std::max(maxDiag, H[(size_t)a * nFree + a]);
  const double ridges[4] = {0.0, maxDiag * 1.e-9, maxDiag * 1.e-6, maxDiag * 1.e-3};
  gsl_matrix *A = gsl_matrix_alloc(nFree, nFree);
  gsl_set_error_handler_off();
  int used = -1;
  for (int t = 0; t < 4; t++) {
    for (int a = 0; a < nFree; a++)
      for (int b = 0; b < nFree; b++)
        gsl_matrix_set(A, a, b, H[(size_t)a * nFree + b] + (a == b ? ridges[t] : 0.0));
    if (gsl_linalg_cholesky_decomp1(A) == 0 && gsl_linalg_cholesky_invert(A) == 0) {
      used = t;
      break;
    }
  }
  if (used < 0) {
    configure().outStream << "Warning: the fit curvature matrix (J^T J) is singular; parameter errors "
                             "and the uncertainty band are unavailable."
                          << std::endl;
    reportInsensitive();
  } else {
    if (used > 0) {
      configure().outStream << "Warning: the fit curvature matrix was near-singular; regularized it "
                               "with a small ridge to obtain the covariance. Uncertainties along the "
                               "poorly-constrained direction(s) are large and approximate."
                            << std::endl;
      reportInsensitive();
      // Strongly correlated pairs (|rho| near 1) are the degenerate combinations.
      struct CorrPair {
        double c;
        int a, b;
      };
      std::vector<CorrPair> pairs;
      for (int a = 0; a < nFree; a++) {
        double va = gsl_matrix_get(A, a, a);
        if (va <= 0.0) continue;
        for (int b = a + 1; b < nFree; b++) {
          double vb = gsl_matrix_get(A, b, b);
          if (vb <= 0.0) continue;
          double c = gsl_matrix_get(A, a, b) / std::sqrt(va * vb);
          if (std::fabs(c) > 0.95) pairs.push_back({c, a, b});
        }
      }
      std::sort(pairs.begin(), pairs.end(),
                [](const CorrPair &x, const CorrPair &y) { return std::fabs(x.c) > std::fabs(y.c); });
      for (size_t k = 0; k < pairs.size() && k < 5; k++)
        configure().outStream << "  Correlated (rho=" << pairs[k].c << "):" << std::endl
                              << "    " << describe(pairs[k].a) << std::endl
                              << "    " << describe(pairs[k].b) << std::endl;
    }
    for (int a = 0; a < nFree; a++) {
      double v = gsl_matrix_get(A, a, a);
      if (v > 0.0) mp.SetError(p2f[a], std::sqrt(v));
    }
    // Name free parameters whose error dwarfs their value (they dominate the band).
    {
      int n = 0;
      for (int a = 0; a < nFree; a++) {
        double v = gsl_matrix_get(A, a, a);
        if (v <= 0.0) continue;
        double val = mp.Value(p2f[a]);
        if (std::sqrt(v) > 10.0 * std::max(std::fabs(val), 1.e-30)) {
          if (!n)
            configure().outStream << "Note: weakly-determined parameter(s) (error > 1000% of value) "
                                     "dominate the uncertainty band. Constrain or fix these to "
                                     "tighten it:"
                                  << std::endl;
          configure().outStream << "    " << describe(a) << std::endl;
          n++;
        }
      }
    }
    // Export the full covariance for cross-section bands.  Columns are the
    // packed free parameters, in the same order as p2f, so we tag them with
    // their parameter identities via a matching index map.
    if (bandCovOut) {
      AZUREParams tp;
      compound()->FillMnParams(tp.GetMinuitParams(), &configure());
      data()->FillMnParams(tp.GetMinuitParams());
      int nMn = tp.GetMinuitParams().Params().size();
      std::vector<bool> fixed(nMn);
      for (int i = 0; i < nMn; i++) fixed[i] = tp.GetMinuitParams().Parameter(i).IsFixed();
      ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), fixed);
      if (pmap.NumPacked() == nFree) {
        // Keep only the R-matrix sub-block: the band is insensitive to norms
        // and energy shifts, so they are dropped from the saved covariance.
        const std::vector<int> rc = RMatrixPackedColumns(pmap);
        const int m = (int)rc.size();
        bandCovOut->cols.resize(m);
        for (int a = 0; a < m; a++) bandCovOut->cols[a] = pmap.Desc(pmap.PackedToFull(rc[a]));
        bandCovOut->M.assign(m, std::vector<double>(m, 0.0));
        for (int a = 0; a < m; a++)
          for (int b = 0; b < m; b++) bandCovOut->M[a][b] = gsl_matrix_get(A, rc[a], rc[b]);
      }
    }
  }
  gsl_matrix_free(A);
}

bool AZURECalc::ResidualsOnly(const vector_r &full, vector_r &residuals) const {
  // Forward-only counterpart of ResidualJacobian (no adjoint): same row order
  // and residual r = (model - data*n)/(cmErr*n) as ComputeResidualJacobian.
  // Pooled working copies keep their applied energy shifts and mapping between calls; a fresh
  // Clone() re-shifted every shifted segment and rebuilt the mapping each time (~1 s here).
  CNuc *lc = GetPooledCNuc();
  EData *ld = GetPooledEData();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());
  if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure());

  residuals.clear();
  for (int i = 1; i <= ld->NumSegments(); i++) {
    ESegment *segment = ld->GetSegment(i);
    if (!segment) continue;
    const double norm = segment->GetNorm();
    for (int pid = 0; pid < segment->NumPoints(); pid++) {
      EPoint *pt = segment->GetPoint(pid + 1);
      if (!pt) continue;
      double th = segment->CalculateTheoreticalCrossSection(pid, lc, configure(), ld);
      double denom = pt->GetCMCrossSectionError() * norm;
      residuals.push_back(denom != 0.0 ? (th - pt->GetCMCrossSection() * norm) / denom : 0.0);
    }
  }
  ReturnPooledCNuc(lc);
  ReturnPooledEData(ld);
  return true;
}

namespace {

// State shared with the GSL least-squares callbacks; owned by RunGSLNonlinear.
struct GSLNlsData {
  const AZURECalc *self;
  int nFree;
  int nData;                    // number of data residual rows
  const std::vector<int> *p2f;  // packed free param -> full-vector index
  const std::vector<double> *lo;
  const std::vector<double> *hi;
  const std::vector<double> *pen_nom;
  const std::vector<double> *pen_inv2;
  const std::vector<int> *penRows;  // free-param indices that carry a prior
  vector_r baseFull;                // full vector holding the fixed params
};

// Project x into [lo,hi] and scatter it into the full parameter vector.
void gslExpandFull(const GSLNlsData *d, const gsl_vector *x, vector_r &full) {
  full = d->baseFull;
  for (int a = 0; a < d->nFree; a++) {
    double v = gsl_vector_get(x, a);
    v = std::min(std::max(v, (*d->lo)[a]), (*d->hi)[a]);
    full[(*d->p2f)[a]] = v;
  }
}

// Residual vector: data rows plus Gaussian-prior rows sqrt(1/sigma^2)*(x-nom).
int gslResidualF(const gsl_vector *x, void *params, gsl_vector *f) {
  const GSLNlsData *d = static_cast<const GSLNlsData *>(params);
  vector_r full;
  gslExpandFull(d, x, full);
  vector_r res;
  if (!d->self->ResidualsOnly(full, res) || (int)res.size() != d->nData) return GSL_EBADFUNC;
  for (int i = 0; i < d->nData; i++) gsl_vector_set(f, i, res[i]);
  for (size_t k = 0; k < d->penRows->size(); k++) {
    int a = (*d->penRows)[k];
    double s = std::sqrt((*d->pen_inv2)[a]);
    gsl_vector_set(f, d->nData + (int)k, s * (full[(*d->p2f)[a]] - (*d->pen_nom)[a]));
  }
  return GSL_SUCCESS;
}

// Residual Jacobian: analytic data-row block plus the constant prior rows.
int gslResidualDf(const gsl_vector *x, void *params, gsl_matrix *Jm) {
  const GSLNlsData *d = static_cast<const GSLNlsData *>(params);
  vector_r full;
  gslExpandFull(d, x, full);
  vector_r r, Jflat;
  std::vector<int> p2f;
  if (!d->self->ResidualJacobian(full, r, Jflat, p2f)) return GSL_EBADFUNC;
  if ((int)p2f.size() != d->nFree || (int)r.size() != d->nData) return GSL_EBADFUNC;
  const int nFree = d->nFree;
  for (int i = 0; i < d->nData; i++)
    for (int a = 0; a < nFree; a++)
      gsl_matrix_set(Jm, i, a, Jflat[(size_t)i * nFree + a]);
  for (size_t k = 0; k < d->penRows->size(); k++) {
    int row = d->nData + (int)k;
    for (int a = 0; a < nFree; a++) gsl_matrix_set(Jm, row, a, 0.0);
    int a = (*d->penRows)[k];
    gsl_matrix_set(Jm, row, a, std::sqrt((*d->pen_inv2)[a]));
  }
  return GSL_SUCCESS;
}

// Per-iteration progress line (chi^2 = ||residual||^2, priors included).
void gslProgress(const size_t iter, void *params,
                 const gsl_multifit_nlinear_workspace *w) {
  const GSLNlsData *d = static_cast<const GSLNlsData *>(params);
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  double chi2 = 0.0;
  gsl_blas_ddot(f, f, &chi2);
  d->self->configure().outStream << "\r\tGSL-LM iteration: " << std::setw(4) << iter
                                 << "  Chi-Squared: " << chi2 << "        ";
  d->self->configure().outStream.flush();

  // Intermediate output, so a long GSL fit also leaves usable files behind
  // before it converges rather than only at the end.
  if (iter > 0 && iter % AZURECalc::kIterationOutputInterval == 0) {
    vector_r full = d->baseFull;
    const gsl_vector *x = gsl_multifit_nlinear_position(w);
    for (int a = 0; a < d->nFree; a++) full[(*d->p2f)[a]] = gsl_vector_get(x, a);
    d->self->WriteIterationOutput(full);
  }
}

}  // namespace

double AZURECalc::RunGSLNonlinear(AZUREParams &params, int maxIter,
                                  BandCovariance *bandCovOut) const {
  ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
  const int nFull = mp.Params().size();
  vector_r full(nFull);
  for (int i = 0; i < nFull; i++) full[i] = mp.Value(i);

  vector_r r0, J0, x;
  std::vector<int> p2f;
  std::vector<double> lo, hi, pen_nom, pen_inv2;
  const int nFree = PrepareFreeParams(full, mp, r0, J0, p2f, x, lo, hi, pen_nom, pen_inv2);
  if (nFree < 0) return -1.0;  // analytic Jacobian unsupported -> caller falls back
  if (nFree == 0) return Chi2Value(full);
  const int nData = (int)r0.size();

  // Gaussian priors on norms / shifts / nuisance params become extra residual rows.
  std::vector<int> penRows;
  for (int a = 0; a < nFree; a++)
    if (pen_inv2[a] > 0.0) penRows.push_back(a);
  const int nRes = nData + (int)penRows.size();

  GSLNlsData d;
  d.self = this;
  d.nFree = nFree;
  d.nData = nData;
  d.p2f = &p2f;
  d.lo = &lo;
  d.hi = &hi;
  d.pen_nom = &pen_nom;
  d.pen_inv2 = &pen_inv2;
  d.penRows = &penRows;
  d.baseFull = full;

  gsl_multifit_nlinear_fdf fdf;
  fdf.f = gslResidualF;
  fdf.df = gslResidualDf;
  fdf.fvv = nullptr;  // second directional derivative via finite differences of f
  fdf.n = nRes;
  fdf.p = nFree;
  fdf.params = &d;

  // Geodesic-accelerated LM in a trust region; the default solver factorizes J
  // directly (QR) rather than forming J^T J.
  gsl_multifit_nlinear_parameters gp = gsl_multifit_nlinear_default_parameters();
  gp.trs = gsl_multifit_nlinear_trs_lmaccel;

  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w = gsl_multifit_nlinear_alloc(T, &gp, nRes, nFree);
  if (!w) return -1.0;

  gsl_vector *x0 = gsl_vector_alloc(nFree);
  for (int a = 0; a < nFree; a++) gsl_vector_set(x0, a, x[a]);

  gsl_set_error_handler_off();
  gsl_multifit_nlinear_init(x0, &fdf, w);

  const double xtol = 1.0e-9, gtol = 1.0e-9, ftol = 1.0e-9;
  int info = 0;
  int status = gsl_multifit_nlinear_driver(maxIter, xtol, gtol, ftol,
                                           gslProgress, &d, &info, w);
  configure().outStream << std::endl;

  std::string reason = "reached the maximum number of iterations";
  if (status == GSL_SUCCESS)
    reason = (info == 1) ? "converged (small step)" : "converged (small gradient)";
  else if (status != GSL_EMAXITER)
    reason = std::string("stopped: ") + gsl_strerror(status);
  configure().outStream << "GSL trust-region solver " << reason << " after "
                        << gsl_multifit_nlinear_niter(w) << " iterations." << std::endl;

  // Read back the (projected) best-fit parameters.
  gsl_vector *xf = gsl_multifit_nlinear_position(w);
  for (int a = 0; a < nFree; a++) {
    double v = gsl_vector_get(xf, a);
    v = std::min(std::max(v, lo[a]), hi[a]);
    full[p2f[a]] = v;
  }
  for (int i = 0; i < nFull; i++) mp.SetValue(i, full[i]);

  gsl_vector_free(x0);
  gsl_multifit_nlinear_free(w);

  double cost = Chi2Value(full);

  // Errors and band covariance from (J^T J + priors)^{-1}, as the LM path does.
  FinalizeLeastSquaresCovariance(full, p2f, pen_inv2, mp, bandCovOut);

  return cost;
}

void AZURECalc::AddNuisanceGradient(const vector_r &p, std::vector<double> &grad) const {
  if (!limitsManager_) return;

  AZUREParams tempParams;
  compound()->FillMnParams(tempParams.GetMinuitParams(), &configure());
  data()->FillMnParams(tempParams.GetMinuitParams());

  std::vector<int> nonFixedToActualIndex;
  for (int i = 0; i < tempParams.GetMinuitParams().Params().size(); i++) {
    if (!tempParams.GetMinuitParams().Parameter(i).IsFixed() ||
        tempParams.GetMinuitParams().Parameter(i).GetName().find("segment") != std::string::npos) {
      nonFixedToActualIndex.push_back(i);
    }
  }

  for (int nf = 0; nf < (int)nonFixedToActualIndex.size() && nf < (int)p.size(); nf++) {
    int actualIndex = nonFixedToActualIndex[nf];
    std::string paramName = tempParams.GetMinuitParams().Parameter(actualIndex).GetName();
    if (paramName.find("norm") != std::string::npos || paramName.find("shift") != std::string::npos) continue;
    if (!limitsManager_->IsNuisanceParameterByIndex(nf)) continue;
    double nominalValue = limitsManager_->GetConvertedNominalValueByIndex(nf);
    double paramError = limitsManager_->GetConvertedErrorByIndex(nf);
    if (paramError > 0.0) {
      // d/d p[nf] of ((p[nf]-nominal)/error)^2.
      grad[nf] += 2.0 * (p[nf] - nominalValue) / (paramError * paramError);
    }
  }
}

double AZURECalc::CalculateNuisanceChiSquared(const vector_r &p) const {
  double nuisanceChiSquared = 0.0;

  // Create temporary AZUREParams to get parameter names
  AZUREParams tempParams;
  compound()->FillMnParams(tempParams.GetMinuitParams(), &configure());
  data()->FillMnParams(tempParams.GetMinuitParams());

  // Build mapping from non-fixed parameter index to actual parameter index
  std::vector<int> nonFixedToActualIndex;
  for (int i = 0; i < tempParams.GetMinuitParams().Params().size(); i++) {
    if (!tempParams.GetMinuitParams().Parameter(i).IsFixed() || tempParams.GetMinuitParams().Parameter(i).GetName().find("segment") != std::string::npos) {
      nonFixedToActualIndex.push_back(i);
    }
  }

  // Check each non-fixed parameter to see if it's marked as nuisance
  for (int nonFixedIndex = 0; nonFixedIndex < nonFixedToActualIndex.size() && nonFixedIndex < p.size(); nonFixedIndex++) {
    int actualIndex = nonFixedToActualIndex[nonFixedIndex];
    std::string paramName = tempParams.GetMinuitParams().Parameter(actualIndex).GetName();

    // If norm or shift in param name, skip
    if (paramName.find("norm") != std::string::npos || paramName.find("shift") != std::string::npos) {
      continue;
    }

    // First check if this parameter is marked as nuisance (fast check)
    if (!limitsManager_->IsNuisanceParameterByIndex(nonFixedIndex)) {
      continue;  // Skip if not a nuisance parameter
    }

    // Only do expensive conversions if parameter is marked as nuisance
    double nominalValue = limitsManager_->GetConvertedNominalValueByIndex(nonFixedIndex);
    double paramError = limitsManager_->GetConvertedErrorByIndex(nonFixedIndex);

    // If we got valid values (non-zero error means this parameter has valid nuisance settings)
    if (paramError > 0.0) {
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
  int pool_size = 4;  // default minimum
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
CNuc *AZURECalc::GetPooledCNuc() const {
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
EData *AZURECalc::GetPooledEData() const {
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
void AZURECalc::ReturnPooledCNuc(CNuc *obj) const {
  if (!obj) return;

  std::lock_guard<std::mutex> lock(pool_mutex_);
  cnuc_pool_.push(std::unique_ptr<CNuc>(obj));
}

/*!
 * Return an EData object to the pool for reuse
 */
void AZURECalc::ReturnPooledEData(EData *obj) const {
  if (!obj) return;

  std::lock_guard<std::mutex> lock(pool_mutex_);
  edata_pool_.push(std::unique_ptr<EData>(obj));
}

/*!
 * Write parameters to file
 */
void AZURECalc::WriteParameters(AZUREParams &params, const Config &configure) const {
  char filename[256];
  snprintf(filename, sizeof(filename), "%sparam.fit", configure.outputdir.c_str());
  std::ofstream out;
  out.open(filename);
  if (out) {
    out.precision(7);
    for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
      out << std::setw(20) << params.GetMinuitParams().GetName(i)
          << std::scientific << std::setw(20) << params.GetMinuitParams().Value(i)
          << std::scientific << std::setw(20) << params.GetMinuitParams().Error(i) << std::endl;
    }
    out.flush();
    out.close();
  } else
    configure.outStream << "Could not save param.fit file." << std::endl;
}