// AZURECalc, the gradient: d chi^2 / d p of the evaluation of AZURECalc.cpp as Minuit2 asks for it (Gradient,
// GradientBy: the T-matrix adjoint for levels and norms, the THM rows, central differences for what the adjoint
// does not cover, the priors, and the AZURE_GRAD_CHECK comparison), the same gradient from the residual Jacobian
// (GradientFromJacobian), and the Gauss-Newton second derivatives Minuit2 seeds its step with (G2, EnableG2).  One
// unit because all of them are the derivative of one function and share the fixed-parameter mask and the
// finite-difference step.  What they keep between calls is in the caller's FitterStates (states_.gradient).
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
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <mutex>
#include <sstream>
#include "ParameterLabel.h"
#include "Config.h"
#include "CompoundSystem.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "Unit.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "AZUREGrad.h"
#include "CovarianceBand.h"
#include "THMMatrixFunc.h"
#include "ThmCotangents.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {


// Central-difference step for parameter x0.  The usual scale*(|x0| + 1) is
// absolute for |x0| < 1, which is right for widths and energies that pass
// through zero but not for a normalization: a yield normalized in absolute
// units (11B+a CASPAR thick-target yields, n ~ 2e-5) moved by 1e-6, i.e. 5%,
// and its finite-difference gradient came out 14x too large.  A norm is
// positive and scale-free, so its step is relative.
double FdStep(ParamKind kind, double x0, double scale) {
  if (kind == ParamKind::Norm && x0 != 0.0) return scale * std::fabs(x0);
  return scale * (std::fabs(x0) + 1.0);
}

// chi^2 contribution of a single segment: its data term plus its own energy-shift nuisance penalty, evaluated at the
// segment's current shift, as AZURECalc::Chi2Value has it, so a finite difference of this over the segment's shift
// equals one of the whole chi^2 when no other segment's contribution depends on this shift (a self-contained segment:
// the eligibility guard in Gradient()).  A folded segment on the fixed nodes is folded again, no node recomputed.

}  // namespace

std::vector<double> AZURECalc::Gradient(const std::vector<double> &p) const { return GradientBy(p, GradientRoute()); }

std::vector<double> AZURECalc::GradientBy(const std::vector<double> &p, const GradientRoute &route) const {
  if (data()->IsFit() && StopRequested()) return std::vector<double>(p.size(), 0.0);  // a stopped fit: see operator()
  AZUREProfile::ScopedReport profileReport(configure().profile, "gradient");
  AZUREProfile::Timer profileTimer(configure().profile, AZUREProfile::kGradient);
  {
    std::lock_guard<std::mutex> lk(states_.gradient.mutex);
    states_.gradient.at = p;
  }
  std::vector<double> grad(p.size(), 0.0);
  const auto tg0 = std::chrono::steady_clock::now();
  // With G2 enabled, one residual Jacobian gives both the gradient (2 J^T r +
  // penalties) and Minuit's G2, which it requests at the same point next.
  if (g2Enabled_ && !configure().environment.gradientByAdjoint && !route.adjointOnly && GradientFromJacobian(p, grad)) {
    const azure2::support::Environment &environment = azure2::support::ProcessEnvironment();
    if (environment.Set(azure2::support::Switch::kGradTiming)) {
      int c = ++configure().diagnostics->jacobianPathTimings;
      if (c <= 5 || c % 50 == 0)
        std::cerr << "[grad-timing] call " << c << " (Jacobian path, gradient + G2): "
                  << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tg0).count() << " ms" << std::endl;
    }
    // AZURE_GRAD_CHECK: on the first call also take the adjoint gradient (which
    // itself runs the sampled finite-difference check) and compare all components.
    if (environment.Set(azure2::support::Switch::kGradCheck) && !configure().diagnostics->jacobianPathChecked.exchange(true)) {
      GradientRoute adjoint;
      adjoint.adjointOnly = true;
      adjoint.jacobianPathGradient = &grad;
      const std::vector<double> ga = GradientBy(p, adjoint);
      AZUREParams fpc;
      compound()->FillMnParams(fpc, state_, inputs_, &configure());
      data()->FillMnParams(fpc);
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
      if (environment.Text(azure2::support::Switch::kGradCheck) == "exit") ExitRun(configure(), 0);
    }
    {
      std::lock_guard<std::mutex> lk(states_.gradient.mutex);
      states_.gradient.value = grad;
    }
    return grad;
  }
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);

  // --- Analytic energy / reduced-width block via the shared adjoint engine. ---
  // Pooled working copies keep their applied energy shifts and mapping between calls; a fresh
  // Clone() re-shifted every shifted segment and rebuilt the mapping each time (~1 s here).
  const CompoundSystem *lc = compound();
  std::unique_ptr<LevelState> pooled = GetPooledState();
  LevelState &ls = *pooled;
  WorkingData *working = GetPooledData();
  WorkingCopy *ld = working->data.get();
  azure2::Results &results = working->results;
  lc->FillCompoundFromParams(p, ls);
  ld->FillNormsFromParams(p);
  ld->FillEnergyShiftsFromParams(p, results, lc, &ls, &inputs_, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r *sdp = nullptr;
  if (brune) {
    lc->CalcShiftFunctions(configure(), ls);
    shiftDeriv = BuildShiftDerivTable(lc, ls, configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, std::vector<bool>());
  GradAccum accum;
  accum.Init(lc, ls);

  // chi^2 = sum (fit - data*n)^2 / (err*n)^2 + penalties.  The data-term model
  // cotangent is d(chi2)/d(model) = 2 r / err^2.  The same per-point model also
  // gives the analytic normalization gradient, accumulated per segment here so
  // it needs no extra forward pass:
  //   d(chi2)/d(n_s) = sum_pts [ -2 r data/(err^2 n^2) - 2 r^2/(err^2 n^3) ], r = model/J - data*n (DataResidual).
  std::vector<double> normData(ld->Data().NumSegments() + 1, 0.0);
  // A data covariance couples the points: their cotangent is 2 (M^-1 r)_i / (err n J), which
  // needs every residual of the group first (one forward pass of those points).
  if (ld->Data().HasCovariance()) ld->PrepareCovarianceSeeds(lc, ls, configure(), results);
  FitBarFn fitBarFn = [&](const Segment *seg, int i, int pid, double model) -> double {
    Unit *pt = seg->GetPoint(*ld, pid + 1);
    if (!pt) return 0.0;
    double norm = ld->State(*seg).norm;
    const double dataval = seg->LabValue(pid), labErr = seg->LabError(pid), jac = Segment::Jacobian(ld->State(*seg), 0, pid);
    const azure2::Residual res = azure2::DataResidual(model, jac, dataval, labErr, norm);
    const double r = res.difference, err = res.error;
    if (err == 0.0) return 0.0;
    if (seg->GetCovarianceGroup() >= 0) {
      // chi2 = r^T M^-1 r with r_i = (model/J - d n)/(err n): d chi2/d r_i = 2 z_i, z = M^-1 r;
      // d r_i/d model = 1/(err n J), d r_i/d n = -(model/J)/(err n^2).
      const double z = pt->GetCovarianceSeed(results);
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

  // THM data in the same pass where it can take them (ThmCotangents.h): hoes is their adjoint's context, else null.
  const std::unique_ptr<ThmAdjointContext> hoes = azure2::ThmGradientContext(lc, ls, *ld, configure(), sdp);
  bool eg = AccumulateEGammaGradient(lc, ls, ld, results, configure(), pmap, sdp, fitBarFn, accum, hoes.get());
  const auto tg1 = std::chrono::steady_clock::now();
  if (!eg && azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kGradDebug)) {
    std::cerr << "[grad] analytic energy/gamma adjoint bailed -> full finite "
                 "differences (no speed-up). An unsupported segment/config is "
                 "present (RMC, or a not-yet-handled case)."
              << std::endl;
  }
  // The model data (T9); where their gradient is not available every parameter falls back to finite differences.
  if (eg && !Objective::Engine().ModelDataGradient(ObjectiveContext{lc, ls, inputs_, ld, results, p, configure()}, pmap, sdp, accum))
    eg = false;
  if (eg) {
    AddParkPenaltyGradient(lc, ls, configure(), accum); // Park J > 0 wall
    AddAbundancePenaltyGradient(*ld, configure(), accum);  // abundance priors, sum-to-one
    accum.Scatter(pmap, grad);                      // energies + reduced widths, abundances, screening
    AddNuisanceGradient(p, grad);                   // nuisance penalty
    // The diagonal priors (DiagonalPriors.h) of screening, Doppler temperature and sample thickness on p; density and
    // apparatus priors are in the central differences of the whole chi^2 below, the norms' in the loop after this.
    const DiagonalPriors priors = DiagonalPriorsOf(*ld, pmap);
    for (const DiagonalPrior &e : priors)
      if ((e.source == PriorSource::kScreening || e.source == PriorSource::kDoppler || e.source == PriorSource::kThickness) &&
          e.index < (int)grad.size())
        grad[e.index] += GaussianPriorGradient(p[e.index], e.centre, e.sigma);
    for (int s = 1; s <= ld->Data().NumSegments(); s++) {  // normalizations
      const Segment *seg = ld->Data().GetSegment(s);
      if (!seg || !seg->IsVaryNorm() || seg->IsNormGroupFollower()) continue;  // a follower's norm is fixed
      // THM segments never enter the analytic accumulation (their normData
      // stays 0); their norm gradient -- data term AND penalty -- comes from
      // the THM-only finite differences below.
      if (seg->IsTHM()) continue;
      int idx = pmap.NormIndex(s);
      if (idx < 0 || idx >= (int)grad.size()) continue;
      double g = normData[s];
      // The norm's prior, on the working copy's norm (set from p: FillNormsFromParams).
      for (const DiagonalPrior &e : priors)
        if (e.source == PriorSource::kNorm && e.group == s) g += GaussianPriorGradient(e.value, e.centre, e.sigma);
      grad[idx] = g;
    }
  }

  // Fixed-parameter mask: Minuit ignores the gradient of fixed parameters, so
  // do not waste finite differences on them (large fits fix most energy shifts).
  // The minimizer's flags (FixedMask), not only the file's.
  const std::vector<bool> fixedMask = FixedMask();
  const int nMn = (int)fixedMask.size();

  // THM (HOES) segments are excluded from the T-matrix adjoint; their
  // energy/gamma contribution comes from the HOES adjoint below (their norms,
  // and everything with AZURE2_THM_FD=1, by central differences of the
  // THM-only chi^2).
  bool haveTHM = false;
  for (int s = 1; s <= data()->Data().NumSegments() && !haveTHM; s++) {
    const Segment *seg = data()->Data().GetSegment(s);
    if (seg && seg->IsTHM()) haveTHM = true;
  }

  // Else THM rows by the HOES adjoint (ComputeTHMRows): the THM data term is sum_i r_i^2, the norm profiled as
  // Chi2Value(.., true) does, so its level / coherent-background gradient is 2 J^T r.  THM norms keep the central
  // differences of the THM-only chi^2 below (with its norm penalty); AZURE2_THM_FD=1 for every THM column.
  bool thmAnalytic = eg && hoes;
  if (eg && haveTHM && !hoes && !ThmFiniteDifferencesRequested(configure())) {
    ParamIndexMap tmap = BuildParamIndexMap(lc, ld, fixedMask);
    std::vector<THMRows> thm = ComputeTHMRows(lc, ls, inputs_, ld, results, configure(), p, tmap);
    const int nCols = tmap.NumPacked();
    for (const THMRows &tr : thm)
      for (int col = 0; col < nCols; col++) {
        const int f = tmap.PackedToFull(col);
        const ParamKind kind = tmap.Desc(f).kind;
        if (!IsRMatrixParameter(kind) && kind != ParamKind::ThmCoherent) continue;
        double g = 0.0;
        for (size_t i = 0; i < tr.r.size(); i++) g += tr.r[i] * tr.J[i * nCols + col];
        if (f < (int)grad.size()) grad[f] += 2.0 * g;
      }
    thmAnalytic = true;
  }

  // Segment-local energy-shift gradient is exact only when a shift's effect is
  // confined to its own segments' data term.  Nuisance penalties handled by the
  // limits manager can couple a shift to the global chi^2, so disable the local
  // path when one is active.
  const bool localShiftOk = !HasNuisanceParameters();
  // Components taken by a finite difference of the whole chi2 (Chi2Value) already hold the prior-covariance term;
  // every other one (analytic, or a segment-local shift difference) gets its gradient added below.
  std::vector<char> fdWhole(p.size(), 0);

  // --- Finite differences only for what is left: non-fixed energy shifts, the
  //     THM-segment contribution to energies/gammas/norms, densities and stopping scales (the whole chi^2 with
  //     its priors; analytic with the kernels branch), and (if the analytic path bailed) the energy/gamma and norm blocks too. ---
  for (int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++) {
    if (idx < nMn && fixedMask[idx]) continue;
    ParamKind kind = pmap.Desc(idx).kind;
    // A THM coherent background (cbackground=) moves THM points only: the
    // THM-only chi^2 below differentiates it completely.  An electron-screening
    // potential is analytic (the screening factor of the classic points; a THM
    // segment carries none).
    if (eg && (kind == ParamKind::Screening || kind == ParamKind::Doppler || kind == ParamKind::Thickness))
      continue;  // analytic (dev: Doppler weights, the transmission / yield transform)
    if (eg && (IsRMatrixParameter(kind) || kind == ParamKind::Norm || kind == ParamKind::ThmCoherent)) {
      if (!haveTHM) continue;  // fully analytic
      if (thmAnalytic && kind != ParamKind::Norm) continue;  // THM part added above
      if (kind == ParamKind::Norm) {
        // A norm belongs to exactly one segment: non-THM norms are fully
        // analytic (handled above); THM norms are fully finite-difference.
        const Segment *seg = data()->Data().GetSegment(pmap.Desc(idx).segment);
        if (!seg || !seg->IsTHM()) continue;
      }
      double x0 = p[idx];
      double h = FdStep(kind, x0, 1.0e-6);
      vector_r pp = p;
      pp[idx] = x0 + h;
      vector_r pm = p;
      pm[idx] = x0 - h;
      // += : energies/gammas already carry the analytic non-THM part (a THM
      // norm's grad entry is still 0). Penalty terms inside the THM-only chi^2
      // that do not depend on p[idx] cancel in the central difference.
      grad[idx] += (Chi2Value(pp, true) - Chi2Value(pm, true)) / (2.0 * h);
      continue;
    }
    double x0 = p[idx];
    double h = FdStep(kind, x0, 1.0e-6);

    // Fast path: a free energy shift on self-contained segments only moves those
    // segments' points (its own, plus its shift group's followers for a group's
    // constant shift or a shift-function parameter), so finite-difference just
    // their chi^2 on the already-filled clone instead of re-solving the entire
    // dataset twice; with a step inside the planned nodes (ShiftedFold.h, question 5), else the whole chi^2.
    if ((kind == ParamKind::EnergyShift || kind == ParamKind::EnergyShiftSqrt || kind == ParamKind::EnergyShiftFn) &&
        localShiftOk && pmap.Desc(idx).segment >= 1) {
      const int s = pmap.Desc(idx).segment;
      const ShiftParamTarget target = ShiftTargets(ld, pmap.Desc(idx));
      // Not for a segment of a THM experiment (its chi^2 is not its own), nor of a covariance group (a group's
      // chi2 is not a sum over its segments).
      bool inThmGroup = false, inCovariance = false;
      for (int m : target.segs) {
        inThmGroup = inThmGroup || ld->Data().ThmGroupOf(m) >= 0;
        inCovariance = inCovariance || ld->Data().GetSegment(m)->GetCovarianceGroup() >= 0;
      }
      const auto move = [&](double v) {
        target.Set(ld, v);
        for (int m : target.segs) ld->Data().GetSegment(m)->UpdatePointEnergiesWithShift(*ld, lc, &configure());
      };
      if (const double step = !inCovariance && target.Local(ld) && !inThmGroup
                                  ? azure2::InsideShiftStep(*ld, target.segs, x0, h, move) : 0.0;
          step > 0.0) {
        const double d0 = target.Get(ld);
        const std::vector<double> held = azure2::HeldNorms(*ld, target.segs, lc, ls, configure(), results);
        auto groupChi = [&](double v) {
          move(v);  // every segment the parameter moves is moved before any is calculated
          double c = 0.0;
          for (std::size_t k = 0; k < target.segs.size(); k++)
            c += azure2::SegmentLocalChi2(*ld->Data().GetSegment(target.segs[k]), lc, ls, configure(), *ld, results, held[k]);
          return c;
        };
        double chiP = groupChi(x0 + step);
        double chiM = groupChi(x0 - step);
        move(d0);  // restore base state on the clone
        grad[idx] = (chiP - chiM) / (2.0 * step);

        // Optional self-check: compare against the exact full-dataset finite
        // difference for this shift.  Enable with AZURE_SHIFT_GRAD_CHECK=1 to
        // validate the fast path on real data (one gradient evaluation), then
        // disable it for the actual fit.
        if (azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kShiftGradCheck)) {
          vector_r pp = p;
          pp[idx] = x0 + step;
          vector_r pm = p;
          pm[idx] = x0 - step;
          double gFull = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * step);
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
    fdWhole[idx] = 1;
  }
  // Prior parameter covariance: 2 C^-1 (p - mu) on its parameters.
  if (const PriorCovariance *pc = ld->Data().GetPriorCovariance()) {
    std::vector<double> pg(p.size(), 0.0);
    pc->AddGradient(p, pg);
    for (size_t i = 0; i < grad.size() && i < pg.size(); i++)
      if (!fdWhole[i]) grad[i] += pg[i];
  }

  const auto tg2 = std::chrono::steady_clock::now();
  const azure2::support::Environment &environment = azure2::support::ProcessEnvironment();
  if (environment.Set(azure2::support::Switch::kGradCheck) || environment.Set(azure2::support::Switch::kGradTiming)) {
    int c = ++configure().diagnostics->gradientTimings;
    if (c <= 5 || c % 50 == 0) {
      auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
      int nfd = 0;
      for (int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++) {
        if (idx < nMn && fixedMask[idx]) continue;
        ParamKind kind = pmap.Desc(idx).kind;
        if (!(eg && !haveTHM && (IsRMatrixParameter(kind) || kind == ParamKind::Norm)) &&
            !(eg && (kind == ParamKind::Screening || kind == ParamKind::Doppler || kind == ParamKind::Thickness))) nfd++;
      }
      std::cerr << "[grad-timing] call " << c << ": setup+adjoint " << ms(tg0, tg1) << " ms (analytic " << (eg ? "ok" : "BAILED")
                << "), finite-difference part " << ms(tg1, tg2) << " ms for " << nfd << " parameters" << std::endl;
    }
  }

  // Optional full self-check, AZURE_GRAD_CHECK=1 (or =exit to stop the run after
  // it): on the first call, compare every non-fixed gradient component with a
  // central finite difference of the side-effect-free chi-squared and list the
  // worst disagreements.  Minuit itself never checks (CheckGradient() is false).
  if (environment.Set(azure2::support::Switch::kGradCheck)) {
    const std::string &gc = environment.Text(azure2::support::Switch::kGradCheck);
    if (!configure().diagnostics->gradientChecked.exchange(true)) {
      static const char *const kindName[] = {"energy", "gamma", "norm", "shift", "shift_sqrt", "shift_fn", "abundance", "screening",
                                       "density", "doppler", "thickness", "cbkg", "apparatus", "radiative",
                                       "isospin_shift", "isospin_mix"};
      struct Row { double score, g, fd; int idx; ParamDesc d; double gj, scoreJ; };
      std::vector<Row> rows;
      constexpr int kNumKinds = sizeof(kindName) / sizeof(kindName[0]);  // one per ParamKind
      int nBad[kNumKinds] = {}, nAll[kNumKinds] = {}, nBadJ[kNumKinds] = {};
      const bool haveJ = route.jacobianPathGradient && route.jacobianPathGradient->size() == grad.size();
      const auto tc0 = std::chrono::steady_clock::now();
      const double f0 = Chi2Value(p);
      std::cerr << "[grad-check] one Chi2Value: " << std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count() << " s" << std::endl;
      std::vector<int> freeIdx;
      for (int idx = 0; idx < (int)p.size() && idx < pmap.NumFull(); idx++)
        if (!(idx < nMn && fixedMask[idx])) freeIdx.push_back(idx);
      const int nCheck = static_cast<int>(environment.Integer(azure2::support::Switch::kGradCheckCount, 16));
      const int stride = std::max(1, (int)freeIdx.size() / std::max(1, nCheck));
      std::vector<int> checkIdx;
      for (int q = 0; q < (int)freeIdx.size(); q += stride) checkIdx.push_back(freeIdx[q]);
      // Targeted additions: AZURE_GRAD_CHECK_IDX = "i,j,a-b" (parameter
      // indices) and AZURE_GRAD_CHECK_SEGS = "s,a-b" (1-based segment numbers
      // in data order: their norm and energy-shift parameters).  Fixed
      // parameters are skipped.
      auto parseList = [](const std::string &spec) {
        std::vector<int> out;
        std::stringstream ss(spec);
        std::string item;
        while (std::getline(ss, item, ',')) {
          if (item.empty()) continue;
          size_t dash = item.find('-', 1);
          int lo = std::atoi(item.substr(0, dash).c_str());
          int hi = (dash == std::string::npos) ? lo : std::atoi(item.substr(dash + 1).c_str());
          for (int v = lo; v <= hi; v++) out.push_back(v);
        }
        return out;
      };
      std::vector<int> extra;
      if (environment.Set(azure2::support::Switch::kGradCheckIndices))
        for (int v : parseList(environment.Text(azure2::support::Switch::kGradCheckIndices))) extra.push_back(v);
      if (environment.Set(azure2::support::Switch::kGradCheckSegments))
        for (int sNum : parseList(environment.Text(azure2::support::Switch::kGradCheckSegments))) {
          if (sNum < 1 || sNum > ld->Data().NumSegments()) continue;
          extra.push_back(pmap.NormIndex(sNum));
          extra.push_back(pmap.EnergyShiftIndex(sNum));
          extra.push_back(pmap.EnergyShiftSqrtIndex(sNum));
          for (int j = 0; pmap.EnergyShiftFnIndex(sNum, j) >= 0; j++) extra.push_back(pmap.EnergyShiftFnIndex(sNum, j));
        }
      for (int v : extra)
        if (std::find(freeIdx.begin(), freeIdx.end(), v) != freeIdx.end() &&
            std::find(checkIdx.begin(), checkIdx.end(), v) == checkIdx.end())
          checkIdx.push_back(v);
      for (int idx : checkIdx) {
        // AZURE_GRAD_CHECK_H scales the finite-difference step (default 1e-6
        // relative): a component that moves with it is truncation error of the
        // FD (e.g. near a narrow level), not of the analytic value.
        const double hScale = environment.Number(azure2::support::Switch::kGradCheckStep, 1.0e-6);
        double x0 = p[idx], h = FdStep(pmap.Desc(idx).kind, x0, hScale);
        vector_r pp = p, pm = p;
        pp[idx] = x0 + h;
        pm[idx] = x0 - h;
        double fd = (Chi2Value(pp) - Chi2Value(pm)) / (2.0 * h);
        auto scoreOf = [&](double g) {
          return std::fabs(g - fd) / (std::fabs(fd) + std::fabs(g) + 1.0e-6 * std::fabs(f0) + 1.0);
        };
        double score = scoreOf(grad[idx]);
        double gj = haveJ ? (*route.jacobianPathGradient)[idx] : 0.0, scoreJ = haveJ ? scoreOf(gj) : 0.0;
        ParamDesc d = pmap.Desc(idx);
        int k = (int)d.kind;
        nAll[k]++;
        if (score > 1.0e-2) nBad[k]++;
        if (scoreJ > 1.0e-2) nBadJ[k]++;
        rows.push_back({score, grad[idx], fd, idx, d, gj, scoreJ});
      }
      std::sort(rows.begin(), rows.end(),
                [](const Row &a, const Row &b) { return std::max(a.score, a.scoreJ) > std::max(b.score, b.scoreJ); });
      std::cerr << "[grad-check] chi2 " << f0 << ", analytic block " << (eg ? "used" : "BAILED (all finite differences)")
                << "; components with relative disagreement > 1e-2:";
      for (int k = 0; k < kNumKinds; k++)
        if (nAll[k]) std::cerr << " " << kindName[k] << " " << nBad[k] << "/" << nAll[k];
      if (haveJ) {
        std::cerr << "; Jacobian path:";
        for (int k = 0; k < kNumKinds; k++)
          if (nAll[k]) std::cerr << " " << kindName[k] << " " << nBadJ[k] << "/" << nAll[k];
      }
      std::cerr << std::endl;
      // Full precision, so two runs (e.g. two builds, or two thread counts) can be
      // compared value by value; AZURE_GRAD_CHECK_ROWS lists more than 60.
      const std::streamsize prec = std::cerr.precision(15);
      const int maxRows = static_cast<int>(environment.Integer(azure2::support::Switch::kGradCheckRows, 60));
      for (int r = 0; r < (int)rows.size() && r < maxRows; r++) {
        const Row &w = rows[r];
        std::cerr << "[grad-check] idx " << w.idx << " " << kindName[(int)w.d.kind] << " J" << w.d.jGroup << " L" << w.d.level
                  << " ch" << w.d.channel << " seg" << w.d.segment << "  analytic " << w.g << "  fd " << w.fd
                  << "  rel " << w.score;
        if (haveJ) std::cerr << "  jacobian " << w.gj << "  relJ " << w.scoreJ;
        std::cerr << std::endl;
      }
      std::cerr.precision(prec);
      if (gc == "exit" && !route.jacobianPathGradient) ExitRun(configure(), 0);  // else the caller exits after it
    }
  }

  ReturnPooledState(std::move(pooled));
  ReturnPooledData(working);

  return grad;
}

std::vector<bool> AZURECalc::FixedMask() const {
  AZUREParams tp;
  compound()->FillMnParams(tp, state_, inputs_, &configure());
  data()->FillMnParams(tp);
  const int nMn = tp.GetMinuitParams().Params().size();
  if ((int)fixedMask_.size() == nMn) return fixedMask_;
  std::vector<bool> fixed(nMn);
  for (int i = 0; i < nMn; i++) fixed[i] = tp.GetMinuitParams().Parameter(i).IsFixed();
  return fixed;
}

bool AZURECalc::GradientFromJacobian(const std::vector<double> &p, std::vector<double> &grad) const {
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
  // Diagonal priors (x - x0)^2 / sigma^2 (Objective.h): gradient 2 (x - x0) / sigma^2, curvature 2 / sigma^2.
  ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), std::vector<bool>());
  ForEachDiagonalPrior(data(), pmap, [&](int idx, double nominal, double e) {
    if (idx >= (int)p.size()) return;
    grad[idx] += GaussianPriorGradient(p[idx], nominal, e);
    g2[idx] += GaussianPriorCurvature(e);
  });
  // Park's wall and the abundance sum-to-one term are residual rows (AppendModelRows): in J^T r above.
  states_.gradient.g2At = p;
  states_.gradient.g2 = g2;
  return true;
}

std::vector<double> AZURECalc::G2(const std::vector<double> &p) const {
  if (!states_.gradient.g2.empty() && states_.gradient.g2At == p) return states_.gradient.g2;
  // A stopped fit: a unit curvature, so that MIGRAD's step from the zero gradient is finite and zero.
  if (data()->IsFit() && StopRequested()) return std::vector<double>(p.size(), 1.0);
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
  // Penalty curvatures, (x - x0)^2 / sigma^2 -> 2 / sigma^2: the diagonal priors (Objective.h).  The abundance
  // sum-to-one term and Park's wall are residual rows (AppendModelRows), in the sum above.
  ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), std::vector<bool>());
  ForEachDiagonalPrior(data(), pmap, [&](int idx, double, double e) {
    if (idx < (int)g2.size()) g2[idx] += GaussianPriorCurvature(e);
  });
  if (azure2::support::ProcessEnvironment().Set(azure2::support::Switch::kGradTiming)) {
    int c = ++configure().diagnostics->g2Timings;
    if (c <= 5 || c % 50 == 0)
      std::cerr << "[g2-timing] call " << c << ": " << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count()
                << " ms" << std::endl;
  }
  return g2;
}

bool AZURECalc::EnableG2(const std::vector<double> &p) {
  g2Enabled_ = false;
  if (HasNuisanceParameters()) return false;  // nuisance penalties are not in the Jacobian
  std::vector<double> g2 = G2(p);
  g2Enabled_ = !g2.empty();
  return g2Enabled_;
}
