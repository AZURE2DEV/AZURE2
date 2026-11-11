// The calculation of one point: EPoint::Calculate (the R-, A- or THM-matrix route, the mapped points it fills, the
// probes of a target-effect point), the thin-target and observable-wise integrations of its sub-points, the marking
// of a point whose calculation failed, and the scope that switches point mapping off on a thread.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include "AMatrixFunc.h"
#include "CNuc.h"
#include "Config.h"
#include "EData.h"
#include "EPoint.h"
#include "RMatrixFunc.h"
#include "THMMatrixFunc.h"
#include "TargetEffect.h"

namespace {
thread_local bool tlPointMappingOff = false;
}

ScopedPointMappingOff::ScopedPointMappingOff() : previous_(tlPointMappingOff) { tlPointMappingOff = true; }
ScopedPointMappingOff::~ScopedPointMappingOff() { tlPointMappingOff = previous_; }
bool ScopedPointMappingOff::Active() { return tlPointMappingOff; }

void EPoint::MarkFailed(EData *data, azure2::Results &results, const char *role, const Config &configure,
                        const char *why) {
  // Every fixed quantity NaN and no angular-distribution coefficient, on this point and every point mapped onto it.
  const azure2::MappedPoints mappedPts = this->PlanMapped(data);
  for (int i = 0; i <= mappedPts.size(); i++) {
    EPoint *p = i ? mappedPts[i - 1] : this;
    std::fill(p->Rec(results), p->Rec(results) + azure2::Results::kNumFixed, std::numeric_limits<double>::quiet_NaN());
    if (p->num_ang_slots_) p->Rec(results)[azure2::Results::AngularOffset()] = 0.0;
  }
  static std::atomic<bool> reported(false);
  if (!reported.exchange(true)) {
#pragma omp critical(azure_failed_point)
    configure.outStream << "WARNING: the calculation failed at segment " << this->GetSegmentKey() << " (" << role
                        << "), E_lab " << this->GetLabEnergy() << " MeV, angle " << this->GetLabAngle() << " deg: " << why
                        << "  The point is NaN (further failures are not reported)." << std::endl;
  }
}

void EPoint::Calculate(EData *data, CNuc *theCNuc, const Config &configure, azure2::Results &results, EPoint *parent,
                       int subPointNum) {
  // Test hook (tests/pyazr/failed_point_test.py): AZURE_TEST_FAIL_E=<lab energy, MeV> makes the calculation of a
  // data or test point at that energy throw, AZURE_TEST_FAIL_EXIT=<pair key> only for that exit pair, so the
  // handling of a failed point (azure2::EvaluateModelPoint) can be exercised.
  static const double testFailE = [] { const char *v = std::getenv("AZURE_TEST_FAIL_E"); return v ? std::atof(v) : -1.0; }();
  static const int testFailExit = [] { const char *v = std::getenv("AZURE_TEST_FAIL_EXIT"); return v ? std::atoi(v) : 0; }();
  if (testFailE >= 0.0 && !parent && std::fabs(this->GetLabEnergy() - testFailE) <= 1e-9 * std::max(1.0, testFailE) &&
      (testFailExit == 0 || this->GetExitKey() == testFailExit))
    throw std::runtime_error("injected failure (AZURE_TEST_FAIL_E).");
  // Points mapped onto this one (EData::MapData), unless mapping is switched
  // off on this thread (ScopedPointMappingOff).
  const azure2::MappedPoints mappedPts = ScopedPointMappingOff::Active() ? azure2::MappedPoints() : this->PlanMapped(data);
  const int numMapped = mappedPts.size();
  if (!this->IsTargetEffect() || !data->GetTargetEffect(this->GetTargetEffectNum())->IsSubPointEffect()) {
    // THM (modified R-matrix) segments: reuse the shared interior (A-matrix)
    // but assemble the half-off-shell cross section instead of the standard
    // T-matrix cross section. Kept fully per-segment so a project may mix THM
    // and conventional data against the same levels.
    if (this->IsTHM()) {
      // Reused across the points of this thread, as the A-matrix function
      // below: its per-J-group buffers are sized once, then only refilled
      // (Reset repoints it at this compound and configuration first, so
      // nothing of an earlier session is read).
      static thread_local THMMatrixFunc reusableThmFunc(theCNuc, configure);
      reusableThmFunc.Reset(theCNuc, configure);
      THMMatrixFunc &thmFunc = reusableThmFunc;
      thmFunc.ClearMatrices();
      thmFunc.FillMatrices(this);
      thmFunc.InvertMatrices();
      thmFunc.CalculateTHMCrossSection(this, results);
      if (subPointNum && parent) {
        const azure2::MappedPoints parentMapped = parent->PlanMapped(data);
        for (int i = 1; i <= parentMapped.size(); i++) {
          EPoint *mappedSubPoint = parentMapped[i - 1]->GetSubPoint(*data, subPointNum);
          thmFunc.CalculateTHMCrossSection(mappedSubPoint, results);
        }
      } else {
        const azure2::MappedPoints ownMapped = this->PlanMapped(data);
        for (int i = 1; i <= ownMapped.size(); i++) thmFunc.CalculateTHMCrossSection(ownMapped[i - 1], results);
      }
      return;
    }

    GenMatrixFunc *theMatrixFunc;
    // The A-matrix function object is reused across the points calculated by
    // this thread.  Its per-JGroup buffers and the level-matrix factorization
    // workspace are sized once and then only refilled, which removes the
    // allocation traffic that dominated the old per-point construction.
    // Recursion is not a concern: the branch that recurses into sub-points does
    // not use theMatrixFunc.
    static thread_local AMatrixFunc reusableAMatrixFunc(theCNuc, configure);
    std::unique_ptr<RMatrixFunc> ownedRMatrixFunc;
    if (configure.paramMask & Config::USE_AMATRIX) {
      reusableAMatrixFunc.Reset(theCNuc, configure);
      theMatrixFunc = &reusableAMatrixFunc;
    } else {
      ownedRMatrixFunc.reset(new RMatrixFunc(theCNuc, configure));
      theMatrixFunc = ownedRMatrixFunc.get();
    }
    theMatrixFunc->ClearMatrices();
    theMatrixFunc->FillMatrices(this);
    theMatrixFunc->InvertMatrices();
    theMatrixFunc->CalculateTMatrix(this);
    theMatrixFunc->CalculateCrossSection(this, results);
    // A mapped point may have another exit pair (a component of a SUM, RATIO or total-capture segment, EData::MapData):
    // the solved level matrix is the same, only its T-matrix is built anew.  Mapped points are ordered by exit, so
    // each exit's T-matrix is built once.
    int tExitKey = this->GetExitKey();
    auto mappedCrossSection = [&](EPoint *mapped) {
      if (mapped->GetExitKey() != tExitKey) {
        theMatrixFunc->RecalculateTMatrix(mapped);
        tExitKey = mapped->GetExitKey();
      }
      theMatrixFunc->CalculateCrossSection(mapped, results);
    };
    if (subPointNum && parent) {
      const azure2::MappedPoints parentMapped =
          ScopedPointMappingOff::Active() ? azure2::MappedPoints() : parent->PlanMapped(data);
      for (int i = 1; i <= parentMapped.size(); i++) mappedCrossSection(parentMapped[i - 1]->GetSubPoint(*data, subPointNum));
    } else {
      for (int i = 1; i <= numMapped; i++) mappedCrossSection(mappedPts[i - 1]);
    }
  } else {
    TargetEffect *effect = data->GetTargetEffect(this->GetTargetEffectNum());
    if (effect->IsTargetIntegration() && !effect->IsThickTarget()) {
      this->CalculateThinTarget(*data, theCNuc, configure, results, effect);
      return;
    }
    const double blend = this->GetTargetBlendWeight();
    const double autoTol = effect->GetAutoTolerance();
    const bool hasMapped = numMapped > 0;

    // The bare value of the observable at the central energy is needed both
    // for edge blending and for the automatic per-point decision.  It is
    // obtained from a probe copy that pretends not to carry the effect, so the
    // ordinary single-point machinery above stays the single source of truth.
    double bareValue = 0.0, bareE1 = 0.0, bareE2 = 0.0;
    bool haveBare = false;
    if (blend < 1.0 || autoTol > 0.) {
      EPoint probe(*this);
      azure2::Results probeResults;  // the probe's own record: it must not write the point's
      probe.UsePrivateRecord(probeResults, results);
      probe.SetTargetEffectNum(0);
      probe.DetachFromPlan();
      probe.SetSubPointRange(0, 0);
      probe.Calculate(data, theCNuc, configure, probeResults);
      bareValue = probe.GetFitCrossSection(probeResults);
      bareE1 = probe.GetFitE1CrossSection(probeResults), bareE2 = probe.GetFitE2CrossSection(probeResults);
      haveBare = true;
    }

    // Automatic application: probe the first, central and last sub-point and
    // estimate the size of the effect from the curvature across the window
    // (plus the centroid shift when the target thickness is integrated).  If
    // the effect would change the observable by less than the tolerance, the
    // bare value is used and the remaining sub-points are never evaluated, so
    // any discontinuity this introduces is bounded by the tolerance itself.
    // Points with mapped observables are skipped only when every one of them
    // passes its own estimate against its own bare probe.
    if (autoTol > 0. && haveBare && this->NumSubPoints() >= 3) {
      int iFirst = 1;
      int iLast = this->NumSubPoints();
      int iMid = 1;
      double dBest = 1e100;
      for (int i = 1; i <= this->NumSubPoints(); i++) {
        double d = std::fabs(this->GetSubPoint(*data, i)->GetCMEnergy() - this->GetCMEnergy());
        if (d < dBest) { dBest = d; iMid = i; }
      }
      const int probes[3] = {iFirst, iMid, iLast};
      for (int k = 0; k < 3; k++) {
        EPoint *subPoint = this->GetSubPoint(*data, probes[k]);
        if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE)
          subPoint->RecalcEDependentValues(data, theCNuc, configure);
        if (hasMapped)
          subPoint->Calculate(data, theCNuc, configure, results, this, probes[k]);
        else
          subPoint->Calculate(data, theCNuc, configure, results);
      }
      bool allPass = true;
      for (int m = 0; m <= numMapped && allPass; m++) {
        EPoint *host = (m == 0) ? this : mappedPts[m - 1];
        double sFirst = host->GetSubPoint(*data, iFirst)->GetFitCrossSection(results);
        double sMid = host->GetSubPoint(*data, iMid)->GetFitCrossSection(results);
        double sLast = host->GetSubPoint(*data, iLast)->GetFitCrossSection(results);
        double estimate = std::fabs(0.5 * (sFirst + sLast) - sMid);
        if (effect->IsTargetIntegration()) estimate += 0.5 * std::fabs(sLast - sFirst);
        if (estimate > autoTol * std::max(std::fabs(sMid), 1e-300)) allPass = false;
      }
      if (allPass) {
        // the bare point's sigma and its E1 and E2 parts (what an E1 / E2 line reads)
        this->SetFitCrossSection(results, bareValue);
        this->SetFitE1CrossSection(results, bareE1), this->SetFitE2CrossSection(results, bareE2);
        for (int m = 1; m <= numMapped; m++) {
          EPoint *mappedPoint = mappedPts[m - 1];
          EPoint mappedProbe(*mappedPoint);
          azure2::Results probeResults;
          mappedProbe.UsePrivateRecord(probeResults, results);
          mappedProbe.SetTargetEffectNum(0);
          mappedProbe.DetachFromPlan();
          mappedProbe.SetSubPointRange(0, 0);
          mappedProbe.Calculate(data, theCNuc, configure, probeResults);
          mappedPoint->SetFitCrossSection(results, mappedProbe.GetFitCrossSection(probeResults));
          mappedPoint->SetFitE1CrossSection(results, mappedProbe.GetFitE1CrossSection(probeResults));
          mappedPoint->SetFitE2CrossSection(results, mappedProbe.GetFitE2CrossSection(probeResults));
        }
        return;
      }
    }

    for (int i = 1; i <= this->NumSubPoints(); i++) {
      EPoint *subPoint = this->GetSubPoint(*data, i);
      // Recalculate energy-dependent values for sub-points too if using external capture
      if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
        subPoint->RecalcEDependentValues(data, theCNuc, configure);
      }
      if (numMapped > 0)
        subPoint->Calculate(data, theCNuc, configure, results, this, i);
      else
        subPoint->Calculate(data, theCNuc, configure, results);
    }
    this->IntegrateTargetEffectForObservable(*data, configure, results);
    if (blend < 1.0 && haveBare)
      this->SetFitCrossSection(results, blend * this->GetFitCrossSection(results) + (1.0 - blend) * bareValue);
    for (int i = 1; i <= numMapped; i++) {
      EPoint *mappedPoint = mappedPts[i - 1];
      mappedPoint->IntegrateTargetEffectForObservable(*data, configure, results);
      if (blend < 1.0) {
        // The mapped point is another observable at the same energy: its bare
        // value needs its own probe.
        EPoint mappedProbe(*mappedPoint);
        azure2::Results probeResults;
        mappedProbe.UsePrivateRecord(probeResults, results);
        mappedProbe.SetTargetEffectNum(0);
        mappedProbe.DetachFromPlan();
        mappedProbe.SetSubPointRange(0, 0);
        mappedProbe.Calculate(data, theCNuc, configure, probeResults);
        mappedPoint->SetFitCrossSection(results, blend * mappedPoint->GetFitCrossSection(results) +
                                                     (1.0 - blend) * mappedProbe.GetFitCrossSection(probeResults));
      }
    }
  }
}

/*!
 * The sub-point stage of a thin-target integration.  Unlike the other kernels it falls back to sigma at the
 * depth-averaged energy E_bar = E - D (SetTargetFallbackOnly), never to a plain point at the data energy: the data
 * of a target integration carry the beam energy, and sigma(E_beam) is off by the centroid shift ~Delta/2 (13C+a
 * ND2021r with windows around the narrow levels: chi2 420,221 against 402,111 integrated everywhere, 2026-11-05).
 *  - outside every window: the one to three fallback sub-points, then the kernel (which returns sigma(E_bar));
 *  - at a window edge: the full grid, blended with sigma(E_bar) inside the kernel;
 *  - automatic tolerance: the first and last sub-points and the two around E_bar first; if the chord between
 *    the end nodes stays within the tolerance of sigma(E_bar) (a curvature bound, about three times the error of
 *    sigma(E_bar) for a box), sigma(E_bar) is used and no other sub-point is evaluated.
 */
void EPoint::CalculateThinTarget(EData &data, CNuc *theCNuc, const Config &configure, azure2::Results &results,
                                 TargetEffect *effect) {
  const azure2::MappedPoints mappedPts =
      ScopedPointMappingOff::Active() ? azure2::MappedPoints() : this->PlanMapped(&data);
  const int numMapped = mappedPts.size();
  const bool hasMapped = numMapped > 0;
  const int n = this->NumSubPoints();
  auto calcSub = [&](int i) {
    EPoint *subPoint = this->GetSubPoint(data, i);
    if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) subPoint->RecalcEDependentValues(&data, theCNuc, configure);
    if (hasMapped)
      subPoint->Calculate(&data, theCNuc, configure, results, this, i);
    else
      subPoint->Calculate(&data, theCNuc, configure, results);
  };
  auto integrateAll = [&]() {
    this->IntegrateTargetEffectForObservable(data, configure, results);
    for (int m = 1; m <= numMapped; m++)
      mappedPts[m - 1]->IntegrateTargetEffectForObservable(data, configure, results);
  };
  const double autoTol = effect->GetAutoTolerance();
  if (autoTol > 0. && !target_fallback_only_ && n >= 4) {
    const double eBar = this->GetCMEnergy() - target_fallback_depth_;
    int jb = 1;  // the sub-points jb, jb + 1 bracket E_bar (descending energies)
    for (int j = 1; j < n; j++)
      if (this->GetSubPoint(data, j + 1)->GetCMEnergy() <= eBar) { jb = j; break; }
    if (eBar <= this->GetSubPoint(data, n)->GetCMEnergy()) jb = n - 1;
    const int probes[4] = {1, jb, jb + 1, n};
    std::vector<bool> done(n + 1, false);
    for (int k = 0; k < 4; k++)
      if (!done[probes[k]]) { calcSub(probes[k]); done[probes[k]] = true; }
    const double eF = this->GetSubPoint(data, 1)->GetCMEnergy(), eL = this->GetSubPoint(data, n)->GetCMEnergy();
    const double ea = this->GetSubPoint(data, jb)->GetCMEnergy(), eb = this->GetSubPoint(data, jb + 1)->GetCMEnergy();
    const double tChord = (eF > eL) ? (eF - eBar) / (eF - eL) : 0.0;
    const double tBar = (ea > eb) ? std::min(std::max((ea - eBar) / (ea - eb), 0.0), 1.0) : 0.0;
    bool allPass = true;
    for (int m = 0; m <= numMapped && allPass; m++) {
      EPoint *host = (m == 0) ? this : mappedPts[m - 1];
      const double sF = host->GetSubPoint(data, 1)->GetFitCrossSection(results);
      const double sL = host->GetSubPoint(data, n)->GetFitCrossSection(results);
      const double sBar = (1.0 - tBar) * host->GetSubPoint(data, jb)->GetFitCrossSection(results) +
                          tBar * host->GetSubPoint(data, jb + 1)->GetFitCrossSection(results);
      const double chord = (1.0 - tChord) * sF + tChord * sL;
      if (std::fabs(chord - sBar) > autoTol * std::max(std::fabs(sBar), 1e-300)) allPass = false;
    }
    if (allPass) {
      // sigma(E_bar) from the two bracketing sub-points, through the one integrator (A_y, coefficients and
      // capture components keep their own combinations)
      EPoint::TargetFallbackScope self(this);
      std::vector<std::unique_ptr<EPoint::TargetFallbackScope>> mapped;
      for (int m = 1; m <= numMapped; m++)
        mapped.emplace_back(new EPoint::TargetFallbackScope(mappedPts[m - 1]));
      integrateAll();
      return;
    }
    for (int i = 1; i <= n; i++)
      if (!done[i]) calcSub(i);
    integrateAll();
    return;
  }
  for (int i = 1; i <= n; i++) calcSub(i);
  integrateAll();
}

/*!
 * A segment compared against the E1 or E2 component of an angle-integrated
 * capture cross section (isDiff 5/6) reads that component off the parent
 * point, which a sub-point integration never filled: only the total was
 * integrated.  The kernel is linear in the sub-point values, so each
 * component is integrated with the same combiner by swapping it into the
 * sub-points' cross-section slot, exactly as the analyzing power is handled.
 */
void EPoint::IntegrateTargetEffectComponents(const EData &data, const Config &configure, azure2::Results &results) {
  if (settings_->crossSectionComponent == 0) return;
  const int n = this->NumSubPoints();
  if (n <= 0) return;
  std::vector<double> sigma(n);
  for (int i = 1; i <= n; i++) sigma[i - 1] = this->GetSubPoint(data, i)->GetFitCrossSection(results);
  const double total = this->GetFitCrossSection(results);

  for (int i = 1; i <= n; i++)
    this->GetSubPoint(data, i)->SetFitCrossSection(results, this->GetSubPoint(data, i)->GetFitE1CrossSection(results));
  this->IntegrateTargetEffect(data, configure, results);
  this->SetFitE1CrossSection(results, this->GetFitCrossSection(results));

  for (int i = 1; i <= n; i++)
    this->GetSubPoint(data, i)->SetFitCrossSection(results, this->GetSubPoint(data, i)->GetFitE2CrossSection(results));
  this->IntegrateTargetEffect(data, configure, results);
  this->SetFitE2CrossSection(results, this->GetFitCrossSection(results));

  for (int i = 1; i <= n; i++) this->GetSubPoint(data, i)->SetFitCrossSection(results, sigma[i - 1]);
  this->SetFitCrossSection(results, total);
}

/*!
 * Integrates the observable over the target, which is not the same operation
 * for an analyzing power as for a cross section.
 *
 * A_y is a ratio, so averaging it over the target thickness with equal weights
 * is wrong: the parts of the target where the cross section is large must
 * count for more. The physically meaningful quantity is the ratio of the
 * integrated polarized and unpolarized yields,
 *
 *   <A_y> = Int A_y(E) sigma(E) dE / Int sigma(E) dE
 *
 * which is what a measurement actually returns. Both integrals are done with
 * the existing yield integrator -- once on sigma, once on the product -- so the
 * quadrature, straggling and energy-loss treatment stay identical and there is
 * only one integration scheme to maintain.
 */
void EPoint::IntegrateTargetEffectForObservable(const EData &data, const Config &configure, azure2::Results &results) {
  if (this->IsCoefficient() && !this->GetCoefficientAbsolute()) {
    // A normalised Legendre coefficient is a ratio like A_y: the absolute A_k = a_k A_0
    // and A_0 are each folded with the integrator (linear in the cross section), then
    // divided -- averaging the ratio itself would be wrong.  An absolute A_k is a
    // cross-section-like quantity and takes the plain path (its slot holds a_k A_0).
    const int n = this->NumSubPoints();
    std::vector<double> ak(n), a0(n);
    for (int i = 1; i <= n; i++) {
      ak[i - 1] = this->GetSubPoint(data, i)->GetFitCrossSection(results);
      a0[i - 1] = this->GetSubPoint(data, i)->GetAngDistA0(results);
    }
    for (int i = 1; i <= n; i++) this->GetSubPoint(data, i)->SetFitCrossSection(results, a0[i - 1]);
    this->IntegrateTargetEffect(data, configure, results);
    const double denominator = this->GetFitCrossSection(results);
    for (int i = 1; i <= n; i++) this->GetSubPoint(data, i)->SetFitCrossSection(results, ak[i - 1] * a0[i - 1]);
    this->IntegrateTargetEffect(data, configure, results);
    const double numerator = this->GetFitCrossSection(results);
    for (int i = 1; i <= n; i++) this->GetSubPoint(data, i)->SetFitCrossSection(results, ak[i - 1]);
    this->SetFitCrossSection(results, std::fabs(denominator) > 0.0 ? numerator / denominator : 0.0);
    return;
  }
  if (!this->IsAnalyzingPower()) {
    this->IntegrateTargetEffect(data, configure, results);
    this->IntegrateTargetEffectComponents(data, configure, results);
    return;
  }

  const int n = this->NumSubPoints();
  std::vector<double> sigma(n);
  for (int i = 1; i <= n; i++) sigma[i - 1] = this->GetSubPoint(data, i)->GetFitCrossSection(results);

  this->IntegrateTargetEffect(data, configure, results);
  const double denominator = this->GetFitCrossSection(results);

  for (int i = 1; i <= n; i++)
    this->GetSubPoint(data, i)->SetFitCrossSection(results, sigma[i - 1] * this->GetSubPoint(data, i)->GetAnalyzingPower(results));
  this->IntegrateTargetEffect(data, configure, results);
  const double numerator = this->GetFitCrossSection(results);

  for (int i = 1; i <= n; i++) this->GetSubPoint(data, i)->SetFitCrossSection(results, sigma[i - 1]);

  this->SetFitCrossSection(results, std::fabs(denominator) > 0.0 ? numerator / denominator : 0.0);
}
