// The calculation of one point outside the energy pass (Unit::Calculate: a host with the points mapped onto it, or a
// point alone), the folds of a point's sub-points after they are calculated (Unit::FoldSubPoints: the probes, the
// automatic tolerance, the thin-target fallback, the observable-wise integrations), the marking of a point whose
// calculation failed, and the scope that switches point mapping off on a thread.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include "CNuc.h"
#include "CalculationPlan.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "Unit.h"
#include "TargetEffect.h"

namespace {
thread_local bool tlPointMappingOff = false;
}

ScopedPointMappingOff::ScopedPointMappingOff() : previous_(tlPointMappingOff) { tlPointMappingOff = true; }
ScopedPointMappingOff::~ScopedPointMappingOff() { tlPointMappingOff = previous_; }
bool ScopedPointMappingOff::Active() { return tlPointMappingOff; }

void Unit::MarkFailed(WorkingCopy *data, azure2::Results &results, const char *role, const Config &configure,
                        const char *why) {
  // Every fixed quantity NaN and no angular-distribution coefficient, on this point and every point mapped onto it.
  const azure2::MappedPoints mappedPts = this->PlanMapped(data);
  for (int i = 0; i <= mappedPts.size(); i++) {
    Unit *p = i ? mappedPts[i - 1] : this;
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

bool Unit::InjectedFailure() const {
  // Test hook (tests/pyazr/failed_point_test.py): AZURE_TEST_FAIL_E=<lab energy, MeV> makes the calculation of a
  // data or test point at that energy fail, AZURE_TEST_FAIL_EXIT=<pair key> only for that exit pair, so the
  // handling of a failed point (azure2::EvaluateModelPoint, WorkingCopy::CalculatePointsInParallel) can be exercised.
  static const double testFailE = [] { const char *v = std::getenv("AZURE_TEST_FAIL_E"); return v ? std::atof(v) : -1.0; }();
  static const int testFailExit = [] { const char *v = std::getenv("AZURE_TEST_FAIL_EXIT"); return v ? std::atoi(v) : 0; }();
  return testFailE >= 0.0 && std::fabs(this->GetLabEnergy() - testFailE) <= 1e-9 * std::max(1.0, testFailE) &&
         (testFailExit == 0 || this->GetExitKey() == testFailExit);
}

void Unit::Calculate(WorkingCopy *data, CNuc *theCNuc, const Config &configure, azure2::Results &results) {
  if (this->InjectedFailure()) throw std::runtime_error("injected failure (AZURE_TEST_FAIL_E).");
  // Points mapped onto this one (WorkingCopy::MapData): they share its energy, and each of its sub-points' energies,
  // unless mapping is switched off on this thread (ScopedPointMappingOff).
  const azure2::MappedPoints mappedPts = ScopedPointMappingOff::Active() ? azure2::MappedPoints() : this->PlanMapped(data);
  const int numMapped = mappedPts.size();
  std::vector<Unit *> units(1 + numMapped);
  if (!this->IsTargetEffect() || !azure2::HasSubPointFold(*data, *this)) {
    units[0] = this;
    for (int i = 1; i <= numMapped; i++) units[i] = mappedPts[i - 1];
    azure2::SolveAtOneEnergy(theCNuc, configure, results, units.data(), (int)units.size());
    return;
  }
  // One energy per sub-point: this point's sub-point i and sub-point i of each mapped point; then the folds.
  this->PrepareSubPoints(data, theCNuc, configure);
  for (int i = 1; i <= this->NumSubPoints(); i++) {
    units[0] = this->GetSubPoint(*data, i);
    for (int m = 1; m <= numMapped; m++) units[m] = mappedPts[m - 1]->GetSubPoint(*data, i);
    azure2::SolveAtOneEnergy(theCNuc, configure, results, units.data(), (int)units.size());
  }
  this->FoldSubPoints(data, theCNuc, configure, results);
}

void Unit::PrepareSubPoints(WorkingCopy *data, CNuc *theCNuc, const Config &configure) {
  for (int i = 1; i <= this->NumSubPoints(); i++) {
    Unit *subPoint = this->GetSubPoint(*data, i);
    // Recalculate energy-dependent values for sub-points too if using external capture
    if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
      subPoint->RecalcEDependentValues(data, theCNuc, configure);
    }
  }
}

void Unit::FoldSubPoints(WorkingCopy *data, CNuc *theCNuc, const Config &configure, azure2::Results &results) {
  const azure2::MappedPoints mappedPts = ScopedPointMappingOff::Active() ? azure2::MappedPoints() : this->PlanMapped(data);
  const int numMapped = mappedPts.size();
  const TargetEffect *effect = data->Data().GetTargetEffect(this->GetTargetEffectNum());
  if (effect->IsTargetIntegration() && !effect->IsThickTarget()) {
    this->FoldThinTarget(*data, results, configure, effect, mappedPts);
    return;
  }
  const double blend = this->GetTargetBlendWeight();
  const double autoTol = effect->GetAutoTolerance();

  // The bare value of the observable at the central energy is needed both
  // for edge blending and for the automatic per-point decision.  It is
  // obtained from a probe copy that pretends not to carry the effect, so the
  // ordinary single-point machinery stays the single source of truth.
  double bareValue = 0.0, bareE1 = 0.0, bareE2 = 0.0;
  bool haveBare = false;
  if (blend < 1.0 || autoTol > 0.) {
    Unit probe(*this);
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

  // Automatic application: from the first, central and last sub-point,
  // estimate the size of the effect from the curvature across the window
  // (plus the centroid shift when the target thickness is integrated).  If
  // the effect would change the observable by less than the tolerance, the
  // bare value is used, so any discontinuity this introduces is bounded by
  // the tolerance itself.  Points with mapped observables take the bare value
  // only when every one of them passes its own estimate against its own bare
  // probe.  (The sub-points are all calculated in the energy pass; the
  // decision reads three of them.)
  if (autoTol > 0. && haveBare && this->NumSubPoints() >= 3) {
    int iFirst = 1;
    int iLast = this->NumSubPoints();
    int iMid = 1;
    double dBest = 1e100;
    for (int i = 1; i <= this->NumSubPoints(); i++) {
      double d = std::fabs(this->GetSubPoint(*data, i)->GetCMEnergy() - this->GetCMEnergy());
      if (d < dBest) { dBest = d; iMid = i; }
    }
    bool allPass = true;
    for (int m = 0; m <= numMapped && allPass; m++) {
      Unit *host = (m == 0) ? this : mappedPts[m - 1];
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
        Unit *mappedPoint = mappedPts[m - 1];
        Unit mappedProbe(*mappedPoint);
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

  this->IntegrateTargetEffectForObservable(*data, configure, results);
  if (blend < 1.0 && haveBare)
    this->SetFitCrossSection(results, blend * this->GetFitCrossSection(results) + (1.0 - blend) * bareValue);
  for (int i = 1; i <= numMapped; i++) {
    Unit *mappedPoint = mappedPts[i - 1];
    mappedPoint->IntegrateTargetEffectForObservable(*data, configure, results);
    if (blend < 1.0) {
      // The mapped point is another observable at the same energy: its bare
      // value needs its own probe.
      Unit mappedProbe(*mappedPoint);
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

/*!
 * The fold of a thin-target integration, with its sigma(E_bar) fallback.  Unlike the other kernels it falls back to
 * sigma at the depth-averaged energy E_bar = E - D (SetTargetFallbackOnly), never to a plain point at the data energy:
 * the data of a target integration carry the beam energy, and sigma(E_bar) is off by the centroid shift ~Delta/2
 * (13C+a ND2021r with windows around the narrow levels: chi2 420,221 against 402,111 integrated everywhere,
 * 2026-11-05).
 *  - outside every window: the one to three fallback sub-points, then the kernel (which returns sigma(E_bar));
 *  - at a window edge: the full grid, blended with sigma(E_bar) inside the kernel;
 *  - automatic tolerance: from the first and last sub-points and the two around E_bar, if the chord between the end
 *    nodes stays within the tolerance of sigma(E_bar) (a curvature bound, about three times the error of
 *    sigma(E_bar) for a box), sigma(E_bar) is used.
 */
void Unit::FoldThinTarget(WorkingCopy &data, azure2::Results &results, const Config &configure,
                          const TargetEffect *effect, const azure2::MappedPoints &mappedPts) {
  const int numMapped = mappedPts.size();
  const int n = this->NumSubPoints();
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
    const double eF = this->GetSubPoint(data, 1)->GetCMEnergy(), eL = this->GetSubPoint(data, n)->GetCMEnergy();
    const double ea = this->GetSubPoint(data, jb)->GetCMEnergy(), eb = this->GetSubPoint(data, jb + 1)->GetCMEnergy();
    const double tChord = (eF > eL) ? (eF - eBar) / (eF - eL) : 0.0;
    const double tBar = (ea > eb) ? std::min(std::max((ea - eBar) / (ea - eb), 0.0), 1.0) : 0.0;
    bool allPass = true;
    for (int m = 0; m <= numMapped && allPass; m++) {
      Unit *host = (m == 0) ? this : mappedPts[m - 1];
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
      Unit::TargetFallbackScope self(this);
      std::vector<std::unique_ptr<Unit::TargetFallbackScope>> mapped;
      for (int m = 1; m <= numMapped; m++)
        mapped.emplace_back(new Unit::TargetFallbackScope(mappedPts[m - 1]));
      integrateAll();
      return;
    }
  }
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
void Unit::IntegrateTargetEffectComponents(const WorkingCopy &data, const Config &configure, azure2::Results &results) {
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
void Unit::IntegrateTargetEffectForObservable(const WorkingCopy &data, const Config &configure, azure2::Results &results) {
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
