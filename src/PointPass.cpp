// The calculation of one point outside the energy pass (Unit::Calculate), the fold of a point's sub-points after
// they are calculated (Unit::FoldSubPoints, observable by observable), and the marking of a point whose calculation
// failed.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "azure2/support/Environment.h"
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "CNuc.h"
#include "CalculationPlan.h"
#include "Config.h"
#include "PointFold.h"
#include "WorkingCopy.h"
#include "Unit.h"

void Site::MarkFailed(azure2::Results &results, const char *role, const Config &configure, const char *why) {
  // Every fixed quantity NaN; the observables of the angular part NaN (PoisonObservable), so that every point reading them is NaN (ValueFor), and an
  // angular distribution's model record gets none (CopyAngularDistsToModel).
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::fill(Rec(results), Rec(results) + azure2::Results::kNumFixed, nan);
  this->PoisonObservable(results);
  if (configure.FirstTime("failed point")) {
#pragma omp critical(azure_failed_point)
    configure.outStream << "WARNING: the calculation failed at segment " << this->GetSegmentKey() << " (" << role
                        << "), E_lab " << this->GetLabEnergy() << " MeV, angle " << this->GetLabAngle() << " deg: " << why
                        << "  The point is NaN (further failures are not reported)." << std::endl;
  }
}

bool Site::InjectedFailure() const {
  // Test hook (tests/pyazr/failed_point_test.py): AZURE_TEST_FAIL_E=<lab energy, MeV> makes the calculation of a
  // data or test point at that energy fail, AZURE_TEST_FAIL_EXIT=<pair key> only for that exit pair, so the
  // handling of a failed point (azure2::EvaluateModelPoint, WorkingCopy::CalculatePointsInParallel) can be exercised.
  static const double testFailE = azure2::support::ProcessEnvironment().Number(azure2::support::Switch::kTestFailEnergy, -1.0);
  static const int testFailExit = static_cast<int>(azure2::support::ProcessEnvironment().Integer(azure2::support::Switch::kTestFailExit, 0));
  return testFailE >= 0.0 && std::fabs(this->GetLabEnergy() - testFailE) <= 1e-9 * std::max(1.0, testFailE) &&
         (testFailExit == 0 || this->GetExitKey() == testFailExit);
}

void Unit::Calculate(WorkingCopy *data, CNuc *theCNuc, const Config &configure, azure2::Results &results) {
  if (this->InjectedFailure()) throw std::runtime_error("injected failure (AZURE_TEST_FAIL_E).");
  Site *unit = this;
  if (!this->IsTargetEffect() || !azure2::HasSubPointFold(*data, *this)) {
    azure2::SolveAtOneEnergy(theCNuc, data->ThmCoherentValues(), configure, results, &unit, 1);
    return;
  }
  // One energy per sub-point; then the fold.
  this->PrepareSubPoints(data, theCNuc, configure);
  for (int i = 1; i <= this->NumSubPoints(); i++) {
    unit = this->GetSubPoint(*data, i);
    azure2::SolveAtOneEnergy(theCNuc, data->ThmCoherentValues(), configure, results, &unit, 1);
  }
  data->ContractRowSetOf(*this, results);
  this->FoldSubPoints(*data, configure, results);
}

void Unit::PrepareSubPoints(WorkingCopy *data, CNuc *theCNuc, const Config &configure) {
  for (int i = 1; i <= this->NumSubPoints(); i++) {
    if (!data->WritesSlot(*this, i)) continue;  // a slot of the point's group that another point writes
    Site *subPoint = this->GetSubPoint(*data, i);
    // Recalculate energy-dependent values for sub-points too if using external capture
    if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
      subPoint->RecalcEDependentValues(data, theCNuc, configure);
    }
  }
}

/*!
 * A segment compared against the E1 or E2 component of an angle-integrated
 * capture cross section (isDiff 5/6) reads that component off the parent
 * point, which a sub-point integration never filled: only the total was
 * integrated.  The kernel is linear in the sub-point values, so each
 * component is integrated with the same combiner (FoldValues), exactly as the
 * analyzing power is handled.
 */
void Unit::IntegrateTargetEffectComponents(const WorkingCopy &data, const Config &configure, azure2::Results &results) {
  if (settings_->crossSectionComponent == 0) return;
  const int n = this->NumSubPoints();
  if (n <= 0) return;
  std::vector<double> e1(n), e2(n);
  for (int i = 1; i <= n; i++) {
    e1[i - 1] = this->GetSubPoint(data, i)->GetFitE1CrossSection(results);
    e2[i - 1] = this->GetSubPoint(data, i)->GetFitE2CrossSection(results);
  }
  this->SetFitE1CrossSection(results, FoldValues(*this, data, configure, e1));
  this->SetFitE2CrossSection(results, FoldValues(*this, data, configure, e2));
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
 * only one integration scheme to maintain.  The slots' records are read, never
 * written: other points of the group read them too (SlotGroup.h).
 */
bool Unit::FoldsRatio() const {
  return kSiteObservables[(int)Observable()].values == 2 || (IsCoefficient() && !GetCoefficientAbsolute());
}

void Unit::RatioPartsAt(const WorkingCopy &data, int slot, const azure2::Results &results, double &weight,
                        double &element) const {
  const Site &site = *GetSubPoint(data, slot);
  if (kSiteObservables[(int)Observable()].values == 2) {
    weight = SlotValue(site, slot, results);
    element = SlotAnalyzingPower(site, slot, results);
  } else {
    // A normalised Legendre coefficient a_k: the absolute A_k = a_k A_0 and A_0 are folded (an absolute A_k is a
    // cross-section-like quantity and takes the plain path).
    weight = site.GetAngDistA0(results);
    element = site.ValueFor(*this, slot, results);
  }
}

void Unit::RatioWeightRoute(int &component, bool &ayCrossSection) const {
  const bool twoValues = kSiteObservables[(int)Observable()].values == 2;
  component = twoValues ? 0 : -1;
  ayCrossSection = twoValues;
}

void Unit::FoldSubPoints(const WorkingCopy &data, const Config &configure, azure2::Results &results) {
  if (!this->FoldsRatio()) {
    FoldPoint(*this, data, configure, results);
    this->IntegrateTargetEffectComponents(data, configure, results);
    return;
  }
  const int n = this->NumSubPoints();
  std::vector<double> weight(n), product(n);
  for (int i = 1; i <= n; i++) {
    double element = 0.0;
    this->RatioPartsAt(data, i, results, weight[i - 1], element);
    product[i - 1] = weight[i - 1] * element;
  }
  const double denominator = FoldValues(*this, data, configure, weight);
  const double numerator = FoldValues(*this, data, configure, product);
  this->SetFitCrossSection(results, std::fabs(denominator) > 0.0 ? numerator / denominator : 0.0);
}
