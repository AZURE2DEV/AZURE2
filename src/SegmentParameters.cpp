// The fit parameters of a segment: its constant and sqrt(E) energy shifts with their current, last and nominal
// values, the shift and normalization groups it belongs to, its energy-shift function (expression and
// coefficients), and the prior penalty of its shift.
#include <limits>
#include <string>
#include <vector>

#include "Config.h"
#include "Equation.h"
#include "Segment.h"

/*!
 * Returns the nominal energy shift value for the data segment.
 */

double Segment::GetNominalEnergyShift() const {
  return energyShiftNominal_;
}

/*!
 * Returns the energy shift error for the data segment.
 */

double Segment::GetEnergyShiftError() const {
  return energyShiftError_;
}

/*!
 * Returns true if the energy shift is varied in the fit.
 */

bool Segment::IsVaryEnergyShift() const {
  return varyEnergyShift_;
}

/*!
 * Returns the uncertainty of the sqrt(E) coefficient.
 */

double Segment::GetEnergyShiftSqrtError() const {
  return energyShiftSqrtError_;
}

/*!
 * Returns true if the sqrt(E) coefficient is a fit parameter.
 */

std::string Segment::GetShiftGroup() const { return shiftGroup_; }

int Segment::GetShiftLeaderPos() const { return shiftLeaderPos_; }

bool Segment::IsShiftGroupFollower() const { return shiftLeaderPos_ >= 0; }

std::string Segment::GetShiftFnGroup() const { return shiftFnGroup_; }

void Segment::SetShiftGroupFollower(int pos, double shift) {
  shiftLeaderPos_ = pos;
  varyEnergyShift_ = false;  // the group's parameter belongs to its first member
  startEnergyShift_ = shift;
}

std::string Segment::GetNormGroup() const { return normGroup_; }

int Segment::GetNormLeaderPos() const { return normLeaderPos_; }

bool Segment::IsNormGroupFollower() const { return normLeaderPos_ >= 0; }

void Segment::SetNormGroupFollower(int pos, double norm, double nominal) {
  normLeaderPos_ = pos;
  // The group's norm, its nominal value and its penalty belong to the first member.
  // varyNorm_ is kept: a follower with a varied norm keeps its (fixed) parameter, so
  // the parameter layout, and every index after it, does not depend on the groups.
  startNorm_ = norm;
  dataNormNominal_ = nominal;
  dataNormError_ = 0.0;
}

bool Segment::IsVaryEnergyShiftSqrt() const {
  return varyEnergyShiftSqrt_;
}

/*!
 * The complete lab-energy shift of a point whose original lab energy is
 * labEnergy: the constant term plus b*sqrt(E/MeV).  The sqrt term is
 * evaluated at the point's own (unshifted) energy, so the mapping is the
 * same function E -> E + a + b*sqrt(E) for every point and every call.
 */

double Segment::TotalEnergyShift(const SegmentState &state, double labEnergy) const {
  const double shift = (state.energyShiftSqrt == 0.0) ? state.energyShift
                                                      : state.energyShift + state.energyShiftSqrt * sqrt(fabs(labEnergy));
  if (!hasShiftFn_) return shift;
  return shift + state.shiftFn.Evaluate(*shiftFnConfig_, labEnergy);
}

double Segment::ShiftFnTerm(const SegmentState &state, double labEnergy) const {
  return hasShiftFn_ ? state.shiftFn.Evaluate(*shiftFnConfig_, labEnergy) : 0.0;
}

bool Segment::OwnsShiftFn() const { return !shiftFnDeclared_.empty(); }
std::string Segment::GetDeclaredShiftFn() const { return shiftFnDeclared_; }
const std::vector<double> &Segment::GetNominalShiftFnValues() const { return shiftFnNominal_; }
const std::vector<double> &Segment::GetShiftFnErrors() const { return shiftFnError_; }
const std::vector<double> &Segment::GetShiftFnCentres() const { return shiftFnCentre_; }
bool Segment::IsVaryShiftFn(int j) const { return j >= 0 && j < (int)shiftFnVary_.size() && shiftFnVary_[j]; }
bool Segment::HasSqrtShiftBlock() const { return hasSqrtShiftBlock_; }
bool Segment::HasShiftFn() const { return hasShiftFn_; }
std::string Segment::GetShiftFnExpression() const { return hasShiftFn_ ? shiftFnExpression_ : std::string(); }
int Segment::GetShiftFnSourcePos() const { return shiftFnSourcePos_; }
int Segment::GetShiftFnParamStart() const { return shiftFnParamStart_; }
void Segment::SetShiftFnParamStart(int start) { shiftFnParamStart_ = start; }
void Segment::SetShiftFnValue(SegmentState &state, int j, double value) const {
  if (!hasShiftFn_ || j < 0 || j >= (int)state.shiftFnValues.size()) return;
  state.shiftFnValues[j] = value;
  state.shiftFn.SetParameter(j, value, *shiftFnConfig_);
}

void Segment::InstallShiftFn(SegmentState &state, const Segment &owner, int ownerPos, const Config &configure) {
  const int n = (int)owner.shiftFnNominal_.size();
  state.shiftFn = Equation(owner.shiftFnDeclared_, n, configure);
  shiftFnExpression_ = state.shiftFn.GetEquation();
  shiftFnConfig_ = &configure;
  hasShiftFn_ = true;
  numShiftFnParams_ = n;
  shiftFnSourcePos_ = ownerPos;
  state.shiftFnValues.assign(n, 0.0);
  for (int j = 0; j < n; j++) SetShiftFnValue(state, j, owner.shiftFnNominal_[j]);
  // Never "already applied": the first FillEnergyShiftsFromParams re-energises the points.
  state.lastShiftFnValues.assign(n, std::numeric_limits<double>::quiet_NaN());
}

