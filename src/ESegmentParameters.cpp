// The fit parameters of a segment: its constant and sqrt(E) energy shifts with their current, last and nominal
// values, the shift and normalization groups it belongs to, its energy-shift function (expression and
// coefficients), and the prior penalty of its shift.
#include <limits>
#include <string>
#include <vector>

#include "Config.h"
#include "Equation.h"
#include "ESegment.h"

/*!
 * Returns the nominal energy shift value for the data segment.
 */

double ESegment::GetNominalEnergyShift() const {
  return energyShiftNominal_;
}

/*!
 * Returns the energy shift error for the data segment.
 */

double ESegment::GetEnergyShiftError() const {
  return energyShiftError_;
}

/*!
 * Returns true if the energy shift is varied in the fit.
 */

bool ESegment::IsVaryEnergyShift() const {
  return varyEnergyShift_;
}

/*!
 * Returns the uncertainty of the sqrt(E) coefficient.
 */

double ESegment::GetEnergyShiftSqrtError() const {
  return energyShiftSqrtError_;
}

/*!
 * Returns true if the sqrt(E) coefficient is a fit parameter.
 */

std::string ESegment::GetShiftGroup() const { return shiftGroup_; }

int ESegment::GetShiftLeaderPos() const { return shiftLeaderPos_; }

bool ESegment::IsShiftGroupFollower() const { return shiftLeaderPos_ >= 0; }

std::string ESegment::GetShiftFnGroup() const { return shiftFnGroup_; }

void ESegment::SetShiftGroupFollower(int pos, double shift) {
  shiftLeaderPos_ = pos;
  varyEnergyShift_ = false;  // the group's parameter belongs to its first member
  startEnergyShift_ = shift;
}

std::string ESegment::GetNormGroup() const { return normGroup_; }

int ESegment::GetNormLeaderPos() const { return normLeaderPos_; }

bool ESegment::IsNormGroupFollower() const { return normLeaderPos_ >= 0; }

void ESegment::SetNormGroupFollower(int pos, double norm, double nominal) {
  normLeaderPos_ = pos;
  // The group's norm, its nominal value and its penalty belong to the first member.
  // varyNorm_ is kept: a follower with a varied norm keeps its (fixed) parameter, so
  // the parameter layout, and every index after it, does not depend on the groups.
  startNorm_ = norm;
  dataNormNominal_ = nominal;
  dataNormError_ = 0.0;
}

bool ESegment::IsVaryEnergyShiftSqrt() const {
  return varyEnergyShiftSqrt_;
}

/*!
 * The complete lab-energy shift of a point whose original lab energy is
 * labEnergy: the constant term plus b*sqrt(E/MeV).  The sqrt term is
 * evaluated at the point's own (unshifted) energy, so the mapping is the
 * same function E -> E + a + b*sqrt(E) for every point and every call.
 */

double ESegment::TotalEnergyShift(const SegmentState &state, double labEnergy) const {
  const double shift = (state.energyShiftSqrt == 0.0) ? state.energyShift
                                                      : state.energyShift + state.energyShiftSqrt * sqrt(fabs(labEnergy));
  if (!hasShiftFn_) return shift;
  return shift + state.shiftFn.Evaluate(*shiftFnConfig_, labEnergy);
}

double ESegment::ShiftFnTerm(const SegmentState &state, double labEnergy) const {
  return hasShiftFn_ ? state.shiftFn.Evaluate(*shiftFnConfig_, labEnergy) : 0.0;
}

bool ESegment::OwnsShiftFn() const { return !shiftFnDeclared_.empty(); }
std::string ESegment::GetDeclaredShiftFn() const { return shiftFnDeclared_; }
const std::vector<double> &ESegment::GetNominalShiftFnValues() const { return shiftFnNominal_; }
const std::vector<double> &ESegment::GetShiftFnErrors() const { return shiftFnError_; }
bool ESegment::IsVaryShiftFn(int j) const { return j >= 0 && j < (int)shiftFnVary_.size() && shiftFnVary_[j]; }
bool ESegment::HasSqrtShiftBlock() const { return hasSqrtShiftBlock_; }
bool ESegment::HasShiftFn() const { return hasShiftFn_; }
std::string ESegment::GetShiftFnExpression() const { return hasShiftFn_ ? shiftFnExpression_ : std::string(); }
int ESegment::GetShiftFnSourcePos() const { return shiftFnSourcePos_; }
int ESegment::GetShiftFnParamStart() const { return shiftFnParamStart_; }
void ESegment::SetShiftFnParamStart(int start) { shiftFnParamStart_ = start; }
void ESegment::SetShiftFnValue(SegmentState &state, int j, double value) const {
  if (!hasShiftFn_ || j < 0 || j >= (int)state.shiftFnValues.size()) return;
  state.shiftFnValues[j] = value;
  state.shiftFn.SetParameter(j, value, *shiftFnConfig_);
}

void ESegment::InstallShiftFn(SegmentState &state, const ESegment &owner, int ownerPos, const Config &configure) {
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

// Every energy-shift term is penalised from ZERO (no shift): (dE/sigma)^2, (b/sigma)^2, (a_j/sigma)^2.  The value
// on the segment line is only the starting value.  Before 2026-10-25 the penalty was measured from that value, so
// a refit seeded from an earlier fit's shifts carried those shifts as its prior centre (DeBoer: "It should measure
// the penalty from 0").  The normalization penalty is unchanged (measured from the nominal norm).  On thm an
// explicit "segment_N_energy_shift prior_centre c" row of <parameterSettings> moves the constant term's centre
// to c (GetNominalEnergyShift); without one it is 0, as on dev.
void ESegment::AddShiftPenalty(const SegmentState &state, double &chi) const {
  // A prior the prior covariance replaces is not added (objective.md D3).
  if (varyEnergyShift_ && energyShiftError_ != 0.0 && !shiftPriorReplaced_)
    chi += pow((state.energyShift - energyShiftNominal_) / energyShiftError_, 2.0);  // nominal: 0 or the prior_centre row
  if (varyEnergyShiftSqrt_ && energyShiftSqrtError_ != 0.0 && !shiftSqrtPriorReplaced_) chi += pow(state.energyShiftSqrt / energyShiftSqrtError_, 2.0);
  // The shift function's parameters are penalised once, on the segment that owns them.
  if (OwnsShiftFn() && hasShiftFn_)
    for (int j = 0; j < (int)state.shiftFnValues.size(); j++)
      if (shiftFnVary_[j] && shiftFnError_[j] != 0.0 && !IsShiftFnPriorReplaced(j)) chi += pow(state.shiftFnValues[j] / shiftFnError_[j], 2.0);
}

void ESegment::ShiftPenaltyParts(const SegmentState &state, double &constant, double &sqrtTerm, double &function) const {
  constant = sqrtTerm = function = 0.0;
  if (varyEnergyShift_ && energyShiftError_ != 0.0 && !shiftPriorReplaced_)
    constant = pow((state.energyShift - energyShiftNominal_) / energyShiftError_, 2.0);
  if (varyEnergyShiftSqrt_ && energyShiftSqrtError_ != 0.0 && !shiftSqrtPriorReplaced_)
    sqrtTerm = pow(state.energyShiftSqrt / energyShiftSqrtError_, 2.0);
  if (OwnsShiftFn() && hasShiftFn_)
    for (int j = 0; j < (int)state.shiftFnValues.size(); j++)
      if (shiftFnVary_[j] && shiftFnError_[j] != 0.0 && !IsShiftFnPriorReplaced(j))
        function += pow(state.shiftFnValues[j] / shiftFnError_[j], 2.0);
}

