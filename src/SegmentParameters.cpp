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

// Every energy-shift term is penalised from the centre of its prior, ((dE - c)/sigma)^2, ((b - c_b)/sigma)^2,
// ((a_j - c_j)/sigma)^2, with the centre as written and 0 where none is written (an old line; DeBoer: "It should
// measure the penalty from 0").  The value on the line is only the start (docs/upgrade/parameter_kinds.md section 4).
// An old line's "segment_N_energy_shift prior_centre c" row of <parameterSettings> moves the constant term's centre
// to c (GetNominalEnergyShift).
void Segment::AddShiftPenalty(const SegmentState &state, double &chi) const {
  // A prior the prior covariance replaces is not added (objective.md D3).
  if (varyEnergyShift_ && energyShiftError_ != 0.0 && !shiftPriorReplaced_)
    chi += pow((state.energyShift - energyShiftNominal_) / energyShiftError_, 2.0);  // nominal: 0 or the prior_centre row
  if (varyEnergyShiftSqrt_ && energyShiftSqrtError_ != 0.0 && !shiftSqrtPriorReplaced_)
    chi += pow((state.energyShiftSqrt - energyShiftSqrtCentre_) / energyShiftSqrtError_, 2.0);
  // The shift function's parameters are penalised once, on the segment that owns them.
  if (OwnsShiftFn() && hasShiftFn_)
    for (int j = 0; j < (int)state.shiftFnValues.size(); j++)
      if (shiftFnVary_[j] && shiftFnError_[j] != 0.0 && !IsShiftFnPriorReplaced(j))
        chi += pow((state.shiftFnValues[j] - shiftFnCentre_[j]) / shiftFnError_[j], 2.0);
}

void Segment::ShiftPenaltyParts(const SegmentState &state, double &constant, double &sqrtTerm, double &function) const {
  constant = sqrtTerm = function = 0.0;
  if (varyEnergyShift_ && energyShiftError_ != 0.0 && !shiftPriorReplaced_)
    constant = pow((state.energyShift - energyShiftNominal_) / energyShiftError_, 2.0);
  if (varyEnergyShiftSqrt_ && energyShiftSqrtError_ != 0.0 && !shiftSqrtPriorReplaced_)
    sqrtTerm = pow((state.energyShiftSqrt - energyShiftSqrtCentre_) / energyShiftSqrtError_, 2.0);
  if (OwnsShiftFn() && hasShiftFn_)
    for (int j = 0; j < (int)state.shiftFnValues.size(); j++)
      if (shiftFnVary_[j] && shiftFnError_[j] != 0.0 && !IsShiftFnPriorReplaced(j))
        function += pow((state.shiftFnValues[j] - shiftFnCentre_[j]) / shiftFnError_[j], 2.0);
}

