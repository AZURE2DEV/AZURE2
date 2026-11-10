#include "AMatrixFunc.h"
#include "AdaptiveIntegrationGrid.h"
#include "ALevel.h"
#include "JGroup.h"
#include "PPair.h"
#include <memory>
#include "Profile.h"
#include "AngCoeff.h"
#include "CNuc.h"
#include "Config.h"
#include "ChannelFunc.h"
#include "CoulFunc.h"
#include "Kinematics.h"
#include "ExperimentAdapter.h"
#include "azure2/experiment/Frames.h"
#include "azure2/experiment/Measurement.h"
#include "ECIntegral.h"
#include "ECAmplitudeCache.h"
#include "ECTable.h"
#include "EData.h"
#include "ESegment.h"
#include "RMatrixFunc.h"
#include "ScatteringLength.h"
#include "ShftFunc.h"
#include "ThmFunc.h"
#include "ThmLineshape.h"
#include "THMMatrixFunc.h"
#include "TargetEffect.h"
#include "IntegratedFermiFunc.h"
#include "AChannel.h"
#include "ThmDistortion.h"
#include "ThmDwVertex.h"
#include "ThmAngular.h"
#include <atomic>
#include <cstdint>
#include <iostream>
#include <map>
#include <cmath>
#include <sstream>
#include <assert.h>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <limits>
#include <atomic>
#include <stdexcept>
#include <cstdlib>
#include <algorithm>

/*!
 * This constructor is used if the data point is to be created from a line in a data file.
 * A pointer to the parent segment is passed to the constructor for the intialization of
 * the EPoint object.
 */

EPoint::EPoint(const azure2::experiment::MeasurementRows &rows, int row, ESegment *parent) : rows_(&rows), row_(row) {
  entrance_key_ = parent->GetEntranceKey();
  exit_key_ = parent->GetExitKey();
  cm_angle_ = rows.angle[row];
  lab_angle_ = rows.angle[row];
  nominal_energy_ = 0.0;  // a row reads its energy from the rows
  bin_low_cm_ = rows.binLow[row];
  bin_high_cm_ = rows.binHigh[row];
  cm_energy_ = rows.energy[row];
  lab_energy_ = rows.energy[row];
  excitation_energy_ = rows.energy[row];
  crossSectionComponent_ = parent ? parent->GetCrossSectionComponent() : 0;
  segment_key_ = parent->GetSegmentKey();
  geofactor_ = 0.;
  sfactorconv_ = 0.;
  is_differential_ = parent->IsDifferential();
  is_phase_ = parent->IsPhase();
  is_thm_ = parent->IsTHM();
  thm_weight_ = parent->GetThmWeight();
  thm_lineshape_ = parent->GetThmLineshape();
  thm_window_ = parent->GetThmSpectatorWindow();
  thm_distortion_ = parent->GetThmDistortion();
  thm_dw_ = parent->GetThmDwVertex();
  thm_angle_ = parent->GetThmAngleWindow();
  thm_coherent_ = parent->GetThmCoherent();
  is_ang_dist_ = parent->IsAngularDist() || parent->IsCoefficients();
  max_ang_dist_order_ = parent->GetMaxAngDistOrder();
  j_value_ = parent->GetJ();
  phase_s_value_ = parent->GetPhaseS();
  l_value_ = parent->GetL();
  isUPOS_ = parent->IsUPOS();
  secondaryDecayL_ = parent->GetSecondaryDecayL();
  Ic_ = parent->GetIc();
  delta_ = parent->GetDelta();
  is_mapped_ = false;
  targetEffectNum_ = 0;
  parentData_ = NULL;
  stoppingPower_ = 0.0;
}

/*!
 * This constructor is used if the data point is to be created with no experimental data.
 * A pointer to the parent segment is passed to the constructor for the intialization of
 * the EPoint object.
 */

EPoint::EPoint(double angle, double energy, ESegment *parent) {
  entrance_key_ = parent->GetEntranceKey();
  exit_key_ = parent->GetExitKey();
  lab_angle_ = angle;
  cm_angle_ = angle;
  nominal_energy_ = energy;
  lab_energy_ = energy;
  cm_energy_ = energy;
  excitation_energy_ = energy;
  crossSectionComponent_ = parent ? parent->GetCrossSectionComponent() : 0;
  segment_key_ = parent->GetSegmentKey();
  geofactor_ = 0.;
  sfactorconv_ = 0.;
  is_differential_ = parent->IsDifferential();
  is_phase_ = parent->IsPhase();
  is_thm_ = parent->IsTHM();
  thm_weight_ = parent->GetThmWeight();
  thm_lineshape_ = parent->GetThmLineshape();
  thm_window_ = parent->GetThmSpectatorWindow();
  thm_distortion_ = parent->GetThmDistortion();
  thm_dw_ = parent->GetThmDwVertex();
  thm_angle_ = parent->GetThmAngleWindow();
  thm_coherent_ = parent->GetThmCoherent();
  is_ang_dist_ = parent->IsAngularDist() || parent->IsCoefficients();
  max_ang_dist_order_ = parent->GetMaxAngDistOrder();
  j_value_ = parent->GetJ();
  phase_s_value_ = parent->GetPhaseS();
  l_value_ = parent->GetL();
  isUPOS_ = parent->IsUPOS();
  secondaryDecayL_ = parent->GetSecondaryDecayL();
  Ic_ = parent->GetIc();
  delta_ = parent->GetDelta();
  is_mapped_ = false;
  targetEffectNum_ = 0;
  parentData_ = NULL;
  stoppingPower_ = 0.0;
}

/*!
 * This constructor is used if the data point is to be created with no experimental data
 * and no parent segment.  Such a point is used if AZURE is to be called as an energy-
 * dependent function for reaction rate and target effect calculations with dynamic
 * integration. All options must be set manually.
 */

EPoint::EPoint(double angle, double energy, int entranceKey,
               int exitKey, bool isDifferential, bool isPhase, bool isAngularDist, double jValue, int lValue, int maxAngDistOrder) {
  entrance_key_ = entranceKey;
  exit_key_ = exitKey;
  lab_angle_ = angle;
  cm_angle_ = angle;
  nominal_energy_ = energy;
  lab_energy_ = energy;
  cm_energy_ = energy;
  excitation_energy_ = energy;
  crossSectionComponent_ = 0;
  geofactor_ = 0.;
  sfactorconv_ = 0.;
  is_differential_ = isDifferential;
  is_phase_ = isPhase;
  is_thm_ = false;
  is_ang_dist_ = isAngularDist;
  max_ang_dist_order_ = maxAngDistOrder;
  j_value_ = jValue;
  l_value_ = lValue;
  isUPOS_ = false;
  secondaryDecayL_ = 0;
  Ic_ = 0.0;
  delta_ = 0.0;
  is_mapped_ = false;
  targetEffectNum_ = 0;
  parentData_ = NULL;
  stoppingPower_ = 0.0;
}

/*!
 * Returns true if the point is differential cross section, otherwise returns false.
 */

bool EPoint::IsDifferential() const {
  return is_differential_;
}

/*!
 * Returns true if the point is phase shift, otherwise returns false.
 */

bool EPoint::IsPhase() const {
  return is_phase_;
}

/*!
 * Returns true if the point belongs to a THM (modified R-matrix) segment.
 */

bool EPoint::IsTHM() const {
  return is_thm_;
}

/*!
 * Returns true if the point is angular distribution, otherwise returns false.
 */

bool EPoint::IsAngularDist() const {
  return is_ang_dist_;
}

/*!
 * Returns true if the point is a mapped point, otherwise returns false.  Mapping in
 * AZURE is performed so calculations are not redundantly performed for like energies.
 * Energy dependent quantities are calculated only once for a given energy, and then
 * copied to mapped points.
 */

bool EPoint::IsMapped() const {
  return is_mapped_;
}

/*!
 * This function retruns true if the point has a corresponding TargetEffect object, otherwise it returns false.
 */

bool EPoint::IsTargetEffect() const {
  if (GetTargetEffectNum() != 0)
    return true;
  else
    return false;
}

/*!
 * Returns true if this is an Unobserved Primary, Observed Secondary (UPOS) point.
 */

bool EPoint::IsUPOS() const {
  return isUPOS_;
}

/*!
 * Returns the secondary decay angular momentum for a UPOS point.
 */

int EPoint::GetSecondaryDecayL() const {
  return secondaryDecayL_;
}

/*!
 * Returns the final state spin for a UPOS point.
 */

double EPoint::GetIc() const {
  return Ic_;
}

/*!
 * Returns the multipole mixing ratio for a UPOS point.
 */

double EPoint::GetDelta() const {
  return delta_;
}

/*!
 * Returns the entrance particle pair key of the data point.  The key need not be the same
 * as the position of the pair in the PPair vector.  Pair keys are used in the setup files
 * of AZURE.
 */

int EPoint::GetEntranceKey() const {
  return entrance_key_;
}

/*!
 * Returns the exit particle pair key of the data point. The key need not be the same
 * as the position of the pair in the PPair vector.  Pair keys are used in the setup files
 * of AZURE.
 */

int EPoint::GetExitKey() const {
  return exit_key_;
}

/*!
 * Returns the orbital angular momentum value for the point.  Only applies if the point is
 * phase shift.
 */

int EPoint::GetL() const {
  return l_value_;
}

/*!
 * Returns the total number of sub-points contained within the present objet.  The sub-points are used to calculate the target effect integrals.
 */

int EPoint::NumSubPoints() const {
  return integrationPoints_.size();
}

/*!
 * Returns the position of the corresponding TargetEffect object in the parent EData object.
 */

int EPoint::GetTargetEffectNum() const {
  return targetEffectNum_;
}

/*!
 * Returns the maximum polynomial order of the point is angular distribution.
 */

int EPoint::GetMaxAngDistOrder() const {
  return max_ang_dist_order_;
}

/*!
 * Returns the angle of the data point in the laboratory frame.
 */

double EPoint::GetLabAngle() const {
  return lab_angle_;
}

/*!
 * Returns the angle of the data point in the center of mass frame.
 */

double EPoint::GetCMAngle() const {
  return cm_angle_;
}

/*!
 * Returns the energy of the point in the laboratory frame.
 */

double EPoint::GetLabEnergy() const {
  return lab_energy_;
}

/*!
 * Returns the energy of the point in the center of mass frame.
 */

double EPoint::GetCMEnergy() const {
  return cm_energy_;
}

/*!
 * Returns the energy of the point in compound excitation energy.
 */

double EPoint::GetExcitationEnergy() const {
  return excitation_energy_;
}

double EPoint::GetOriginalEnergy() const {
  return IsRow() ? rows_->energy[row_] : nominal_energy_;
}

/*!
 * The measured value and its uncertainty as read (lab frame for lab data), from the point's row; 0 and 0.1 for a
 * point that is no row.
 */
double EPoint::GetLabCrossSection() const {
  return IsRow() ? rows_->value[row_] : 0.0;
}

double EPoint::GetLabCrossSectionError() const {
  return IsRow() ? rows_->error[row_] : 0.1;
}

int EPoint::GetFileRecord() const {
  return IsRow() ? rows_->record[row_] : 0;
}

/*!
 * An angular-coefficient row or a sub-point of one (the angle column is the Legendre order k): its order and whether
 * the value is A_k.
 */
bool EPoint::IsCoefficient() const {
  return rows_ && rows_->angleIsLegendreOrder;
}

int EPoint::GetCoefficientOrder() const {
  return IsCoefficient() ? (int)std::round(rows_->angle[row_]) : -1;
}

bool EPoint::GetCoefficientAbsolute() const {
  return IsCoefficient() && rows_->absolute;
}

/*!
 * Returns the total spin value for the data point.  Only applies if the
 * point is phase shift.
 */

double EPoint::GetJ() const {
  return j_value_;
}

/*!
 * For target integration to fit yield curves, the stopping power (or, rather, stopping cross section) is calculated and stored for each sub-point.  This function returns the precalculated value.
 */

double EPoint::GetStoppingPower() const {
  return stoppingPower_;
}

/*!
 * Returns the energy loss of the beam in the target for the current EPoint object.
 */

double EPoint::GetTargetThickness() const {
  return targetThickness_;
}

/*!
 * If a point is mapped, returns an EnergyMap structure containing the point to which it is mapped.
 */

EnergyMap EPoint::GetMap() const {
  EnergyMap thisMap;
  if (this->IsMapped()) {
    thisMap.segment = energy_map_.segment;
    thisMap.point = energy_map_.point;
  } else {
    thisMap.segment = 0;
    thisMap.point = 0;
  }
  return thisMap;
}

/*!
 * Initializes a data point to be used in a calculation.  Initilization is done before the fitting process
 * to calculate all energy dependent quantities that do no rely on the R-Matrix fit parameters.
 */

void EPoint::Initialize(CNuc *compound, const Config &configure) {
  this->CalcEDependentValues(compound, configure);
  if (this->IsDifferential()) {
    // Get TargetEffect with Q-coefficients if applicable
    TargetEffect *effect = NULL;
    EData *parentData = this->GetParentData();
    if (parentData && this->IsTargetEffect()) {
      TargetEffect *te = parentData->GetTargetEffect(this->GetTargetEffectNum());
      if (te && te->IsQCoefficients()) {
        effect = te;
      }
    }
    this->CalcLegendreP(configure.maxLOrder, compound, effect);
  }
  this->CalcCoulombAmplitude(compound);
  if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) this->CalculateECAmplitudes(compound, configure);
}

// ---------------------------------------------------------------- lab -> c.m.
// The formulas are azure2::experiment's (experiment/include/azure2/experiment/Frames.h, docs/upgrade/data_model.md
// step S6); these functions read the point's lab values and its pairs, call them, and store the results in the point
// with the state rules they always had (the maximum lab angle and the reason kept from the last map that set them; a
// relativistic Jacobian kept at its load-time value where an energy shift moved the point beyond the kinematic
// range).  ESegment::ConvertPointToCM is the one sequence that calls them.

void EPoint::NoteAngleLimit(const azure2::experiment::AngleLimit &lim) const {
  if (lim.maxSet) labAngleMax_ = lim.max;
  if (!lim.why.empty()) labAngleWhy_ = lim.why;
}

void EPoint::SetFrame(const azure2::experiment::FrameRow &row) {
  cm_energy_ = row.eCM;
  excitation_energy_ = row.eX;
  bin_low_cm_ = row.binLowCM;
  bin_high_cm_ = row.binHighCM;
  cm_angle_ = row.thetaCM;
}

/*!
 * Calculates the excitation energy from the c.m. energy: E_x = E_cm + S + E* of the pair.
 */
void EPoint::ConvertExcitationEnergy(PPair *pPair) {
  excitation_energy_ = azure2::experiment::ExcitationEnergy(cm_energy_, PartitionOf(pPair));
}

/*!
 * Relativistic lab -> c.m. transformation of a detected particle (relativistic=1): particle 1 of \p exitPair (of
 * \p entrancePair for elastic scattering) at the point's lab angle, at its c.m. energy
 * (azure2::experiment::RelativisticLabToCM).  Where the lab angle is reached from two c.m. angles the root nearest
 * \p nrCMAngle is taken.  Returns false if the angle is not reached (GetLabAngleMax, GetLabAngleWhy).
 */
bool EPoint::RelativisticLabToCM(PPair *entrancePair, PPair *exitPair, double nrCMAngle, double &cmAngle,
                                 double &factor) const {
  azure2::experiment::AngleLimit lim;
  const bool ok = azure2::experiment::RelativisticLabToCM(PartitionOf(entrancePair), PartitionOf(exitPair),
                                                          this->GetCMEnergy(), this->GetLabAngle(), nrCMAngle, cmAngle,
                                                          factor, lim);
  this->NoteAngleLimit(lim);
  return ok;
}

/*!
 * dOmega_lab/dOmega_cm of a particle at this point, without storing it (the RATIO frame factors,
 * ESegment::RatioTermFrameFactor).  Relativistic=1: the relativistic map's factor where it reaches the lab angle,
 * otherwise -- unlike the frame of the data, which keeps the load-time factor -- the non-relativistic one.
 */
double EPoint::CalculateCrossSectionConversionFactor(PPair *entrancePair, PPair *exitPair) {
  double conversionFactor;
  double relAngle;
  if (Kinematics::Relativistic() && this->RelativisticLabToCM(entrancePair, exitPair, this->GetCMAngle(), relAngle, conversionFactor)) {
    return conversionFactor;
  }
  return azure2::experiment::ParticleJacobian(PartitionOf(entrancePair), PartitionOf(exitPair), this->GetLabEnergy(),
                                              this->GetLabAngle(), this->GetCMAngle());
}

/*!
 * The lab-to-c.m. conversion factor for gamma rays (the RATIO path of azure2::EvaluateModelPoint).
 */
double EPoint::CalculateCrossSectionGammaConversionFactor(PPair *entrancePair) {
  return azure2::experiment::PhotonJacobian(PartitionOf(entrancePair), this->GetLabEnergy(), this->GetLabAngle());
}

/*!
 * Sets the lab energy to the given value;
 */

void EPoint::SetLabEnergy(double energy) {
  lab_energy_ = energy;
}

/*!
 * Sets the CM energy to the given value;
 */

void EPoint::SetCMEnergy(double energy) {
  cm_energy_ = energy;
}

/*!
 * Sets the lab angle to the given value;
 */

void EPoint::SetLabAngle(double angle) {
  lab_angle_ = angle;
}

/*!
 * Sets the CM angle to the given value;
 */

void EPoint::SetCMAngle(double angle) {
  cm_angle_ = angle;
}

/*!
 * Sets the exit key to the given value;
 */

void EPoint::SetExitKey(int key) {
  exit_key_ = key;
}

/*!
 * Sets the entrance key for the data point
 */
void EPoint::SetEntranceKey(int key) {
  entrance_key_ = key;
}

/*!
 * If a point is mapped, sets the internal attribute indicating which point it is mapped to.
 */

void EPoint::SetMap(int segmentNum, int pointNum) {
  is_mapped_ = true;
  energy_map_.segment = segmentNum;
  energy_map_.point = pointNum;
}

/*!
 * Clears the mapping status of the point, forcing it to be recalculated.
 * Used when energies change (e.g., after energy shifts during fitting) and
 * mapped points are no longer at identical energies.
 */

void EPoint::ClearMapping() {
  is_mapped_ = false;
  own_e_dependent_ = false;
  energy_map_.segment = 0;
  energy_map_.point = 0;
}


/*!
 * Sets the position of the corresponding TargetEffect object in the parent EData object.
 */

void EPoint::SetTargetEffectNum(int targetEffectNum) {
  targetEffectNum_ = targetEffectNum;
}

/// Sets the line shape of the point and of its sub-points.
void EPoint::SetThmLineshape(const ThmLineshape *l) {
  thm_lineshape_ = l;
  for (EPoint &sub : integrationPoints_) sub.thm_lineshape_ = l;
}

void EPoint::SetThmSpectatorWindow(const ThmSpectatorWindow *w) {
  thm_window_ = w;
  for (EPoint &sub : integrationPoints_) sub.thm_window_ = w;
}

void EPoint::SetThmDistortion(const ThmDistortion *d) {
  thm_distortion_ = d;
  for (EPoint &sub : integrationPoints_) sub.thm_distortion_ = d;
}

void EPoint::SetThmDwVertex(const ThmDwVertex *v) {
  thm_dw_ = v;
  for (EPoint &sub : integrationPoints_) sub.thm_dw_ = v;
}

void EPoint::SetThmCoherent(const ThmCoherentBackground *c) {
  thm_coherent_ = c;
  for (EPoint &sub : integrationPoints_) sub.thm_coherent_ = c;
}

void EPoint::SetThmAngleWindow(const ThmAngleWindow *w) {
  thm_angle_ = w;
  for (EPoint &sub : integrationPoints_) sub.thm_angle_ = w;
}

/*!
 * Adds a sub-point to the current point object for target effect integration.
 */

void EPoint::AddSubPoint(EPoint subPoint) {
  // A sub-point is a copy made before the observable was stamped on the parent,
  // so it has to inherit it here. Without this an analyzing-power segment that
  // also carries target effects is integrated as though it were a cross
  // section, and A_y is never computed at all.
  subPoint.is_analyzing_power_ = this->is_analyzing_power_;
  subPoint.is_total_ = this->is_total_;
  subPoint.is_sub_point_ = true;
  subPoint.screen_ue_ = this->screen_ue_;
  subPoint.screen_param_ = this->screen_param_;
  integrationPoints_.push_back(subPoint);
}

/*!
 * This function sets an internal pointer to the parent EData object.
 */

void EPoint::SetParentData(EData *parentData) {
  parentData_ = parentData;
}

/*!
 * This function sets the stopping cross section (effective) for the current
 * EPoint object.  Used only for sub-point involoved in target effect integration.
 */

void EPoint::SetStoppingPower(double stoppingPower) {
  stoppingPower_ = stoppingPower;
}

/*!
 * This functions sets the energy loss of the beam in the target for the current EPoint object.
 */

void EPoint::SetTargetThickness(double targetThickness) {
  targetThickness_ = targetThickness;
}

/*!
 * Sets the angular distribution coefficients.
 */

void EPoint::SetAngularDists(azure2::Results &results, vector_r dists) {
  if (dists.size() > num_ang_slots_)
    throw std::logic_error("EPoint::SetAngularDists: " + std::to_string(dists.size()) + " coefficients, " +
                           std::to_string(num_ang_slots_) + " result slots (segment " + std::to_string(segment_key_) +
                           ")");
  double *r = Rec(results) + azure2::Results::AngularOffset();
  r[0] = (double)dists.size();
  std::copy(dists.begin(), dists.end(), r + 1);
}

/*!
 * The abundance-weighted terms written by the last evaluation: as many as were written (none for a point that is
 * not the point of a sample segment).  A record poisoned and never written in this evaluation gives its capacity
 * of NaN, so the read shows (AZURE_POISON_RESULTS).
 */

std::vector<double> EPoint::GetNaturalTerms(const azure2::Results &results) const {
  if (!num_natural_slots_) return {};
  const double *r = Rec(results) + azure2::Results::NaturalOffset(num_ang_slots_);
  const int n = std::isnan(r[0]) ? (int)num_natural_slots_ : (int)r[0];
  return std::vector<double>(r + 1, r + 1 + n);
}

void EPoint::SetNaturalTerms(azure2::Results &results, const std::vector<double> &terms) {
  if (terms.size() > num_natural_slots_)
    throw std::logic_error("EPoint::SetNaturalTerms: " + std::to_string(terms.size()) + " terms, " +
                           std::to_string(num_natural_slots_) + " result slots (segment " +
                           std::to_string(segment_key_) + ")");
  double *r = Rec(results) + azure2::Results::NaturalOffset(num_ang_slots_);
  r[0] = (double)terms.size();
  std::copy(terms.begin(), terms.end(), r + 1);
}

void EPoint::SetResultRecord(std::size_t slot, int nAng, int nNat) {
  slot_ = slot;
  num_ang_slots_ = (std::uint16_t)nAng;
  num_natural_slots_ = (std::uint16_t)nNat;
}

void EPoint::UsePrivateRecord(azure2::Results &scratch, const azure2::Results &values) {
  const int nAng = this->ResultAngularCapacity();
  const int nNat = num_natural_slots_;
  std::vector<double> next(azure2::Results::RecordSize(nAng, nNat), 0.0);
  // the values of the record the copy was made from: what the copy read before (a copy of a point copied its
  // result members)
  azure2::CopyRecord(Rec(values), num_ang_slots_, num_natural_slots_, next.data(), nAng, nNat);
  scratch.slot.swap(next);
  this->SetResultRecord(0, nAng, nNat);
}

void EPoint::UsePrivateRecord(azure2::Results &scratch) {
  const int nAng = this->ResultAngularCapacity();
  const int nNat = num_natural_slots_;
  scratch.slot.assign(azure2::Results::RecordSize(nAng, nNat), 0.0);
  // the caller's buffer of a new evaluation (a reaction-rate point): NaN under AZURE_POISON_RESULTS, as every other
  if (azure2::Results::PoisonRequested()) scratch.Poison();
  this->SetResultRecord(0, nAng, nNat);
}

void EPoint::NoResultRecord(const azure2::Results &r) const {
  throw std::logic_error(
      "EPoint: a result was read or written on a point without a record in the buffer it was given (segment " +
      std::to_string(segment_key_) + ", E_lab " + std::to_string(lab_energy_) + " MeV, offset " +
      (slot_ == kNoRecord ? std::string("none") : std::to_string(slot_)) + ", buffer of " +
      std::to_string(r.slot.size()) +
      " slots): the data set's layout (EData::LayoutResults), its EData::BeginEvaluation with this buffer, or the "
      "caller's UsePrivateRecord is missing");
}

/*!
 * Returns a pointer to the parent EData object.
 */

azure2::MappedPoints EPoint::PlanMapped() const {
  if (plan_task_ < 0 || !parentData_) return azure2::MappedPoints();
  const azure2::CalculationTask &task = parentData_->Plan().Task(plan_task_);
  // only the host itself: a copy of it (a probe) or a point mapped onto it calculates itself alone
  if (parentData_->PlanPoint(task.host) != this) return azure2::MappedPoints();
  return azure2::MappedPoints(parentData_, &task.mapped);
}

EData *EPoint::GetParentData() const {
  return parentData_;
}


/*!
 * Returns a pointer to the specified sub-point in the current EPoint object.
 */

EPoint *EPoint::GetSubPoint(int subPoint) {
  EPoint *tempPoint;
  if (subPoint <= integrationPoints_.size())
    tempPoint = &integrationPoints_[subPoint - 1];
  else
    tempPoint = NULL;
  return tempPoint;
}

/*!
 * Returns a reference to the vector of EPoint objects containing the subpoints used in target
 * effect integration.
 */

std::vector<EPoint> &EPoint::GetSubPoints() {
  return integrationPoints_;
}


