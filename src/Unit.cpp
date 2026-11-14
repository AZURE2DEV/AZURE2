#include "AMatrixFunc.h"
#include "AdaptiveIntegrationGrid.h"
#include "ALevel.h"
#include "JGroup.h"
#include "PPair.h"
#include "PolarizationFunc.h"
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
#include "WorkingCopy.h"
#include "Segment.h"
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
 * the Unit object.
 */

Unit::Unit(const azure2::experiment::MeasurementRows &rows, int row, const Segment *parent, int term)
{
  settings_ = parent->TermSettingsOf(term).get();
  row_ = row;
  cm_angle_ = rows.angle[row];
  lab_angle_ = rows.angle[row];
  nominal_energy_ = rows.energy[row];
  if (rows.angleIsLegendreOrder) legendre_order_ = (int)std::round(rows.angle[row]);
  bin_low_cm_ = rows.binLow[row];
  bin_high_cm_ = rows.binHigh[row];
  cm_energy_ = rows.energy[row];
  lab_energy_ = rows.energy[row];
  excitation_energy_ = rows.energy[row];
  segment_key_ = parent->GetSegmentKey();
  geofactor_ = 0.;
  sfactorconv_ = 0.;
  targetEffectNum_ = 0;
  stoppingPower_ = 0.0;
}

/*!
 * A point that is no row (a sub-point of a target effect), with the settings of term `term` of the segment.
 */

Site::Site(double angle, double energy, const Segment *parent, int term)
    : settings_(parent->TermSettingsOf(term).get()) {
  lab_angle_ = angle;
  cm_angle_ = angle;
  nominal_energy_ = energy;
  lab_energy_ = energy;
  cm_energy_ = energy;
  excitation_energy_ = energy;
  segment_key_ = parent->TermKey(term);
  geofactor_ = 0.;
  sfactorconv_ = 0.;
  stoppingPower_ = 0.0;
}

Unit::Unit(double angle, double energy, const Segment *parent, int term) : Site(angle, energy, parent, term) {
  targetEffectNum_ = 0;
}

/*!
 * This constructor is used if the data point is to be created with no experimental data
 * and no parent segment.  Such a point is used if AZURE is to be called as an energy-
 * dependent function for reaction rate and target effect calculations with dynamic
 * integration. All options must be set manually; the point holds settings of its own.
 */

Unit::Unit(double angle, double energy, int entranceKey,
               int exitKey, bool isDifferential, bool isPhase, bool isAngularDist, double jValue, int lValue, int maxAngDistOrder) {
  auto settings = std::make_shared<TermSettings>();
  settings->entranceKey = entranceKey;
  settings->exitKey = exitKey;
  settings->differential = isDifferential;
  settings->phase = isPhase;
  settings->angularDist = isAngularDist;
  settings->maxAngDistOrder = maxAngDistOrder;
  settings->j = jValue;
  settings->l = lValue;
  own_settings_ = std::move(settings);
  settings_ = own_settings_.get();
  lab_angle_ = angle;
  cm_angle_ = angle;
  nominal_energy_ = energy;
  lab_energy_ = energy;
  cm_energy_ = energy;
  excitation_energy_ = energy;
  geofactor_ = 0.;
  sfactorconv_ = 0.;
  targetEffectNum_ = 0;
  stoppingPower_ = 0.0;
}

/*!
 * Returns true if the point is differential cross section, otherwise returns false.
 */

bool Site::IsDifferential() const {
  return settings_->differential;
}

/*!
 * Returns true if the point is phase shift, otherwise returns false.
 */

bool Site::IsPhase() const {
  return settings_->phase;
}

/*!
 * Returns true if the point is angular distribution, otherwise returns false.
 */

bool Site::IsAngularDist() const {
  return settings_->angularDist;
}

/*!
 * This function retruns true if the point has a corresponding TargetEffect object, otherwise it returns false.
 */

bool Unit::IsTargetEffect() const {
  if (GetTargetEffectNum() != 0)
    return true;
  else
    return false;
}

/*!
 * Returns true if this is an Unobserved Primary, Observed Secondary (UPOS) point.
 */

bool Site::IsUPOS() const {
  return settings_->upos;
}

/*!
 * Returns the secondary decay angular momentum for a UPOS point.
 */

int Site::GetSecondaryDecayL() const {
  return settings_->secondaryDecayL;
}

/*!
 * Returns the final state spin for a UPOS point.
 */

double Site::GetIc() const {
  return settings_->ic;
}

/*!
 * Returns the multipole mixing ratio for a UPOS point.
 */

double Site::GetDelta() const {
  return settings_->delta;
}

/*!
 * Returns the entrance particle pair key of the data point.  The key need not be the same
 * as the position of the pair in the PPair vector.  Pair keys are used in the setup files
 * of AZURE.
 */

int Site::GetEntranceKey() const {
  return settings_->entranceKey;
}

/*!
 * Returns the exit particle pair key of the data point. The key need not be the same
 * as the position of the pair in the PPair vector.  Pair keys are used in the setup files
 * of AZURE.
 */

int Site::GetExitKey() const {
  return settings_->exitKey;
}

/*!
 * Returns the orbital angular momentum value for the point.  Only applies if the point is
 * phase shift.
 */

int Site::GetL() const {
  return settings_->l;
}

/*!
 * Returns the position of the corresponding TargetEffect object in the parent WorkingCopy object.
 */

int Unit::GetTargetEffectNum() const {
  return targetEffectNum_;
}

/*!
 * Returns the maximum polynomial order of the point is angular distribution.
 */

int Site::GetMaxAngDistOrder() const {
  return settings_->maxAngDistOrder;
}

/*!
 * Returns the angle of the data point in the laboratory frame.
 */

double Site::GetLabAngle() const {
  return lab_angle_;
}

/*!
 * Returns the angle of the data point in the center of mass frame.
 */

double Site::GetCMAngle() const {
  return cm_angle_;
}

/*!
 * Returns the energy of the point in the laboratory frame.
 */

double Site::GetLabEnergy() const {
  return lab_energy_;
}

/*!
 * Returns the energy of the point in the center of mass frame.
 */

double Site::GetCMEnergy() const {
  return cm_energy_;
}

/*!
 * Returns the energy of the point in compound excitation energy.
 */

double Site::GetExcitationEnergy() const {
  return excitation_energy_;
}

double Site::GetOriginalEnergy() const {
  return nominal_energy_;
}

/*!
 * An angular-coefficient row or a sub-point of one (the angle column is the Legendre order k): its order and whether
 * the value is A_k.
 */
bool Site::IsCoefficient() const {
  return legendre_order_ >= 0;
}

int Site::GetCoefficientOrder() const {
  return legendre_order_;
}

bool Site::GetCoefficientAbsolute() const {
  return IsCoefficient() && settings_->coefficientAbsolute;
}

/*!
 * Returns the total spin value for the data point.  Only applies if the
 * point is phase shift.
 */

double Site::GetJ() const {
  return settings_->j;
}

/*!
 * For target integration to fit yield curves, the stopping power (or, rather, stopping cross section) is calculated and stored for each sub-point.  This function returns the precalculated value.
 */

double Site::GetStoppingPower() const {
  return stoppingPower_;
}

/*!
 * Returns the energy loss of the beam in the target for the current Unit object.
 */

double Site::GetTargetThickness() const {
  return targetThickness_;
}

/*!
 * Initializes a data point to be used in a calculation.  Initilization is done before the fitting process
 * to calculate all energy dependent quantities that do no rely on the R-Matrix fit parameters.
 */

void Unit::Initialize(WorkingCopy *data, const CNuc *compound, const Config &configure) {
  this->CalcEDependentValues(data, compound, configure);
  if (this->IsDifferential()) {
    // Get TargetEffect with Q-coefficients if applicable
    const TargetEffect *effect = NULL;
    if (data && this->IsTargetEffect()) {
      const TargetEffect *te = data->Data().GetTargetEffect(this->GetTargetEffectNum());
      if (te && te->IsQCoefficients()) {
        effect = te;
      }
    }
    this->CalcLegendreP(data, configure.maxLOrder, compound, effect);
  }
  this->CalcCoulombAmplitude(data, compound);
  if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) this->CalculateECAmplitudes(data, compound, configure);
}

// ---------------------------------------------------------------- lab -> c.m.
// The formulas are azure2::experiment's (experiment/include/azure2/experiment/Frames.h, docs/upgrade/data_model.md
// step S6); these functions read the point's lab values and its pairs, call them, and store the results in the point
// with the state rules they always had (a relativistic Jacobian kept at its load-time value where an energy shift moved
// the point beyond the kinematic range).  Segment::ConvertRowToCM is the one sequence that calls them; it hands the
// maximum lab angle and the reason of a failed map to its caller, for the message.

void Unit::SetFrame(const azure2::experiment::FrameRow &row) {
  cm_energy_ = row.eCM;
  excitation_energy_ = row.eX;
  bin_low_cm_ = row.binLowCM;
  bin_high_cm_ = row.binHighCM;
  cm_angle_ = row.thetaCM;
  frame_jacobian_ = row.jacobian;
}

/*!
 * Calculates the excitation energy from the c.m. energy: E_x = E_cm + S + E* of the pair.
 */
void Site::ConvertExcitationEnergy(const PPair *pPair) {
  excitation_energy_ = azure2::experiment::ExcitationEnergy(cm_energy_, PartitionOf(pPair));
}

/*!
 * Relativistic lab -> c.m. transformation of a detected particle (relativistic=1): particle 1 of \p exitPair (of
 * \p entrancePair for elastic scattering) at the point's lab angle, at its c.m. energy
 * (azure2::experiment::RelativisticLabToCM).  Where the lab angle is reached from two c.m. angles the root nearest
 * \p nrCMAngle is taken.  Returns false if the angle is not reached.
 */
bool Site::RelativisticLabToCM(const PPair *entrancePair, const PPair *exitPair, double nrCMAngle, double &cmAngle,
                                 double &factor) const {
  azure2::experiment::AngleLimit lim;
  const bool ok = azure2::experiment::RelativisticLabToCM(PartitionOf(entrancePair), PartitionOf(exitPair),
                                                          this->GetCMEnergy(), this->GetLabAngle(), nrCMAngle, cmAngle,
                                                          factor, lim);
  return ok;
}

/*!
 * dOmega_lab/dOmega_cm of a particle at this point, without storing it (the RATIO frame factors,
 * Segment::RatioTermFrameFactor).  Relativistic=1: the relativistic map's factor where it reaches the lab angle,
 * otherwise -- unlike the frame of the data, which keeps the load-time factor -- the non-relativistic one.
 */
double Site::CalculateCrossSectionConversionFactor(const PPair *entrancePair, const PPair *exitPair) {
  double conversionFactor;
  double relAngle;
  if (entrancePair->FunctionSettings().relativistic && this->RelativisticLabToCM(entrancePair, exitPair, this->GetCMAngle(), relAngle, conversionFactor)) {
    return conversionFactor;
  }
  return azure2::experiment::ParticleJacobian(PartitionOf(entrancePair), PartitionOf(exitPair), this->GetLabEnergy(),
                                              this->GetLabAngle(), this->GetCMAngle());
}

/*!
 * The lab-to-c.m. conversion factor for gamma rays (the RATIO path of azure2::EvaluateModelPoint).
 */
double Site::CalculateCrossSectionGammaConversionFactor(const PPair *entrancePair) {
  return azure2::experiment::PhotonJacobian(PartitionOf(entrancePair), this->GetLabEnergy(), this->GetLabAngle());
}

/*!
 * Sets the lab energy to the given value;
 */

void Site::SetLabEnergy(double energy) {
  lab_energy_ = energy;
}

/*!
 * Sets the CM energy to the given value;
 */

void Site::SetCMEnergy(double energy) {
  cm_energy_ = energy;
}

/*!
 * Sets the lab angle to the given value;
 */

void Site::SetLabAngle(double angle) {
  lab_angle_ = angle;
}

/*!
 * Sets the CM angle to the given value;
 */

void Site::SetCMAngle(double angle) {
  cm_angle_ = angle;
}

/*!
 * Sets the exit key to the given value;
 */


/*!
 * Sets the entrance key for the data point
 */


/*!
 * Sets the position of the corresponding TargetEffect object in the parent WorkingCopy object.
 */

void Unit::SetTargetEffectNum(int targetEffectNum) {
  targetEffectNum_ = targetEffectNum;
}

/// Sets the line shape of the point and of its sub-points.






void Unit::MakeSubPoint(Site &subPoint) const {
  // A sub-point shares its term's settings (the observable, A_y or total included) with its parent.
  subPoint.TakeSlotSettingsOf(*this);
}

void Site::TakeSlotSettingsOf(const Site &point) {
  is_sub_point_ = true;
  screen_ue_ = point.screen_ue_;
  screen_param_ = point.screen_param_;
}

/*!
 * This function sets the stopping cross section (effective) for the current
 * Unit object.  Used only for sub-point involoved in target effect integration.
 */

void Site::SetStoppingPower(double stoppingPower) {
  stoppingPower_ = stoppingPower;
}

/*!
 * This functions sets the energy loss of the beam in the target for the current Unit object.
 */

void Site::SetTargetThickness(double targetThickness) {
  targetThickness_ = targetThickness;
}

void Unit::CopyAngularDistsToModel(azure2::Results &results) const {
  if (!num_model_ang_) return;
  double *to = results.model.data() + ModelRecord(results) + 1;
  const double held = num_ang_slots_ ? Rec(results)[azure2::Results::AngularOffset()] : 0.0;
  if (std::isnan(held)) {  // poisoned, never written: shows as such (GetNumModelAngularDists)
    to[0] = held;
    return;
  }
  if (!(held > 0.0) || std::isnan(GetAngDistA0(results))) {  // none, or a failed point (MarkFailed)
    to[0] = 0.0;
    return;
  }
  const int count = std::min<int>(num_model_ang_, Config::maxLOrder + 1);
  to[0] = (double)count;
  for (int k = 0; k < count; k++) to[1 + k] = CoefficientOfOrder(*this, results, k) / GetAngDistA0(results);
}

/*!
 * The abundance-weighted terms written by the last evaluation: as many as were written (none for a point that is
 * not the point of a sample segment).  A record poisoned and never written in this evaluation gives its capacity
 * of NaN, so the read shows (AZURE_POISON_RESULTS).
 */

std::vector<double> Site::GetNaturalTerms(const azure2::Results &results) const {
  if (!num_natural_slots_) return {};
  const double *r = Rec(results) + azure2::Results::NaturalOffset(num_ang_slots_);
  const int n = std::isnan(r[0]) ? (int)num_natural_slots_ : (int)r[0];
  return std::vector<double>(r + 1, r + 1 + n);
}

void Site::SetNaturalTerms(azure2::Results &results, const std::vector<double> &terms) {
  if (terms.size() > num_natural_slots_)
    throw std::logic_error("Unit::SetNaturalTerms: " + std::to_string(terms.size()) + " terms, " +
                           std::to_string(num_natural_slots_) + " result slots (segment " +
                           std::to_string(segment_key_) + ")");
  double *r = Rec(results) + azure2::Results::NaturalOffset(num_ang_slots_);
  r[0] = (double)terms.size();
  std::copy(terms.begin(), terms.end(), r + 1);
}

void Site::SetResultRecord(std::size_t slot, int nAng, int nNat) {
  slot_ = slot;
  num_ang_slots_ = (std::uint16_t)nAng;
  num_natural_slots_ = (std::uint16_t)nNat;
}

void Site::UsePrivateRecord(azure2::Results &scratch) {
  const int nAng = this->ResultAngularCapacity(RecordCapacity{Config::maxLOrder + 1, 0});
  const int nNat = num_natural_slots_;
  scratch.slot.assign(azure2::Results::RecordSize(nAng, nNat), 0.0);
  // the caller's buffer of a new evaluation (a reaction-rate point): NaN under AZURE_POISON_RESULTS, as every other
  if (azure2::Results::PoisonRequested()) scratch.Poison();
  this->SetResultRecord(0, nAng, nNat);
}

void Site::NoResultRecord(const azure2::Results &r) const {
  throw std::logic_error(
      "Unit: a result was read or written on a point without a record in the buffer it was given (segment " +
      std::to_string(segment_key_) + ", E_lab " + std::to_string(lab_energy_) + " MeV, offset " +
      (slot_ == kNoRecord ? std::string("none") : std::to_string(slot_)) + ", buffer of " +
      std::to_string(r.slot.size()) +
      " slots): the data set's layout (WorkingCopy::LayoutResults), its WorkingCopy::BeginEvaluation with this buffer, or the "
      "caller's UsePrivateRecord is missing");
}

void Unit::CalcLegendreP(WorkingCopy *data, int maxL, const CNuc *theCNuc, const TargetEffect *targetEffect) {
  Site::CalcLegendreP(data, maxL, theCNuc, targetEffect);
  for (int i = 1; i <= NumSubPoints(); i++)
    if (data->WritesSlot(*this, i)) GetSubPoint(*data, i)->CalcLegendreP(data, maxL, theCNuc, targetEffect);
}

azure2::RowSet Unit::RowSetOver(const WorkingCopy &data, int first, int count, int maxL, const CNuc *theCNuc,
                                const TargetEffect *targetEffect, const azure2::experiment::FrameRule *labFrame) const {
  azure2::RowSet set;
  set.kind = Observable();
  set.first = first;
  set.count = count;
  set.valuesPerSlot = kSiteObservables[(int)set.kind].values;
  set.labFrame = labFrame != nullptr;
  set.perSlot = set.labFrame || (maxL >= 0 && RowDependsOnEnergy(theCNuc, targetEffect));
  const int positions = set.perSlot ? count : 1;
  if (!set.labFrame) {
    set.angles.assign(1, GetCMAngle());
  } else {
    // a lab angle: theta_cm(E_n; theta_lab) and J(E_n; theta_lab) at each slot's energy (energy_grid.md 4.3); NaN at
    // a slot where the lab angle is beyond the kinematic maximum (a margin node of a shift's reach, say): a fold that
    // reads it is a failed point (FoldPoint), one that does not is unaffected
    set.frame = *labFrame;
    set.thetaLab = GetLabAngle();
    set.angles.resize((std::size_t)count);
    set.jacobian.resize((std::size_t)count);
    for (int k = 0; k < count; k++) {
      azure2::experiment::AngleLimit lim;
      if (!azure2::experiment::LabAngleAt(set.frame, data.SlotAt(first + k).GetCMEnergy(), set.thetaLab,
                                          set.angles[k], set.jacobian[k], lim))
        set.angles[k] = set.jacobian[k] = std::numeric_limits<double>::quiet_NaN();
    }
  }
  if (!set.perSlot) {
    set.rowLength = legendreP_->size();
    set.rows = *legendreP_;
  } else {
    set.rowLength = (std::size_t)std::max(maxL + 1, 0);
    set.rows.resize(set.rowLength * (std::size_t)count);
    for (int k = 0; k < count; k++)
      LegendreRowAt(data.SlotAt(first + k).GetCMEnergy(), set.AngleAt(k), maxL, theCNuc, targetEffect,
                    set.rows.data() + (std::size_t)k * set.rowLength);
  }
  if (set.kind == SiteObservable::kAnalyzingPower) {
    // the factors of A_y at the readers' angle: Y_l'^mu (particle exits) and P_k^1 (capture)
    const PPair *in = theCNuc->GetPair(theCNuc->GetPairNumFromKey(GetEntranceKey()));
    const PPair *out = theCNuc->GetPair(theCNuc->GetPairNumFromKey(GetExitKey()));
    set.muMax = Polarization::MaxProjection(in->GetJ(1), in->GetJ(2), out->GetJ(1), out->GetJ(2));
    set.p1.assign(set.rowLength * (std::size_t)positions, 0.0);
    vector_c table;
    for (int k = 0; k < positions; k++) {
      const double theta = set.AngleAt(k) * pi / 180.0;
      Polarization::HarmonicTable(theta, Config::maxLOrder, set.muMax, table);
      set.harmonicsLength = table.size();
      set.harmonics.insert(set.harmonics.end(), table.begin(), table.end());
      double *p1 = set.p1.data() + (std::size_t)k * set.rowLength;
      for (std::size_t l = 1; l < set.rowLength; l++) p1[l] = AngCoeff::LegendreP1((int)l, std::cos(theta));
    }
  }
  if (set.kind == SiteObservable::kElastic ||
      (set.kind == SiteObservable::kAnalyzingPower && GetEntranceKey() == GetExitKey())) {
    // f_C depends on the energy: at the readers' angle and the energies of the slots
    const PPair *pair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(GetEntranceKey()));
    set.coulomb.resize(2 * (std::size_t)count);
    bool spinDependent = false;
    for (int k = 0; k < count; k++)
      CoulombAmplitudesAt(pair, data.SlotAt(first + k).GetCMEnergy(), set.AngleAt(k), set.coulomb[2 * (std::size_t)k],
                          set.coulomb[2 * (std::size_t)k + 1], spinDependent);
  }
  return set;
}

void Unit::CalcEDependentValues(WorkingCopy *data, const CNuc *theCNuc, const Config &configure, const ChannelRowSource *rows) {
  Site::CalcEDependentValues(data, theCNuc, configure, rows);
  for (int i = 1; i <= NumSubPoints(); i++)
    if (data->WritesSlot(*this, i)) GetSubPoint(*data, i)->CalcEDependentValues(data, theCNuc, configure, rows);
}

void Unit::ChannelKeys(WorkingCopy *data, const CNuc *theCNuc, std::vector<ChannelTable::Key> &keys) {
  Site::ChannelKeys(data, theCNuc, keys);
  for (int i = 1; i <= NumSubPoints(); i++)
    if (data->WritesSlot(*this, i)) GetSubPoint(*data, i)->ChannelKeys(data, theCNuc, keys);
}

void Unit::CalcCoulombAmplitude(WorkingCopy *data, const CNuc *theCNuc) {
  Site::CalcCoulombAmplitude(data, theCNuc);
  for (int i = 1; i <= NumSubPoints(); i++)
    if (data->WritesSlot(*this, i)) GetSubPoint(*data, i)->CalcCoulombAmplitude(data, theCNuc);
}

void Unit::CalculateECAmplitudes(WorkingCopy *data, const CNuc *theCNuc, const Config &configure, const ChannelTable *table) {
  Site::CalculateECAmplitudes(data, theCNuc, configure, table);
  for (int i = 1; i <= NumSubPoints(); i++)
    if (data->WritesSlot(*this, i)) GetSubPoint(*data, i)->CalculateECAmplitudes(data, theCNuc, configure, table);
}

// The slots are emptied before any is refilled: a slot keeps the storage of its rows (Site::ClearEDependentRows).
// An energy shift moves each slot and recomputes it on its own (Segment::UpdatePointEnergyWithShift), so it calls
// this on the point's own site only.
void Unit::RecalcEDependentValues(WorkingCopy *data, const CNuc *theCNuc, const Config &configure, const ChannelRowSource *rows) {
  if (EDependentCurrent()) return;
  for (int i = 1; i <= NumSubPoints(); i++)
    if (data->WritesSlot(*this, i)) GetSubPoint(*data, i)->ClearEDependentRows();
  Site::RecalcEDependentValues(data, theCNuc, configure, rows);
  for (int i = 1; i <= NumSubPoints(); i++) {
    if (!data->WritesSlot(*this, i)) continue;
    Site &slot = *GetSubPoint(*data, i);
    slot.CalcEDependentValues(data, theCNuc, configure, rows);
    slot.CalcCoulombAmplitude(data, theCNuc);
  }
}
