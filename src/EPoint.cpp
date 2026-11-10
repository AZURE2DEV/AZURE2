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
#include "DataLine.h"
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
#include "NuclearPotentialManager.h"
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

EPoint::EPoint(DataLine dataLine, ESegment *parent) {
  entrance_key_ = parent->GetEntranceKey();
  exit_key_ = parent->GetExitKey();
  cm_angle_ = dataLine.angle();
  lab_angle_ = dataLine.angle();
  original_energy_ = dataLine.energy();
  if (dataLine.numExtra() >= 2) {
    bin_low_lab_ = dataLine.extra(0);
    bin_high_lab_ = dataLine.extra(1);
    bin_low_cm_ = bin_low_lab_;
    bin_high_cm_ = bin_high_lab_;
  }
  double shiftedEnergy = dataLine.energy();  //+ parent->GetEnergyShift();
  // Don't allow energies below 0.005 MeV
  if (shiftedEnergy < 0.005) {
    shiftedEnergy = dataLine.energy();
  }
  cm_energy_ = dataLine.energy();
  lab_energy_ = shiftedEnergy;
  excitation_energy_ = dataLine.energy();
  crossSectionComponent_ = parent ? parent->GetCrossSectionComponent() : 0;
  segment_key_ = parent->GetSegmentKey();
  cm_crosssection_ = dataLine.crossSection();
  cm_dcrosssection_ = dataLine.error();
  lab_crosssection_ = dataLine.crossSection();
  lab_dcrosssection_ = dataLine.error();
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
  crossSectionKinFactor_ = 1.0;
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
  original_energy_ = energy;
  lab_energy_ = energy;
  cm_energy_ = energy;
  excitation_energy_ = energy;
  crossSectionComponent_ = parent ? parent->GetCrossSectionComponent() : 0;
  segment_key_ = parent->GetSegmentKey();
  cm_crosssection_ = 0.;
  cm_dcrosssection_ = 0.1;
  lab_crosssection_ = 0.;
  lab_dcrosssection_ = 0.1;
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
  crossSectionKinFactor_ = 1.0;
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
  original_energy_ = energy;
  lab_energy_ = energy;
  cm_energy_ = energy;
  excitation_energy_ = energy;
  crossSectionComponent_ = 0;
  cm_crosssection_ = 0.;
  cm_dcrosssection_ = 0.1;
  lab_crosssection_ = 0.;
  lab_dcrosssection_ = 0.1;
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
  return original_energy_;
}

/*!
 * Returns the experimental cross section in the laboratory frame.
 */

double EPoint::GetLabCrossSection() const {
  return lab_crosssection_;
}

/*!
 * Returms the experimental cross section in the center of mass frame.
 */

double EPoint::GetCMCrossSection() const {
  return cm_crosssection_;
}

/*!
 * Returns the experimental uncertainty in the laboratory frame.
 */

double EPoint::GetLabCrossSectionError() const {
  return lab_dcrosssection_;
}

/*!
 * Returns the experimental uncertainty in the center of mass frame.
 */

double EPoint::GetCMCrossSectionError() const {
  return cm_dcrosssection_;
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
 * Returns the kinematic factor to convert cross section from lab to cm angle;
 */

double EPoint::GetCrossSectionKinFactor() const {
  return crossSectionKinFactor_;
}

static complex ECPathwayAmplitude(CNuc *, const Config &, PPair *, int, double, KGroup *, ECMGroup *, double,
                                  const ChannelTable *table = nullptr);

/*!
 * Returns the external capture amplitude for a given external reaction pathway
 * specified by positions in the KGroup and subsequent ECMGroup vectors.
 */

complex EPoint::GetECAmplitude(int kGroupNum, int ecMGroupNum) const {
  return ec_amplitudes_[kGroupNum - 1][ecMGroupNum - 1];
}

/*!
 * Returns the external capture amplitude at the point's current energy, which
 * an energy shift may have moved away from the energy the amplitude was
 * computed at.
 *
 * With the energy tables (Config::useECTable, the default) a moved point is
 * evaluated at its new energy, from the same tables as every other amplitude,
 * and the result is kept until the energy moves again.  The final level's
 * energy is the one the point's own amplitudes were computed with (recorded at
 * the first shifted evaluation), so shifted and unshifted segments stay
 * consistent even where an MCMC moves that level.
 *
 * With --no-ec-table the amplitude is interpolated linearly between the
 * energies of the segment's own points (ECAmplitudeCache) and held flat beyond
 * its ends, as before the tables.
 */
complex EPoint::GetECAmplitudeWithShift(int kGroupNum, int ecMGroupNum, CNuc *theCNuc, const Config &configure) const {
  // Get current energy including shifts
  double currentEnergy = this->GetCMEnergy();

  if (configure.useECTable) {
    const size_t k = kGroupNum - 1, m = ecMGroupNum - 1;
    if (k < ec_energies_.size() && m < ec_energies_[k].size() && ec_energies_[k][m] == currentEnergy)
      return ec_amplitudes_[k][m];
    // Sized from the requested indices: a point without stored amplitudes (a sub-point or a fit's
    // working copy) has an empty ec_amplitudes_, and sizing from it indexed past the end.
    if (ec_shifted_.size() <= k) ec_shifted_.resize(k + 1);
    if (ec_shifted_[k].size() <= m) ec_shifted_[k].resize(m + 1);
    ShiftedECAmplitude &memo = ec_shifted_[k][m];
    if (memo.valid && memo.energy == currentEnergy) return memo.amplitude;

    PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
    int ir = theCNuc->GetPairNumFromKey(this->GetExitKey());
    KGroup *theKGroup = entrancePair->GetDecay(ir)->GetKGroup(kGroupNum);
    ECMGroup *theECMGroup = theKGroup->GetECMGroup(ecMGroupNum);
    int j = theECMGroup->GetJGroupNum();
    if (!memo.hasLevelEnergy) {
      memo.levelEnergy = theCNuc->GetJGroup(j)->GetLevel(theECMGroup->GetLevelNum())->GetECEnergy();
      memo.hasLevelEnergy = true;
    }
    memo.amplitude = ECPathwayAmplitude(theCNuc, configure, entrancePair, j, memo.levelEnergy, theKGroup,
                                        theECMGroup, currentEnergy);
    memo.energy = currentEnergy;
    memo.valid = true;
    return memo.amplitude;
  }

  // Use ECAmplitudeCache for interpolation at any energy
  if (g_ecAmplitudeCache) {
    ECAmplitudeCache::AmplitudeKey key;
    key.kGroupNum = kGroupNum;
    key.ecMGroupNum = ecMGroupNum;
    key.entranceKey = this->GetEntranceKey();
    key.exitKey = this->GetExitKey();
    // Use segment key directly - same as used in EData::CalculateECAmplitudes
    key.segmentKey = segment_key_;

    // Simple interpolation from cache
    complex interpolatedAmplitude = g_ecAmplitudeCache->GetInterpolatedAmplitude(key, currentEnergy);
    return interpolatedAmplitude;
  }

  // Fallback to original cached amplitude if no global cache available
  return GetECAmplitude(kGroupNum, ecMGroupNum);
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

void EPoint::ApplyAngleLimit(const azure2::experiment::AngleLimit &lim) const {
  if (lim.maxSet) labAngleMax_ = lim.max;
  if (!lim.why.empty()) labAngleWhy_ = lim.why;
}

/*!
 * Calculates center of mass energy.  When a data point is initialized, the same energy is copied into
 * the attributes for center of mass and lab energy.  If this function is called, the center of mass energy
 * attribute is overwritten with the value calculated from the lab energy attribute.
 */

void EPoint::ConvertLabEnergy(PPair *pPair) {
  // E_x = E_cm + S + E* of the pair (audit 01 B3, 2026-10-28)
  const azure2::experiment::BeamEnergiesCM e = azure2::experiment::BeamToCM(
      PartitionOf(pPair), this->GetLabEnergy(), bin_low_lab_, bin_high_lab_, Kinematics::Relativistic());
  cm_energy_ = e.eCM;
  excitation_energy_ = e.eX;
  bin_low_cm_ = e.binLowCM;
  bin_high_cm_ = e.binHighCM;
}

/*!
 * Calculates excitation energy. When a data point is initialized, only the excitation energy is calculated
 * if the data is already in the center of mass frame. Else it is calculated by ConvertLabEnergy.
 */

void EPoint::ConvertExcitationEnergy(PPair *pPair) {
  excitation_energy_ = azure2::experiment::ExcitationEnergy(cm_energy_, PartitionOf(pPair));
}

/*!
 * Calculates the total decay energy from the light particle decay energy, assuming the parent nucleus was at
 * rest when it decayed.  When a data point is initialized, the same energy is copied into
 * the attributes for center of mass and lab energy.  If this function is called, the center of mass energy
 * attribute is overwritten with the value calculated from the lab energy attribute.
 */

void EPoint::ConvertDecayEnergy(PPair *pPair) {
  const azure2::experiment::Partition p = PartitionOf(pPair);
  cm_energy_ = azure2::experiment::DecayToCM(p, this->GetLabEnergy(), Kinematics::Relativistic());
  excitation_energy_ = azure2::experiment::ExcitationEnergy(cm_energy_, p);
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
  this->ApplyAngleLimit(lim);
  return ok;
}

/*!
 * Calculates center of mass angles for elastic scattering (azure2::experiment::ElasticLabToCM; relativistic=1:
 * RelativisticLabToCM, with the non-relativistic angle as the branch).  False beyond the kinematic maximum.
 */
bool EPoint::ConvertLabAngle(PPair *pPair) {
  azure2::experiment::AngleLimit lim;
  double nrAngle;
  const bool ok = azure2::experiment::ElasticLabToCM(PartitionOf(pPair), this->GetLabAngle(), nrAngle, lim);
  this->ApplyAngleLimit(lim);
  if (!ok) return false;
  cm_angle_ = nrAngle;
  double cmAngle, factor;
  if (Kinematics::Relativistic()) {
    if (!this->RelativisticLabToCM(pPair, pPair, cm_angle_, cmAngle, factor)) return false;  // why: set by the map
    cm_angle_ = cmAngle;
  }
  return true;
}

/*!
 * Calculates center of mass angles.  When a data point is initialized, the same angle is copied into
 * the attributes for center of mass and lab angles.  If this function is called, the center of mass angle
 * attribute is overwritten with the value calculated from the lab angle attribute.  This version of the
 * overloaded function is for non-elastic particle channels.
 */

bool EPoint::ConvertLabAngle(PPair *entrancePair, PPair *exitPair, const Config &configure) {
  // relativistic=1: the relativistic map (RelativisticLabToCM) decides whether the lab angle is reached and gives
  // the c.m. angle; the Marion-Young angle below only picks the branch.  Until 2026-10-29 its non-relativistic
  // limits decided, so a point between the two maxima was converted with the NON-relativistic angle and Jacobian
  // without a word (7Li(p,n) at 1.90 MeV, 44.9 deg: c.m. angle off by 0.3 deg, Jacobian by 9.5 %).
  if (Kinematics::Relativistic()) {
    const bool nrOk = this->ConvertLabAngleNR(entrancePair, exitPair);
    double relAngle, relFactor;
    // no non-relativistic angle (beyond its maximum or threshold): the forward root
    if (!this->RelativisticLabToCM(entrancePair, exitPair, nrOk ? cm_angle_ : 0.0, relAngle, relFactor)) return false;
    cm_angle_ = relAngle;
    return true;
  }
  if (!this->ConvertLabAngleNR(entrancePair, exitPair)) return false;
  return true;
}

/*!
 * The non-relativistic (Marion-Young) lab -> c.m. angle of a reaction exit (azure2::experiment::ReactionLabToCM);
 * false below the exit threshold or beyond the maximum lab angle (labAngleWhy_, labAngleMax_), cm_angle_ otherwise.
 */
bool EPoint::ConvertLabAngleNR(PPair *entrancePair, PPair *exitPair) {
  azure2::experiment::AngleLimit lim;
  double cmAngle;
  const bool ok = azure2::experiment::ReactionLabToCM(PartitionOf(entrancePair), PartitionOf(exitPair),
                                                      this->GetLabEnergy(), this->GetLabAngle(), cmAngle, lim);
  this->ApplyAngleLimit(lim);
  if (!ok) return false;
  cm_angle_ = cmAngle;
  return true;
}

/*!
 * Calculates center of mass angles for gamma rays: the aberration of a photon emitted with the centre-of-mass
 * velocity beta, cos(theta_cm) = (cos(theta_lab) - beta) / (1 - beta cos(theta_lab)).
 */
void EPoint::ConvertLabAngleGammas(PPair *entrancePair) {
  const double beta = azure2::experiment::CMVelocity(PartitionOf(entrancePair), this->GetLabEnergy());
  cm_angle_ = azure2::experiment::PhotonLabToCM(beta, this->GetLabAngle());
}

/*!
 * Calculates center of mass cross sections.  When a data point is initialized, the same cross section is copied into
 * the attributes for center of mass and lab cross section.  If this function is called, the center of mass cross section
 * attribute is overwritten with the value calculated from the lab cross section attribute.
 */

void EPoint::ConvertCrossSection(PPair *entrancePair, PPair *exitPair) {
  double conversionFactor;
  double relAngle;
  if (Kinematics::Relativistic()) {
    // Unreached only for a point an energy shift moved past the kinematic limit mid-fit (the angle converter
    // warned): the factor stays the one the point was loaded with (KeepLoadFrame), not a non-relativistic one.
    if (this->RelativisticLabToCM(entrancePair, exitPair, this->GetCMAngle(), relAngle, conversionFactor))
      crossSectionKinFactor_ = conversionFactor;
    cm_crosssection_ = this->GetLabCrossSection() * crossSectionKinFactor_;
    cm_dcrosssection_ = this->GetLabCrossSectionError() * crossSectionKinFactor_;
    return;
  }
  conversionFactor = azure2::experiment::ParticleJacobian(PartitionOf(entrancePair), PartitionOf(exitPair),
                                                          this->GetLabEnergy(), this->GetLabAngle(), this->GetCMAngle());
  cm_crosssection_ = this->GetLabCrossSection() * conversionFactor;
  cm_dcrosssection_ = this->GetLabCrossSectionError() * conversionFactor;
  crossSectionKinFactor_ = conversionFactor;
}

/*!
 * dOmega_lab/dOmega_cm of a particle at this point, without storing it (the RATIO frame factors,
 * ESegment::RatioTermFrameFactor).  Relativistic=1: the relativistic map's factor where it reaches the lab angle,
 * otherwise -- unlike ConvertCrossSection, which keeps the load-time factor -- the non-relativistic one.
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
 * Calculates center of mass cross sections for gamma rays.
 */

void EPoint::ConvertCrossSectionGammas(PPair *entrancePair) {
  // dsigma/dOmega_cm = dsigma/dOmega_lab * dOmega_lab/dOmega_cm, dOmega_lab/dOmega_cm = (1 - beta cos theta_lab)^2 /
  // (1 - beta^2) (azure2::experiment::PhotonJacobian).  (Before 2026-10-23 this was (1 - beta)^2/(1 - beta cos
  // theta_lab), wrong at first order in beta: an overall -2 beta and a P1-shaped 3 beta cos(theta) -- -1.9 % at
  // 90 deg and -4.8 % at 167 deg for 12C(a,g) at 3 MeV.)
  double conversionFactor = this->CalculateCrossSectionGammaConversionFactor(entrancePair);
  cm_crosssection_ = this->GetLabCrossSection() * conversionFactor;
  cm_dcrosssection_ = this->GetLabCrossSectionError() * conversionFactor;
  // the output files take the fitted c.m. cross section back to the lab with this factor (EData::WriteOutputFiles)
  crossSectionKinFactor_ = conversionFactor;
}

/*!
 * The lab-to-c.m. conversion factor for gamma rays (the one ConvertCrossSectionGammas stores; also the
 * RATIO path of azure2::EvaluateModelPoint).
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
 * The external-capture amplitude of one pathway (KGroup, ECMGroup of the decay
 * to the point's exit pair) at c.m. energy cmEnergy, for the final level of
 * J-group j at energy levelEnergy: the entrance Coulomb and hard-sphere phases,
 * the photon-energy factor and the radial integral (from ECTable, or directly
 * with --no-ec-table).  CalculateECAmplitudes calls it at the point's energy;
 * GetECAmplitudeWithShift at a shifted one.
 */
static complex ECPathwayAmplitude(CNuc *theCNuc, const Config &configure, PPair *entrancePair, int j,
                                  double levelEnergy, KGroup *theKGroup, ECMGroup *theECMGroup,
                                  double cmEnergy,
                                  const ChannelTable *table) {
  // A closed entrance channel -- at or below threshold, or so close above it
  // that P = 0 in double precision (ChannelFunc) -- captures nothing: the
  // amplitude carries the entrance Coulomb wave, whose routines return NaN or
  // garbage there (e.g. a Gaussian-folding sub-point at E_cm <= 0).
  if (ChannelFunc(entrancePair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC)).IsClosed(cmEnergy))
    return complex(0.0, 0.0);
  double inEnergy = cmEnergy + entrancePair->GetSepE() + entrancePair->GetExE();
  // entrance Phase Calculations;
  CoulFunc theCoulombFunction(entrancePair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
  // Every pathway of this point with the same l asks for the same exp(-i phi) at a_c (CoulFunc::SetRadialScratch).
  if (!(configure.paramMask & Config::CALCULATE_REACTION_RATE)) theCoulombFunction.SetRadialScratch(true);

  double eta = Kinematics::Channel(entrancePair).Eta(cmEnergy);
  complex expCP(1.0, 0.0);
  for (int ll = 1; ll <= theECMGroup->GetL(); ll++)
    expCP *= complex((double)ll / sqrt(pow((double)ll, 2.0) + pow(eta, 2.0)),
                     eta / sqrt(pow((double)ll, 2.0) + pow(eta, 2.0)));
  // exp(-i phi) at a_c from the point's channel table where it holds the open channel, at E_c = (E - S) - e
  // (the table's expression; within one ulp of cmEnergy); otherwise computed.
  int entranceNum = 0;
  for (int p = 1; p <= theCNuc->NumPairs() && !entranceNum; p++)
    if (theCNuc->GetPair(p) == entrancePair) entranceNum = p;
  const ChannelTable::Row *row = table ? table->Find(entranceNum, theECMGroup->GetL(),
                                                     (inEnergy - entrancePair->GetSepE()) - entrancePair->GetExE())
                                       : nullptr;
  complex expHSP = (row && row->P > 0.0) ? row->expHSP
                                         : theCoulombFunction.HardSpherePhase(theECMGroup->GetL(), entrancePair->GetChRad(), cmEnergy);

  double sqrtGammaPene = pow((inEnergy - levelEnergy) / hbarc, theECMGroup->GetMult() + 0.5);

  AChannel *theFinalChannel = theCNuc->GetJGroup(j)->GetChannel(theECMGroup->GetFinalChannel());
  PPair *theFinalPair = theCNuc->GetPair(theFinalChannel->GetPairNum());
  int theInitialLValue;
  double theInitialSValue;
  if (theECMGroup->IsChannelCapture()) {
    MGroup *theChanCapMGroup = entrancePair->GetDecay(theECMGroup->GetChanCapDecay())->GetKGroup(theECMGroup->GetChanCapKGroup())->GetMGroup(theECMGroup->GetChanCapMGroup());
    theInitialLValue = theCNuc->GetJGroup(theChanCapMGroup->GetJNum())->GetChannel(theChanCapMGroup->GetChpNum())->GetL();
    theInitialSValue = theCNuc->GetJGroup(theChanCapMGroup->GetJNum())->GetChannel(theChanCapMGroup->GetChpNum())->GetS();
  } else {
    theInitialLValue = theECMGroup->GetL();
    theInitialSValue = theKGroup->GetS();
  }

  complex integrals;
  if (configure.useECTable) {
    integrals = ECTable::Integral(theFinalPair, configure, theInitialLValue, theFinalChannel->GetL(),
                                  theInitialSValue, theFinalChannel->GetS(),
                                  theECMGroup->GetJ(), theCNuc->GetJGroup(j)->GetJ(),
                                  theECMGroup->GetMult(), theECMGroup->GetRadType(),
                                  inEnergy, levelEnergy, theECMGroup->IsChannelCapture(), table,
                                  theFinalChannel->GetPairNum());
  } else {
    ECIntegral theECIntegral(theFinalPair, configure);
    integrals = theECIntegral(theInitialLValue, theFinalChannel->GetL(),
                              theInitialSValue, theFinalChannel->GetS(),
                              theECMGroup->GetJ(), theCNuc->GetJGroup(j)->GetJ(),
                              theECMGroup->GetMult(), theECMGroup->GetRadType(),
                              inEnergy, levelEnergy,
                              theECMGroup->IsChannelCapture());
  }

  // calculate the total radial integral
  return expCP * expHSP * sqrtGammaPene * integrals;
}

/*!
 * Calculates the external capture amplitudes for the data point.
 * The amplitudes are calculated for all reaction pathways with corresponding
 * entrance and exit pairs.
 */

void EPoint::CalculateECAmplitudes(CNuc *theCNuc, const Config &configure, const ChannelTable *table) {
  int aa = theCNuc->GetPairNumFromKey(this->GetEntranceKey());
  if (theCNuc->GetPair(aa)->GetPType() == 20) return;
  if (theCNuc->GetPair(aa)->IsEntrance()) {
    PPair *entrancePair = theCNuc->GetPair(aa);
    for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
      for (int la = 1; la <= theCNuc->GetJGroup(j)->NumLevels(); la++) {
        if (theCNuc->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
          ALevel *ecLevel = theCNuc->GetJGroup(j)->GetLevel(la);
          int ir = theCNuc->GetPairNumFromKey(this->GetExitKey());
          if (ecLevel->GetECPairNum() == ir) {
            for (int k = 1; k <= entrancePair->GetDecay(ir)->NumKGroups(); k++) {
              KGroup *theKGroup = entrancePair->GetDecay(ir)->GetKGroup(k);
              for (int ecm = 1; ecm <= theKGroup->NumECMGroups(); ecm++) {
                ECMGroup *theECMGroup = theKGroup->GetECMGroup(ecm);
                complex ecAmplitude = ECPathwayAmplitude(theCNuc, configure, entrancePair, j, ecLevel->GetECEnergy(),
                                                         theKGroup, theECMGroup, this->GetCMEnergy(), table);
                this->AddECAmplitude(k, ecm, ecAmplitude);
              }
            }
          }
        }
      }
    }
  }
  for (int i = 1; i <= this->NumSubPoints(); i++) {
    this->GetSubPoint(i)->CalculateECAmplitudes(theCNuc, configure, table);
  }
}

/*!
 * Adds an external capture amplitude with reference to a specified reaction pathway.
 */

void EPoint::AddECAmplitude(int kGroupNum, int ecMGroupNum, complex ecAmplitude, double energy) {
  vector_c d_amp;
  vector_r d_energy;
  while (kGroupNum > ec_amplitudes_.size()) {
    ec_amplitudes_.push_back(d_amp);
    ec_energies_.push_back(d_energy);
  }
  ec_amplitudes_[kGroupNum - 1].push_back(ecAmplitude);
  ec_energies_[kGroupNum - 1].push_back(energy);
  assert(ecMGroupNum == ec_amplitudes_[kGroupNum - 1].size());
  assert(ec_amplitudes_[kGroupNum - 1].size() == ec_energies_[kGroupNum - 1].size());

  // Add to shared cache for interpolation across different EPoints.  The energy
  // tables evaluate shifted points directly (GetECAmplitudeWithShift), so the
  // cache is only filled without them.
  if (g_ecAmplitudeCache && !(g_config && g_config->useECTable)) {
    ECAmplitudeCache::AmplitudeKey key;
    key.kGroupNum = kGroupNum;
    key.ecMGroupNum = ecMGroupNum;
    key.entranceKey = this->GetEntranceKey();
    key.exitKey = this->GetExitKey();
    key.segmentKey = segment_key_;
    g_ecAmplitudeCache->AddAmplitude(key, energy, ecAmplitude);

    // A non-finite amplitude (never cached, see ECAmplitudeCache::AddAmplitude) is replaced by the value
    // interpolated from the neighbouring cached energies.
    if (!std::isfinite(std::real(ecAmplitude)) || !std::isfinite(std::imag(ecAmplitude))) {
      std::cout << "WARNING: Replacing EC Amplitude with interpolated value from cache for kGroup "
                << kGroupNum << ", ecMGroup " << ecMGroupNum
                << " at energy " << energy << " MeV." << std::endl;
      complex interpAmplitude = g_ecAmplitudeCache->GetInterpolatedAmplitude(key, energy);
      std::cout << "         Original Amplitude: " << ecAmplitude
                << ", Interpolated Amplitude: " << interpAmplitude << std::endl;
      ec_amplitudes_[kGroupNum - 1].back() = interpAmplitude;
    }
  }
}

void EPoint::AddECAmplitude(int kGroupNum, int ecMGroupNum, complex ecAmplitude) {
  // Default behavior - use current CM energy
  double currentEnergy = this->GetCMEnergy();
  AddECAmplitude(kGroupNum, ecMGroupNum, ecAmplitude, currentEnergy);
}

/*!
 * Clears all EC amplitudes and energies. Used when cloning to avoid duplicate accumulation.
 */
void EPoint::ClearECAmplitudes() {
  ec_amplitudes_.clear();
  ec_energies_.clear();
  ec_shifted_.clear();
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


