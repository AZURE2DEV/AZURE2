#include "AMatrixFunc.h"
#include "Profile.h"
#include "AngCoeff.h"
#include "CNuc.h"
#include "Config.h"
#include "CoulFunc.h"
#include "Kinematics.h"
#include "DataLine.h"
#include "ECIntegral.h"
#include "ECAmplitudeCache.h"
#include "ECTable.h"
#include "EData.h"
#include "ESegment.h"
#include "RMatrixFunc.h"
#include "ScatteringLength.h"
#include "ShftFunc.h"
#include "TargetEffect.h"
#include "IntegratedFermiFunc.h"
#include <iostream>
#include <assert.h>
#include <fstream>
#include <memory>
#include <mutex>

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
  fitcrosssection_ = 0.;
  fitE1crosssection_ = 0.;
  fitE2crosssection_ = 0.;
  sfactorconv_ = 0.;
  is_differential_ = parent->IsDifferential();
  is_phase_ = parent->IsPhase();
  is_ang_dist_ = parent->IsAngularDist() || parent->IsCoefficients();
  max_ang_dist_order_ = parent->GetMaxAngDistOrder();
  j_value_ = parent->GetJ();
  l_value_ = parent->GetL();
  isUPOS_ = parent->IsUPOS();
  secondaryDecayL_ = parent->GetSecondaryDecayL();
  Ic_ = parent->GetIc();
  delta_ = parent->GetDelta();
  is_mapped_ = false;
  targetEffectNum_ = 0;
  parentData_ = NULL;
  stoppingPower_ = 0.0;
  angleKinFactor_ = 1.0;
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
  fitcrosssection_ = 0.;
  fitE1crosssection_ = 0.;
  fitE2crosssection_ = 0.;
  sfactorconv_ = 0.;
  is_differential_ = parent->IsDifferential();
  is_phase_ = parent->IsPhase();
  is_ang_dist_ = parent->IsAngularDist() || parent->IsCoefficients();
  max_ang_dist_order_ = parent->GetMaxAngDistOrder();
  j_value_ = parent->GetJ();
  l_value_ = parent->GetL();
  isUPOS_ = parent->IsUPOS();
  secondaryDecayL_ = parent->GetSecondaryDecayL();
  Ic_ = parent->GetIc();
  delta_ = parent->GetDelta();
  is_mapped_ = false;
  targetEffectNum_ = 0;
  parentData_ = NULL;
  stoppingPower_ = 0.0;
  angleKinFactor_ = 1.0;
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
  fitcrosssection_ = 0.;
  fitE1crosssection_ = 0.;
  fitE2crosssection_ = 0.;
  sfactorconv_ = 0.;
  is_differential_ = isDifferential;
  is_phase_ = isPhase;
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
 * The maximum order of the Legendre polynomials stored in the point object.
 */

int EPoint::GetMaxLOrder() const {
  return legendreP_.size() - 1;
}

/*!
 * Returns the orbital angular momentum value for the point.  Only applies if the point is
 * phase shift.
 */

int EPoint::GetL() const {
  return l_value_;
}

/*!
 * Returns the number of points mapped to the current point.
 */

int EPoint::NumLocalMappedPoints() const {
  return local_mapped_points_.size();
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
 * Return the number of angular distribution coefficients in the vector.
 */

int EPoint::GetNumAngularDists() const {
  return angularDists_.size();
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
 * Returns the Legendre polynomial specified by an order.
 */

double EPoint::GetLegendreP(int lOrder) const {
  return legendreP_[lOrder];
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
 * Returns the geometrical cross section factor \f$ \frac{\pi}{k^2} \f$.
 */

double EPoint::GetGeometricalFactor() const {
  return geofactor_;
}

/*!
 * Returns the calculated AZURE cross section.
 */

double EPoint::GetFitCrossSection() const {
  return fitcrosssection_;
}

double EPoint::GetFitE1CrossSection() const {
  return fitE1crosssection_;
}

double EPoint::GetFitE2CrossSection() const {
  return fitE2crosssection_;
}

/*!
 * Returns the multiplicative conversion factor from cross section to s-factor.
 */

double EPoint::GetSFactorConversion() const {
  return sfactorconv_;
}

/*!
 * Returns the square root of the penetrability for a channel specified by
 * positions in the JGroup and subsequent AChannel vectors.
 */

double EPoint::GetSqrtPenetrability(int jGroupNum, int channelNum) const {
  return penetrabilities_[jGroupNum - 1][channelNum - 1];
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
 * Returns the angular distribution coefficient corresponding to the given order;
 */

double EPoint::GetAngularDist(int order) const {
  return angularDists_[order];
}

/*!
 * Returns the kinematic factor to convert angle from lab to cm angle;
 */

double EPoint::GetAngleKinFactor() const {
  return angleKinFactor_;
}

/*!
 * Returns the kinematic factor to convert cross section from lab to cm angle;
 */

double EPoint::GetCrossSectionKinFactor() const {
  return crossSectionKinFactor_;
}

/*!
 * Returns the \f$ L_o \f$ diagonal matrix element for a channel specified
 * by positions in the JGroup and subsequent AChannel vectors.
 */

complex EPoint::GetLoElement(int jGroupNum, int channelNum) const {
  return lo_elements_[jGroupNum - 1][channelNum - 1];
}

/*!
 * Returns the factor \f$ \exp(\omega_c) \f$ where \f$ \omega_c \f$ is the Coulomb
 * phase shift.  The  channel us specified
 * by positions in the JGroup and subsequent AChannel vectors.
 */

complex EPoint::GetExpCoulombPhase(int jGroupNum, int channelNum) const {
  return coulombphase_[jGroupNum - 1][channelNum - 1];
}

/*!
 * Returns the factor \f$ \exp(\delta_c) \f$ where \f$ \delta_c \f$ is the hard sphere
 * phase shift.  The  channel us specified
 * by positions in the JGroup and subsequent AChannel vectors.
 */

complex EPoint::GetExpHardSpherePhase(int jGroupNum, int channelNum) const {
  return hardspherephase_[jGroupNum - 1][channelNum - 1];
}

/*!
 * Returns the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

complex EPoint::GetCoulombAmplitude() const {
  return coulombamplitude_;
}

complex EPoint::GetCoulombAmplitude(double s) const {
  // Every pair but an identical one with spin returns the stored amplitude
  // untouched, so the spin-0 and distinguishable results cannot move.
  if (!spinDependentCoulomb_) return coulombamplitude_;
  const int twoS = (int)std::lround(2.0 * s);
  // Channel spins of two identical particles are integers.
  const double sign = ((twoS / 2) % 2 == 0) ? 1.0 : -1.0;
  return coulombdirect_ + sign * coulombexchange_;
}

static complex ECPathwayAmplitude(CNuc *, const Config &, PPair *, int, double, KGroup *, ECMGroup *, double);

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
      memo.levelEnergy = theCNuc->GetJGroup(j)->GetLevel(theECMGroup->GetLevelNum())->GetE();
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

/*!
 * Calculates lab frame value from center of mass value.
 */

double EPoint::ConvertCMValue(double value, PPair *pPair) {
  value = value * (pPair->GetM(1) + pPair->GetM(2)) / (pPair->GetM(2));
  return value;
}

/*!
 * Calculates center of mass energy for a given lab energy.
 */

double EPoint::ConvertLabValue(double value, PPair *pPair) {
  value = value * (pPair->GetM(2)) / (pPair->GetM(1) + pPair->GetM(2));
  return value;
}

/*!
 * Calculates center of mass energy.  When a data point is initialized, the same energy is copied into
 * the attributes for center of mass and lab energy.  If this function is called, the center of mass energy
 * attribute is overwritten with the value calculated from the lab energy attribute.
 */

void EPoint::ConvertLabEnergy(PPair *pPair) {
  if (Kinematics::Relativistic()) {
    // T = sqrt(s) - m1 - m2 with s = (m1 + m2)^2 + 2 m2 E_lab (Kinematics.h)
    const double m1 = pPair->GetM(1) * uconv, m2 = pPair->GetM(2) * uconv;
    cm_energy_ = Kinematics::CMFromLab(m1, m2, this->GetLabEnergy());
    excitation_energy_ = cm_energy_ + pPair->GetSepE();
    bin_low_cm_ = Kinematics::CMFromLab(m1, m2, bin_low_lab_);
    bin_high_cm_ = Kinematics::CMFromLab(m1, m2, bin_high_lab_);
    return;
  }
  double factor = (pPair->GetM(2)) / (pPair->GetM(1) + pPair->GetM(2));
  cm_energy_ = this->GetLabEnergy() * factor;
  excitation_energy_ = cm_energy_ + pPair->GetSepE();
  bin_low_cm_ = bin_low_lab_ * factor;
  bin_high_cm_ = bin_high_lab_ * factor;
}

/*!
 * Calculates excitation energy. When a data point is initialized, only the excitation energy is calculated
 * if the data is already in the center of mass frame. Else it is calculated by ConvertLabEnergy.
 */

void EPoint::ConvertExcitationEnergy(PPair *pPair) {
  excitation_energy_ = cm_energy_ + pPair->GetSepE();
}

/*!
 * Calculates the total decay energy from the light particle decay energy, assuming the parent nucleus was at
 * rest when it decayed.  When a data point is initialized, the same energy is copied into
 * the attributes for center of mass and lab energy.  If this function is called, the center of mass energy
 * attribute is overwritten with the value calculated from the lab energy attribute.
 */

void EPoint::ConvertDecayEnergy(PPair *pPair) {
  if (Kinematics::Relativistic()) {
    // two-body decay at rest: equal momenta, T = E1 + sqrt(m2^2 + p^2) - m2
    cm_energy_ = Kinematics::CMFromDecay(pPair->GetM(1) * uconv, pPair->GetM(2) * uconv, this->GetLabEnergy());
    excitation_energy_ = cm_energy_ + pPair->GetSepE();
    return;
  }
  cm_energy_ = this->GetLabEnergy() /
      (pPair->GetM(2)) *
      (pPair->GetM(1) + pPair->GetM(2));
  excitation_energy_ = cm_energy_ + pPair->GetSepE();
}

/*!
 * Calculates center of mass angles.  When a data point is initialized, the same angle is copied into
 * the attributes for center of mass and lab angles.  If this function is called, the center of mass angle
 * attribute is overwritten with the value calculated from the lab angle attribute.  This version of the
 * overloaded function is for scattering channels.
 */

/*!
 * Relativistic lab -> c.m. transformation of a detected particle (relativistic=1, Kinematics.h): particle 1 of
 * \p exitPair (of \p entrancePair for elastic scattering) at lab angle theta_l.  With the c.m. boost
 * beta = P_lab / E_lab, gamma = E_lab / sqrt(s) and g = beta / beta_3* (beta_3* the detected particle's c.m.
 * velocity, from T_f = E_x - S_f), tan theta_l = sin theta* / (gamma (cos theta* + g)), and
 * dOmega_lab / dOmega_cm = gamma (1 + g cos theta*) / (gamma^2 (g + cos theta*)^2 + sin^2 theta*)^(3/2).
 * Where the lab angle is reached from two c.m. angles (g > 1) the root nearest \p nrCMAngle is taken, i.e. the
 * non-relativistic branch.  Returns false if the angle is not reached.
 */
bool EPoint::RelativisticLabToCM(PPair *entrancePair, PPair *exitPair, double nrCMAngle, double &cmAngle,
                                 double &factor) const {
  const double m1 = entrancePair->GetM(1) * uconv, m2 = entrancePair->GetM(2) * uconv;
  const double m3 = exitPair->GetM(1) * uconv, m4 = exitPair->GetM(2) * uconv;
  const double ti = this->GetCMEnergy();
  const double tf = ti + entrancePair->GetSepE() + entrancePair->GetExE() - exitPair->GetSepE() - exitPair->GetExE();
  const double e1 = Kinematics::LabFromCM(m1, m2, ti);
  const double sqrtS = m1 + m2 + ti, eTot = e1 + m1 + m2;
  const double pLab = sqrt(e1 * (e1 + 2.0 * m1));
  const double beta = pLab / eTot, gam = eTot / sqrtS;
  const double pf = Kinematics::K(m3, m4, tf) * hbarc, e3 = sqrt(m3 * m3 + pf * pf);
  const double g = beta / (pf / e3);
  const double thl = this->GetLabAngle() * pi / 180.;
  double x;  // cos theta*
  if (fabs(sin(thl)) < 1e-12) {
    x = cos(thl) > 0 ? 1.0 : -1.0;
    if (g > 1.0) x = (fabs(nrCMAngle) < 90.0) ? 1.0 : -1.0;
  } else {
    const double t2 = pow(tan(thl), 2.0), gt2 = gam * gam * t2;
    const double disc = 1.0 + gt2 * (1.0 - g * g);
    if (disc < 0.0) return false;
    const double r1 = (-gt2 * g + sqrt(disc)) / (1.0 + gt2), r2 = (-gt2 * g - sqrt(disc)) / (1.0 + gt2);
    // a root is valid if cos theta* + g has the sign of cos theta_l
    const bool v1 = (r1 + g) * cos(thl) > 0.0 && fabs(r1) <= 1.0, v2 = (r2 + g) * cos(thl) > 0.0 && fabs(r2) <= 1.0;
    const double nrX = cos(nrCMAngle * pi / 180.);
    if (v1 && v2) x = (fabs(r1 - nrX) <= fabs(r2 - nrX)) ? r1 : r2;
    else if (v1) x = r1;
    else if (v2) x = r2;
    else return false;
  }
  x = std::max(-1.0, std::min(1.0, x));
  cmAngle = acos(x) * 180. / pi;
  factor = gam * (1.0 + g * x) / pow(gam * gam * (g + x) * (g + x) + 1.0 - x * x, 1.5);
  return true;
}

void EPoint::ConvertLabAngle(PPair *pPair) {
  cm_angle_ = this->GetLabAngle() + 180. / pi * asin(pPair->GetM(1) / pPair->GetM(2) * sin(pi / 180. * this->GetLabAngle()));
  double cmAngle, factor;
  if (Kinematics::Relativistic() && this->RelativisticLabToCM(pPair, pPair, cm_angle_, cmAngle, factor)) cm_angle_ = cmAngle;
}

/*!
 * Calculates center of mass angles.  When a data point is initialized, the same angle is copied into
 * the attributes for center of mass and lab angles.  If this function is called, the center of mass angle
 * attribute is overwritten with the value calculated from the lab angle attribute.  This version of the
 * overloaded function is for non-elastic particle channels.
 */

void EPoint::ConvertLabAngle(PPair *entrancePair, PPair *exitPair, const Config &configure) {
  double qValue = entrancePair->GetSepE() + entrancePair->GetExE() - exitPair->GetSepE() - exitPair->GetExE();
  double a13 = (entrancePair->GetM(1) * exitPair->GetM(1)) * this->GetLabEnergy() / (this->GetLabEnergy() + qValue) /
      (entrancePair->GetM(1) + entrancePair->GetM(2)) / (exitPair->GetM(1) + exitPair->GetM(2));
  double a24 = (entrancePair->GetM(2) * exitPair->GetM(2)) * (1 + entrancePair->GetM(1) / entrancePair->GetM(2) * qValue / (this->GetLabEnergy() + qValue)) /
      (entrancePair->GetM(1) + entrancePair->GetM(2)) / (exitPair->GetM(1) + exitPair->GetM(2));

  if (a13 > a24) {
    double thetaMax = asin(sqrt(a24 / a13)) * 180. / pi;
    if (thetaMax < this->GetLabAngle()) configure.outStream << std::endl
                                                            << "Lab Angle (" << this->GetLabAngle()
                                                            << " degrees) is not kinematically possible.  Maximum angle is "
                                                            << thetaMax << " degrees." << std::endl;
    assert(thetaMax >= this->GetLabAngle());
  }

  double E3PerEt = a13 * pow(cos(this->GetLabAngle() * pi / 180.) + sqrt(a24 / a13 - pow(sin(this->GetLabAngle() * pi / 180.), 2.0)), 2.0);
  double tempE3PerEt = a13 * pow(cos((this->GetLabAngle() + 0.001) * pi / 180.) + sqrt(a24 / a13 - pow(sin((this->GetLabAngle() + 0.001) * pi / 180.), 2.0)), 2.0);
  double slope = sqrt(tempE3PerEt / a24) * sin((this->GetLabAngle() + 0.001) * pi / 180.) - sqrt(E3PerEt / a24) * sin(this->GetLabAngle() * pi / 180.);
  bool switchDomain = false;
  if (slope < 0.) switchDomain = true;

  cm_angle_ = 180. / pi * asin(sqrt(E3PerEt / a24) * sin(this->GetLabAngle() * pi / 180.));
  if (switchDomain) cm_angle_ = 180. - cm_angle_;
  double relAngle, relFactor;
  if (Kinematics::Relativistic() && this->RelativisticLabToCM(entrancePair, exitPair, cm_angle_, relAngle, relFactor))
    cm_angle_ = relAngle;
  // define angle kin factor to convert back to lab frame
  angleKinFactor_ = this->GetLabAngle() / cm_angle_;
  this->SetAngleKinFactor(angleKinFactor_);
}

/*!
 * Calculates center of mass angles as in ConvertLabAngle but uses the reletivistic equations of Iliadis C.37 and C.38.
 * When a data point is initialized, the same angle is copied into
 * the attributes for center of mass and lab angles.  If this function is called, the center of mass angle
 * attribute is overwritten with the value calculated from the lab angle attribute.  This version of the
 * overloaded function is for non-elastic particle channels.
 */

void EPoint::ConvertCMAngle(PPair *entrancePair, PPair *exitPair, const Config &configure) {
  double qValue = entrancePair->GetSepE() + entrancePair->GetExE() - exitPair->GetSepE() - exitPair->GetExE();
  double gamma = sqrt(entrancePair->GetM(1) * exitPair->GetM(1) * this->GetLabEnergy() / (exitPair->GetM(1) * (exitPair->GetM(1) + exitPair->GetM(2)) * qValue + exitPair->GetM(1) * (exitPair->GetM(1) + exitPair->GetM(2) - entrancePair->GetM(1)) * this->GetLabEnergy()));
  double temp_lab_angle = 180. / pi * acos((gamma + cos(this->GetCMAngle() * pi / 180.)) / sqrt(1 + gamma * gamma + 2. * gamma * cos(this->GetCMAngle() * pi / 180.)));
  //  std::cout<<temp_lab_angle<<std::endl;
}

/*!
 * Velocity (in units of c) of the centre of mass of a projectile of kinetic energy E_lab on a target at rest:
 * beta = p c / (E_total) = sqrt(E_lab (E_lab + 2 m1 c^2)) / ((m1 + m2) c^2 + E_lab).
 */
static double CMVelocity(double labEnergy, PPair *entrancePair) {
  return sqrt(labEnergy * (labEnergy + 2.0 * entrancePair->GetM(1) * uconv)) /
      ((entrancePair->GetM(1) + entrancePair->GetM(2)) * uconv + labEnergy);
}

/*!
 * Calculates center of mass angles for gamma rays: the aberration of a photon emitted with the centre-of-mass
 * velocity beta, cos(theta_cm) = (cos(theta_lab) - beta) / (1 - beta cos(theta_lab)).
 */
void EPoint::ConvertLabAngleGammas(PPair *entrancePair) {
  double beta = CMVelocity(this->GetLabEnergy(), entrancePair);
  cm_angle_ = 180. / pi * acos((cos(this->GetLabAngle() * pi / 180.) - beta) / (1 - beta * cos(this->GetLabAngle() * pi / 180.)));
  // as for particles: the factor that takes the angle back to the lab for the output files
  angleKinFactor_ = (cm_angle_ > 0.) ? this->GetLabAngle() / cm_angle_ : 1.0;
  this->SetAngleKinFactor(angleKinFactor_);
}

/*!
 * Calculates center of mass cross sections.  When a data point is initialized, the same cross section is copied into
 * the attributes for center of mass and lab cross section.  If this function is called, the center of mass cross section
 * attribute is overwritten with the value calculated from the lab cross section attribute.
 */

void EPoint::ConvertCrossSection(PPair *entrancePair, PPair *exitPair) {
  double conversionFactor;
  double relAngle;
  if (Kinematics::Relativistic() && this->RelativisticLabToCM(entrancePair, exitPair, this->GetCMAngle(), relAngle, conversionFactor)) {
    cm_crosssection_ = this->GetLabCrossSection() * conversionFactor;
    cm_dcrosssection_ = this->GetLabCrossSectionError() * conversionFactor;
    crossSectionKinFactor_ = conversionFactor;
    this->SetCrossSectionKinFactor(crossSectionKinFactor_);
    return;
  }
  if (this->GetLabAngle() == 0.0 || this->GetLabAngle() == 180.0) {
    double m1 = entrancePair->GetM(1);
    double m2 = entrancePair->GetM(2);
    double m3 = exitPair->GetM(1);
    double m4 = exitPair->GetM(2);
    double e1 = this->GetLabEnergy();
    double qValue = entrancePair->GetSepE() + entrancePair->GetExE() - exitPair->GetSepE() - exitPair->GetExE();
    double et = e1 + qValue;
    double a = m1 * m4 * e1 / (m1 + m2) / (m3 + m4) / et;
    double b = m1 * m3 * e1 / (m1 + m2) / (m3 + m4) / et;
    double c = m2 * m3 / (m1 + m2) / (m3 + m4) * (1 + m1 * qValue / m2 / et);
    double d = m2 * m4 / (m1 + m2) / (m3 + m4) * (1 + m1 * qValue / m2 / et);
    double e3et = b + d + 2 * pow(a * c, 0.5) * cos(pi / 180. * this->GetCMAngle());
    conversionFactor = pow(a * c, 0.5) * pow(d / b - pow(sin(this->GetLabAngle() * pi / 180.), 2.0), 0.5) / e3et;
  } else {
    conversionFactor = pow(sin(pi / 180. * this->GetLabAngle()) / sin(pi / 180. * this->GetCMAngle()), 2.0) * cos(pi / 180. * (this->GetCMAngle() - this->GetLabAngle()));
  }
  cm_crosssection_ = this->GetLabCrossSection() * conversionFactor;
  cm_dcrosssection_ = this->GetLabCrossSectionError() * conversionFactor;
  crossSectionKinFactor_ = conversionFactor;
  this->SetCrossSectionKinFactor(crossSectionKinFactor_);
}

/* Calculates the conversion factor */

double EPoint::CalculateCrossSectionConversionFactor(PPair *entrancePair, PPair *exitPair) {
  double conversionFactor;
  double relAngle;
  if (Kinematics::Relativistic() && this->RelativisticLabToCM(entrancePair, exitPair, this->GetCMAngle(), relAngle, conversionFactor)) {
    return conversionFactor;
  }
  if (this->GetLabAngle() == 0.0 || this->GetLabAngle() == 180.0) {
    double m1 = entrancePair->GetM(1);
    double m2 = entrancePair->GetM(2);
    double m3 = exitPair->GetM(1);
    double m4 = exitPair->GetM(2);
    double e1 = this->GetLabEnergy();
    double qValue = entrancePair->GetSepE() + entrancePair->GetExE() - exitPair->GetSepE() - exitPair->GetExE();
    double et = e1 + qValue;
    double a = m1 * m4 * e1 / (m1 + m2) / (m3 + m4) / et;
    double b = m1 * m3 * e1 / (m1 + m2) / (m3 + m4) / et;
    double c = m2 * m3 / (m1 + m2) / (m3 + m4) * (1 + m1 * qValue / m2 / et);
    double d = m2 * m4 / (m1 + m2) / (m3 + m4) * (1 + m1 * qValue / m2 / et);
    double e3et = b + d + 2 * pow(a * c, 0.5) * cos(pi / 180. * this->GetCMAngle());
    conversionFactor = pow(a * c, 0.5) * pow(d / b - pow(sin(this->GetLabAngle() * pi / 180.), 2.0), 0.5) / e3et;
  } else {
    conversionFactor = pow(sin(pi / 180. * this->GetLabAngle()) / sin(pi / 180. * this->GetCMAngle()), 2.0) * cos(pi / 180. * (this->GetCMAngle() - this->GetLabAngle()));
  }
  return conversionFactor;
}

/*!
 * Calculates center of mass cross sections for gamma rays.
 */

void EPoint::ConvertCrossSectionGammas(PPair *entrancePair) {
  // dsigma/dOmega_cm = dsigma/dOmega_lab * dOmega_lab/dOmega_cm.  For a photon, from the aberration above,
  // dOmega_lab/dOmega_cm = d cos(theta_lab)/d cos(theta_cm) = (1 - beta^2)/(1 + beta cos theta_cm)^2
  //                      = (1 - beta cos theta_lab)^2 / (1 - beta^2).
  // (Before 2026-10-23 this was (1 - beta)^2/(1 - beta cos theta_lab), wrong at first order in beta: an overall
  // -2 beta and a P1-shaped 3 beta cos(theta) -- -1.9 % at 90 deg and -4.8 % at 167 deg for 12C(a,g) at 3 MeV.)
  double conversionFactor = this->CalculateCrossSectionGammaConversionFactor(entrancePair);
  cm_crosssection_ = this->GetLabCrossSection() * conversionFactor;
  cm_dcrosssection_ = this->GetLabCrossSectionError() * conversionFactor;
  // the output files take the fitted c.m. cross section back to the lab with this factor (EData::WriteOutputFiles)
  crossSectionKinFactor_ = conversionFactor;
  this->SetCrossSectionKinFactor(crossSectionKinFactor_);
}

/*!
 * The lab-to-c.m. conversion factor for gamma rays (the one ConvertCrossSectionGammas stores; also the
 * RATIO path of ESegment::CalculateTheoreticalCrossSection).  Until 2026-10-26 this copy still had the
 * pre-2026-10-23 formula (1 - beta)^2 / (1 - beta cos theta_lab); a gamma RATIO of one entrance pair was
 * unaffected (the same wrong factor in numerator and denominator).
 */

double EPoint::CalculateCrossSectionGammaConversionFactor(PPair *entrancePair) {
  double beta = CMVelocity(this->GetLabEnergy(), entrancePair);
  return pow(1.0 - beta * cos(this->GetLabAngle() * pi / 180.0), 2.0) / (1.0 - beta * beta);
}

/*!
 * Adds a Legendre polynomial to the vector.  Polynomials are added to the vector in the order as L=0,1,2,...
 */

void EPoint::AddLegendreP(double polynomial) {
  legendreP_.push_back(polynomial);
}

/*!
 * Clears all Legendre polynomials stored in the point.
 * This is needed when recalculating polynomials after energy shifts.
 */

void EPoint::ClearLegendrePolynomials() {
  legendreP_.clear();
}

/*!
 * Sets the geometrical cross section factor \f$ \frac{\pi}{k^2} \f$.
 */

void EPoint::SetGeometricalFactor(double geoFactor) {
  geofactor_ = geoFactor;
}

/*!
 * Sets the calculated AZURE cross section.
 */

void EPoint::SetFitCrossSection(double crossSection) {
  fitcrosssection_ = crossSection;
}

void EPoint::SetFitE1CrossSection(double e1crossSection) {
  fitE1crosssection_ = e1crossSection;
}

void EPoint::SetFitE2CrossSection(double e2crossSection) {
  fitE2crosssection_ = e2crossSection;
}

/*!
 * Sets the multiplicative conversion from cross section to s-factor.
 */

void EPoint::SetSFactorConversion(double conversion) {
  sfactorconv_ = conversion;
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
 * Sets the excitation energy to the given value;
 */

void EPoint::SetExcitationEnergy(double energy) {
  excitation_energy_ = energy;
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
 * Sets the kinematic factor to convert angle from lab to cm angle;
 */

void EPoint::SetAngleKinFactor(double anglekinfactor) {
  angleKinFactor_ = anglekinfactor;
}

/*!
 * Sets the kinematic factor to convert cross section from lab to cm angle;
 */

void EPoint::SetCrossSectionKinFactor(double crosssectionkinfactor) {
  crossSectionKinFactor_ = crosssectionkinfactor;
}

/*!
 * Gauss-Legendre nodes and weights on [-1, 1] (Newton iteration on P_n), built once.
 */
namespace {
struct GaussLegendreRule {
  std::vector<double> x, w;
  explicit GaussLegendreRule(int n) : x(n), w(n) {
    for (int i = 0; i < (n + 1) / 2; i++) {
      double z = cos(pi * (i + 0.75) / (n + 0.5)), dp = 1.0;
      for (int it = 0; it < 100; it++) {
        double p0 = 1.0, p1 = z;
        for (int k = 2; k <= n; k++) {
          const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;
          p0 = p1;
          p1 = p2;
        }
        dp = n * (z * p1 - p0) / (z * z - 1.0);
        const double dz = p1 / dp;
        z -= dz;
        if (fabs(dz) < 1e-15) break;
      }
      x[i] = -z;
      x[n - 1 - i] = z;
      w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
  }
};

const GaussLegendreRule &QTransformRule() {
  static const GaussLegendreRule rule(96);
  return rule;
}

/// P_0(x) .. P_n(x) into p (size n + 1).
void LegendreSeries(int n, double x, double *p) {
  p[0] = 1.0;
  if (n >= 1) p[1] = x;
  for (int l = 2; l <= n; l++) p[l] = (2.0 * l - 1.0) / l * x * p[l - 1] - (l - 1.0) / l * p[l - 2];  // as before 2026-10-25
}
}  // namespace

/*!
 * Lab-frame motion of the detected particle (particle 1 of the exit pair) for the finite-geometry attenuation
 * coefficients: with cos theta* the c.m. angle, cos theta_lab = gb (x + g) / sqrt(gb^2 (x + g)^2 + 1 - x^2) and
 * dOmega_lab/dOmega_cm = gb (1 + g x) / (gb^2 (g + x)^2 + 1 - x^2)^(3/2).  Non-relativistic: gb = 1 and
 * g^2 = m1 m3 E / (m2 m4 (E + Q)) (E the c.m. energy; the same g as ConvertLabAngle's a13/a24); relativistic=1:
 * the c.m. boost gamma and beta / beta_3* of RelativisticLabToCM.  Returns false where the transform is not used.
 */
bool EPoint::QLabBoost(PPair *entrancePair, PPair *exitPair, double &g, double &gb) const {
  const double ti = this->GetCMEnergy();
  const double tf = ti + entrancePair->GetSepE() + entrancePair->GetExE() - exitPair->GetSepE() - exitPair->GetExE();
  if (ti <= 0.0 || tf <= 0.0) return false;
  if (Kinematics::Relativistic()) {
    const double m1 = entrancePair->GetM(1) * uconv, m2 = entrancePair->GetM(2) * uconv;
    const double m3 = exitPair->GetM(1) * uconv, m4 = exitPair->GetM(2) * uconv;
    const double e1 = Kinematics::LabFromCM(m1, m2, ti);
    const double sqrtS = m1 + m2 + ti, eTot = e1 + m1 + m2;
    const double beta = sqrt(e1 * (e1 + 2.0 * m1)) / eTot;
    const double pf = Kinematics::K(m3, m4, tf) * hbarc, e3 = sqrt(m3 * m3 + pf * pf);
    g = beta / (pf / e3);
    gb = eTot / sqrtS;
  } else {
    g = sqrt(entrancePair->GetM(1) * exitPair->GetM(1) * ti / (entrancePair->GetM(2) * exitPair->GetM(2) * tf));
    gb = 1.0;
  }
  return std::isfinite(g) && g < 1.0;
}

/*!
 * Calculates Legendre polynomials up to a maximum order.  The polynomials are added, in order, to a vector.
 *
 * With attenuation coefficients Q_K (finite detector geometry, M. E. Rose, Phys. Rev. 91, 610 (1953)) the stored
 * factors become the coefficients of the R-matrix's c.m. Legendre coefficients c_L in the measured quantity:
 * - photon exits, UPOS (the detected radiation is the secondary gamma ray), analyzing powers, non-particle entrances
 *   and charged-particle elastic scattering: P_L(x) Q_L, Q as entered (for elastic scattering the Coulomb amplitude is
 *   not a finite Legendre series, so the lab transform below does not apply; Q on such data stays an approximation);
 * - a particle exit otherwise: the detector averages the LAB angular distribution, so Q_K multiplies its LAB Legendre
 *   coefficients b_K = sum_L M_KL c_L, M_KL = (2K+1)/2 int dx P_K(x_lab(x)) P_L(x) (x the c.m. cosine; the lab and c.m.
 *   solid-angle elements cancel).  Back at the point's c.m. angle, with J = dOmega_lab/dOmega_cm,
 *     factor_L = P_L(x) + J sum_{K < nQ} (Q_K - 1) M_KL P_K(x_lab),
 *   exact for Q_K = 1 beyond the coefficients given; Q_K = 1 for every K gives P_L(x).
 * Before 2026-10-25 the particle-exit case multiplied P_L by B = U Q, U an O(gamma^2) lab -> c.m. coefficient matrix
 * evaluated at gamma^2 instead of gamma: B is not the attenuation of the lab distribution (it changed the result even
 * with every Q_K = 1), a non-particle entrance gave every factor 0, and more than maxL coefficients read past the end.
 */

void EPoint::CalcLegendreP(int maxL, CNuc *theCNuc, TargetEffect *targetEffect) {
  if (maxL >= 0) {
    std::vector<double> p(maxL + 1);
    const double x = cos(this->GetCMAngle() * pi / 180.0);
    LegendreSeries(maxL, x, p.data());
    const int numQ = targetEffect ? targetEffect->NumQCoefficients() : 0;
    if (numQ > 0) {
      PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
      PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetExitKey()));
      const bool elastic = this->GetEntranceKey() == this->GetExitKey();
      const bool charged = entrancePair->GetZ(1) != 0 && entrancePair->GetZ(2) != 0;
      double g = 0.0, gb = 1.0;
      const bool labTransform = entrancePair->GetPType() == 0 && exitPair->GetPType() == 0 && !this->IsUPOS() &&
          !this->IsAnalyzingPower() && !(elastic && charged) && this->QLabBoost(entrancePair, exitPair, g, gb);
      if (!labTransform) {
        for (int l = 0; l <= maxL && l < numQ; l++) p[l] *= targetEffect->GetQCoefficient(l);
      } else {
        const GaussLegendreRule &rule = QTransformRule();
        const int n = static_cast<int>(rule.x.size());
        auto toLab = [g, gb](double c) { return gb * (c + g) / sqrt(gb * gb * (c + g) * (c + g) + 1.0 - c * c); };
        const double jac = gb * (1.0 + g * x) / pow(gb * gb * (g + x) * (g + x) + 1.0 - x * x, 1.5);
        std::vector<double> pc((maxL + 1) * n), pl(numQ * n), pK(numQ);
        for (int i = 0; i < n; i++) {
          LegendreSeries(maxL, rule.x[i], &pc[i * (maxL + 1)]);
          LegendreSeries(numQ - 1, toLab(rule.x[i]), &pl[i * numQ]);
        }
        LegendreSeries(numQ - 1, toLab(x), pK.data());
        std::vector<double> add(maxL + 1, 0.0);
        for (int K = 0; K < numQ; K++) {
          const double dq = targetEffect->GetQCoefficient(K) - 1.0;
          if (dq == 0.0) continue;
          const double f = jac * dq * pK[K] * (2.0 * K + 1.0) / 2.0;
          for (int i = 0; i < n; i++) {
            const double wk = f * rule.w[i] * pl[i * numQ + K];
            for (int l = 0; l <= maxL; l++) add[l] += wk * pc[i * (maxL + 1) + l];
          }
        }
        for (int l = 0; l <= maxL; l++) p[l] += add[l];
      }
    }
    for (int l = 0; l <= maxL; l++) this->AddLegendreP(p[l]);
  }
  for (int i = 1; i <= this->NumSubPoints(); i++) {
    this->GetSubPoint(i)->CalcLegendreP(maxL, theCNuc, targetEffect);
  }
}


/*!
 * This function calculates several energy dependent quantities simultaniously.
 * This includes the geometrical cross section, the s-factor conversion, the \f$ L_o \f$ matrix
 * elements, the square root of the penetrability, and the exponentials of the Coulomb phase shifts
 * and hard sphere phase shfits.
 */

void EPoint::CalcEDependentValues(CNuc *theCNuc, const Config &configure, bool includeSubPoints) {
  AZUREProfile::Timer profileTimer(AZUREProfile::kEDependent);
  PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
  PPair *exitPair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetExitKey()));

  double inEnergy;
  double geofactor;
  double sfactorconv;
  if (theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()))->GetPType() == 20) {
    inEnergy = this->GetCMEnergy() + exitPair->GetSepE() + exitPair->GetExE();
    geofactor = 1.;
    sfactorconv = 1.;
  } else {
    inEnergy = this->GetCMEnergy() + entrancePair->GetSepE() + entrancePair->GetExE();
    geofactor = Kinematics::Channel(entrancePair).PiOverK2(this->GetCMEnergy());
    sfactorconv = this->GetCMEnergy() * exp(2 * pi * sqrt(uconv / 2.) * fstruc * entrancePair->GetZ(1) * entrancePair->GetZ(2) * sqrt(entrancePair->GetRedMass() / this->GetCMEnergy()));
  }
  this->SetGeometricalFactor(geofactor);
  this->SetSFactorConversion(sfactorconv);
  // Electron screening: f = exp(pi eta U_e / E), so d ln f / d U_e = pi eta / E (charged particle pairs only).
  screen_k_ = 0.0;
  if (entrancePair->GetPType() == 0 && entrancePair->GetZ(1) * entrancePair->GetZ(2) != 0 && this->GetCMEnergy() > 0.0) {
    const double eta = Kinematics::Channel(entrancePair).Eta(this->GetCMEnergy());
    screen_k_ = pi * eta / this->GetCMEnergy();
  }

  // Only the groups of the entrance pair's compound nucleus take part in this
  // point: a group of another nucleus in the file has thresholds on a
  // different energy axis, and its rows stay empty (every consumer skips it).
  const int nucleus = entrancePair->GetNucleus();
  for (int j = 1; j <= theCNuc->NumJGroups(); j++) {
    if (theCNuc->GetJGroup(j)->IsInRMatrix() && theCNuc->GetJGroup(j)->GetNucleus() == nucleus) {
      JGroup *theJGroup = theCNuc->GetJGroup(j);
      for (int ch = 1; ch <= theJGroup->NumChannels(); ch++) {
        AChannel *theChannel = theJGroup->GetChannel(ch);
        PPair *thePair = theCNuc->GetPair(theChannel->GetPairNum());
        int lValue = theChannel->GetL();
        double localEnergy = inEnergy - thePair->GetSepE() - thePair->GetExE();
        if (thePair->GetPType() == 0) {
          if (fabs(localEnergy) <= ScatteringLength::kThresholdTolerance) {
            // Exactly at the channel's threshold (the entrance channel of a scattering length):
            // the E -> 0 limit, P = 0 and S_l(0) from the zero-energy Coulomb solutions, in
            // place of the Whittaker log-derivative at zero binding, which is ill-conditioned.
            complex loElement(ScatteringLength::ZeroEnergyShift(thePair, lValue) - theChannel->GetBoundaryCondition(), 0.0);
            this->AddLoElement(j, ch, loElement);
            this->AddSqrtPenetrability(j, ch, 0.0);
            this->AddExpCoulombPhase(j, ch, 1.0);
            this->AddExpHardSpherePhase(j, ch, 1.0);
          } else if (localEnergy < 0.0) {
            ShftFunc theShiftFunction(thePair);
            double localShift = theShiftFunction(lValue, inEnergy);
            double boundary = theChannel->GetBoundaryCondition();
            complex loElement(localShift - boundary, 0.0);
            this->AddLoElement(j, ch, loElement);
            this->AddSqrtPenetrability(j, ch, 0.0);
            this->AddExpCoulombPhase(j, ch, 1.0);
            this->AddExpHardSpherePhase(j, ch, 1.0);
            //            if(ch==1){
            //            std::ofstream coul_check;
            //            coul_check.open("coul_check_neg.chk",std::ios::app);
            //            coul_check << this->GetCMEnergy() << "  " << lValue << "  " << eta << "  " << coul.F << "  " << coul.dF << "  " << coul.G << "  " << coul.dG << std::endl;
            //            coul_check << this->GetEntranceKey() << "," << this->GetCMEnergy() << "," << lValue << "," << 0.0 << "," << localShift << std::endl;
            //            coul_check.close();
            //	    }
          } else {
            CoulFunc theCoulombFunction(thePair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));
            double radius = thePair->GetChRad();
            double localPene = theCoulombFunction.Penetrability(lValue, radius, localEnergy);
            double localShift = theCoulombFunction.PEShift(lValue, radius, localEnergy);
            double boundary = theChannel->GetBoundaryCondition();
            complex loElement(localShift - boundary, localPene);
            const Kinematics::Channel kin(thePair);
            double eta = kin.Eta(localEnergy);
            double rho = kin.Rho(radius, localEnergy);
            complex expCP(1.0, 0.0);
            for (int ll = 1; ll <= theChannel->GetL(); ll++)
              expCP *= complex((double)ll / sqrt(pow(eta, 2.0) + pow((double)ll, 2.0)),
                               eta / sqrt(pow(eta, 2.0) + pow((double)ll, 2.0)));
            complex expHSP = theCoulombFunction.HardSpherePhase(lValue, radius, localEnergy);
            this->AddLoElement(j, ch, loElement);
            this->AddSqrtPenetrability(j, ch, sqrt(localPene));
            this->AddExpCoulombPhase(j, ch, expCP);
            this->AddExpHardSpherePhase(j, ch, expHSP);
            //            if(ch==1){
            //            std::ofstream coul_check;
            //            coul_check.open("coul_check_pos.chk",std::ios::app);
            //            coul_check << this->GetLabEnergy() << "  " << lValue << "  " << rho << "  " << eta << "  "
            //                       << coul.F << "  " << coul.dF << "  " << coul.G << "  " << coul.dG << std::endl;
            //            coul_check << this->GetEntranceKey() << "," << this->GetCMEnergy() << "," << lValue << "," << expHSP << "," << localShift << std::endl;
            //            coul_check.close();
            //            }
          }
        } else if (thePair->GetPType() == 10) {
          complex loElement = complex(0.0, 0.0);
          this->AddLoElement(j, ch, loElement);
          double sqrtPene = (configure.paramMask & Config::USE_RMC_FORMALISM) ? 1. : pow(localEnergy / hbarc, (double)lValue + 0.5);
          this->AddSqrtPenetrability(j, ch, sqrtPene);
          this->AddExpCoulombPhase(j, ch, 1.0);
          this->AddExpHardSpherePhase(j, ch, 1.0);
        } else if (thePair->GetPType() == 20) {
          complex loElement = complex(0.0, 0.0);
          this->AddLoElement(j, ch, loElement);
          IntegratedFermiFunc fermiFunc(thePair->GetZ(1));
          double endPointE = thePair->GetSepE() - inEnergy;
          double sqrtPene = (1. + endPointE / 0.510998903 <= 1.) ? 0. : sqrt(fermiFunc(1. + endPointE / 0.510998903, exitPair->GetZ(1) + exitPair->GetZ(2), thePair->GetChRad()));
          this->AddSqrtPenetrability(j, ch, sqrtPene);
          this->AddExpCoulombPhase(j, ch, 1.0);
          this->AddExpHardSpherePhase(j, ch, 1.0);
        }
      }
    }
  }
  if (includeSubPoints) {
    for (int i = 1; i <= this->NumSubPoints(); i++) {
      this->GetSubPoint(i)->CalcEDependentValues(theCNuc, configure);
    }
  }
  for (int i = 1; i <= this->NumLocalMappedPoints(); i++) {
    EPoint *mappedPoint = this->GetLocalMappedPoint(i);
    if (mappedPoint->OwnsEDependentValues()) continue;
    mappedPoint->geofactor_ = geofactor_;
    mappedPoint->sfactorconv_ = sfactorconv_;
    mappedPoint->screen_k_ = screen_k_;
    mappedPoint->lo_elements_ = lo_elements_;
    mappedPoint->penetrabilities_ = penetrabilities_;
    mappedPoint->coulombphase_ = coulombphase_;
    mappedPoint->hardspherephase_ = hardspherephase_;
    if (!includeSubPoints) continue;
    for (int ii = 1; ii <= this->NumSubPoints(); ii++) {
      EPoint *subMappedPoint = mappedPoint->GetSubPoint(ii);
      subMappedPoint->geofactor_ = this->GetSubPoint(ii)->geofactor_;
      subMappedPoint->sfactorconv_ = this->GetSubPoint(ii)->sfactorconv_;
      subMappedPoint->screen_k_ = this->GetSubPoint(ii)->screen_k_;
      subMappedPoint->lo_elements_ = this->GetSubPoint(ii)->lo_elements_;
      subMappedPoint->penetrabilities_ = this->GetSubPoint(ii)->penetrabilities_;
      subMappedPoint->coulombphase_ = this->GetSubPoint(ii)->coulombphase_;
      subMappedPoint->hardspherephase_ = this->GetSubPoint(ii)->hardspherephase_;
    }
  }
  eDependentEnergy_ = this->GetCMEnergy();
  eDependentValid_ = true;
}

/*!
 * Empties the per-JGroup rows of the energy-dependent arrays, keeping their storage.
 */
void EPoint::ClearEDependentRows() {
  for (auto &row : lo_elements_) row.clear();
  for (auto &row : penetrabilities_) row.clear();
  for (auto &row : coulombphase_) row.clear();
  for (auto &row : hardspherephase_) row.clear();
}

/*!
 * Recalculates energy-dependent values using the current (possibly shifted) energy.
 * This is needed when energy shifts are applied after initialization.
 */
void EPoint::RecalcEDependentValues(CNuc *theCNuc, const Config &configure, bool includeSubPoints) {
  // Nothing to do while the energy is the one the values were computed at.
  // Every sub-point of a target-effect point used to pay for this on every
  // evaluation whenever external capture was in use (Coulomb and Whittaker
  // functions for each channel, serial): 40 s per evaluation for two
  // convolved 16N beta-delayed alpha spectra on the 12C+alpha model.
  if (eDependentValid_ && this->GetCMEnergy() == eDependentEnergy_) return;
  // Clear existing energy-dependent values first.  The per-JGroup rows are
  // emptied in place rather than destroyed, so CalcEDependentValues refills the
  // same storage instead of freeing and re-allocating every row of every point
  // on every energy-shift step.  That churn, done on the pooled data sets by
  // whichever OpenMP thread happens to evaluate them, fragmented the malloc
  // arenas and made resident memory creep up for the whole fit.  The Add*
  // methods index the rows by J-group, so the values stored are unchanged.
  this->ClearEDependentRows();
  // CalcEDependentValues also refills the sub-points (appending to their rows),
  // so empty theirs too: otherwise every call doubled their rows, and the
  // doubled storage would now be kept.  The callers re-energise and recalculate
  // each sub-point afterwards in any case.
  //
  // ESegment::UpdatePointEnergiesWithShift passes includeSubPoints = false: it
  // moves each sub-point to its shifted energy and recalculates it right
  // after, so recomputing the sub-points here, still at their OLD energies,
  // only doubled the cost of every shift step (all Coulomb functions of every
  // sub-point, twice, serially).
  if (includeSubPoints)
    for (int i = 1; i <= this->NumSubPoints(); i++) this->GetSubPoint(i)->ClearEDependentRows();

  // Recalculate with current energy
  this->CalcEDependentValues(theCNuc, configure, includeSubPoints);

  // Also recalculate Coulomb amplitude with current energy
  this->CalcCoulombAmplitude(theCNuc, includeSubPoints);
}

/*!
 * Copies the energy-dependent values (geometrical factor, S-factor conversion,
 * L_o elements, penetrabilities, Coulomb and hard-sphere phases) from another
 * point at exactly the same c.m. energy and with the same entrance pair --
 * a host of point mapping (EData::MapData) or its sub-point.  Those values do
 * not depend on the angle; the Coulomb amplitude does, and is not copied.
 * Returns false, and copies nothing, if the energies differ.
 */
bool EPoint::CopyEDependentValuesFrom(const EPoint &other) {
  if (other.GetCMEnergy() != this->GetCMEnergy() || other.GetEntranceKey() != this->GetEntranceKey() ||
      other.GetExitKey() != this->GetExitKey() || !other.eDependentValid_ ||
      other.eDependentEnergy_ != other.GetCMEnergy())
    return false;
  geofactor_ = other.geofactor_;
  sfactorconv_ = other.sfactorconv_;
  screen_k_ = other.screen_k_;
  lo_elements_ = other.lo_elements_;
  penetrabilities_ = other.penetrabilities_;
  coulombphase_ = other.coulombphase_;
  hardspherephase_ = other.hardspherephase_;
  eDependentEnergy_ = other.eDependentEnergy_;
  eDependentValid_ = true;
  return true;
}

/*!
 * Adds an \f$ L_o \f$ matrix element with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void EPoint::AddLoElement(int jGroupNum, int channelNum, complex loElement) {
  vector_c d;
  while (jGroupNum > lo_elements_.size()) lo_elements_.push_back(d);
  lo_elements_[jGroupNum - 1].push_back(loElement);
  assert(channelNum = lo_elements_[jGroupNum - 1].size());
}

/*!
 * Adds a square root of penetrability  with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void EPoint::AddSqrtPenetrability(int jGroupNum, int channelNum, double sqrtPene) {
  vector_r d;
  while (jGroupNum > penetrabilities_.size()) penetrabilities_.push_back(d);
  penetrabilities_[jGroupNum - 1].push_back(sqrtPene);
  assert(channelNum = penetrabilities_[jGroupNum - 1].size());
}

/*!
 * Adds an exponential of the Coulomb phase shift with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void EPoint::AddExpCoulombPhase(int jGroupNum, int channelNum, complex expShift) {
  vector_c d;
  while (jGroupNum > coulombphase_.size()) coulombphase_.push_back(d);
  coulombphase_[jGroupNum - 1].push_back(expShift);
  assert(channelNum = coulombphase_[jGroupNum - 1].size());
}

/*!
 * Adds an exponential of the hard sphere phase shift with reference to positions in the JGroup and subsequent
 * AChannel vectors.
 */

void EPoint::AddExpHardSpherePhase(int jGroupNum, int channelNum, complex expShift) {
  vector_c d;
  while (jGroupNum > hardspherephase_.size()) hardspherephase_.push_back(d);
  hardspherephase_[jGroupNum - 1].push_back(expShift);
  assert(channelNum = hardspherephase_[jGroupNum - 1].size());
}

/*!
 * Calculates the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

void EPoint::CalcCoulombAmplitude(CNuc *theCNuc, bool includeSubPoints) {
  if (this->GetEntranceKey() == this->GetExitKey()) {
    PPair *entrancePair = theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
    int z1 = entrancePair->GetZ(1);
    int z2 = entrancePair->GetZ(2);
    double redmass = entrancePair->GetRedMass();
    double energy = this->GetCMEnergy();
    double angle = this->GetCMAngle();
    double eta = Kinematics::Channel(entrancePair).rel ? Kinematics::Channel(entrancePair).Eta(energy)
                                                       : sqrt(uconv / 2.) * fstruc * z1 * z2 * sqrt(redmass / energy);
    double cal = (1.0 / (2.0 * sqrt(pi))) * eta * (1.0 / pow(sin(angle * pi / 360.0), 2.));
    double cex = 2.0 * eta * log(sin(angle * pi / 360.0));
    complex calpha(cal * cos(cex), -cal * sin(cex));
    // Identical-particle Mott amplitude: the exchange term f_C(pi - theta).
    // Exchanging two identical particles of spin j multiplies a channel-spin
    // state |s v> by (-1)^(2j-s) and the wave function must pick up
    // (-1)^(2j), so the exchange term enters channel spin s with (-1)^s --
    // not with one boson/fermion sign for all s. Only for j = 0 (s = 0 alone)
    // do the two coincide, and there the symmetrized amplitude is stored as
    // before and reproduces the spin-0 Mott cross section. For j != 0 the
    // direct and exchange terms are kept apart and combined per channel spin
    // by GetCoulombAmplitude(s).
    spinDependentCoulomb_ = false;
    if (entrancePair->IsIdentical()) {
      double cal_p = (1.0 / (2.0 * sqrt(pi))) * eta * (1.0 / pow(cos(angle * pi / 360.0), 2.));
      double cex_p = 2.0 * eta * log(cos(angle * pi / 360.0));
      complex calpha_p(cal_p * cos(cex_p), -cal_p * sin(cex_p));
      if (entrancePair->HasSpinDependentExchange()) {
        spinDependentCoulomb_ = true;
        coulombdirect_ = calpha;
        coulombexchange_ = calpha_p;
      } else {
        calpha = calpha + double(entrancePair->GetIdenticalSign()) * calpha_p;
      }
    }
    this->SetCoulombAmplitude(calpha);
  } else {
    spinDependentCoulomb_ = false;
    this->SetCoulombAmplitude(complex(0., 0.));
  }
  if (!includeSubPoints) return;
  for (int i = 1; i <= this->NumSubPoints(); i++) {
    this->GetSubPoint(i)->CalcCoulombAmplitude(theCNuc);
  }
}

/*!
 * Sets the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

void EPoint::SetCoulombAmplitude(complex amplitude) {
  coulombamplitude_ = amplitude;
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
                                  double cmEnergy) {
  double inEnergy = cmEnergy + entrancePair->GetSepE() + entrancePair->GetExE();
  // entrance Phase Calculations;
  CoulFunc theCoulombFunction(entrancePair, !!(configure.paramMask & Config::USE_GSL_COULOMB_FUNC));

  double eta = Kinematics::Channel(entrancePair).Eta(cmEnergy);
  complex expCP(1.0, 0.0);
  for (int ll = 1; ll <= theECMGroup->GetL(); ll++)
    expCP *= complex((double)ll / sqrt(pow((double)ll, 2.0) + pow(eta, 2.0)),
                     eta / sqrt(pow((double)ll, 2.0) + pow(eta, 2.0)));
  complex expHSP = theCoulombFunction.HardSpherePhase(theECMGroup->GetL(), entrancePair->GetChRad(), cmEnergy);

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
                                  inEnergy, levelEnergy, theECMGroup->IsChannelCapture());
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

void EPoint::CalculateECAmplitudes(CNuc *theCNuc, const Config &configure) {
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
                complex ecAmplitude = ECPathwayAmplitude(theCNuc, configure, entrancePair, j, ecLevel->GetE(),
                                                         theKGroup, theECMGroup, this->GetCMEnergy());
                this->AddECAmplitude(k, ecm, ecAmplitude);
              }
            }
          }
        }
      }
    }
  }
  for (int i = 1; i <= this->NumSubPoints(); i++) {
    this->GetSubPoint(i)->CalculateECAmplitudes(theCNuc, configure);
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
 * Calculates the cross section for a data point based on the fit parameters in the compound nucleus.
 */

namespace {
thread_local bool tlPointMappingOff = false;
}

ScopedPointMappingOff::ScopedPointMappingOff() : previous_(tlPointMappingOff) { tlPointMappingOff = true; }
ScopedPointMappingOff::~ScopedPointMappingOff() { tlPointMappingOff = previous_; }
bool ScopedPointMappingOff::Active() { return tlPointMappingOff; }

void EPoint::Calculate(CNuc *theCNuc, const Config &configure, EPoint *parent, int subPointNum) {
  // Points mapped onto this one (EData::MapData), unless mapping is switched
  // off on this thread (ScopedPointMappingOff).
  const int numMapped = ScopedPointMappingOff::Active() ? 0 : this->NumLocalMappedPoints();
  if (!this->IsTargetEffect() ||
      !this->GetParentData()->GetTargetEffect(this->GetTargetEffectNum())->IsSubPointEffect()) {
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
    theMatrixFunc->CalculateCrossSection(this);
    // A mapped point may have another exit pair (a component of a SUM, RATIO or total-capture segment, EData::MapData):
    // the solved level matrix is the same, only its T-matrix is built anew.  Mapped points are ordered by exit, so
    // each exit's T-matrix is built once.
    int tExitKey = this->GetExitKey();
    auto mappedCrossSection = [&](EPoint *mapped) {
      if (mapped->GetExitKey() != tExitKey) {
        theMatrixFunc->RecalculateTMatrix(mapped);
        tExitKey = mapped->GetExitKey();
      }
      theMatrixFunc->CalculateCrossSection(mapped);
    };
    if (subPointNum && parent) {
      const int numParentMapped = ScopedPointMappingOff::Active() ? 0 : parent->NumLocalMappedPoints();
      for (int i = 1; i <= numParentMapped; i++)
        mappedCrossSection(parent->GetLocalMappedPoint(i)->GetSubPoint(subPointNum));
    } else {
      for (int i = 1; i <= numMapped; i++) mappedCrossSection(this->GetLocalMappedPoint(i));
    }
  } else {
    TargetEffect *effect = this->GetParentData()->GetTargetEffect(this->GetTargetEffectNum());
    const double blend = this->GetTargetBlendWeight();
    const double autoTol = effect->GetAutoTolerance();
    const bool hasMapped = numMapped > 0;

    // The bare value of the observable at the central energy is needed both
    // for edge blending and for the automatic per-point decision.  It is
    // obtained from a probe copy that pretends not to carry the effect, so the
    // ordinary single-point machinery above stays the single source of truth.
    double bareValue = 0.0;
    bool haveBare = false;
    if (blend < 1.0 || autoTol > 0.) {
      EPoint probe(*this);
      probe.SetTargetEffectNum(0);
      probe.ClearLocalMappedPoints();
      probe.integrationPoints_.clear();
      probe.Calculate(theCNuc, configure);
      bareValue = probe.GetFitCrossSection();
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
        double d = std::fabs(this->GetSubPoint(i)->GetCMEnergy() - this->GetCMEnergy());
        if (d < dBest) { dBest = d; iMid = i; }
      }
      const int probes[3] = {iFirst, iMid, iLast};
      for (int k = 0; k < 3; k++) {
        EPoint *subPoint = this->GetSubPoint(probes[k]);
        if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE)
          subPoint->RecalcEDependentValues(theCNuc, configure);
        if (hasMapped)
          subPoint->Calculate(theCNuc, configure, this, probes[k]);
        else
          subPoint->Calculate(theCNuc, configure);
      }
      bool allPass = true;
      for (int m = 0; m <= numMapped && allPass; m++) {
        EPoint *host = (m == 0) ? this : this->GetLocalMappedPoint(m);
        double sFirst = host->GetSubPoint(iFirst)->GetFitCrossSection();
        double sMid = host->GetSubPoint(iMid)->GetFitCrossSection();
        double sLast = host->GetSubPoint(iLast)->GetFitCrossSection();
        double estimate = std::fabs(0.5 * (sFirst + sLast) - sMid);
        if (effect->IsTargetIntegration()) estimate += 0.5 * std::fabs(sLast - sFirst);
        if (estimate > autoTol * std::max(std::fabs(sMid), 1e-300)) allPass = false;
      }
      if (allPass) {
        this->SetFitCrossSection(bareValue);
        for (int m = 1; m <= numMapped; m++) {
          EPoint *mappedPoint = this->GetLocalMappedPoint(m);
          EPoint mappedProbe(*mappedPoint);
          mappedProbe.SetTargetEffectNum(0);
          mappedProbe.ClearLocalMappedPoints();
          mappedProbe.integrationPoints_.clear();
          mappedProbe.Calculate(theCNuc, configure);
          mappedPoint->SetFitCrossSection(mappedProbe.GetFitCrossSection());
        }
        return;
      }
    }

    for (int i = 1; i <= this->NumSubPoints(); i++) {
      EPoint *subPoint = this->GetSubPoint(i);
      // Recalculate energy-dependent values for sub-points too if using external capture
      if (configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
        subPoint->RecalcEDependentValues(theCNuc, configure);
      }
      if (numMapped > 0)
        subPoint->Calculate(theCNuc, configure, this, i);
      else
        subPoint->Calculate(theCNuc, configure);
    }
    this->IntegrateTargetEffectForObservable(configure);
    if (blend < 1.0 && haveBare)
      this->SetFitCrossSection(blend * this->GetFitCrossSection() + (1.0 - blend) * bareValue);
    for (int i = 1; i <= numMapped; i++) {
      EPoint *mappedPoint = this->GetLocalMappedPoint(i);
      mappedPoint->IntegrateTargetEffectForObservable(configure);
      if (blend < 1.0) {
        // The mapped point is another observable at the same energy: its bare
        // value needs its own probe.
        EPoint mappedProbe(*mappedPoint);
        mappedProbe.SetTargetEffectNum(0);
        mappedProbe.ClearLocalMappedPoints();
        mappedProbe.integrationPoints_.clear();
        mappedProbe.Calculate(theCNuc, configure);
        mappedPoint->SetFitCrossSection(blend * mappedPoint->GetFitCrossSection() +
                                        (1.0 - blend) * mappedProbe.GetFitCrossSection());
      }
    }
  }
}

/*!
 * A segment compared against the E1 or E2 component of an angle-integrated
 * capture cross section (isDiff 5/6) reads that component off the parent
 * point, which a sub-point integration never filled: only the total was
 * integrated.  The kernel is linear in the sub-point values, so each
 * component is integrated with the same combiner by swapping it into the
 * sub-points' cross-section slot, exactly as the analyzing power is handled.
 */
void EPoint::IntegrateTargetEffectComponents(const Config &configure) {
  if (crossSectionComponent_ == 0) return;
  const int n = this->NumSubPoints();
  if (n <= 0) return;
  std::vector<double> sigma(n);
  for (int i = 1; i <= n; i++) sigma[i - 1] = this->GetSubPoint(i)->GetFitCrossSection();
  const double total = this->GetFitCrossSection();

  for (int i = 1; i <= n; i++) this->GetSubPoint(i)->SetFitCrossSection(this->GetSubPoint(i)->GetFitE1CrossSection());
  this->IntegrateTargetEffect(configure);
  this->SetFitE1CrossSection(this->GetFitCrossSection());

  for (int i = 1; i <= n; i++) this->GetSubPoint(i)->SetFitCrossSection(this->GetSubPoint(i)->GetFitE2CrossSection());
  this->IntegrateTargetEffect(configure);
  this->SetFitE2CrossSection(this->GetFitCrossSection());

  for (int i = 1; i <= n; i++) this->GetSubPoint(i)->SetFitCrossSection(sigma[i - 1]);
  this->SetFitCrossSection(total);
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
void EPoint::IntegrateTargetEffectForObservable(const Config &configure) {
  if (this->IsCoefficient() && !this->GetCoefficientAbsolute()) {
    // A normalised Legendre coefficient is a ratio like A_y: the absolute A_k = a_k A_0
    // and A_0 are each folded with the integrator (linear in the cross section), then
    // divided -- averaging the ratio itself would be wrong.  An absolute A_k is a
    // cross-section-like quantity and takes the plain path (its slot holds a_k A_0).
    const int n = this->NumSubPoints();
    std::vector<double> ak(n), a0(n);
    for (int i = 1; i <= n; i++) {
      ak[i - 1] = this->GetSubPoint(i)->GetFitCrossSection();
      a0[i - 1] = this->GetSubPoint(i)->GetAngDistA0();
    }
    for (int i = 1; i <= n; i++) this->GetSubPoint(i)->SetFitCrossSection(a0[i - 1]);
    this->IntegrateTargetEffect(configure);
    const double denominator = this->GetFitCrossSection();
    for (int i = 1; i <= n; i++) this->GetSubPoint(i)->SetFitCrossSection(ak[i - 1] * a0[i - 1]);
    this->IntegrateTargetEffect(configure);
    const double numerator = this->GetFitCrossSection();
    for (int i = 1; i <= n; i++) this->GetSubPoint(i)->SetFitCrossSection(ak[i - 1]);
    this->SetFitCrossSection(std::fabs(denominator) > 0.0 ? numerator / denominator : 0.0);
    return;
  }
  if (!this->IsAnalyzingPower()) {
    this->IntegrateTargetEffect(configure);
    this->IntegrateTargetEffectComponents(configure);
    return;
  }

  const int n = this->NumSubPoints();
  std::vector<double> sigma(n);
  for (int i = 1; i <= n; i++) sigma[i - 1] = this->GetSubPoint(i)->GetFitCrossSection();

  this->IntegrateTargetEffect(configure);
  const double denominator = this->GetFitCrossSection();

  for (int i = 1; i <= n; i++)
    this->GetSubPoint(i)->SetFitCrossSection(sigma[i - 1] *
                                             this->GetSubPoint(i)->GetAnalyzingPower());
  this->IntegrateTargetEffect(configure);
  const double numerator = this->GetFitCrossSection();

  for (int i = 1; i <= n; i++) this->GetSubPoint(i)->SetFitCrossSection(sigma[i - 1]);

  this->SetFitCrossSection(std::fabs(denominator) > 0.0 ? numerator / denominator : 0.0);
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
 * If a point is mapped to the current point, a pointer to the mapped point is added to a vector.
 */

void EPoint::AddLocalMappedPoint(EPoint *point) {
  local_mapped_points_.push_back(point);
}

/*!
 * Clears vector containing pointers to points mapped to the current point.
 */

void EPoint::ClearLocalMappedPoints() {
  local_mapped_points_.clear();
}

/*!
 * Sets the position of the corresponding TargetEffect object in the parent EData object.
 */

void EPoint::SetTargetEffectNum(int targetEffectNum) {
  targetEffectNum_ = targetEffectNum;
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

namespace {
// 2-point Gauss-Legendre over the sub-point grid of `point` (energies descending from sub-point 1),
// restricted to the c.m. window [lo, hi]: the cross section (through `rd`) and, for a yield, the
// stopping power are interpolated linearly between the two sub-points that bracket each node.
// Accumulates num += integral of wNum(E) sigma(E) [/ (stopping(E) 1e24)] and den += integral of
// wDen(E), so a kernel average is num / den and a yield is num alone. The one quadrature every
// experimental effect uses; before 2026-10 it was written out separately in five branches.
template <class V, class Read, class WNum, class WDen>
void GridQuadrature(EPoint *point, Read rd, double lo, double hi, bool divideByStopping, WNum wNum, WDen wDen,
                    V &num, double &den) {
  const double x1 = -1.0 / sqrt(3.0), x2 = 1.0 / sqrt(3.0);
  const int n = point->NumSubPoints();
  for (int j = 1; j < n; j++) {
    const double Ea = point->GetSubPoint(j)->GetCMEnergy();
    const double Eb = point->GetSubPoint(j + 1)->GetCMEnergy();
    if (Eb >= hi || Ea <= lo) continue;
    const double a = std::min(Ea, hi), b = std::max(Eb, lo);
    const double mid = 0.5 * (a + b), half = 0.5 * (a - b);
    if (half < 1.0e-12) continue;
    for (int gp = 0; gp < 2; gp++) {
      const double e = mid + half * (gp == 0 ? x1 : x2);
      const double frac = (Ea - e) / (Ea - Eb);   // 0 at sub-point j, 1 at sub-point j + 1
      V sigma = (1.0 - frac) * rd(j) + frac * rd(j + 1);
      double f = 1.0;
      if (divideByStopping) {
        double stopping = (1.0 - frac) * point->GetSubPoint(j)->GetStoppingPower() + frac * point->GetSubPoint(j + 1)->GetStoppingPower();
        f = 1.0 / stopping / 1.0e24;
      }
      num += (half * wNum(e) * f) * sigma;
      den += half * wDen(e);
    }
  }
}
}  // namespace

/*!
 * Integrates the sub-points of this point into its fitted value for the experimental effect it carries.
 *
 * Gaussian beam spread (a constant sigma or an expression, both through TargetEffect::SigmaCM):
 *   Y(E0) = int K(E'; E0) sigma(E') dE' / int K(E'; E0) dE',  K a Gaussian of width sigma_b(E0) about E0,
 * normalised by its own integral over the sub-point window, so a constant cross section is reproduced
 * exactly however the window is truncated (before 2026-10 the constant-sigma form was not normalised and
 * came out 0.27 % low, and the expression form took its width at E' rather than at the point's energy).
 *
 * Target integration:  Y(E0) = (1/N) int_{E0-T}^{E0} <sigma>(E_d) / epsilon(E_d) dE_d  over the depth, where
 * <sigma>(E_d) is sigma itself, or its Gaussian average about E_d of width sqrt(sigma_b(E_d)^2 +
 * sigma_straggling(E0 - E_d)^2) when there is a beam spread and/or straggling, normalised as above.
 *
 * UDR (time-of-flight) and beam-profile kernels follow; both are normalised averages as well.
 * 2-point Gauss-Legendre on every sub-point interval throughout (GridQuadrature), exact for cubics and
 * indifferent to the adaptive grid's spacing.
 */

template <class V, class Read>
V EPoint::IntegrateTargetEffectImpl(const Config &configure, Read rd) {
  V yield = 0.0;
  TargetEffect *targetEffect = this->GetParentData()->GetTargetEffect(this->GetTargetEffectNum());
  const int numPoints = this->NumSubPoints();
  if (numPoints < 2) return V(0.0);
  const double eHigh = this->GetSubPoint(1)->GetCMEnergy();
  const double eLow = this->GetSubPoint(numPoints)->GetCMEnergy();
  const double centroid = this->GetCMEnergy();
  const double range = targetEffect->convolutionRange;
  auto one = [](double) { return 1.0; };
  // sigma(E) / (epsilon(E) 1e24) interpolated at one energy, for a depth node with no spread to average over
  auto yieldAt = [&](double e) {
    for (int j = 1; j < numPoints; j++) {
      double Ea = this->GetSubPoint(j)->GetCMEnergy(), Eb = this->GetSubPoint(j + 1)->GetCMEnergy();
      if (e <= Ea && e >= Eb) {
        double frac = (Ea - e) / (Ea - Eb);
        V sigma = (1.0 - frac) * rd(j) + frac * rd(j + 1);
        double stopping = (1.0 - frac) * this->GetSubPoint(j)->GetStoppingPower() + frac * this->GetSubPoint(j + 1)->GetStoppingPower();
        return (1.0 / stopping / 1.0e24) * sigma;
      }
    }
    return V(0.0);
  };
  // Areal density (atoms/cm^2) of target traversed from the grid's top node down to c.m. energy e, int dE / epsilon
  // over the sub-points (trapezoid in 1/epsilon; epsilon in c.m. units as the sub-points carry it, and dE_cm /
  // eps_cm = dE_lab / eps_lab).  For Bohr straggling.  A node with a non-positive stopping power contributes nothing.
  std::vector<double> arealCum;
  auto invEps = [](double eps) { return (eps > 0.0) ? 1.0 / eps : 0.0; };
  if (targetEffect->IsBohrStraggling()) {
    arealCum.assign(numPoints + 1, 0.0);
    for (int j = 1; j < numPoints; j++) {
      double Ea = this->GetSubPoint(j)->GetCMEnergy(), Eb = this->GetSubPoint(j + 1)->GetCMEnergy();
      double ia = invEps(this->GetSubPoint(j)->GetStoppingPower()), ib = invEps(this->GetSubPoint(j + 1)->GetStoppingPower());
      arealCum[j + 1] = arealCum[j] + 0.5 * (ia + ib) * (Ea - Eb);
    }
  }
  auto arealTo = [&](double e) {
    if (arealCum.empty() || e >= eHigh) return 0.0;
    for (int j = 1; j < numPoints; j++) {
      double Ea = this->GetSubPoint(j)->GetCMEnergy(), Eb = this->GetSubPoint(j + 1)->GetCMEnergy();
      if (e <= Ea && e >= Eb) {
        double frac = (Ea - e) / (Ea - Eb);
        double ia = invEps(this->GetSubPoint(j)->GetStoppingPower()), ib = invEps(this->GetSubPoint(j + 1)->GetStoppingPower());
        double ie = (1.0 - frac) * ia + frac * ib;
        return arealCum[j] + 0.5 * (ia + ie) * (Ea - e);
      }
    }
    return arealCum[numPoints];
  };

  if (targetEffect->IsTargetIntegration()) {
    const double surface = centroid;
    const double back = centroid - this->GetTargetThickness();
    V outer = 0.0;
    if (!targetEffect->IsGaussian() && !targetEffect->IsStraggling()) {
      // The grid was built to span the target (EData::GridExtent), so it is integrated whole, as it
      // always was: clipping it to [E0 - T, E0] here would drop a point whose segment energy shift
      // moved E0 away from the grid's own centre (13N, segment 1, lowest point).
      double unused = 0.0;
      GridQuadrature<V>(this, rd, eLow, eHigh, true, one, one, outer, unused);
    } else {
      const double x1 = -1.0 / sqrt(3.0), x2 = 1.0 / sqrt(3.0);
      for (int i = 1; i < numPoints; i++) {
        double Ea = this->GetSubPoint(i)->GetCMEnergy(), Eb = this->GetSubPoint(i + 1)->GetCMEnergy();
        if (Eb >= surface || Ea <= back) continue;
        double a = std::min(Ea, surface), b = std::max(Eb, back);
        double mid = 0.5 * (a + b), half = 0.5 * (a - b);
        if (half < 1.0e-10) continue;
        for (int gp = 0; gp < 2; gp++) {
          double eDepth = mid + half * (gp == 0 ? x1 : x2);
          double beam = targetEffect->IsGaussian() ? targetEffect->SigmaCM(eDepth, configure, cm_factor_) : 0.0;
          if (!(beam > 0.)) beam = 0.0;
          double straggle = targetEffect->IsBohrStraggling()
                                ? targetEffect->BohrStragglingSigma(arealTo(eDepth) - arealTo(surface), cm_factor_, bohr_z1_, bohr_z2_)
                                : targetEffect->StragglingSigma(surface - eDepth);
          double sigmaEff = std::sqrt(beam * beam + straggle * straggle);
          V inner = 0.0;
          if (sigmaEff > 1.0e-10) {
            double den = 0.0;
            auto k = [&](double e) { return TargetEffect::GaussianKernel(e, eDepth, sigmaEff); };
            GridQuadrature<V>(this, rd, std::max(eLow, eDepth - range * sigmaEff), std::min(eHigh, eDepth + range * sigmaEff), true, k, k, inner, den);
            inner = (den > 0.0) ? inner / den : V(0.0);
          } else {
            inner = yieldAt(eDepth);
          }
          outer += half * inner;
        }
      }
    }
    // a thick target (density 0) is the yield per incident particle, not divided by a density
    yield = targetEffect->IsThickTarget() ? outer : outer / (targetEffect->GetDensity() * 1.0e-24);
  } else if (targetEffect->IsGaussian()) {
    const double s = targetEffect->SigmaCM(centroid, configure, cm_factor_);
    if (!(s > 0.)) {
      // no spread: the bare value at the point, interpolated on the grid
      for (int j = 1; j < numPoints; j++) {
        double Ea = this->GetSubPoint(j)->GetCMEnergy(), Eb = this->GetSubPoint(j + 1)->GetCMEnergy();
        if (centroid <= Ea && centroid >= Eb) {
          double frac = (Ea - centroid) / (Ea - Eb);
          return (1.0 - frac) * rd(j) + frac * rd(j + 1);
        }
      }
      return V(0.0);
    }
    auto k = [&](double e) { return TargetEffect::GaussianKernel(e, centroid, s); };
    V num = 0.0;
    double den = 0.0;
    GridQuadrature<V>(this, rd, eLow, eHigh, false, k, k, num, den);
    yield = (den > 0.0) ? num / den : V(0.0);
  } else if (targetEffect->IsUdr()) {
    // User-defined numerical resolution function (SAMMY UDR).  The kernel
    // R(tau) is the density of time-of-flight delays at this point's energy;
    // a sub-point at true lab energy E' is seen at the delay
    //   tau(E') = t(E0) - t(E'),
    // so the fitted value is  sum_i w_i sigma_i / sum_i w_i  with
    //   w_i = int R(tau) phi_i(tau) dtau,
    // phi_i the hat function of sub-point i on the delay axis: the exact
    // integral of the piecewise-linear kernel against the linearly
    // interpolated cross section, as in SAMMY's Udr_Resb.  Each product is a
    // quadratic on every piece between breakpoints of either grid, so
    // Simpson's rule evaluates it exactly.
    const UdrKernel *kernel = this->GetUdrKernel();
    int numPoints = this->NumSubPoints();
    if (!kernel || kernel->tau.size() < 2 || numPoints < 2) {
      return V(0.0);
    }
    const double mass = this->GetUdrProjectileMass();
    const double path = targetEffect->GetUdrFlightPath();
    const double cmPerLab = (lab_energy_ != 0.0) ? cm_energy_ / lab_energy_ : 1.0;
    const double t0 = UdrKinematics::TimeOfFlight(lab_energy_, mass, path);
    // Sub-points run from high to low energy, i.e. from large to small delay;
    // collect them in ascending delay.
    std::vector<double> tauSub(numPoints);
    std::vector<V> sigmaSub(numPoints);
    for (int i = 1; i <= numPoints; i++) {
      const EPoint *sp = this->GetSubPoint(i);
      tauSub[numPoints - i] = t0 - UdrKinematics::TimeOfFlight(sp->GetCMEnergy() / cmPerLab, mass, path);
      sigmaSub[numPoints - i] = rd(i);
    }
    std::vector<double> w(numPoints, 0.0);
    const std::vector<double> &T = kernel->tau;
    const std::vector<double> &R = kernel->r;
    size_t j = 0;  // sub-point interval [tauSub[j], tauSub[j+1]]
    for (size_t k = 0; k + 1 < T.size(); k++) {
      double ta = T[k], tb = T[k + 1];
      if (tb <= tauSub.front() || ta >= tauSub.back()) continue;
      while (j + 1 < static_cast<size_t>(numPoints) && tauSub[j + 1] <= ta) j++;
      double lo = std::max(ta, tauSub.front());
      double hiAll = std::min(tb, tauSub.back());
      size_t jj = j;
      while (lo < hiAll && jj + 1 < static_cast<size_t>(numPoints)) {
        double hi = std::min(hiAll, tauSub[jj + 1]);
        if (hi > lo) {
          double sa = tauSub[jj], sb = tauSub[jj + 1], ds = sb - sa;
          auto rOf = [&](double x) { return R[k] + (R[k + 1] - R[k]) * (x - ta) / (tb - ta); };
          double mid = 0.5 * (lo + hi);
          double fl = rOf(lo), fm = rOf(mid), fh = rOf(hi);
          double h = (hi - lo) / 6.0;
          // hat of jj:  (sb - x)/ds ;  hat of jj+1:  (x - sa)/ds
          w[jj] += h * (fl * (sb - lo) + 4.0 * fm * (sb - mid) + fh * (sb - hi)) / ds;
          w[jj + 1] += h * (fl * (lo - sa) + 4.0 * fm * (mid - sa) + fh * (hi - sa)) / ds;
        }
        lo = hi;
        if (lo >= tauSub[jj + 1]) jj++;
      }
    }
    V numerator = 0.0;
    double denominator = 0.0;
    for (int i = 0; i < numPoints; i++) {
      numerator += w[i] * sigmaSub[i];
      denominator += w[i];
    }
    yield = (denominator > 0.0) ? numerator / denominator : V(0.0);
  } else if (targetEffect->IsBeamProfile()) {
    // Beam-profile kernel (photodissociation in a broad, skewed gamma beam
    // with an event-by-event reconstructed energy, cf. Haverson 2026 App. A):
    //
    //   Y = sum_i G(E_i) W(E_i) D(E_i) sigma(E_i) / sum_i G(E_i) W(E_i),
    //
    // G  the absolute beam energy profile (a sum of skewed Gaussians, not
    //    centred on the data point),
    // W  the fraction of the Gaussian detector resolution (sigma s) that
    //    lands inside the point's energy window [a, b], shifted with the
    //    segment's energy shift; 1 when the point has no window,
    // D  the detailed-balance factor sigma(gamma,a)/sigma(a,gamma) relative
    //    to its value at the point's own (unshifted) energy; 1 unless
    //    requested.  The measurement spreads its luminosity over G*W and
    //    averages sigma(gamma,a) = f_db * sigma(a,gamma); the datum is that
    //    average divided by f_db at the point's energy.  So D weights the
    //    numerator only: in the denominator too it would replace f_db(E0) by
    //    its beam average and cancel its own normalisation (a 1-5 % error
    //    when the kernel's centroid is tens of keV off the point energy).
    const double s = cm_factor_ * targetEffect->GetBeamTpcSigma();  // the line's lab sigma, this pair's frame
    const bool window = this->HasBinWindow();
    // The energy shift moved the point (and the sub-point grid); the window
    // is stored unshifted, so shift it by the same c.m. amount.
    double cmPerLab = (lab_energy_ != 0.0) ? cm_energy_ / lab_energy_ : 1.0;
    double originalCM = original_energy_ * cmPerLab;
    double delta = cm_energy_ - originalCM;
    double a = bin_low_cm_ + delta;
    double b = bin_high_cm_ + delta;
    const bool photo = targetEffect->IsBeamPhotodissociation() && photo_mcn_ > 0.0;
    // g(E) = E / E_gamma(E)^2 is the energy dependence of the detailed-balance
    // factor mu c^2 E / E_gamma^2, with the photon energy including the
    // recoil of the compound nucleus: E_gamma^2/(2 M) - E_gamma + (E + Q) = 0.
    auto dbWeight = [&](double e) {
      double sum = e + photo_qgamma_;
      double arg = 1.0 - 2.0 * sum / photo_mcn_;
      double egamma = (arg > 0.0) ? photo_mcn_ * (1.0 - std::sqrt(arg)) : sum;
      return (egamma > 0.0) ? e / (egamma * egamma) : 0.0;
    };
    double dbNorm = photo ? dbWeight(originalCM) : 1.0;
    auto kernel = [&](double e) {
      double k = targetEffect->BeamProfileWeight(e / cm_factor_);  // the profile is lab; e is c.m.
      if (window) {
        if (s > 0.0)
          k *= 0.5 * (std::erf((b - e) / (s * std::sqrt(2.0))) - std::erf((a - e) / (s * std::sqrt(2.0))));
        else
          k *= (e >= a && e <= b) ? 1.0 : 0.0;
      }
      return k;
    };
    auto weighted = [&](double e) {
      return kernel(e) * ((photo && dbNorm > 0.0) ? dbWeight(e) / dbNorm : 1.0);
    };
    V numerator = 0.0;
    double denominator = 0.0;
    GridQuadrature<V>(this, rd, eLow, eHigh, false, weighted, kernel, numerator, denominator);
    yield = (denominator > 0.0) ? numerator / denominator : V(0.0);
  }
  return yield;
}
void EPoint::IntegrateTargetEffect(const Config &configure) {
  this->SetFitCrossSection(this->IntegrateTargetEffectImpl<double>(
      configure, [this](int j) { return this->GetSubPoint(j)->GetFitCrossSection(); }));
}

bool EPoint::TargetEffectWeights(const Config &configure, std::vector<double> &c) {
  const int n = this->NumSubPoints();
  c.assign(n > 0 ? n : 0, 0.0);
  TELinearForm y = this->IntegrateTargetEffectImpl<TELinearForm>(
      configure, [](int j) { return TELinearForm::Unit(j); });
  if (!y.Valid()) return false;
  for (const auto &t : y.Terms()) {
    if (t.first < 1 || t.first > n) return false;
    c[t.first - 1] += t.second;
  }
  return true;
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

void EPoint::SetAngularDists(vector_r dists) {
  angularDists_.clear();
  angularDists_ = dists;
}

/*!
 * Returns a pointer to the parent EData object.
 */

EData *EPoint::GetParentData() const {
  return parentData_;
}

/*!
 * Returns a pointer to a point mapped to the current point specified by a position in the mapped point vector.
 */

EPoint *EPoint::GetLocalMappedPoint(int mappedPointNum) const {
  return local_mapped_points_[mappedPointNum - 1];
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

/*!
 * Returns a reference to the vector of pointers to mapped EPoint objects.
 */

std::vector<EPoint *> &EPoint::GetMappedPoints() {
  return local_mapped_points_;
}
