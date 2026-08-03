#include "AMatrixFunc.h"
#include "AngCoeff.h"
#include "CNuc.h"
#include "Config.h"
#include "CoulFunc.h"
#include "DataLine.h"
#include "ECIntegral.h"
#include "ECAmplitudeCache.h"
#include "EData.h"
#include "ESegment.h"
#include "RMatrixFunc.h"
#include "ShftFunc.h"
#include "TargetEffect.h"
#include "Straggling.h"
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
  entrance_key_=parent->GetEntranceKey();
  exit_key_=parent->GetExitKey();
  cm_angle_=dataLine.angle();
  lab_angle_=dataLine.angle();
  original_energy_=dataLine.energy();
  double shiftedEnergy = dataLine.energy(); //+ parent->GetEnergyShift();
  // Don't allow energies below 0.005 MeV
  if(shiftedEnergy < 0.005) {
    shiftedEnergy = dataLine.energy();
  }
  cm_energy_=dataLine.energy();
  lab_energy_=shiftedEnergy;
  excitation_energy_=dataLine.energy();
  parentSegment_=parent;
  segment_key_=parent->GetSegmentKey();
  cm_crosssection_=dataLine.crossSection();
  cm_dcrosssection_=dataLine.error();
  lab_crosssection_=dataLine.crossSection();
  lab_dcrosssection_=dataLine.error();
  geofactor_=0.;
  fitcrosssection_=0.;
  fitE1crosssection_=0.;
  fitE2crosssection_=0.; 
  sfactorconv_=0.;
  is_differential_=parent->IsDifferential();
  is_phase_=parent->IsPhase();
  is_ang_dist_=parent->IsAngularDist();
  max_ang_dist_order_=parent->GetMaxAngDistOrder();
  j_value_=parent->GetJ();
  l_value_=parent->GetL();
  isUPOS_=parent->IsUPOS();
  secondaryDecayL_=parent->GetSecondaryDecayL();
  Ic_=parent->GetIc();
  delta_=parent->GetDelta();
  is_mapped_=false;
  targetEffectNum_=0;
  parentData_=NULL;
  stoppingPower_=0.0;
  angleKinFactor_=1.0;
  crossSectionKinFactor_=1.0;
}

/*!
 * This constructor is used if the data point is to be created with no experimental data.
 * A pointer to the parent segment is passed to the constructor for the intialization of
 * the EPoint object.
 */

EPoint::EPoint(double angle, double energy, ESegment* parent) {
  entrance_key_=parent->GetEntranceKey();
  exit_key_=parent->GetExitKey();
  lab_angle_=angle;
  cm_angle_=angle;
  original_energy_=energy;
  lab_energy_=energy;
  cm_energy_=energy;
  excitation_energy_=energy;
  parentSegment_=parent;
  segment_key_=parent->GetSegmentKey();
  cm_crosssection_=0.;
  cm_dcrosssection_=0.1;
  lab_crosssection_=0.;
  lab_dcrosssection_=0.1;
  geofactor_=0.;
  fitcrosssection_=0.;
  fitE1crosssection_=0.;
  fitE2crosssection_=0.; 
  sfactorconv_=0.;
  is_differential_=parent->IsDifferential();
  is_phase_=parent->IsPhase();
  is_ang_dist_=parent->IsAngularDist();
  max_ang_dist_order_=parent->GetMaxAngDistOrder();
  j_value_=parent->GetJ();
  l_value_=parent->GetL();
  isUPOS_=parent->IsUPOS();
  secondaryDecayL_=parent->GetSecondaryDecayL();
  Ic_=parent->GetIc();
  delta_=parent->GetDelta();
  is_mapped_=false;
  targetEffectNum_=0;
  parentData_=NULL;
  stoppingPower_=0.0;
  angleKinFactor_=1.0;
  crossSectionKinFactor_=1.0;
}

/*!
 * This constructor is used if the data point is to be created with no experimental data
 * and no parent segment.  Such a point is used if AZURE is to be called as an energy-
 * dependent function for reaction rate and target effect calculations with dynamic
 * integration. All options must be set manually.
 */

EPoint::EPoint(double angle, double energy, int entranceKey, 
	       int exitKey, bool isDifferential, bool isPhase, bool isAngularDist, double jValue, int lValue, int maxAngDistOrder) {
  entrance_key_=entranceKey;
  exit_key_=exitKey;
  lab_angle_=angle;
  cm_angle_=angle;
  original_energy_=energy;
  lab_energy_=energy;
  cm_energy_=energy;
  excitation_energy_=energy;
  parentSegment_=nullptr;
  cm_crosssection_=0.;
  cm_dcrosssection_=0.1;
  lab_crosssection_=0.;
  lab_dcrosssection_=0.1;
  geofactor_=0.;
  fitcrosssection_=0.;
  fitE1crosssection_=0.;
  fitE2crosssection_=0.; 
  sfactorconv_=0.;
  is_differential_=isDifferential;
  is_phase_=isPhase;
  is_ang_dist_=isAngularDist;
  max_ang_dist_order_=maxAngDistOrder;
  j_value_=jValue;
  l_value_=lValue;
  isUPOS_=false;
  secondaryDecayL_=0;
  Ic_=0.0;
  delta_=0.0;
  is_mapped_=false;
  targetEffectNum_=0;
  parentData_=NULL;
  stoppingPower_=0.0;
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
  if(GetTargetEffectNum()!=0) return true;
  else return false;
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
  return legendreP_.size()-1;
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
  return penetrabilities_[jGroupNum-1][channelNum-1];
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
  return lo_elements_[jGroupNum-1][channelNum-1];
}

/*!
 * Returns the factor \f$ \exp(\omega_c) \f$ where \f$ \omega_c \f$ is the Coulomb 
 * phase shift.  The  channel us specified
 * by positions in the JGroup and subsequent AChannel vectors.
 */

complex EPoint::GetExpCoulombPhase(int jGroupNum, int channelNum) const {
  return coulombphase_[jGroupNum-1][channelNum-1];
}

/*!
 * Returns the factor \f$ \exp(\delta_c) \f$ where \f$ \delta_c \f$ is the hard sphere 
 * phase shift.  The  channel us specified
 * by positions in the JGroup and subsequent AChannel vectors.
 */

complex EPoint::GetExpHardSpherePhase(int jGroupNum, int channelNum) const {
  return hardspherephase_[jGroupNum-1][channelNum-1];
}

/*!
 * Returns the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

complex EPoint::GetCoulombAmplitude() const {
  return coulombamplitude_;
}

/*!
 * Returns the external capture amplitude for a given external reaction pathway
 * specified by positions in the KGroup and subsequent ECMGroup vectors.
 */

complex EPoint::GetECAmplitude(int kGroupNum, int ecMGroupNum) const {
  return ec_amplitudes_[kGroupNum-1][ecMGroupNum-1];
}

/*!
 * Returns the external capture amplitude with current energy (including shifts)
 * This method uses the ECAmplitudeCache to get interpolated values.
 */
complex EPoint::GetECAmplitudeWithShift(int kGroupNum, int ecMGroupNum, CNuc *theCNuc, const Config& configure) const {
  // Check if we have valid indices
  //if(kGroupNum < 1 || kGroupNum > ec_amplitudes_.size()) {
  //  return complex(0.0, 0.0);
  //}
  //if(ecMGroupNum < 1 || ecMGroupNum > ec_amplitudes_[kGroupNum-1].size()) {
  //  return complex(0.0, 0.0);
  //}
  
  // Get current energy including shifts
  double currentEnergy = this->GetCMEnergy();
  
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
  if(this->IsMapped()) {
    thisMap.segment=energy_map_.segment;
    thisMap.point=energy_map_.point;
  }
  else {
    thisMap.segment=0;
    thisMap.point=0;
  }
  return thisMap;
}

/*!
 * Initializes a data point to be used in a calculation.  Initilization is done before the fitting process
 * to calculate all energy dependent quantities that do no rely on the R-Matrix fit parameters.
 */

void EPoint::Initialize(CNuc *compound,const Config &configure) {
  this->CalcEDependentValues(compound,configure);
  if(this->IsDifferential()) {
    // Get TargetEffect with Q-coefficients if applicable
    TargetEffect* effect = NULL;
    EData* parentData = this->GetParentData();
    if(parentData && this->IsTargetEffect()) {
      TargetEffect* te = parentData->GetTargetEffect(this->GetTargetEffectNum());
      if(te && te->IsQCoefficients()) {
        effect = te;
      }
    }
    this->CalcLegendreP(configure.maxLOrder,compound,effect);
  }
  this->CalcCoulombAmplitude(compound);
  if(configure.paramMask & Config::USE_EXTERNAL_CAPTURE) this->CalculateECAmplitudes(compound,configure);
}

/*!
 * Calculates lab frame value from center of mass value.
 */

double EPoint::ConvertCMValue(double value, PPair *pPair) {
  value=value*(pPair->GetM(1)+pPair->GetM(2))/(pPair->GetM(2));
  return value;
}

/*!
 * Calculates center of mass energy for a given lab energy.
 */

double EPoint::ConvertLabValue(double value, PPair *pPair) {
  value=value*(pPair->GetM(2))/(pPair->GetM(1)+pPair->GetM(2));
  return value;
}

/*!
 * Calculates center of mass energy.  When a data point is initialized, the same energy is copied into
 * the attributes for center of mass and lab energy.  If this function is called, the center of mass energy
 * attribute is overwritten with the value calculated from the lab energy attribute.
 */

void EPoint::ConvertLabEnergy(PPair *pPair) {
  cm_energy_=this->GetLabEnergy()*
    (pPair->GetM(2))/
    (pPair->GetM(1)+pPair->GetM(2));
  excitation_energy_=cm_energy_+pPair->GetSepE();
}

/*!
 * Calculates excitation energy. When a data point is initialized, only the excitation energy is calculated
 * if the data is already in the center of mass frame. Else it is calculated by ConvertLabEnergy.  
 */

void EPoint::ConvertExcitationEnergy(PPair *pPair) {
  excitation_energy_=cm_energy_+pPair->GetSepE();
}

/*!
 * Calculates the total decay energy from the light particle decay energy, assuming the parent nucleus was at
 * rest when it decayed.  When a data point is initialized, the same energy is copied into
 * the attributes for center of mass and lab energy.  If this function is called, the center of mass energy
 * attribute is overwritten with the value calculated from the lab energy attribute.
 */

void EPoint::ConvertDecayEnergy(PPair *pPair) {
  cm_energy_=this->GetLabEnergy()/
    (pPair->GetM(2))*
    (pPair->GetM(1)+pPair->GetM(2));
  excitation_energy_=cm_energy_+pPair->GetSepE();
}

/*!
 * Calculates center of mass angles.  When a data point is initialized, the same angle is copied into
 * the attributes for center of mass and lab angles.  If this function is called, the center of mass angle
 * attribute is overwritten with the value calculated from the lab angle attribute.  This version of the
 * overloaded function is for scattering channels.
 */

void EPoint::ConvertLabAngle(PPair *pPair) {
  cm_angle_=this->GetLabAngle()+180./pi*asin(pPair->GetM(1)/pPair->GetM(2)*sin(pi/180.*this->GetLabAngle()));
}

/*!
 * Calculates center of mass angles.  When a data point is initialized, the same angle is copied into
 * the attributes for center of mass and lab angles.  If this function is called, the center of mass angle
 * attribute is overwritten with the value calculated from the lab angle attribute.  This version of the
 * overloaded function is for non-elastic particle channels.
 */

void EPoint::ConvertLabAngle(PPair *entrancePair, PPair *exitPair, const Config& configure) {
  double qValue=entrancePair->GetSepE()+entrancePair->GetExE()-exitPair->GetSepE()-exitPair->GetExE();
  double a13=(entrancePair->GetM(1)*exitPair->GetM(1))*this->GetLabEnergy()/(this->GetLabEnergy()+qValue)/
    (entrancePair->GetM(1)+entrancePair->GetM(2))/(exitPair->GetM(1)+exitPair->GetM(2));
  double a24=(entrancePair->GetM(2)*exitPair->GetM(2))*(1+entrancePair->GetM(1)/entrancePair->GetM(2)*qValue/(this->GetLabEnergy()+qValue))/
    (entrancePair->GetM(1)+entrancePair->GetM(2))/(exitPair->GetM(1)+exitPair->GetM(2));

  if(a13>a24) {
    double thetaMax=asin(sqrt(a24/a13))*180./pi;
    if(thetaMax<this->GetLabAngle()) configure.outStream << std::endl << "Lab Angle (" << this->GetLabAngle() 
					       << " degrees) is not kinematically possible.  Maximum angle is " 
					       << thetaMax << " degrees." << std::endl;
    assert(thetaMax>=this->GetLabAngle());
  }

  double E3PerEt=a13*pow(cos(this->GetLabAngle()*pi/180.)+sqrt(a24/a13-pow(sin(this->GetLabAngle()*pi/180.),2.0)),2.0);
  double tempE3PerEt = a13*pow(cos((this->GetLabAngle()+0.001)*pi/180.)+sqrt(a24/a13-pow(sin((this->GetLabAngle()+0.001)*pi/180.),2.0)),2.0);
  double slope = sqrt(tempE3PerEt/a24)*sin((this->GetLabAngle()+0.001)*pi/180.)-sqrt(E3PerEt/a24)*sin(this->GetLabAngle()*pi/180.);
  bool switchDomain=false;
  if(slope<0.) switchDomain=true;
  
  cm_angle_=180./pi*asin(sqrt(E3PerEt/a24)*sin(this->GetLabAngle()*pi/180.));
  if(switchDomain) cm_angle_ = 180.- cm_angle_;
  //define angle kin factor to convert back to lab frame
  angleKinFactor_=this->GetLabAngle()/cm_angle_;
  this->SetAngleKinFactor(angleKinFactor_);
}

/*!
 * Calculates center of mass angles as in ConvertLabAngle but uses the reletivistic equations of Iliadis C.37 and C.38.  
 * When a data point is initialized, the same angle is copied into
 * the attributes for center of mass and lab angles.  If this function is called, the center of mass angle
 * attribute is overwritten with the value calculated from the lab angle attribute.  This version of the
 * overloaded function is for non-elastic particle channels.
 */

void EPoint::ConvertCMAngle(PPair *entrancePair, PPair *exitPair, const Config& configure) {
  double qValue=entrancePair->GetSepE()+entrancePair->GetExE()-exitPair->GetSepE()-exitPair->GetExE();
  double gamma=sqrt(entrancePair->GetM(1)*exitPair->GetM(1)*this->GetLabEnergy()/(exitPair->GetM(1)*
   (exitPair->GetM(1)+exitPair->GetM(2))*qValue+exitPair->GetM(1)*(exitPair->GetM(1)+exitPair->GetM(2)
   -entrancePair->GetM(1))*this->GetLabEnergy()));
  double temp_lab_angle=180./pi*acos((gamma+cos(this->GetCMAngle()*pi/180.))/sqrt(1+gamma*gamma+2.*gamma*cos(this->GetCMAngle()*pi/180.)));
//  std::cout<<temp_lab_angle<<std::endl;
}

/*!
 * Calculates center of mass angles for gamma rays for the relativistic correction.
 */
void EPoint::ConvertLabAngleGammas(PPair *entrancePair) {
  double beta = sqrt(this->GetLabEnergy()*(this->GetLabEnergy()+2.0*entrancePair->GetM(1)*uconv))/
    ((entrancePair->GetM(1)+entrancePair->GetM(2))*uconv + this->GetLabEnergy());
  cm_angle_ = 180./pi*acos((cos(this->GetLabAngle()*pi/180.)-beta)/(1-beta*cos(this->GetLabAngle()*pi/180.)));
}

/*!
 * Calculates center of mass cross sections.  When a data point is initialized, the same cross section is copied into
 * the attributes for center of mass and lab cross section.  If this function is called, the center of mass cross section
 * attribute is overwritten with the value calculated from the lab cross section attribute.
 */

void EPoint::ConvertCrossSection(PPair *entrancePair, PPair *exitPair) {
  double conversionFactor;
  if(this->GetLabAngle()==0.0 || this->GetLabAngle()==180.0) {
    double m1=entrancePair->GetM(1);
    double m2=entrancePair->GetM(2);
    double m3=exitPair->GetM(1);
    double m4=exitPair->GetM(2);
    double e1=this->GetLabEnergy();
    double qValue=entrancePair->GetSepE()+entrancePair->GetExE()-exitPair->GetSepE()-exitPair->GetExE();
    double et=e1+qValue;
    double a=m1*m4*e1/(m1+m2)/(m3+m4)/et;
    double b=m1*m3*e1/(m1+m2)/(m3+m4)/et;
    double c=m2*m3/(m1+m2)/(m3+m4) * (1+m1*qValue/m2/et);
    double d=m2*m4/(m1+m2)/(m3+m4) * (1+m1*qValue/m2/et);
    double e3et=b+d+2*pow(a*c,0.5)*cos(pi/180.*this->GetCMAngle());
    conversionFactor=pow(a*c,0.5)*pow(d/b-pow(sin(this->GetLabAngle()*pi/180.),2.0),0.5)/e3et;
  }
  else {
    conversionFactor=pow(sin(pi/180.*this->GetLabAngle())/sin(pi/180.*this->GetCMAngle()),2.0)*cos(pi/180.*(this->GetCMAngle()-this->GetLabAngle()));
  }
  cm_crosssection_=this->GetLabCrossSection()*conversionFactor;
  cm_dcrosssection_=this->GetLabCrossSectionError()*conversionFactor;
  crossSectionKinFactor_=conversionFactor;
  this->SetCrossSectionKinFactor(crossSectionKinFactor_);
}

/* Calculates the conversion factor */

double EPoint::CalculateCrossSectionConversionFactor(PPair *entrancePair, PPair *exitPair) {
  double conversionFactor;
  if(this->GetLabAngle()==0.0 || this->GetLabAngle()==180.0) {
    double m1=entrancePair->GetM(1);
    double m2=entrancePair->GetM(2);
    double m3=exitPair->GetM(1);
    double m4=exitPair->GetM(2);
    double e1=this->GetLabEnergy();
    double qValue=entrancePair->GetSepE()+entrancePair->GetExE()-exitPair->GetSepE()-exitPair->GetExE();
    double et=e1+qValue;
    double a=m1*m4*e1/(m1+m2)/(m3+m4)/et;
    double b=m1*m3*e1/(m1+m2)/(m3+m4)/et;
    double c=m2*m3/(m1+m2)/(m3+m4) * (1+m1*qValue/m2/et);
    double d=m2*m4/(m1+m2)/(m3+m4) * (1+m1*qValue/m2/et);
    double e3et=b+d+2*pow(a*c,0.5)*cos(pi/180.*this->GetCMAngle());
    conversionFactor=pow(a*c,0.5)*pow(d/b-pow(sin(this->GetLabAngle()*pi/180.),2.0),0.5)/e3et;
  }
  else {
    conversionFactor=pow(sin(pi/180.*this->GetLabAngle())/sin(pi/180.*this->GetCMAngle()),2.0)*cos(pi/180.*(this->GetCMAngle()-this->GetLabAngle()));
  }
  return conversionFactor;
}

/*!
 * Calculates center of mass cross sections for gamma rays.
 */

void EPoint::ConvertCrossSectionGammas(PPair *entrancePair) {

  double conversionFactor;

  double beta = sqrt(this->GetLabEnergy()*(this->GetLabEnergy()+2.0*entrancePair->GetM(1)*uconv))/
    ((entrancePair->GetM(1)+entrancePair->GetM(2))*uconv + this->GetLabEnergy());
  
  conversionFactor = pow(1.0-beta,2.0)/(1-beta*cos(this->GetLabAngle()*pi/180.0));
  cm_crosssection_=this->GetLabCrossSection()*conversionFactor;
  cm_dcrosssection_=this->GetLabCrossSectionError()*conversionFactor;
}

/* Calculates the conversion factor for gamma rays */

double EPoint::CalculateCrossSectionGammaConversionFactor(PPair *entrancePair) {

  double conversionFactor;

  double beta = sqrt(this->GetLabEnergy()*(this->GetLabEnergy()+2.0*entrancePair->GetM(1)*uconv))/
    ((entrancePair->GetM(1)+entrancePair->GetM(2))*uconv + this->GetLabEnergy());
  
  conversionFactor = pow(1.0-beta,2.0)/(1-beta*cos(this->GetLabAngle()*pi/180.0));
  return conversionFactor;
}

/*!
 * Calculates Legendre polynomial coefficient reference frame correction up to a maximum order.  The factors are added, in order, to a vector.
 */
/*
void EPoint::ConvertLabQCoeff(int maxL, TargetEffect* targetEffect, PPair *entrancePair, PPair *exitPair) {
  double m1=entrancePair->GetM(1);
  double m2=entrancePair->GetM(2);
  double m3=exitPair->GetM(1);
  double m4=exitPair->GetM(2);
  double e1=this->GetLabEnergy();
  double qValue=entrancePair->GetSepE()+entrancePair->GetExE()-exitPair->GetSepE()-exitPair->GetExE();
  double gamma = m1*m3/m2/m4*e1/(e1+qValue);
  for (int L=0;L<=maxL;lOrder++){
    double U=0;
    for (int L_p=0;L_p<=maxL;L_p++){
      if(L==L_p) U = 1.0-gamma*gamma*L_p*(L_p+1.0)*(2.0*L_p*L_p+2.0*L_p-3.0)/2.0/(2.0*L_p-1)/(2.0*L_p+3.0);
      else if (L_p == L+1) U = 1.0*gamma*L_p*(L_p-1.0)/(2.0*L_p-1.0);
      else if (L_p == L-1) U = -1.0*gamma*(L_p+1.0)*(L_p+2.0)/(2.0*L_p+3.0);
      else if (L_p == L-2) U = gamma*gamma*(L_p+1.0)*(L_p+2.0)*(L_p+3.0)*(L_p+3.0)/2.0/(2.0*L_p+3.0)/(2.0*L_p+5.0);
      else if (L_p == L+2) U = gamma*gamma*L_p*(L_p-1)*(L_p-2)*(L_p-2)/2/(2*L_p-1)/(2*L_p-3);
    }
    U = U*(2.0*L+1.0)/(2.0*L_p+1.0);
    std::cout<<U<<std::endl;
  }
}
*/

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
  geofactor_=geoFactor;
}

/*!
 * Sets the calculated AZURE cross section.
 */

void EPoint::SetFitCrossSection(double crossSection) {
  fitcrosssection_=crossSection;
}

void EPoint::SetFitE1CrossSection(double e1crossSection) {
  fitE1crosssection_= e1crossSection;
}

void EPoint::SetFitE2CrossSection(double e2crossSection) {
  fitE2crosssection_= e2crossSection;
}

/*!
 * Sets the multiplicative conversion from cross section to s-factor.
 */

void EPoint::SetSFactorConversion(double conversion) {
  sfactorconv_=conversion;
}

/*!
 * Sets the lab energy to the given value;
 */

void EPoint::SetLabEnergy(double energy) {
  lab_energy_=energy;
}

/*!
 * Sets the CM energy to the given value;
 */

void EPoint::SetCMEnergy(double energy) {
  cm_energy_=energy;
}

/*!
 * Sets the excitation energy to the given value;
 */

void EPoint::SetExcitationEnergy(double energy) {
  excitation_energy_=energy;
}

/*!
 * Sets the lab angle to the given value;
 */

void EPoint::SetLabAngle(double angle) {
  lab_angle_=angle;
}

/*!
 * Sets the CM angle to the given value;
 */

void EPoint::SetCMAngle(double angle) {
  cm_angle_=angle;
}

/*!
 * Sets the exit key to the given value;
 */

void EPoint::SetExitKey(int key) {
  exit_key_=key;
}

/*!
 * Sets the entrance key for the data point
 */
void EPoint::SetEntranceKey(int key) {
  entrance_key_=key;
}

/*!
 * Sets the kinematic factor to convert angle from lab to cm angle;
 */

void EPoint::SetAngleKinFactor(double anglekinfactor) {
  angleKinFactor_=anglekinfactor;
}

/*!
 * Sets the kinematic factor to convert cross section from lab to cm angle;
 */

void EPoint::SetCrossSectionKinFactor(double crosssectionkinfactor) {
  crossSectionKinFactor_=crosssectionkinfactor;
}

/*!
 * Calculates Legendre polynomials up to a maximum order.  The polynomials are added, in order, to a vector.
 */

void EPoint::CalcLegendreP(int maxL, CNuc *theCNuc, TargetEffect* targetEffect) {
  std::vector<double> Qcm(maxL, 0.0);
  //double Qcm[maxL]={1.0};
  if(targetEffect && targetEffect->NumQCoefficients()>0 && theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()))->GetPType()==0){
    PPair *entrancePair=theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
    PPair *exitPair=theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetExitKey()));
    double m1=entrancePair->GetM(1);
    double m2=entrancePair->GetM(2);
    double m3=exitPair->GetM(1);
    double m4=exitPair->GetM(2);
    double e1=this->GetLabEnergy();
    double qValue=entrancePair->GetSepE()+entrancePair->GetExE()-exitPair->GetSepE()-exitPair->GetExE();
    double gamma = m1*m3/m2/m4*e1/(e1+qValue);
    int numQ = targetEffect->NumQCoefficients();
    std::vector<std::vector<double> > U(numQ, std::vector<double>(numQ, 0.0));
    //double U[numQ][numQ]={0.0};
    for (int lOrder=0;lOrder<numQ;lOrder++){    
      for (int lOrder_p=0;lOrder_p<numQ;lOrder_p++){
        double tempU=0.0;
        if(lOrder_p==lOrder) tempU = 1.0-gamma*gamma*lOrder_p*(lOrder_p+1.0)*(2.0*lOrder_p*lOrder_p+2.0*lOrder_p-3.0)/2.0/(2.0*lOrder_p-1)/(2.0*lOrder_p+3.0);
        else if (lOrder_p == lOrder+1) tempU = 1.0*gamma*lOrder_p*(lOrder_p-1.0)/(2.0*lOrder_p-1.0);
        else if (lOrder_p == lOrder-1) tempU = -1.0*gamma*(lOrder_p+1.0)*(lOrder_p+2.0)/(2.0*lOrder_p+3.0);
        else if (lOrder_p == lOrder-2) tempU = gamma*gamma*(lOrder_p+1.0)*(lOrder_p+2.0)*(lOrder_p+3.0)*(lOrder_p+3.0)/2.0/(2.0*lOrder_p+3.0)/(2.0*lOrder_p+5.0);
        else if (lOrder_p == lOrder+2) tempU = gamma*gamma*lOrder_p*(lOrder_p-1.0)*(lOrder_p-2.0)*(lOrder_p-2.0)/2.0/(2.0*lOrder_p-1.0)/(2.0*lOrder_p-3.0);
        else tempU = 0.0;
        U[lOrder][lOrder_p] = tempU*(2.0*lOrder+1.0)/(2.0*lOrder_p+1.0);
//        std::cout<<"test"<<std::endl;              
      }    
    }
    std::vector<double> B(numQ, 0.0);
    //double B[numQ] = {0.0};
    for (int lOrder=0;lOrder<numQ;lOrder++){    
      for (int lOrder_p=0;lOrder_p<numQ;lOrder_p++){
        B[lOrder] += U[lOrder][lOrder_p]*targetEffect->GetQCoefficient(lOrder_p);
      }
//      std::cout<<"B = " << B[lOrder] <<", Q = "<<targetEffect->GetQCoefficient(lOrder)<<" numQ = "<<numQ<<std::endl;
    }
    for (int lOrder=0;lOrder<maxL;lOrder++){
      if(lOrder<numQ){
        Qcm[lOrder] = B[lOrder];
        //std::cout<<Qcm[lOrder]<<std::endl;
      }
      else {
        Qcm[lOrder] = 1.;
        //std::cout<<Qcm[lOrder]<<std::endl;
      }    
    } 
  }
   
  double x=cos(this->GetCMAngle()*pi/180.0);
  if(maxL>=0) {
    if(targetEffect && targetEffect->NumQCoefficients()>0)
      this->AddLegendreP(Qcm[0]);
    else this->AddLegendreP(1.0);
    double polyMinusTwo=1.0;
    if(maxL>=1) {
      if(targetEffect && targetEffect->NumQCoefficients()>1)
	this->AddLegendreP(x*Qcm[1]);
      else this->AddLegendreP(x);
      double polyMinusOne=x;
      if(maxL>=2) {
	for(int lOrder=2;lOrder<=maxL;lOrder++) {
	  double poly=(2.0*lOrder-1.0)/lOrder*x*polyMinusOne-
	    (lOrder-1.0)/lOrder*polyMinusTwo;
	  if(targetEffect && targetEffect->NumQCoefficients()>lOrder)
	    this->AddLegendreP(poly*Qcm[lOrder]);
	  else this->AddLegendreP(poly);
	  polyMinusTwo=polyMinusOne;
	  polyMinusOne=poly;
	}
      }
    }
  }
  for(int i=1;i<=this->NumSubPoints();i++) {
    this->GetSubPoint(i)->CalcLegendreP(maxL, theCNuc, targetEffect);
  }
}



/*!
 * This function calculates several energy dependent quantities simultaniously.
 * This includes the geometrical cross section, the s-factor conversion, the \f$ L_o \f$ matrix
 * elements, the square root of the penetrability, and the exponentials of the Coulomb phase shifts
 * and hard sphere phase shfits.
 */

void EPoint::CalcEDependentValues(CNuc *theCNuc, const Config& configure) {
  PPair *entrancePair=theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
  PPair *exitPair=theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetExitKey()));

  double inEnergy;
  double geofactor;
  double sfactorconv;
  if(theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()))->GetPType()==20)  {
    inEnergy=this->GetCMEnergy()+exitPair->GetSepE()+exitPair->GetExE();
    geofactor=1.;
    sfactorconv=1.;
  } else {
    inEnergy=this->GetCMEnergy()+entrancePair->GetSepE()+entrancePair->GetExE();
    geofactor=pi*pow(hbarc,2.)/(2*entrancePair->GetRedMass()*uconv*this->GetCMEnergy());
    sfactorconv=this->GetCMEnergy()*exp(2*pi*sqrt(uconv/2.)*fstruc*entrancePair->GetZ(1)*
					entrancePair->GetZ(2)*sqrt(entrancePair->GetRedMass()
								   /this->GetCMEnergy()));
  }
  this->SetGeometricalFactor(geofactor);
  this->SetSFactorConversion(sfactorconv);

  for(int j=1;j<=theCNuc->NumJGroups();j++) {
    if(theCNuc->GetJGroup(j)->IsInRMatrix()) {
      JGroup *theJGroup=theCNuc->GetJGroup(j);
      for(int ch=1;ch<=theJGroup->NumChannels();ch++) {
	AChannel *theChannel=theJGroup->GetChannel(ch);
	PPair *thePair=theCNuc->GetPair(theChannel->GetPairNum());
	int lValue=theChannel->GetL();
	double localEnergy=inEnergy-thePair->GetSepE()-thePair->GetExE();
	if(thePair->GetPType()==0) {
	  if(localEnergy<=0.0) {
	    ShftFunc theShiftFunction(thePair);
	    double localShift=theShiftFunction(lValue,inEnergy);
	    double boundary=theChannel->GetBoundaryCondition();
	    complex loElement(localShift-boundary,0.0);
	    this->AddLoElement(j,ch,loElement);
	    this->AddSqrtPenetrability(j,ch,0.0);
	    this->AddExpCoulombPhase(j,ch,1.0);
	    this->AddExpHardSpherePhase(j,ch,1.0);
//            if(ch==1){
//            std::ofstream coul_check;
//            coul_check.open("coul_check_neg.chk",std::ios::app);
//            coul_check << this->GetCMEnergy() << "  " << lValue << "  " << eta << "  " << coul.F << "  " << coul.dF << "  " << coul.G << "  " << coul.dG << std::endl;
//            coul_check << this->GetEntranceKey() << "," << this->GetCMEnergy() << "," << lValue << "," << 0.0 << "," << localShift << std::endl;
//            coul_check.close();
//	    }
	  } else {
	    CoulFunc theCoulombFunction(thePair,!!(configure.paramMask&Config::USE_GSL_COULOMB_FUNC));
	    double radius=thePair->GetChRad();
	    double localPene=theCoulombFunction.Penetrability(lValue,radius,localEnergy);
	    double localShift=theCoulombFunction.PEShift(lValue,radius,localEnergy);
	    double boundary=theChannel->GetBoundaryCondition();
	    complex loElement(localShift-boundary,localPene);
	    double redmass=thePair->GetRedMass();
	    double eta=sqrt(uconv/2.)*fstruc*thePair->GetZ(1)*thePair->GetZ(2)*
	      sqrt(redmass/localEnergy);
            double rho = sqrt(2.*uconv)/hbarc*radius*sqrt(redmass*localEnergy);
	    complex expCP(1.0,0.0);
	    for(int ll=1;ll<=theChannel->GetL();ll++) 
	      expCP*=complex((double)ll/sqrt(pow(eta,2.0)+pow((double)ll,2.0)),
			     eta/sqrt(pow(eta,2.0)+pow((double)ll,2.0)));
	    struct CoulWaves 
	      coul=theCoulombFunction(lValue,radius,localEnergy);
	    complex expHSP(coul.G/sqrt(pow(coul.F,2.0)+pow(coul.G,2.0)),
			   -coul.F/sqrt(pow(coul.F,2.0)+pow(coul.G,2.0)));
	    this->AddLoElement(j,ch,loElement);
	    this->AddSqrtPenetrability(j,ch,sqrt(localPene));
	    this->AddExpCoulombPhase(j,ch,expCP);
	    this->AddExpHardSpherePhase(j,ch,expHSP);
//            if(ch==1){
//            std::ofstream coul_check;
//            coul_check.open("coul_check_pos.chk",std::ios::app);
//            coul_check << this->GetLabEnergy() << "  " << lValue << "  " << rho << "  " << eta << "  "
//                       << coul.F << "  " << coul.dF << "  " << coul.G << "  " << coul.dG << std::endl;
//            coul_check << this->GetEntranceKey() << "," << this->GetCMEnergy() << "," << lValue << "," << expHSP << "," << localShift << std::endl;
//            coul_check.close();
//            }
            
	  }
	} else if(thePair->GetPType()==10){
	  complex loElement = complex(0.0,0.0);
	  this->AddLoElement(j,ch,loElement);
	  double sqrtPene = (configure.paramMask & Config::USE_RMC_FORMALISM) ? 1. : pow(localEnergy/hbarc, (double) lValue+0.5);
	  this->AddSqrtPenetrability(j,ch,sqrtPene);
	  this->AddExpCoulombPhase(j,ch,1.0);
	  this->AddExpHardSpherePhase(j,ch,1.0);
	} else if(thePair->GetPType()==20){
	  complex loElement = complex(0.0,0.0);
	  this->AddLoElement(j,ch,loElement);
	  IntegratedFermiFunc fermiFunc(thePair->GetZ(1));
	  double endPointE = thePair->GetSepE()-inEnergy;
	  double sqrtPene = (1.+endPointE/0.510998903<=1.) ? 0. : sqrt(fermiFunc(1.+endPointE/0.510998903,exitPair->GetZ(1)+exitPair->GetZ(2),thePair->GetChRad()));
	  this->AddSqrtPenetrability(j,ch,sqrtPene);
	  this->AddExpCoulombPhase(j,ch,1.0);
	  this->AddExpHardSpherePhase(j,ch,1.0);
	}
      }
    }
  }
  for(int i=1;i<=this->NumSubPoints();i++) {
    this->GetSubPoint(i)->CalcEDependentValues(theCNuc,configure);
  }
  for(int i=1;i<=this->NumLocalMappedPoints();i++) {
    EPoint *mappedPoint=this->GetLocalMappedPoint(i);
    mappedPoint->geofactor_=geofactor_;
    mappedPoint->sfactorconv_=sfactorconv_;
    mappedPoint->lo_elements_=lo_elements_;
    mappedPoint->penetrabilities_=penetrabilities_;
    mappedPoint->coulombphase_=coulombphase_;
    mappedPoint->hardspherephase_=hardspherephase_;
    for(int ii=1;ii<=this->NumSubPoints();ii++) {
      EPoint *subMappedPoint=mappedPoint->GetSubPoint(ii);
      subMappedPoint->geofactor_=this->GetSubPoint(ii)->geofactor_;
      subMappedPoint->sfactorconv_=this->GetSubPoint(ii)->sfactorconv_;
      subMappedPoint->lo_elements_=this->GetSubPoint(ii)->lo_elements_;
      subMappedPoint->penetrabilities_=this->GetSubPoint(ii)->penetrabilities_;
      subMappedPoint->coulombphase_=this->GetSubPoint(ii)->coulombphase_;
      subMappedPoint->hardspherephase_=this->GetSubPoint(ii)->hardspherephase_;
    }
  }
}

/*!
 * Recalculates energy-dependent values using the current (possibly shifted) energy.
 * This is needed when energy shifts are applied after initialization.
 */
void EPoint::RecalcEDependentValues(CNuc *theCNuc, const Config& configure) {
  // Clear existing energy-dependent values first
  lo_elements_.clear();
  penetrabilities_.clear();
  coulombphase_.clear();
  hardspherephase_.clear();
  
  // Recalculate with current energy
  this->CalcEDependentValues(theCNuc, configure);
  
  // Also recalculate Coulomb amplitude with current energy
  this->CalcCoulombAmplitude(theCNuc);

}

/*!
 * Adds an \f$ L_o \f$ matrix element with reference to positions in the JGroup and subsequent
 * AChannel vectors.  
 */

void EPoint::AddLoElement(int jGroupNum, int channelNum, complex loElement) {
  vector_c d;
  while(jGroupNum>lo_elements_.size()) lo_elements_.push_back(d);
  lo_elements_[jGroupNum-1].push_back(loElement);
  assert(channelNum=lo_elements_[jGroupNum-1].size());
}

/*!
 * Adds a square root of penetrability  with reference to positions in the JGroup and subsequent
 * AChannel vectors.  
 */

void EPoint::AddSqrtPenetrability(int jGroupNum, int channelNum, double sqrtPene) {
  vector_r d;
  while(jGroupNum>penetrabilities_.size()) penetrabilities_.push_back(d);
  penetrabilities_[jGroupNum-1].push_back(sqrtPene);
  assert(channelNum=penetrabilities_[jGroupNum-1].size());
}

/*!
 * Adds an exponential of the Coulomb phase shift with reference to positions in the JGroup and subsequent
 * AChannel vectors.  
 */

void EPoint::AddExpCoulombPhase(int jGroupNum, int channelNum, complex expShift) {
  vector_c d;
  while(jGroupNum>coulombphase_.size()) coulombphase_.push_back(d);
  coulombphase_[jGroupNum-1].push_back(expShift);
  assert(channelNum=coulombphase_[jGroupNum-1].size());
}

/*!
 * Adds an exponential of the hard sphere phase shift with reference to positions in the JGroup and subsequent
 * AChannel vectors.  
 */

void EPoint::AddExpHardSpherePhase(int jGroupNum, int channelNum, complex expShift) {
  vector_c d;
  while(jGroupNum>hardspherephase_.size()) hardspherephase_.push_back(d);
  hardspherephase_[jGroupNum-1].push_back(expShift);
  assert(channelNum=hardspherephase_[jGroupNum-1].size());
}

/*!
 * Calculates the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

void EPoint::CalcCoulombAmplitude(CNuc *theCNuc) {
  if(this->GetEntranceKey()==this->GetExitKey()) {
    PPair *entrancePair=theCNuc->GetPair(theCNuc->GetPairNumFromKey(this->GetEntranceKey()));
    int z1=entrancePair->GetZ(1);
    int z2=entrancePair->GetZ(2);
    double redmass=entrancePair->GetRedMass();
    double energy=this->GetCMEnergy();
    double angle=this->GetCMAngle();
    double eta=sqrt(uconv/2.)*fstruc*z1*z2*
      sqrt(redmass/energy);
    double cal=(1.0/(2.0*sqrt(pi)))*eta*(1.0/pow(sin(angle*pi/360.0),2.));
    double cex=2.0*eta*log(sin(angle*pi/360.0));
    complex calpha(cal*cos(cex),-cal*sin(cex));
    // Identical-particle Mott amplitude: add f_C(pi - theta) with the
    // boson/fermion sign. For two identical 0+ bosons the sign is +1 and
    // the resulting |C|^2 reproduces the Mott Coulomb cross section.
    if(entrancePair->IsIdentical()) {
      double cal_p=(1.0/(2.0*sqrt(pi)))*eta*(1.0/pow(cos(angle*pi/360.0),2.));
      double cex_p=2.0*eta*log(cos(angle*pi/360.0));
      complex calpha_p(cal_p*cos(cex_p),-cal_p*sin(cex_p));
      calpha = calpha + double(entrancePair->GetIdenticalSign())*calpha_p;
    }
    this->SetCoulombAmplitude(calpha);
  } else this->SetCoulombAmplitude(complex(0.,0.));
  for(int i=1;i<=this->NumSubPoints();i++) {
    this->GetSubPoint(i)->CalcCoulombAmplitude(theCNuc);
  }
}

/*!
 * Sets the Coulomb amplitude \f$ C_\alpha \f$ for the data point.
 */

void EPoint::SetCoulombAmplitude(complex amplitude) {
  coulombamplitude_=amplitude;
}

/*!
 * Calculates the external capture amplitudes for the data point.
 * The amplitudes are calculated for all reaction pathways with corresponding
 * entrance and exit pairs.
 */

void EPoint::CalculateECAmplitudes(CNuc *theCNuc, const Config& configure) {
  int aa=theCNuc->GetPairNumFromKey(this->GetEntranceKey());
  if(theCNuc->GetPair(aa)->GetPType()==20) return;
  if(theCNuc->GetPair(aa)->IsEntrance()) {
    PPair *entrancePair=theCNuc->GetPair(aa);
    for(int j=1;j<=theCNuc->NumJGroups();j++) {
      for(int la=1;la<=theCNuc->GetJGroup(j)->NumLevels();la++) {
	if(theCNuc->GetJGroup(j)->GetLevel(la)->IsECLevel()) {
	  ALevel *ecLevel = theCNuc->GetJGroup(j)->GetLevel(la);
	  int ir=theCNuc->GetPairNumFromKey(this->GetExitKey());
	  if(ecLevel->GetECPairNum()==ir) {
	    double inEnergy=this->GetCMEnergy()+entrancePair->GetSepE()+entrancePair->GetExE();
	    for(int k=1;k<=entrancePair->GetDecay(ir)->NumKGroups();k++) {
	      KGroup *theKGroup=entrancePair->GetDecay(ir)->GetKGroup(k);
	      for(int ecm=1;ecm<=theKGroup->NumECMGroups();ecm++) {
		ECMGroup *theECMGroup=theKGroup->GetECMGroup(ecm);
		//entrance Phase Calculations;
		CoulFunc theCoulombFunction(entrancePair,!!(configure.paramMask&Config::USE_GSL_COULOMB_FUNC));
		struct CoulWaves 
		  coul=theCoulombFunction(theECMGroup->GetL(),entrancePair->GetChRad(),
					  this->GetCMEnergy());		
		double eta=sqrt(uconv/2.)*fstruc*entrancePair->GetZ(1)*entrancePair->GetZ(2)*
		  sqrt(entrancePair->GetRedMass()/this->GetCMEnergy());
		complex expCP(1.0,0.0);
		for(int ll=1;ll<=theECMGroup->GetL();ll++) 
		  expCP*=complex((double)ll/sqrt(pow((double)ll,2.0)+pow(eta,2.0)),
				 eta/sqrt(pow((double)ll,2.0)+pow(eta,2.0)));
		complex expHSP(coul.G/sqrt(pow(coul.F,2.0)+pow(coul.G,2.0)),
			       -coul.F/sqrt(pow(coul.F,2.0)+pow(coul.G,2.0)));
		
		double levelEnergy=ecLevel->GetE();
		double sqrtGammaPene=pow((inEnergy-levelEnergy)/hbarc,theECMGroup->GetMult()+0.5);
		
		//Initialize variables
		AChannel *theFinalChannel = theCNuc->GetJGroup(j)->GetChannel(theECMGroup->GetFinalChannel());
		PPair *theFinalPair=theCNuc->GetPair(theFinalChannel->GetPairNum());		
		int theInitialLValue;
		double theInitialSValue;
		if(theECMGroup->IsChannelCapture()) {
		  MGroup *theChanCapMGroup=entrancePair->GetDecay(theECMGroup->GetChanCapDecay())->
		    GetKGroup(theECMGroup->GetChanCapKGroup())->GetMGroup(theECMGroup->GetChanCapMGroup());
		  theInitialLValue=theCNuc->GetJGroup(theChanCapMGroup->GetJNum())->
		    GetChannel(theChanCapMGroup->GetChpNum())->GetL();
		  theInitialSValue=theCNuc->GetJGroup(theChanCapMGroup->GetJNum())->
		    GetChannel(theChanCapMGroup->GetChpNum())->GetS();
		} else {
		  theInitialLValue=theECMGroup->GetL();
		  theInitialSValue=theKGroup->GetS();
		}	    
		
		ECIntegral theECIntegral(theFinalPair,configure);
		complex integrals = theECIntegral(theInitialLValue, theFinalChannel->GetL(), 
						  theInitialSValue, theFinalChannel->GetS(),
						  theECMGroup->GetJ(), theCNuc->GetJGroup(j)->GetJ(),
						  theECMGroup->GetMult(), theECMGroup->GetRadType(),
						  inEnergy, levelEnergy,
						  theECMGroup->IsChannelCapture());
		
		//calculate the total radial integral
		complex ecAmplitude=expCP*expHSP*sqrtGammaPene*integrals;
		this->AddECAmplitude(k,ecm,ecAmplitude);
	      }
	    }
	  }
	}
      }
    }
  }
  for(int i=1;i<=this->NumSubPoints();i++) {
    this->GetSubPoint(i)->CalculateECAmplitudes(theCNuc,configure);
  }
}

/*!
 * Adds an external capture amplitude with reference to a specified reaction pathway.
 */

void EPoint::AddECAmplitude(int kGroupNum, int ecMGroupNum, complex ecAmplitude, double energy) {
  vector_c d_amp;
  vector_r d_energy;
  while(kGroupNum>ec_amplitudes_.size()) {
    ec_amplitudes_.push_back(d_amp);
    ec_energies_.push_back(d_energy);
  }
  ec_amplitudes_[kGroupNum-1].push_back(ecAmplitude);
  ec_energies_[kGroupNum-1].push_back(energy);
  assert(ecMGroupNum==ec_amplitudes_[kGroupNum-1].size());
  assert(ec_amplitudes_[kGroupNum-1].size() == ec_energies_[kGroupNum-1].size());
  
  // Add to shared cache for interpolation across different EPoints
  if (g_ecAmplitudeCache) {
    ECAmplitudeCache::AmplitudeKey key;
    key.kGroupNum = kGroupNum;
    key.ecMGroupNum = ecMGroupNum;
    key.entranceKey = this->GetEntranceKey();
    key.exitKey = this->GetExitKey();
    key.segmentKey = segment_key_;
    g_ecAmplitudeCache->AddAmplitude(key, energy, ecAmplitude);

    // Now get the interpolated amplitude from the cache to replace the stored value if ecAmplitudes is 0 or NaN
    if(std::isnan(std::real(ecAmplitude)) || std::isnan(std::imag(ecAmplitude))) {
        std::cout << "WARNING: Replacing EC Amplitude with interpolated value from cache for kGroup " 
                  << kGroupNum << ", ecMGroup " << ecMGroupNum 
                  << " at energy " << energy << " MeV." << std::endl;
      complex interpAmplitude = g_ecAmplitudeCache->GetInterpolatedAmplitude(key, energy);
      std::cout << "         Original Amplitude: " << ecAmplitude 
                << ", Interpolated Amplitude: " << interpAmplitude << std::endl;
      ec_amplitudes_[kGroupNum-1].back() = interpAmplitude;
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
}

/*!
 * Calculates the cross section for a data point based on the fit parameters in the compound nucleus.
 */

void EPoint::Calculate(CNuc* theCNuc,const Config &configure, EPoint *parent, int subPointNum) {

  if(!this->IsTargetEffect()||
     (!this->GetParentData()->GetTargetEffect(this->GetTargetEffectNum())->IsConvolution()&&
      !this->GetParentData()->GetTargetEffect(this->GetTargetEffectNum())->IsTargetIntegration()&&
      !this->GetParentData()->GetTargetEffect(this->GetTargetEffectNum())->IsConvCoefficients())) {
    
    GenMatrixFunc *theMatrixFunc;
    // The A-matrix function object is reused across the points calculated by
    // this thread.  Its per-JGroup buffers and the level-matrix factorization
    // workspace are sized once and then only refilled, which removes the
    // allocation traffic that dominated the old per-point construction.
    // Recursion is not a concern: the branch that recurses into sub-points does
    // not use theMatrixFunc.
    static thread_local AMatrixFunc reusableAMatrixFunc(theCNuc,configure);
    std::unique_ptr<RMatrixFunc> ownedRMatrixFunc;
    if(configure.paramMask & Config::USE_AMATRIX) {
      reusableAMatrixFunc.Reset(theCNuc,configure);
      theMatrixFunc=&reusableAMatrixFunc;
    } else {
      ownedRMatrixFunc.reset(new RMatrixFunc(theCNuc,configure));
      theMatrixFunc=ownedRMatrixFunc.get();
    }
    theMatrixFunc->ClearMatrices();
    theMatrixFunc->FillMatrices(this);
    theMatrixFunc->InvertMatrices();
    theMatrixFunc->CalculateTMatrix(this);
    theMatrixFunc->CalculateCrossSection(this);
    if(subPointNum&&parent) {
      for(int i=1;i<=parent->NumLocalMappedPoints();i++) {
	EPoint *mappedSubPoint = parent->GetLocalMappedPoint(i)->
	  GetSubPoint(subPointNum);
	theMatrixFunc->CalculateCrossSection(mappedSubPoint);
      }
    } else {
      for(int i=1;i<=this->NumLocalMappedPoints();i++) {
	EPoint *mappedPoint = this->GetLocalMappedPoint(i);
	theMatrixFunc->CalculateCrossSection(mappedPoint);
      }
    }
  }
  else {
    for(int i = 1; i<=this->NumSubPoints();i++) {
      EPoint *subPoint=this->GetSubPoint(i);
      // Recalculate energy-dependent values for sub-points too if using external capture
      if(configure.paramMask & Config::USE_EXTERNAL_CAPTURE) {
        subPoint->RecalcEDependentValues(theCNuc, configure);
      }
      if(this->NumLocalMappedPoints()>0)
	    subPoint->Calculate(theCNuc,configure,this,i);
      else subPoint->Calculate(theCNuc,configure);
    }
    this->IntegrateTargetEffectForObservable(configure);
    for(int i=1;i<=this->NumLocalMappedPoints();i++)
      this->GetLocalMappedPoint(i)->IntegrateTargetEffectForObservable(configure);
  }
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
void EPoint::IntegrateTargetEffectForObservable(const Config& configure) {
  if(!this->IsAnalyzingPower()) {
    this->IntegrateTargetEffect(configure);
    return;
  }

  const int n = this->NumSubPoints();
  std::vector<double> sigma(n);
  for(int i=1;i<=n;i++) sigma[i-1] = this->GetSubPoint(i)->GetFitCrossSection();

  this->IntegrateTargetEffect(configure);
  const double denominator = this->GetFitCrossSection();

  for(int i=1;i<=n;i++)
    this->GetSubPoint(i)->SetFitCrossSection(sigma[i-1] *
                                             this->GetSubPoint(i)->GetAnalyzingPower());
  this->IntegrateTargetEffect(configure);
  const double numerator = this->GetFitCrossSection();

  for(int i=1;i<=n;i++) this->GetSubPoint(i)->SetFitCrossSection(sigma[i-1]);

  this->SetFitCrossSection(std::fabs(denominator) > 0.0 ? numerator/denominator : 0.0);
}

/*!
 * If a point is mapped, sets the internal attribute indicating which point it is mapped to.
 */

void EPoint::SetMap(int segmentNum, int pointNum) {
  is_mapped_=true;
  energy_map_.segment=segmentNum;
  energy_map_.point=pointNum;
}

/*!
 * Clears the mapping status of the point, forcing it to be recalculated.
 * Used when energies change (e.g., after energy shifts during fitting) and
 * mapped points are no longer at identical energies.
 */

void EPoint::ClearMapping() {
  is_mapped_=false;
  energy_map_.segment=0;
  energy_map_.point=0;
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
  targetEffectNum_=targetEffectNum;
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
  subPoint.is_sub_point_ = true;
  integrationPoints_.push_back(subPoint);
}

/*!
 * This function is called to integrate the vector of sub-points to determine
 * the yield considering a given target effect.
 *
 * Integration method:
 * - Gaussian quadrature for non-uniform grids
 * - Uses 2-point Gauss-Legendre quadrature on each interval
 * - Third-order accurate, exact for cubic polynomials
 * - Works correctly with variable energy spacing from adaptive grid
 * - More robust than Simpson's for non-uniform grids
 * - Original nested Simpson's still used for convolution+targetIntegration
 */

void EPoint::IntegrateTargetEffect(const Config& configure) {
  double yield=0.0;
  TargetEffect *targetEffect=this->GetParentData()->GetTargetEffect(this->GetTargetEffectNum());

  // All integration now uses Gaussian quadrature which works for non-uniform (adaptive) grids
  if(targetEffect->IsConvolution()&&targetEffect->IsTargetIntegration()) {
    // Gaussian quadrature for convolution + target integration
    // Uses 2-point Gauss-Legendre for both outer (depth) and inner (convolution) integrals

    int numPoints = this->NumSubPoints();
    if(numPoints < 2) {
      this->SetFitCrossSection(0.0);
      return;
    }

    // Get beam sigma and straggling parameters
    double beamSigma = targetEffect->GetSigma();
    bool useStraggling = targetEffect->IsStraggling();
    double stragglingCoeff = targetEffect->GetStragglingCoefficient();
    double convRange = targetEffect->convolutionRange;

    // Gauss-Legendre 2-point weights and nodes
    const double w1 = 1.0;
    const double w2 = 1.0;
    const double x1 = -1.0/sqrt(3.0);
    const double x2 =  1.0/sqrt(3.0);

    // Target depth integration range: from surface to back of target
    double surfaceEnergy = this->GetCMEnergy();
    double backEnergy = surfaceEnergy - this->GetTargetThickness();

    // Available energy range from subpoints
    double minEnergy = this->GetSubPoint(numPoints)->GetCMEnergy();
    double maxEnergy = this->GetSubPoint(1)->GetCMEnergy();

    double outerIntegral = 0.0;

    // Outer integral: over target depth using Gaussian quadrature
    // Loop over all subpoint intervals and integrate those within target range
    for(int i = 1; i < numPoints; i++) {
      double Ea_outer = this->GetSubPoint(i)->GetCMEnergy();
      double Eb_outer = this->GetSubPoint(i+1)->GetCMEnergy();

      // Skip intervals completely outside target range [backEnergy, surfaceEnergy]
      // Note: Ea_outer > Eb_outer (energy decreases with index)
      if(Eb_outer > surfaceEnergy || Ea_outer < backEnergy) continue;

      // Clip interval to target range
      double Ea_clipped = std::min(Ea_outer, surfaceEnergy);
      double Eb_clipped = std::max(Eb_outer, backEnergy);

      double mid_outer = (Ea_clipped + Eb_clipped) / 2.0;
      double half_range_outer = (Ea_clipped - Eb_clipped) / 2.0;

      if(half_range_outer < 1.0e-10) continue;  // Skip tiny intervals

      // Evaluate at two Gauss points
      for(int gp = 0; gp < 2; gp++) {
        double xi = (gp == 0) ? x1 : x2;
        double wi = (gp == 0) ? w1 : w2;

        double E_depth = mid_outer + half_range_outer * xi;

        // Calculate energy loss and effective sigma at this depth
        // Use distance from surface (surfaceEnergy) for straggling calculation
        double deltaE = surfaceEnergy - E_depth;
        double deltaE_keV = deltaE * 1000.0;

        double effectiveSigma = beamSigma;
        if(useStraggling && deltaE_keV > 0) {
          double stragglingSigma_keV = stragglingCoeff * std::sqrt(deltaE_keV);
          double stragglingSigma = stragglingSigma_keV / 1000.0;
          effectiveSigma = std::sqrt(beamSigma * beamSigma + stragglingSigma * stragglingSigma);
        }

        // Inner integral: convolution at this depth
        // Integration range: E_depth ± convRange * effectiveSigma
        double convLow = E_depth - convRange * effectiveSigma;
        double convHigh = E_depth + convRange * effectiveSigma;

        // Clamp to available energy range
        if(convLow < minEnergy) convLow = minEnergy;
        if(convHigh > maxEnergy) convHigh = maxEnergy;

        double innerIntegral = 0.0;

        // Inner integral: convolution over energy distribution
        for(int j = 1; j < numPoints; j++) {
          double Ea_inner = this->GetSubPoint(j)->GetCMEnergy();
          double Eb_inner = this->GetSubPoint(j+1)->GetCMEnergy();

          // Skip intervals completely outside convolution range
          if(Eb_inner > convHigh || Ea_inner < convLow) continue;

          // Clip interval to convolution range
          double Ea_inner_clipped = std::min(Ea_inner, convHigh);
          double Eb_inner_clipped = std::max(Eb_inner, convLow);

          double mid_inner = (Ea_inner_clipped + Eb_inner_clipped) / 2.0;
          double half_range_inner = (Ea_inner_clipped - Eb_inner_clipped) / 2.0;

          if(half_range_inner < 1.0e-10) continue;  // Skip tiny intervals

          // Inner Gauss points
          for(int gp_inner = 0; gp_inner < 2; gp_inner++) {
            double xi_inner = (gp_inner == 0) ? x1 : x2;
            double wi_inner = (gp_inner == 0) ? w1 : w2;

            double E_conv = mid_inner + half_range_inner * xi_inner;

            // Interpolate cross-section and stopping power at E_conv
            // Use the original (unclipped) interval for interpolation
            double alpha = (Ea_inner - E_conv) / (Ea_inner - Eb_inner);
            double sigma = (1.0 - alpha) * this->GetSubPoint(j)->GetFitCrossSection() +
                          alpha * this->GetSubPoint(j+1)->GetFitCrossSection();
            double stopping = (1.0 - alpha) * this->GetSubPoint(j)->GetStoppingPower() +
                             alpha * this->GetSubPoint(j+1)->GetStoppingPower();

            // Convolution factor (Gaussian centered at E_depth)
            double convFactor = std::pow(2.0*pi, -0.5) / effectiveSigma *
                               std::exp(-std::pow(E_conv - E_depth, 2.0) / (2.0 * effectiveSigma * effectiveSigma));

            double integrand = sigma / stopping / 1e24 * convFactor;
            innerIntegral += std::abs(half_range_inner) * wi_inner * integrand;
          }
        }

        outerIntegral += std::abs(half_range_outer) * wi * innerIntegral;
      }
    }

    yield = outerIntegral / (targetEffect->GetDensity() * 1.0e-24);
  } 
  else if(targetEffect->IsConvolution()) {
    // Gaussian quadrature for non-uniform grids
    // 2-point Gauss-Legendre: exact for cubic polynomials
    double integral=0.0;
    double centroid=this->GetCMEnergy();
    int numPoints = this->NumSubPoints();

    // Gauss-Legendre 2-point weights and nodes (on [-1,1])
    const double w1 = 1.0;
    const double w2 = 1.0;
    const double x1 = -1.0/sqrt(3.0);  // ≈ -0.577350
    const double x2 =  1.0/sqrt(3.0);  // ≈  0.577350

    // Integrate over each interval [E_i, E_{i+1}]
    for(int i=0; i<numPoints-1; i++) {
      double Ea = this->GetSubPoint(i+1)->GetCMEnergy();
      double Eb = this->GetSubPoint(i+2)->GetCMEnergy();

      // Transform from [Ea,Eb] to standard interval [-1,1]
      double mid = (Ea + Eb) / 2.0;
      double half_range = (Ea - Eb) / 2.0;

      // Gauss points in physical domain
      double E_gauss1 = mid + half_range * x1;
      double E_gauss2 = mid + half_range * x2;

      // Evaluate integrand at Gauss points (linear interpolation from grid points)
      // f1 at E_gauss1
      double alpha1 = (Eb - E_gauss1) / (Eb - Ea);
      double sigma1 = alpha1 * this->GetSubPoint(i+1)->GetFitCrossSection() +
                      (1.0 - alpha1) * this->GetSubPoint(i+2)->GetFitCrossSection();
      double f1 = sigma1 * targetEffect->GetConvolutionFactor(E_gauss1, centroid);

      // f2 at E_gauss2
      double alpha2 = (Eb - E_gauss2) / (Eb - Ea);
      double sigma2 = alpha2 * this->GetSubPoint(i+1)->GetFitCrossSection() +
                      (1.0 - alpha2) * this->GetSubPoint(i+2)->GetFitCrossSection();
      double f2 = sigma2 * targetEffect->GetConvolutionFactor(E_gauss2, centroid);

      // Gaussian quadrature formula
      integral += std::abs(half_range) * (w1 * f1 + w2 * f2);
    }

    yield=integral;
  } 
  else if(targetEffect->IsTargetIntegration()) {
    // Target integration with optional straggling
    int numPoints = this->NumSubPoints();
    if(numPoints < 2) {
      this->SetFitCrossSection(0.0);
      return;
    }

    bool useStraggling = targetEffect->IsStraggling();

    if(useStraggling) {
      // Target integration with straggling: convolve with straggling-induced energy spread
      // Similar to convolution+targetIntegration but sigma is only from straggling (no beam sigma)

      double stragglingCoeff = targetEffect->GetStragglingCoefficient();
      double convRange = targetEffect->convolutionRange;

      // Target depth integration range: from surface to back of target
      double surfaceEnergy = this->GetCMEnergy();
      double backEnergy = surfaceEnergy - this->GetTargetThickness();

      // Energy at the highest subpoint (used for straggling calculation)
      double initialEnergy = this->GetSubPoint(1)->GetCMEnergy();

      // Available energy range from subpoints
      double minEnergy = this->GetSubPoint(numPoints)->GetCMEnergy();
      double maxEnergy = this->GetSubPoint(1)->GetCMEnergy();

      // Gauss-Legendre 2-point weights and nodes
      const double w1 = 1.0;
      const double w2 = 1.0;
      const double x1 = -1.0/sqrt(3.0);
      const double x2 =  1.0/sqrt(3.0);

      double outerIntegral = 0.0;

      // Outer integral: over target depth using Gaussian quadrature
      // Loop over all subpoint intervals and integrate those within target range
      for(int i = 1; i < numPoints; i++) {
        double Ea_outer = this->GetSubPoint(i)->GetCMEnergy();
        double Eb_outer = this->GetSubPoint(i+1)->GetCMEnergy();

        // Skip intervals completely outside target range [backEnergy, surfaceEnergy]
        // Note: Ea_outer > Eb_outer (energy decreases with index)
        if(Eb_outer > surfaceEnergy || Ea_outer < backEnergy) continue;

        // Clip interval to target range
        double Ea_clipped = std::min(Ea_outer, surfaceEnergy);
        double Eb_clipped = std::max(Eb_outer, backEnergy);

        double mid_outer = (Ea_clipped + Eb_clipped) / 2.0;
        double half_range_outer = (Ea_clipped - Eb_clipped) / 2.0;

        if(half_range_outer < 1.0e-10) continue;

        // Evaluate at two Gauss points
        for(int gp = 0; gp < 2; gp++) {
          double xi = (gp == 0) ? x1 : x2;
          double wi = (gp == 0) ? w1 : w2;

          double E_depth = mid_outer + half_range_outer * xi;

          // Calculate straggling sigma at this depth
          // Use distance from surface (surfaceEnergy), not from first subpoint
          double deltaE = surfaceEnergy - E_depth;
          double deltaE_keV = deltaE * 1000.0;

          double stragglingSigma = 0.0;
          if(deltaE_keV > 0) {
            double stragglingSigma_keV = stragglingCoeff * std::sqrt(deltaE_keV);
            stragglingSigma = stragglingSigma_keV / 1000.0;
          }

          double innerIntegral = 0.0;

          if(stragglingSigma > 1.0e-10) {
            // Convolve with straggling distribution
            double convLow = E_depth - convRange * stragglingSigma;
            double convHigh = E_depth + convRange * stragglingSigma;
            if(convLow < minEnergy) convLow = minEnergy;
            if(convHigh > maxEnergy) convHigh = maxEnergy;

            for(int j = 1; j < numPoints; j++) {
              double Ea_inner = this->GetSubPoint(j)->GetCMEnergy();
              double Eb_inner = this->GetSubPoint(j+1)->GetCMEnergy();

              if(Eb_inner > convHigh || Ea_inner < convLow) continue;

              double Ea_inner_clipped = std::min(Ea_inner, convHigh);
              double Eb_inner_clipped = std::max(Eb_inner, convLow);

              double mid_inner = (Ea_inner_clipped + Eb_inner_clipped) / 2.0;
              double half_range_inner = (Ea_inner_clipped - Eb_inner_clipped) / 2.0;

              if(half_range_inner < 1.0e-10) continue;

              for(int gp_inner = 0; gp_inner < 2; gp_inner++) {
                double xi_inner = (gp_inner == 0) ? x1 : x2;
                double wi_inner = (gp_inner == 0) ? w1 : w2;

                double E_conv = mid_inner + half_range_inner * xi_inner;

                // Interpolate cross-section and stopping power
                double alpha = (Ea_inner - E_conv) / (Ea_inner - Eb_inner);
                double sigma = (1.0 - alpha) * this->GetSubPoint(j)->GetFitCrossSection() +
                              alpha * this->GetSubPoint(j+1)->GetFitCrossSection();
                double stopping = (1.0 - alpha) * this->GetSubPoint(j)->GetStoppingPower() +
                                 alpha * this->GetSubPoint(j+1)->GetStoppingPower();

                // Straggling convolution factor
                double convFactor = std::pow(2.0*pi, -0.5) / stragglingSigma *
                                   std::exp(-std::pow(E_conv - E_depth, 2.0) / (2.0 * stragglingSigma * stragglingSigma));

                double integrand = sigma / stopping / 1e24 * convFactor;
                innerIntegral += std::abs(half_range_inner) * wi_inner * integrand;
              }
            }
          } else {
            // No straggling at surface - just use the cross-section directly
            // Interpolate at E_depth
            for(int j = 1; j < numPoints; j++) {
              double Ea = this->GetSubPoint(j)->GetCMEnergy();
              double Eb = this->GetSubPoint(j+1)->GetCMEnergy();
              if(E_depth <= Ea && E_depth >= Eb) {
                double alpha = (Ea - E_depth) / (Ea - Eb);
                double sigma = (1.0 - alpha) * this->GetSubPoint(j)->GetFitCrossSection() +
                              alpha * this->GetSubPoint(j+1)->GetFitCrossSection();
                double stopping = (1.0 - alpha) * this->GetSubPoint(j)->GetStoppingPower() +
                                 alpha * this->GetSubPoint(j+1)->GetStoppingPower();
                innerIntegral = sigma / stopping / 1e24;
                break;
              }
            }
          }

          outerIntegral += std::abs(half_range_outer) * wi * innerIntegral;
        }
      }

      yield = outerIntegral / (targetEffect->GetDensity() * 1.0e-24);
    } else {
      // Standard target integration without straggling
      // Gaussian quadrature for non-uniform grids
      double integral=0.0;

      // Gauss-Legendre 2-point weights and nodes (on [-1,1])
      const double w1 = 1.0;
      const double w2 = 1.0;
      const double x1 = -1.0/sqrt(3.0);
      const double x2 =  1.0/sqrt(3.0);

      // Integrate over each interval [E_i, E_{i+1}]
      for(int i=0; i<numPoints-1; i++) {
        double Ea = this->GetSubPoint(i+1)->GetCMEnergy();
        double Eb = this->GetSubPoint(i+2)->GetCMEnergy();

        // Transform from [Ea,Eb] to standard interval [-1,1]
        double mid = (Ea + Eb) / 2.0;
        double half_range = (Ea - Eb) / 2.0;

        // Gauss points in physical domain
        double E_gauss1 = mid + half_range * x1;
        double E_gauss2 = mid + half_range * x2;

        // Evaluate integrand at Gauss points (linear interpolation from grid points)
        // f1 at E_gauss1
        double alpha1 = (Eb - E_gauss1) / (Eb - Ea);
        double sigma1 = alpha1 * this->GetSubPoint(i+1)->GetFitCrossSection() +
                        (1.0 - alpha1) * this->GetSubPoint(i+2)->GetFitCrossSection();
        double stopping1 = alpha1 * this->GetSubPoint(i+1)->GetStoppingPower() +
                           (1.0 - alpha1) * this->GetSubPoint(i+2)->GetStoppingPower();
        double f1 = sigma1 / stopping1 / 1e24;

        // f2 at E_gauss2
        double alpha2 = (Eb - E_gauss2) / (Eb - Ea);
        double sigma2 = alpha2 * this->GetSubPoint(i+1)->GetFitCrossSection() +
                        (1.0 - alpha2) * this->GetSubPoint(i+2)->GetFitCrossSection();
        double stopping2 = alpha2 * this->GetSubPoint(i+1)->GetStoppingPower() +
                           (1.0 - alpha2) * this->GetSubPoint(i+2)->GetStoppingPower();
        double f2 = sigma2 / stopping2 / 1e24;

        // Gaussian quadrature formula
        integral += std::abs(half_range) * (w1 * f1 + w2 * f2);
      }

      yield=integral/(targetEffect->GetDensity()*1.E-24);
    }
  }
  else if(targetEffect->IsConvCoefficients()) {
    // Gaussian quadrature for both integrals
    // 2-point Gauss-Legendre: exact for cubic polynomials
    double integral=0.0;
    double integralC=0.0;
    double centroid=this->GetCMEnergy();
    int numPoints = this->NumSubPoints();

    // Gauss-Legendre 2-point weights and nodes (on [-1,1])
    const double w1 = 1.0;
    const double w2 = 1.0;
    const double x1 = -1.0/sqrt(3.0);  // ≈ -0.577350
    const double x2 =  1.0/sqrt(3.0);  // ≈  0.577350

    // Integrate over each interval [E_i, E_{i+1}]
    for(int i=0; i<numPoints-1; i++) {
      double Ea = this->GetSubPoint(i+1)->GetCMEnergy();
      double Eb = this->GetSubPoint(i+2)->GetCMEnergy();

      // Transform from [Ea,Eb] to standard interval [-1,1]
      double mid = (Ea + Eb) / 2.0;
      double half_range = (Ea - Eb) / 2.0;

      // Gauss points in physical domain
      double E_gauss1 = mid + half_range * x1;
      double E_gauss2 = mid + half_range * x2;

      // Evaluate integrand at Gauss points (linear interpolation from grid points)
      // Point 1
      double alpha1 = (Eb - E_gauss1) / (Eb - Ea);
      double sigma1 = alpha1 * this->GetSubPoint(i+1)->GetFitCrossSection() +
                      (1.0 - alpha1) * this->GetSubPoint(i+2)->GetFitCrossSection();
      double conv1 = targetEffect->CalculateConvolutionFactor(E_gauss1, centroid, configure);
      double f1 = sigma1 * conv1;

      // Point 2
      double alpha2 = (Eb - E_gauss2) / (Eb - Ea);
      double sigma2 = alpha2 * this->GetSubPoint(i+1)->GetFitCrossSection() +
                      (1.0 - alpha2) * this->GetSubPoint(i+2)->GetFitCrossSection();
      double conv2 = targetEffect->CalculateConvolutionFactor(E_gauss2, centroid, configure);
      double f2 = sigma2 * conv2;

      // Gaussian quadrature formula for both integrals
      integral += std::abs(half_range) * (w1 * f1 + w2 * f2);
      integralC += std::abs(half_range) * (w1 * conv1 + w2 * conv2);
    }

    yield=integral/integralC;
  }
  this->SetFitCrossSection(yield);
}

/*!
 * This function sets an internal pointer to the parent EData object.
 */

void EPoint::SetParentData(EData* parentData) {
  parentData_=parentData;
}

/*!
 * This function sets the stopping cross section (effective) for the current
 * EPoint object.  Used only for sub-point involoved in target effect integration.
 */

void EPoint::SetStoppingPower(double stoppingPower) {
 stoppingPower_=stoppingPower;
}

/*!
 * This functions sets the energy loss of the beam in the target for the current EPoint object.
 */

void EPoint::SetTargetThickness(double targetThickness) {
  targetThickness_=targetThickness;
}

/*!
 * Sets the angular distribution coefficients.
 */

void EPoint::SetAngularDists(vector_r dists) {
  angularDists_.clear();
  angularDists_=dists;
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

EPoint* EPoint::GetLocalMappedPoint(int mappedPointNum) const {
  return local_mapped_points_[mappedPointNum-1];
}

/*!
 * Returns a pointer to the specified sub-point in the current EPoint object.
 */

EPoint* EPoint::GetSubPoint(int subPoint) {
  EPoint *tempPoint;
  if(subPoint<=integrationPoints_.size()) tempPoint=&integrationPoints_[subPoint-1];
  else tempPoint= NULL;
  return tempPoint;
}

/*!
 * Returns a reference to the vector of EPoint objects containing the subpoints used in target
 * effect integration.
 */

std::vector<EPoint>& EPoint::GetSubPoints() {
  return integrationPoints_;
}

/*!
 * Returns a reference to the vector of pointers to mapped EPoint objects.
 */

std::vector<EPoint*>& EPoint::GetMappedPoints() {
  return local_mapped_points_;
}
