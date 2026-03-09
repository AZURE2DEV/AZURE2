#include "TargetEffect.h"
#include <sstream>
#include <iostream>

/*!
 * Constructor reads directly from an std::ifstream pointing to the target
 * effect input file.  If a valid target effect is read, a TargetEffect object is
 * created.
 */

TargetEffect::TargetEffect(std::istream &stream,const Config& configure) {
  int isActive;
  std::string segmentList;
  int numIntegrationPoints;
  int isConvolution;
  double sigma;
  int isTargetIntegration;
  double density;
  std::string stoppingPowerEq;
  int numParameters;
  vector_r parameters; 
  int isQCoefficients;
  int numQCoefficients;
  vector_r qCoefficients;
  stream >> isActive >> segmentList >> numIntegrationPoints >> isConvolution
	 >> sigma >> isTargetIntegration >> density >> stoppingPowerEq 
	 >> numParameters;
  if(!stream.eof()) {
    for(int i=0;i<numParameters;i++) {
      double tempParameter;
      stream >> tempParameter;
      parameters.push_back(tempParameter);
    }
    stream >> isQCoefficients >> numQCoefficients;
    for(int i=0;i<numQCoefficients;i++) {
      double tempQCoefficient;
      stream >> tempQCoefficient;
      qCoefficients.push_back(tempQCoefficient);
    }
    isQCoefficients_ = (isQCoefficients==1) ? true : false;
    qCoefficients_=qCoefficients;
    size_t found=0;
    while(found!=std::string::npos) {
      found=segmentList.find('\"');
      if(found!=std::string::npos) segmentList.erase(found,1);
    }
    found=0;
    while(found!=std::string::npos) {
      found=stoppingPowerEq.find('\"');
      if(found!=std::string::npos) stoppingPowerEq.erase(found,1);
    }
    if(isActive==1) isActive_=true;
    else isActive_=false;
    segmentsList_=segmentList;
    numIntegrationPoints_=numIntegrationPoints;
    if(isConvolution==1) isConvolution_=true;
    else isConvolution_=false;
    sigma_=sigma;
    if(isTargetIntegration==1) isTargetIntegration_=true;
    else isTargetIntegration_=false;
    density_=density;
    if(isTargetIntegration_) {
      stoppingPowerEq_.Initialize(stoppingPowerEq,numParameters,configure);
      for(int i=0;i<numParameters;i++) {
	stoppingPowerEq_.SetParameter(i,parameters[i],configure);
      }
    }
  }
}

/*!
 * Returns true if the target effect was marked as active in the target effects 
 * input file, otherwise returns false.
 */

bool TargetEffect::IsActive() const {
  return isActive_;
}

/*!
 * Returns true if the target effect contains Gaussian beam convolution,
 * otherwise returns false.
 */

bool TargetEffect::IsConvolution() const {
  return isConvolution_;
}

/*!
 * Returns true if the target effect contains target integration, otherwise
 * returns false.
 */

bool TargetEffect::IsTargetIntegration() const {
  return isTargetIntegration_;
}

/*!
 * Returns true if the target effect contains attenuation coefficients, otherwise
 * returns false.
 */

bool TargetEffect::IsQCoefficients() const {
  return isQCoefficients_;
}

/*!
 * Returns the number of sub-points specified for the target effect in 
 * the input file.
 */

int TargetEffect::NumSubPoints() const {
  return numIntegrationPoints_;
}

/*!
 * Returns the number of attenuation coefficients for the target effect in 
 * the input file.
 */

int TargetEffect::NumQCoefficients() const {
  return qCoefficients_.size();
}

/*!
 * Returns the sigma of the Guassian for beam convolution.
 */

double TargetEffect::GetSigma() const {
  return sigma_;
}

/*!
 * Returns the density of the target in atoms/cm^2.  Only needed for 
 * target integration, not Gaussian beam convolution.
 */

double TargetEffect::GetDensity() const {
  return density_;
}

/*!
 * Calculates the Target thickness from the stopping cross section and
 * the target density as a function of energy.
 */

double TargetEffect::TargetThickness(double energy, const Config& configure)  {
  return this->GetStoppingPowerEq()->Evaluate(configure,energy)*this->GetDensity();
}

/*!
 * Returns the attenuation coefficients for the given order specified in by the target effect.
 */

double TargetEffect::GetQCoefficient(int order)  const {
  return (qCoefficients_.size()>order) ? qCoefficients_[order] : 1.;
}

/*!
 * Sets the convolution sigma to a new value.
 */

void TargetEffect::SetSigma(double sigma) {
  sigma_=sigma;
}

/*!
 * Sets the number of sub-points for the TargetEffect object.
 */

void TargetEffect::SetNumSubPoints(int numPoints) {
  numIntegrationPoints_=numPoints;
}

/*!
 * Parses and returns a vector of integers corresponding to the 
 * segment list specified as a string.  The segments list contains the
 * segments for which the target effect is applicable.
 */

std::vector<int> TargetEffect::GetSegmentsList() const {
  std::vector<int> tempList;
  int i=0;
  int lastSegNum=0;
  bool inclusive=false;
  while(i<segmentsList_.length()) {
    if(segmentsList_[i]>='0'&&segmentsList_[i]<='9') {
      std::string tempString;
      while(segmentsList_[i]!=','&&segmentsList_[i]!='-'&&
	    i<segmentsList_.length()) {
	tempString+=segmentsList_[i];
	i++;
      }
      std::istringstream stm;
      stm.str(tempString);
      int tempSegNum;stm>>tempSegNum;
      if(inclusive==true) for(int j=lastSegNum+1;j<=tempSegNum;j++) 
			    tempList.push_back(j);
      else tempList.push_back(tempSegNum);
      lastSegNum=tempSegNum;
    }
    if(segmentsList_[i]=='-') inclusive=true;
    else inclusive =false;
    i++;
  }      
  return tempList;
}

/*!
 * Returns the Equation object corresponding to the parametrized stopping 
 * cross section.
 */

Equation *TargetEffect::GetStoppingPowerEq() {
  Equation *tempEquation;
  tempEquation=&stoppingPowerEq_;
  return tempEquation;
}

/*!
 * Returns the multiplicative convolution factor for evaluation of the integrand
 * of a target effect.  
 */

double TargetEffect::GetConvolutionFactor(double energy, double centroid) const {
  double sigma=this->GetSigma();
  return pow(2.*pi,-0.5)/sigma*exp(-pow(energy-centroid,2.0)/2.0/pow(sigma,2.0));
}
