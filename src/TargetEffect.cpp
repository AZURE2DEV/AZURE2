#include "TargetEffect.h"
#include <sstream>
#include <iostream>
// #include "Straggling.h"  // ERYA straggling integration - commented for now

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
  int isConvCoefficients;
  int numConvCoefficients;
  vector_r convCoefficients;
  std::string convolutionEq;

  stream >> isActive >> segmentList >> numIntegrationPoints >> isConvolution >> sigma >> isTargetIntegration >> density >> stoppingPowerEq >> numParameters;
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

    stream >> isConvCoefficients >> convolutionEq >> numConvCoefficients;
    for(int i=0;i<numConvCoefficients;i++) {
      double tempConvCoefficient;
      stream >> tempConvCoefficient;
      convCoefficients.push_back(tempConvCoefficient);
    }
    isConvCoefficients_ = (isConvCoefficients==1) ? true : false;
    convCoefficients_=convCoefficients;

    // ERYA straggling integration initialization (commented for now)
    /*
    // Read straggling flag and target element for ERYA integration
    int isStraggling;
    int targetElement;
    stream >> isStraggling >> targetElement;
    
    isStraggling_ = (isStraggling == 1);
    targetElement_ = targetElement;
    
    // Initialize straggling calculator if needed
    if (isStraggling_ && isTargetIntegration_) {
        stragglingCalculator_ = std::make_unique<Straggling>(targetElement_, 1000.0); // 1 MeV default
        stragglingCalculator_->setOxideLayer(true, 0.15); // Default oxide layer
    }
    */

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
    found=0;
    while(found!=std::string::npos) {
      found=convolutionEq.find('\"');
      if(found!=std::string::npos) convolutionEq.erase(found,1);
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
    if(isConvCoefficients_) {
      convolutionEq_.Initialize(convolutionEq,numConvCoefficients,configure);
      for(int i=0;i<numConvCoefficients;i++) {
        convolutionEq_.SetParameter(i,convCoefficients[i],configure);
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
 * Returns true if the target effect contains convolution coefficients, otherwise
 * returns false.
 */

bool TargetEffect::IsConvCoefficients() const {
  return isConvCoefficients_;
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
 * Returns the number of convolution coefficients for the target effect in 
 * the input file.
 */

int TargetEffect::NumConvCoefficients() const {
  return convCoefficients_.size();
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
 * Returns the convolution coefficients for the given order specified in by the target effect.
 */

double TargetEffect::GetConvCoefficient(int order) const {
  return (convCoefficients_.size()>order) ? convCoefficients_[order] : 1.;
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
 * Returns the Equation object corresponding to the parametrized convolution.
 */

Equation *TargetEffect::GetConvolutionEq() {
  Equation *tempEquation;
  tempEquation=&convolutionEq_;
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

/*!
 * Calculates the sigma of the Gaussian beam convolution as a function of energy.
 */

double TargetEffect::CalculateSigma(double energy, const Config& configure) {
  double sigma=0;
  // Check if convolution equation is defined
  if(this->GetConvolutionEq()->GetEquation()!="") {
    sigma=this->GetConvolutionEq()->Evaluate(configure,energy);
  }
  else sigma=this->GetSigma();
  return sigma;
}

/*!
 * Returns the multiplicative convolution factor for evaluation of the integrand
 * of a target effect with energy dependent sigma.  
 */

double TargetEffect::CalculateConvolutionFactor(double energy, double centroid, const Config& configure) {
  double sigma=this->CalculateSigma(energy,configure);
  return pow(2.*pi,-0.5)/sigma*exp(-pow(energy-centroid,2.0)/2.0/pow(sigma,2.0));
}

// ERYA straggling integration methods (commented for now)
/*
bool TargetEffect::IsStraggling() const {
  return isStraggling_;
}

int TargetEffect::GetTargetElement() const {
  return targetElement_;
}

void TargetEffect::SetTargetElement(int element) {
  targetElement_ = element;
  
  // Reinitialize straggling calculator if needed
  if (isStraggling_ && isTargetIntegration_) {
    stragglingCalculator_ = std::make_unique<Straggling>(element, 1000.0);
    stragglingCalculator_->setOxideLayer(true, 0.15);
  }
}

double TargetEffect::CalculateDoubleConvolution(double energy, double centroid, const Config& configure) {
  // This method implements the double convolution approach from ERYA notebook
  // It combines Gaussian beam convolution with energy straggling
  
  if (!isStraggling_ || !stragglingCalculator_) {
    // Fall back to regular convolution if straggling is not enabled
    return CalculateConvolutionFactor(energy, centroid, configure);
  }
  
  // Calculate energy loss in target
  double deltaE = centroid - energy;
  if (deltaE < 0.0) deltaE = 0.0;
  
  // Get beam sigma
  double beamSigma = CalculateSigma(energy, configure);
  
  // Calculate double convolution using ERYA methodology
  std::vector<double> resultEnergy, resultIntensity;
  stragglingCalculator_->calculateDoubleConvolution(centroid, beamSigma, deltaE, 
                                                  resultEnergy, resultIntensity);
  
  // Find the intensity at the requested energy point
  // This is a simplified lookup - full implementation would use interpolation
  double convolutionFactor = 0.0;
  double minDistance = 1e10;
  
  for (size_t i = 0; i < resultEnergy.size(); ++i) {
    double distance = std::abs(resultEnergy[i] - energy);
    if (distance < minDistance) {
      minDistance = distance;
      convolutionFactor = resultIntensity[i];
    }
  }
  
  return convolutionFactor;
}
*/
