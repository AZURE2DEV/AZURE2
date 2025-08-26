#include "AZURECalc.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include <iostream>
#include <iomanip>

double AZURECalc::operator()(const vector_r& p) const {

  int thisIteration=data()->Iterations();
  data()->Iterate();
  bool isFit=data()->IsFit();

  CNuc * localCompound = NULL;
  EData *localData = NULL;
  if(isFit) {
    localCompound = compound()->Clone();
    localData = data()->Clone();
  } else {
    localCompound = compound();
    localData = data();
  }

  // Dump all parameters values
  for (size_t i = 0; i < p.size(); ++i) {
    std::cout << "Parameter " << i << ": " << p[i] << std::endl;
  }

  //Fill Compound Nucleus From Minuit Parameters
  localCompound->FillCompoundFromParams(p);
  localData->FillNormsFromParams(p);
  localData->FillEnergyShiftsFromParams(p,localData,localCompound,&configure());
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());
  
  //loop over segments and points
  double chiSquared=0.0;
  double segmentChiSquared=0.0;
  ESegmentIterator firstSumIterator = localData->GetSegments().end();
  ESegmentIterator lastSumIterator = localData->GetSegments().end();
  for(EDataIterator data=localData->begin();data!=localData->end();data++) {
    if(data.segment()->GetPoints().begin()==data.point()) {
      segmentChiSquared=0.0;
      if(data.segment()->IsTotalCapture()) {
	firstSumIterator=data.segment();
	lastSumIterator=data.segment()+data.segment()->IsTotalCapture()-1;
      } 
    }
    if(!data.point()->IsMapped()) data.point()->Calculate(localCompound,configure());
    if(firstSumIterator!=localData->GetSegments().end()&&
       data.segment()!=lastSumIterator) continue;
    double fitCrossSection=data.point()->GetFitCrossSection();
    ESegmentIterator thisSegment = data.segment();
    if(data.segment()==lastSumIterator) {
      int pointIndex=data.point()-data.segment()->GetPoints().begin()+1;
      for(ESegmentIterator it=firstSumIterator;it<data.segment();it++) 
	fitCrossSection+=it->GetPoint(pointIndex)->GetFitCrossSection();
      thisSegment = firstSumIterator;
    }
    double dataNorm=thisSegment->GetNorm();
    double CrossSection=data.point()->GetCMCrossSection()*dataNorm;
    double CrossSectionError=data.point()->GetCMCrossSectionError()*dataNorm;
    double chi=(fitCrossSection-CrossSection)/CrossSectionError;
    double pointChiSquared=pow(chi,2.0);
    segmentChiSquared+=pointChiSquared;
    if(data.segment()->GetPoints().end()-1==data.point()) {
      if(!isFit) thisSegment->SetSegmentChiSquared(segmentChiSquared);
      if(data.segment()==lastSumIterator) {
	firstSumIterator=localData->GetSegments().end();
	lastSumIterator=localData->GetSegments().end();
      }
      double dataNormNominal=thisSegment->GetNominalNorm();
      double dataNormError=dataNormNominal/100.*thisSegment->GetNormError();
      if(dataNormError!=0.)
	segmentChiSquared += pow((dataNorm-dataNormNominal)/dataNormError,2.0);
      
      // Add energy shift chi-squared contribution if energy shift is varied
      if(thisSegment->IsVaryEnergyShift()) {
        double energyShift=thisSegment->GetEnergyShift();
        double energyShiftNominal=thisSegment->GetNominalEnergyShift();
        double energyShiftError=thisSegment->GetEnergyShiftError();
        if(energyShiftError!=0.)
          segmentChiSquared += pow((energyShift-energyShiftNominal)/energyShiftError,2.0);
      }
      
      chiSquared+=segmentChiSquared;
    }
  }

  std::cout << "Chi2: " << chiSquared << std::endl;

  // Add nuisance parameter chi-squared contributions
  if(limitsManager_) {
    chiSquared += CalculateNuisanceChiSquared(p);
  }

  if(!localData->IsErrorAnalysis()&&thisIteration!=0) {
    if(thisIteration%100==0) configure().outStream
			       << "\r\tIteration: " << std::setw(6) << thisIteration
			       << " Chi-Squared: " << chiSquared;  configure().outStream.flush();

    if(thisIteration%1000==0) {
      localData->WriteOutputFiles(configure(),isFit);
      localCompound->TransformOut(configure());
      localCompound->PrintTransformParams(configure());
    }
  }
  if(isFit) {
    delete localCompound;
    delete localData;
  }
  if(configure().stopFlag&&isFit) return 0.;
  else return chiSquared;
}

double AZURECalc::CalculateNuisanceChiSquared(const vector_r& p) const {
  double nuisanceChiSquared = 0.0;
  
  // Create temporary AZUREParams to get parameter names
  AZUREParams tempParams;
  compound()->FillMnParams(tempParams.GetMinuitParams());
  data()->FillMnParams(tempParams.GetMinuitParams());
  
  // Check each parameter to see if it's marked as nuisance
  for(int i = 0; i < tempParams.GetMinuitParams().Params().size() && i < p.size(); i++) {
    std::string paramName = tempParams.GetMinuitParams().Parameter(i).GetName();

    // If norm or shift in param name, skip
    if(paramName.find("norm") != std::string::npos || paramName.find("shift") != std::string::npos) {
      continue;
    }
    
    // Check if this parameter is marked as nuisance in ParameterLimitsManager
    if(limitsManager_->IsNuisanceParameter(paramName)) {
      double paramValue = p[i];
      // Use converted nominal value (physical to reduced for width parameters)
      double nominalValue = limitsManager_->GetConvertedNominalValue(paramName);
      // Use converted error (physical to reduced for width parameters)
      double paramError = limitsManager_->GetConvertedError(paramName);
      
      if(paramError > 0.0) {
        double deviation = (paramValue - nominalValue) / paramError;
        nuisanceChiSquared += deviation * deviation;
      }
    }
  }
  
  return nuisanceChiSquared;
}
