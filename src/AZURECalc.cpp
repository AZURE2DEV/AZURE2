#include "AZURECalc.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "GSLException.h"
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

  //Fill Compound Nucleus From Minuit Parameters
  localCompound->FillCompoundFromParams(p);
  localData->FillNormsFromParams(p);
  localData->FillEnergyShiftsFromParams(p,localData,localCompound,&configure());
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  // Process segments with components - use new integrated calculation method
  double chiSquared=0.0;
  for(int i = 1; i <= localData->NumSegments(); i++) {
    ESegment* segment = localData->GetSegment(i);
    if(segment) {
      // Recalculate points using the new combined calculation method
      for(int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        double theoreticalValue = segment->CalculateTheoreticalCrossSection(pointIdx, localCompound, configure(), localData);
        EPoint* point = segment->GetPoint(pointIdx + 1);
        if(point) {
          point->SetFitCrossSection(theoreticalValue);
        }
      }
      
      // Recalculate chi-squared for this segment with components
      double segmentChiSquared = 0.0;
      for(int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
        EPoint* point = segment->GetPoint(pointIdx + 1);
        if(point) {
          double residual = point->GetFitCrossSection() - point->GetCMCrossSection() * segment->GetNorm();
          double error = point->GetCMCrossSectionError() * segment->GetNorm();
          if(error != 0.0) {
            segmentChiSquared += (residual * residual) / (error * error);
          }
        }
      }

      // Add normalization chi-squared contribution
      double dataNorm=segment->GetNorm();
      double dataNormNominal=segment->GetNominalNorm();
      double dataNormError=dataNormNominal/100.*segment->GetNormError();
      if(dataNormError!=0.){
	      chiSquared += pow((dataNorm-dataNormNominal)/dataNormError,2.0);
      }

      // Add energy shift chi-squared contribution
      if(segment->IsVaryEnergyShift()) {
        double energyShift=segment->GetEnergyShift();
        double energyShiftNominal=segment->GetNominalEnergyShift();
        double energyShiftError=segment->GetEnergyShiftError();
        if(energyShiftError!=0.){
          chiSquared += pow((energyShift-energyShiftNominal)/energyShiftError,2.0);
        }
      }

      segment->SetSegmentChiSquared(segmentChiSquared);
      chiSquared += segmentChiSquared;
    }
  }

  // Add nuisance parameter chi-squared contributions
  if(limitsManager_) {
    chiSquared += CalculateNuisanceChiSquared(p);
  }

  if(!localData->IsErrorAnalysis()&&thisIteration!=0) {
    if(thisIteration%10==0) configure().outStream
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
