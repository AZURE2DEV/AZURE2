#include "AZURECalc.h"
#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "ParameterLimitsManager.h"
#include "AZUREParams.h"
#include "GSLException.h"
#include <iostream>
#include <iomanip>
#include <thread>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

double AZURECalc::operator()(const vector_r& p) const {

  int thisIteration=data()->Iterations();
  data()->Iterate();
  bool isFit=data()->IsFit();

  CNuc *localCompound = NULL;
  EData *localData = NULL;
  if(isFit) {

    // New multithreading with object pools
    if (!pools_initialized_) {
      InitializePools();
    }
    
    // Get objects from pool
    localCompound = GetPooledCNuc();
    localData = GetPooledEData();

    // Old multithreading
    //localCompound = compound()->Clone();
    //localData = data()->Clone();
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

    if(thisIteration%100==0) {
      AZUREParams params;
      localCompound->FillMnParams(params.GetMinuitParams(), &configure());
      localData->FillMnParams(params.GetMinuitParams());
      WriteParameters(params,configure());
      localData->WriteOutputFiles(configure(),isFit);
      localCompound->TransformOut(configure());
      localCompound->PrintTransformParams(configure());
    }
  }
  if(isFit) {

    // New multithreading with object pools
    ReturnPooledCNuc(localCompound);
    ReturnPooledEData(localData);
    
    // Old multithreading
    //delete localCompound;
    //delete localData;
  }

  // Make a check if chiSquared is NaN
  if(std::isnan(chiSquared)) {
      // In that case return infinite since MINUIT2 can have issues with NaN values
      return std::numeric_limits<double>::infinity();
  }

  if(configure().stopFlag&&isFit) return 0.;
  else return chiSquared;
}

double AZURECalc::CalculateNuisanceChiSquared(const vector_r& p) const {
  double nuisanceChiSquared = 0.0;
  
  // Create temporary AZUREParams to get parameter names
  AZUREParams tempParams;
  compound()->FillMnParams(tempParams.GetMinuitParams(), &configure());
  data()->FillMnParams(tempParams.GetMinuitParams());
  
  // Build mapping from non-fixed parameter index to actual parameter index
  std::vector<int> nonFixedToActualIndex;
  for(int i = 0; i < tempParams.GetMinuitParams().Params().size(); i++) {
    if(!tempParams.GetMinuitParams().Parameter(i).IsFixed() || tempParams.GetMinuitParams().Parameter(i).GetName().find("segment") != std::string::npos) {
      nonFixedToActualIndex.push_back(i);
    }
  }
  
  // Check each non-fixed parameter to see if it's marked as nuisance
  for(int nonFixedIndex = 0; nonFixedIndex < nonFixedToActualIndex.size() && nonFixedIndex < p.size(); nonFixedIndex++) {
    int actualIndex = nonFixedToActualIndex[nonFixedIndex];
    std::string paramName = tempParams.GetMinuitParams().Parameter(actualIndex).GetName();

    // If norm or shift in param name, skip
    if(paramName.find("norm") != std::string::npos || paramName.find("shift") != std::string::npos) {
      continue;
    }
    
    // First check if this parameter is marked as nuisance (fast check)
    if(!limitsManager_->IsNuisanceParameterByIndex(nonFixedIndex)) {
      continue; // Skip if not a nuisance parameter
    }
    
    // Only do expensive conversions if parameter is marked as nuisance
    double nominalValue = limitsManager_->GetConvertedNominalValueByIndex(nonFixedIndex);
    double paramError = limitsManager_->GetConvertedErrorByIndex(nonFixedIndex);
    
    // If we got valid values (non-zero error means this parameter has valid nuisance settings)
    if(paramError > 0.0) {
      double paramValue = p[nonFixedIndex];      
      double deviation = (paramValue - nominalValue) / paramError;
      nuisanceChiSquared += deviation * deviation;
    }
  }
  
  return nuisanceChiSquared;
}

/*!
 * Initialize object pools with pre-allocated CNuc and EData objects
 */
void AZURECalc::InitializePools() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);
  if (pools_initialized_) return;
  
  // Calculate pool size based on OpenMP threads (fixes interaction with OpenMP)
  int pool_size = 4; // default minimum
#ifdef _OPENMP
  pool_size = std::max(4, omp_get_max_threads());
#else
  pool_size = std::max(4, static_cast<int>(std::thread::hardware_concurrency()));
#endif
  
  // Pre-allocate CNuc objects by cloning once
  for (int i = 0; i < pool_size; ++i) {
    cnuc_pool_.push(std::unique_ptr<CNuc>(compound()->Clone()));
  }
  
  // Pre-allocate EData objects by cloning once
  for (int i = 0; i < pool_size; ++i) {
    edata_pool_.push(std::unique_ptr<EData>(data()->Clone()));
  }
  
  pools_initialized_ = true;
}

/*!
 * Get a CNuc object from the pool, creating new if pool is empty
 */
CNuc* AZURECalc::GetPooledCNuc() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);
  
  if (!cnuc_pool_.empty()) {
    auto obj = std::move(cnuc_pool_.top());
    cnuc_pool_.pop();
    return obj.release();
  }
  
  // Fallback: create new if pool is empty (shouldn't happen often)
  return compound_->Clone();
}

/*!
 * Get an EData object from the pool, creating new if pool is empty
 */
EData* AZURECalc::GetPooledEData() const {
  std::lock_guard<std::mutex> lock(pool_mutex_);
  
  if (!edata_pool_.empty()) {
    auto obj = std::move(edata_pool_.top());
    edata_pool_.pop();
    return obj.release();
  }
  
  // Fallback: create new if pool is empty (shouldn't happen often)
  return data_->Clone();
}

/*!
 * Return a CNuc object to the pool for reuse
 */
void AZURECalc::ReturnPooledCNuc(CNuc* obj) const {
  if (!obj) return;
  
  std::lock_guard<std::mutex> lock(pool_mutex_);
  cnuc_pool_.push(std::unique_ptr<CNuc>(obj));
}

/*!
 * Return an EData object to the pool for reuse  
 */
void AZURECalc::ReturnPooledEData(EData* obj) const {
  if (!obj) return;
  
  std::lock_guard<std::mutex> lock(pool_mutex_);
  edata_pool_.push(std::unique_ptr<EData>(obj));
}

/*!
 * Write parameters to file
 */
void AZURECalc::WriteParameters(AZUREParams& params, const Config& configure) const {
  char filename[256];
  sprintf(filename,"%sparam.fit",configure.outputdir.c_str());
  std::ofstream out;
  out.open(filename);
  if(out) {
    out.precision(7);
    for(int i=0;i<params.GetMinuitParams().Params().size();i++) {
      out << std::setw(20) << params.GetMinuitParams().GetName(i)
	  << std::scientific << std::setw(20) <<  params.GetMinuitParams().Value(i)
	  << std::scientific << std::setw(20) <<  params.GetMinuitParams().Error(i) << std::endl;
    }
    out.flush();
    out.close();
  } else configure.outStream << "Could not save param.fit file." << std::endl;
}