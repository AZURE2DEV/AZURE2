#include "AZUREAPI.h"
#include "AZUREParams.h"

#include "GSLException.h"
#include "CoulFuncCache.h"
#include "ECAmplitudeCache.h"

#include "Config.h"
#include "CNuc.h"
#include "EData.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <new>
#include <cmath>

bool AZUREAPI::Initialize( ){

  // Initialize caches for performance
  InitializeCoulFuncCache();
  InitializeECAmplitudeCache();

  configure().paramMask |= Config::USE_EXTERNAL_CAPTURE;

  std::string file;
  if( configure().paramMask & Config::CALCULATE_WITH_DATA ) file = configure().outputdir + "intEC.dat";
  else file = configure().outputdir + "intEC.extrap";

  std::ifstream in(file.c_str());
  if( !in ) configure().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  else configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  configure().integralsfile=file;

  // Initialize EC Integral caching system

  std::string cacheFile;
  if (configure().paramMask & Config::CALCULATE_WITH_DATA) {
    cacheFile = configure().outputdir + "intEC_cache.dat";
  } else {
    cacheFile = configure().outputdir + "intEC_cache.extrap";
  }

  // FIXME: It crashes on Linux (but fine on Mac)
  //if( compound_ != nullptr ) delete compound_;
  //if( data_ != nullptr ) delete data_;

  data_ = new EData( );
  compound_ = new CNuc( );

  //configure().outStream << "Filling Compound Nucleus..." << std::endl;
  if(compound()->Fill(configure())==-1) {
    //configure().outStream << "Could not fill compound nucleus from file." << std::endl;
    return -1;
  } else if(compound()->NumPairs()==0 || compound()->NumJGroups()==0) {
    //configure().outStream << "No nuclear data exists. Calculation not possible." << std::endl; 
    return -1;
  } 
  if((configure().screenCheckMask|configure().fileCheckMask) & 
     Config::CHECK_COMPOUND_NUCLEUS) compound()->PrintNuc(configure());

  if(!(configure().paramMask & Config::CALCULATE_REACTION_RATE)) {
    //Fill the data object from the segments and data file
    //  Compound object is passed to the function for pair key verification and
    //  center of mass conversions, s-factor conversions, etc.
    //configure().outStream << "Filling Data Structures..." << std::endl;
    if(configure().paramMask & Config::CALCULATE_WITH_DATA) {
      if(data()->Fill(configure(),compound())==-1) {
	//configure().outStream << "Could not fill data object from file." << std::endl;
	return -1;
      } else if(data()->NumSegments()==0) {
	//configure().outStream << "There is no data provided." << std::endl;
	return -1;
      }
    } else {
      if(data()->MakePoints(configure(),compound())==-1) {
	//configure().outStream << "Could not fill data object from file." << std::endl;
	return -1;
      } else if(data()->NumSegments()==0) {
	//configure().outStream << "Extrapolation segments produce no data." << std::endl;
	return -1;
      }
    } 
    if((configure().fileCheckMask|configure().screenCheckMask) & Config::CHECK_DATA)
      data()->PrintData(configure());
  } else {
    if(!compound()->IsPairKey(configure().rateParams.entrancePair)||!compound()->IsPairKey(configure().rateParams.exitPair)) {
      //configure().outStream << "Reaction rate pairs do not exist in compound nucleus." << std::endl;
      return -1;
    } else {
      compound()->GetPair(compound()->GetPairNumFromKey(configure().rateParams.entrancePair))->SetEntrance();
    }
  }

  //Initialize compound nucleus object
  try {
    compound()->Initialize(configure());
  } catch (GSLException e) {
    configure().outStream << e.what() << std::endl;
    configure().outStream << std::endl
			  << "Calculation was aborted." << std::endl;
    return -1;
  }

  UpdateParameters( );

  if(data()->Initialize(compound(),configure())==-1) return -1;

  return 0;
  
}

bool AZUREAPI::UpdateParameters( ) {

  all_.clear( );
  all_rwa_.clear( );
  names_.clear( );
  fixed_.clear( );
  values_.clear( );
  values_rwa_.clear( );

  AZUREParams params;
  compound()->FillMnParams(params.GetMinuitParams(), &configure());
  data()->FillMnParams(params.GetMinuitParams());

  compound()->FillCompoundFromParams(params.GetMinuitParams( ).Params( ));

  compound()->CalcShiftFunctions( configure() );
  compound()->TransformOut( configure() );

  for(int i = 0; i < params.GetMinuitParams().Params().size(); i++){
    //if( !params.GetMinuitParams().Parameter(i).IsFixed( ) ){
    names_.push_back( params.GetMinuitParams().Parameter(i).GetName() );
    //}
    all_rwa_.push_back( params.GetMinuitParams().Parameter(i).Value() );
    fixed_.push_back( params.GetMinuitParams().Parameter(i).IsFixed() );
    if( !fixed_.back() ) {
      values_rwa_.push_back( params.GetMinuitParams().Parameter(i).Value() );
    }

  }

  all_ = compound()->GetTransformParams( configure() );

  // Norms and shifts are missing from all_ and need to be added using the all_rwa_ values
  for( int i = all_.size(); i < all_rwa_.size(); ++i ){
    all_.push_back( all_rwa_[i] );
  }

  for (int i = 0; i < all_.size(); ++i) {
    //std::cout << "all_[" << i << "] = " << all_[i] << std::endl;
    //std::cout << "fixed_[" << i << "] = " << fixed_[i] << std::endl;
    //std::cout << "names_[" << i << "] = " << names_[i] << std::endl;
    if( !fixed_[i] ){
      values_.push_back( all_[i] );
    }
  }

  return true;

}

int AZUREAPI::UpdateSegments(vector_r& p) {

  calculatedConv_.clear( );
  calculatedEnergies_.clear( );
  calculatedAngles_.clear( );
  calculatedSegments_.clear( );
  calculatedSegmentsE1_.clear( );
  calculatedSegmentsE2_.clear( );
  calculatedExcitationEnergies_.clear( );

  int k = 0;
  vector_r params_ = all_;
  for( int i = 0; i < all_.size( ); ++i ){
    if( !fixed_[i] ){
      params_[i] = p[k];
      ++k;
    }
  }

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  localCompound = compound();
  localData = data();

  AZUREParams params;
  localCompound->FillCompoundFromParamsPhysical(params_);
  bool isValid = localCompound->TransformIn( configure( ) );

  if( !isValid ) return 0;

  localCompound->FillMnParams(params.GetMinuitParams());
  localData->FillMnParams(params.GetMinuitParams());
  localData->FillEnergyShiftsFromParams(params_,localData,localCompound,&configure());
  localCompound->FillCompoundFromParams(params.GetMinuitParams( ).Params( ));
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  int newKey  = -1;
  int prevKey = -1;
  int nSegments = 0;

  std::vector<ESegment>& segments = localData->GetSegments( );
  for( int i = 0; i < segments.size( ); ++i ){

    newKey = segments[i].GetSegmentKey( );
    if( prevKey == newKey ) continue;
    prevKey = newKey; ++nSegments;

    std::vector<EPoint>& data = segments[i].GetPoints();

    std::vector<double> cross, crossE1, crossE2, energies, angles, conv, excitationEnergies;

    // Handle component segments using the new integrated calculation method
    if (segments[i].HasComponents()) {
      for( int k = 0; k < data.size( ); ++k ){
        // Use the new component-aware calculation method
        double theoreticalValue = segments[i].CalculateTheoreticalCrossSection(k, localCompound, configure(), localData);

        // Update the point's fit cross section with the combined result
        data[k].SetFitCrossSection(theoreticalValue);

        cross.push_back( theoreticalValue );
        crossE1.push_back( data[k].GetFitE1CrossSection() );
        crossE2.push_back( data[k].GetFitE2CrossSection() );
        angles.push_back( data[k].GetCMAngle() );
        energies.push_back( data[k].GetCMEnergy( ) );
        conv.push_back( data[k].GetSFactorConversion() );
        excitationEnergies.push_back( data[k].GetExcitationEnergy( ) );
      }
    } else {
      // Regular segment calculation (existing logic)
      for( int k = 0; k < data.size( ); ++k ){

        if(!data[k].IsMapped()) {
          try {
            data[k].Calculate(localCompound,configure());
          } catch (GSLException& e) {
            // Skip this point if GSL calculation fails
            continue;
          }
        }

        cross.push_back( data[k].GetFitCrossSection() );
        crossE1.push_back( data[k].GetFitE1CrossSection() );
        crossE2.push_back( data[k].GetFitE2CrossSection() );
        angles.push_back( data[k].GetCMAngle() );
        energies.push_back( data[k].GetCMEnergy( ) );
        conv.push_back( data[k].GetSFactorConversion() );
        excitationEnergies.push_back( data[k].GetExcitationEnergy( ) );

      }
    }

    calculatedConv_.push_back( conv );
    calculatedSegments_.push_back( cross );
    calculatedSegmentsE1_.push_back( crossE1 );
    calculatedSegmentsE2_.push_back( crossE2 );
    calculatedEnergies_.push_back( energies );
    calculatedAngles_.push_back( angles );
    calculatedExcitationEnergies_.push_back( excitationEnergies );

  }

  return calculatedSegments_.size( );

}

int AZUREAPI::UpdateSegmentsRWA(vector_r& p) {

  calculatedConv_.clear( );
  calculatedEnergies_.clear( );
  calculatedAngles_.clear( );
  calculatedSegments_.clear( );
  calculatedSegmentsE1_.clear( );
  calculatedSegmentsE2_.clear( );
  calculatedExcitationEnergies_.clear( );

  int k = 0;
  vector_r params_ = all_rwa_;
  for( int i = 0; i < all_rwa_.size( ); ++i ){
    if( !fixed_[i] ){
      params_[i] = p[k];
      ++k;
    }
  }

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  localCompound = compound();
  localData = data();

  AZUREParams params;
  localCompound->FillCompoundFromParams(params_);
  localData->FillEnergyShiftsFromParams(params_,localData,localCompound,&configure());
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  int newKey  = -1;
  int prevKey = -1;
  int nSegments = 0;

  std::vector<ESegment>& segments = localData->GetSegments( );
  for( int i = 0; i < segments.size( ); ++i ){

    newKey = segments[i].GetSegmentKey( );
    if( prevKey == newKey ) continue;
    prevKey = newKey; ++nSegments;

    std::vector<EPoint>& data = segments[i].GetPoints();

    std::vector<double> cross, crossE1, crossE2, energies, angles, conv, excitationEnergies;

    // Handle component segments using the new integrated calculation method
    if (segments[i].HasComponents()) {
      for( int k = 0; k < data.size( ); ++k ){
        // Use the new component-aware calculation method
        double theoreticalValue = segments[i].CalculateTheoreticalCrossSection(k, localCompound, configure(), localData);

        // Update the point's fit cross section with the combined result
        data[k].SetFitCrossSection(theoreticalValue);

        cross.push_back( theoreticalValue );
        crossE1.push_back( data[k].GetFitE1CrossSection() );
        crossE2.push_back( data[k].GetFitE2CrossSection() );
        angles.push_back( data[k].GetCMAngle() );
        energies.push_back( data[k].GetCMEnergy( ) );
        conv.push_back( data[k].GetSFactorConversion() );
        excitationEnergies.push_back( data[k].GetExcitationEnergy( ) );
      }
    } else {
      // Regular segment calculation (existing logic)
      for( int k = 0; k < data.size( ); ++k ){

        if(!data[k].IsMapped()) {
          try {
            data[k].Calculate(localCompound,configure());
          } catch (GSLException& e) {
            // Skip this point if GSL calculation fails
            continue;
          }
        }

        cross.push_back( data[k].GetFitCrossSection() );
        crossE1.push_back( data[k].GetFitE1CrossSection() );
        crossE2.push_back( data[k].GetFitE2CrossSection() );
        angles.push_back( data[k].GetCMAngle() );
        energies.push_back( data[k].GetCMEnergy( ) );
        conv.push_back( data[k].GetSFactorConversion() );
        excitationEnergies.push_back( data[k].GetExcitationEnergy( ) );

      }
    }

    calculatedConv_.push_back( conv );
    calculatedSegments_.push_back( cross );
    calculatedSegmentsE1_.push_back( crossE1 );
    calculatedSegmentsE2_.push_back( crossE2 );
    calculatedEnergies_.push_back( energies );
    calculatedAngles_.push_back( angles );
    calculatedExcitationEnergies_.push_back( excitationEnergies );

  }

  return calculatedSegments_.size( );

}

int AZUREAPI::UpdateSegmentsAllRWA(vector_r& p) {

  calculatedConv_.clear( );
  calculatedEnergies_.clear( );
  calculatedAngles_.clear( );
  calculatedSegments_.clear( );
  calculatedSegmentsE1_.clear( );
  calculatedSegmentsE2_.clear( );
  calculatedExcitationEnergies_.clear( );

  int k = 0;
  for( int i = 0; i < all_rwa_.size( ); ++i ){
    all_rwa_[i] = p[i];
    if( !fixed_[i] ){
      values_[k] = p[i];
      ++k;
    }
  }
  vector_r params_ = all_rwa_;

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  localCompound = compound();
  localData = data();

  AZUREParams params;
  localCompound->FillCompoundFromParams(params_);
  localData->FillEnergyShiftsFromParams(params_,localData,localCompound,&configure());
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  int newKey  = -1;
  int prevKey = -1;
  int nSegments = 0;

  std::vector<ESegment>& segments = localData->GetSegments( );
  for( int i = 0; i < segments.size( ); ++i ){

    newKey = segments[i].GetSegmentKey( );
    if( prevKey == newKey ) continue;
    prevKey = newKey; ++nSegments;

    std::vector<EPoint>& data = segments[i].GetPoints();

    std::vector<double> cross, crossE1, crossE2, energies, angles, conv, excitationEnergies;

    // Handle component segments using the new integrated calculation method
    if (segments[i].HasComponents()) {
      for( int k = 0; k < data.size( ); ++k ){
        // Use the new component-aware calculation method
        double theoreticalValue = segments[i].CalculateTheoreticalCrossSection(k, localCompound, configure(), localData);

        // Update the point's fit cross section with the combined result
        data[k].SetFitCrossSection(theoreticalValue);

        cross.push_back( theoreticalValue );
        crossE1.push_back( data[k].GetFitE1CrossSection() );
        crossE2.push_back( data[k].GetFitE2CrossSection() );
        angles.push_back( data[k].GetCMAngle() );
        energies.push_back( data[k].GetCMEnergy( ) );
        conv.push_back( data[k].GetSFactorConversion() );
        excitationEnergies.push_back( data[k].GetExcitationEnergy( ) );
      }
    } else {
      // Regular segment calculation (existing logic)
      for( int k = 0; k < data.size( ); ++k ){

        if(!data[k].IsMapped()) {
          try {
            data[k].Calculate(localCompound,configure());
          } catch (GSLException& e) {
            // Skip this point if GSL calculation fails
            continue;
          }
        }

        cross.push_back( data[k].GetFitCrossSection() );
        crossE1.push_back( data[k].GetFitE1CrossSection() );
        crossE2.push_back( data[k].GetFitE2CrossSection() );
        angles.push_back( data[k].GetCMAngle() );
        energies.push_back( data[k].GetCMEnergy( ) );
        conv.push_back( data[k].GetSFactorConversion() );
        excitationEnergies.push_back( data[k].GetExcitationEnergy( ) );

      }
    }

    calculatedConv_.push_back( conv );
    calculatedSegments_.push_back( cross );
    calculatedSegmentsE1_.push_back( crossE1 );
    calculatedSegmentsE2_.push_back( crossE2 );
    calculatedEnergies_.push_back( energies );
    calculatedAngles_.push_back( angles );
    calculatedExcitationEnergies_.push_back( excitationEnergies );

  }

  return calculatedSegments_.size( );

}

// Transform RWA parameters to physical values
vector_r AZUREAPI::TransformRWAParameters(const vector_r& p) const {

  vector_r params = all_rwa_;
  int k = 0;
  for( int i = 0; i < p.size( ); ++i ){
    if( !fixed_[i] ){
      params[i] = p[k];
      ++k;
    }
  }

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  localCompound = compound();
  localData = data();

  localCompound->FillCompoundFromParams(params);

  localCompound->TransformOut( configure() );

  vector_r transformedParams = compound()->GetTransformParams( configure() );

  // Get only non fixed parameters
  vector_r transformed;
  k = 0;
  for( int i = 0; i < transformedParams.size( ); ++i ){
    if( !fixed_[i] ){
      transformed.push_back( transformedParams[i] );
      ++k;
    }
  }

  return transformed;

}

// Transform RWA parameters to physical values
vector_r AZUREAPI::TransformAllRWAParameters(const vector_r& p) const {

  vector_r params = all_rwa_;
  int k = 0;
  for( int i = 0; i < p.size( ); ++i ){
    params[i] = p[k];
    ++k;
  }

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  localCompound = compound();
  localData = data();

  localCompound->FillCompoundFromParams(params);

  localCompound->TransformOut( configure() );

  vector_r transformedParams = compound()->GetTransformParams( configure() );

  // Get only non fixed parameters
  vector_r transformed;
  k = 0;
  for( int i = 0; i < transformedParams.size( ); ++i ){
    if( !fixed_[i] ){
      transformed.push_back( transformedParams[i] );
      ++k;
    }
  }

  return transformed;

}

bool AZUREAPI::CalculateExternalCapture( ){

  configure().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  data()->CalculateECAmplitudes( compound( ), configure( ) );
  configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  return true;

}

int AZUREAPI::UpdateData( ) {

  dataEnergies_.clear( );
  dataAngles_.clear( );
  dataSegments_.clear( );
  dataSegmentsErrors_.clear( );
  dataConv_.clear( );
  dataExcitationEnergies_.clear( );

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  localCompound = compound();
  localData = data();

  int newKey  = -1;
  int prevKey = -1;
  int nSegments = 0;

  std::vector<ESegment>& segments = localData->GetSegments( );
  for( int i = 0; i < segments.size( ); ++i ){

    newKey = segments[i].GetSegmentKey( );
    if( prevKey == newKey ) continue;
    prevKey = newKey; ++nSegments;

    std::vector<EPoint>& data = segments[i].GetPoints();

    std::vector<double> energies, angles, cross, crossErr, conv, excitationEnergies;
    for( int k = 0; k < data.size( ); ++k ){

      energies.push_back( data[k].GetCMEnergy( ) );
      angles.push_back( data[k].GetCMAngle( ) );
      cross.push_back( data[k].GetCMCrossSection() );
      crossErr.push_back( data[k].GetCMCrossSectionError() );
      conv.push_back( data[k].GetSFactorConversion() );
      excitationEnergies.push_back( data[k].GetExcitationEnergy( ) );

    }

    dataEnergies_.push_back( energies );
    dataSegments_.push_back( cross );
    dataAngles_.push_back( angles );
    dataSegmentsErrors_.push_back( crossErr );
    dataConv_.push_back( conv );
    dataExcitationEnergies_.push_back( excitationEnergies );

  }

  return nSegments;

}

void AZUREAPI::UpdateNorms( ) {

  norms_.clear( );
  normsErrors_.clear( );

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  localCompound = compound();
  localData = data();

  int newKey  = -1;
  int prevKey = -1;
  int nSegments = 0;

  std::vector<ESegment>& segments = localData->GetSegments( );
  for( int i = 0; i < segments.size( ); ++i ){
    
    newKey = segments[i].GetSegmentKey( );
    if( prevKey == newKey ) continue;
    prevKey = newKey; ++nSegments;

    double norm = segments[i].GetNominalNorm( );
    double normErr = segments[i].GetNormError( );

    norms_.push_back( norm );
    normsErrors_.push_back( normErr );

  }

}

// Set AZURE2 to calculate data points
void AZUREAPI::SetData( ) { 
  configure().paramMask |= Config::CALCULATE_WITH_DATA; 
}

// Set AZURE2 to calculate extrapolations
void AZUREAPI::SetExtrap( ) { 
  configure().paramMask &= ~Config::CALCULATE_WITH_DATA; 
}

// Set radius to a fixed value
void AZUREAPI::SetRadius( int idx, double r ) {
  
  if( compound_ != nullptr ) delete compound_;
  if( data_ != nullptr ) delete data_;

  compound_ = new CNuc;
  data_     = new EData;

  std::pair<int,double> pair = std::make_pair( idx, r );

  compound()->Fill( configure( ), pair  );
  data()->Fill(configure(),compound());

  configure().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  compound( )->Initialize( configure( ) );
  data( )->Initialize( compound( ), configure( ) );
  configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;

}

double AZUREAPI::CalculateChi2RWA(const vector_r& rwaParams) const {

  int k = 0;
  vector_r params_ = all_rwa_;
  for( int i = 0; i < all_rwa_.size( ); ++i ){
    if( !fixed_[i] ){
      params_[i] = rwaParams[k];
      ++k;
    }
  }

  double chiSquared = 0.0;
  
  CNuc* localCompound = compound()->Clone();
  EData* localData = data()->Clone();
  
  // Fill compound nucleus and data with RWA parameters
  localCompound->FillCompoundFromParams(params_);
  localData->FillNormsFromParams(params_);
  localData->FillEnergyShiftsFromParams(params_, localData, localCompound, &configure());
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());
  
  // Process segments with components - use new integrated calculation method
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

      segment->SetSegmentChiSquared(segmentChiSquared);
      chiSquared += segmentChiSquared;
    }
  }
  
  delete localCompound;
  delete localData;
  
  return chiSquared;
}

double AZUREAPI::CalculateChi2Physical(const vector_r& physicalParams) const {
  int k = 0;
  vector_r params_ = all_;
  for( int i = 0; i < all_.size( ); ++i ){
    if( !fixed_[i] ){
      params_[i] = physicalParams[k];
      ++k;
    }
  }

  double chiSquared = 0.0;

  CNuc* localCompound = NULL;
  EData* localData = NULL;
  localCompound = compound()->Clone();
  localData = data()->Clone();

  AZUREParams params;
  localCompound->FillCompoundFromParamsPhysical(params_);
  bool isValid = localCompound->TransformIn( configure( ) );

  if( !isValid ) return 0;

  localCompound->FillMnParams(params.GetMinuitParams());
  localData->FillMnParams(params.GetMinuitParams());
  localData->FillEnergyShiftsFromParams(params_,localData,localCompound,&configure());
  localCompound->FillCompoundFromParams(params.GetMinuitParams( ).Params( ));
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());
  
  // Process segments with components - use new integrated calculation method
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

      segment->SetSegmentChiSquared(segmentChiSquared);
      chiSquared += segmentChiSquared;
    }
  }
  
  delete localCompound;
  delete localData;

  return chiSquared;
}

double AZUREAPI::CalculateLnLRWA(const vector_r& params) const {

  // ln(2*pi), used in the Gaussian error-normalization term. Defined as a
  // literal so it does not depend on M_PI being available on every platform.
  static const double kLn2Pi = 1.8378770664093453;

  // The input vector packs the non-fixed RWA parameters first, followed by one
  // error-inflation factor per segment. Split it on the number of non-fixed
  // parameters.
  int nRwa = 0;
  for( size_t i = 0; i < fixed_.size( ); ++i ) if( !fixed_[i] ) ++nRwa;

  vector_r rwaParams( params.begin( ), params.begin( ) + nRwa );
  vector_r inflation( params.begin( ) + nRwa, params.end( ) );

  // Map the non-fixed RWA values back into the full parameter vector.
  int k = 0;
  vector_r params_ = all_rwa_;
  for( int i = 0; i < all_rwa_.size( ); ++i ){
    if( !fixed_[i] ){
      params_[i] = rwaParams[k];
      ++k;
    }
  }

  CNuc* localCompound = compound()->Clone();
  EData* localData = data()->Clone();

  // Fill compound nucleus and data with RWA parameters (norms included).
  localCompound->FillCompoundFromParams(params_);
  localData->FillNormsFromParams(params_);
  localData->FillEnergyShiftsFromParams(params_, localData, localCompound, &configure());
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  double lnL = 0.0;

  // Walk the raw segments so every point contributes, but advance the
  // inflation index only when the segment key changes, keeping it aligned with
  // the norms() / UpdateData() ordering (one inflation factor per segment).
  int prevKey = -1;
  int segIdx  = -1;

  for(int i = 1; i <= localData->NumSegments(); i++) {
    ESegment* segment = localData->GetSegment(i);
    if(!segment) continue;

    int key = segment->GetSegmentKey();
    if( key != prevKey ){ ++segIdx; prevKey = key; }

    double f = ( segIdx >= 0 && segIdx < (int)inflation.size( ) ) ? inflation[segIdx] : 0.0;
    double norm = segment->GetNorm();

    for(int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
      double model = segment->CalculateTheoreticalCrossSection(pointIdx, localCompound, configure(), localData);
      EPoint* point = segment->GetPoint(pointIdx + 1);
      if(!point) continue;
      point->SetFitCrossSection(model);

      double residual  = model - point->GetCMCrossSection() * norm;
      double baseError = point->GetCMCrossSectionError() * norm;
      double inflTerm  = f * model;
      double var = baseError * baseError + inflTerm * inflTerm;

      if( var > 0.0 ) {
        lnL += -0.5 * ( residual * residual / var + std::log( var ) + kLn2Pi );
      }
    }
  }

  delete localCompound;
  delete localData;

  return lnL;
}

double AZUREAPI::CalculateLnLCovRWA(const vector_r& params) const {

  // ln(2*pi) literal (see CalculateLnLRWA).
  static const double kLn2Pi = 1.8378770664093453;

  // Split the input into non-fixed RWA parameters and per-segment inflation.
  int nRwa = 0;
  for( size_t i = 0; i < fixed_.size( ); ++i ) if( !fixed_[i] ) ++nRwa;

  vector_r rwaParams( params.begin( ), params.begin( ) + nRwa );
  vector_r inflation( params.begin( ) + nRwa, params.end( ) );

  int k = 0;
  vector_r params_ = all_rwa_;
  for( int i = 0; i < all_rwa_.size( ); ++i ){
    if( !fixed_[i] ){
      params_[i] = rwaParams[k];
      ++k;
    }
  }

  CNuc* localCompound = compound()->Clone();
  EData* localData = data()->Clone();

  localCompound->FillCompoundFromParams(params_);
  localData->FillNormsFromParams(params_);
  localData->FillEnergyShiftsFromParams(params_, localData, localCompound, &configure());
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) localCompound->CalcShiftFunctions(configure());

  double lnL = 0.0;

  int prevKey = -1;
  int segIdx  = -1;

  for(int i = 1; i <= localData->NumSegments(); i++) {
    ESegment* segment = localData->GetSegment(i);
    if(!segment) continue;

    int key = segment->GetSegmentKey();
    if( key != prevKey ){ ++segIdx; prevKey = key; }

    double f = ( segIdx >= 0 && segIdx < (int)inflation.size( ) ) ? inflation[segIdx] : 0.0;
    double norm = segment->GetNorm();

    // The covariance block is diagonal-plus-rank-1:
    //   C = D + v v^T,   D = diag(sv_i),   v_i = f * model_i
    // so its inverse (Sherman-Morrison) and determinant (matrix-determinant
    // lemma) are available in closed form. This needs only the running sums
    //   A  = sum_i r_i^2 / sv_i              (= r^T D^{-1} r)
    //   W  = sum_i model_i r_i / sv_i        (-> v^T D^{-1} r = f * W)
    //   Mm = sum_i model_i^2 / sv_i          (-> v^T D^{-1} v = f^2 * Mm)
    //   lnDetD = sum_i ln(sv_i)              (= ln det D)
    // computed in a single O(n) pass, instead of forming and factorizing the
    // n x n matrix.
    int n = 0;
    double A = 0.0, W = 0.0, Mm = 0.0, lnDetD = 0.0;
    bool degenerate = false;

    for(int pointIdx = 0; pointIdx < segment->NumPoints(); pointIdx++) {
      double m = segment->CalculateTheoreticalCrossSection(pointIdx, localCompound, configure(), localData);
      EPoint* point = segment->GetPoint(pointIdx + 1);
      if(!point) continue;
      point->SetFitCrossSection(m);

      double baseError = point->GetCMCrossSectionError() * norm;
      double sv = baseError * baseError;
      if( sv <= 0.0 ){
        // Zero statistical error makes D (and, for >1 such point, C) singular:
        // the closed form is undefined. Reject rather than divide by zero.
        degenerate = true;
        break;
      }

      double resid = m - point->GetCMCrossSection() * norm;
      double inv = 1.0 / sv;
      A      += resid * resid * inv;
      W      += m * resid * inv;
      Mm     += m * m * inv;
      lnDetD += std::log( sv );
      ++n;
    }

    if( degenerate ){
      lnL = -std::numeric_limits<double>::infinity();
      break;
    }
    if( n == 0 ) continue;

    double f2 = f * f;
    double denom = 1.0 + f2 * Mm;                 // 1 + v^T D^{-1} v
    double chi2  = A - ( f2 * W * W ) / denom;     // r^T C^{-1} r
    double lnDet = lnDetD + std::log( denom );     // ln det C

    lnL += -0.5 * ( chi2 + lnDet + n * kLn2Pi );
  }

  delete localCompound;
  delete localData;

  return lnL;
}

// Function to get indeces of non-fixed normalization parameters
vector_r AZUREAPI::GetNormalizationIndices( ) {

  int k = 0;
  vector_r indices;
  int totalParams = all_rwa_.size();
  for(int i = 0; i < totalParams; ++i) {
    if(!fixed_[i]) {
      // Check if "norm" is in the parameter name
      if(names_[i].find("norm") != std::string::npos) {
        indices.push_back(k);
      }
      ++k;
    }
  }

  return indices;
}

// Function to get indeces of non-fixed energy shift parameters
vector_r AZUREAPI::GetEnergyShiftIndices( ) {

  int k = 0;
  vector_r indices;
  int totalParams = all_rwa_.size();
  for(int i = 0; i < totalParams; ++i) {
    if(!fixed_[i]) {
      // Check if "Eshift" is in the parameter name
      if(names_[i].find("shift") != std::string::npos) {
        indices.push_back(k);
      }
      ++k;
    }
  }

  return indices;
}