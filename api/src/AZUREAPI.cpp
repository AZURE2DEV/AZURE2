#include "AZUREAPI.h"
#include "AZUREParams.h"

#include "GSLException.h"
#include "CoulFuncCache.h"
#include "ECAmplitudeCache.h"

#include "Config.h"
#include "CNuc.h"
#include "EData.h"
#include "ESegment.h"
#include "EPoint.h"
#include "TargetEffect.h"
#include "AMatrixFunc.h"
#include "AZUREGrad.h"
#include "CovarianceBand.h"
#include "CoulFunc.h"
#include "ECIntegral.h"

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
  calculatedAngularDists_.clear( );
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
    // Angular-distribution coefficients, self-describing per point:
    // a count followed by that many Legendre coefficients. Points that
    // are not part of an angular-distribution segment contribute a zero.
    std::vector<double> angularDists;

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
        angularDists.push_back( (double)data[k].GetNumAngularDists( ) );
        for( int a = 0; a < data[k].GetNumAngularDists( ); ++a )
          angularDists.push_back( data[k].GetAngularDist( a ) );
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
        angularDists.push_back( (double)data[k].GetNumAngularDists( ) );
        for( int a = 0; a < data[k].GetNumAngularDists( ); ++a )
          angularDists.push_back( data[k].GetAngularDist( a ) );

      }
    }

    calculatedConv_.push_back( conv );
    calculatedSegments_.push_back( cross );
    calculatedSegmentsE1_.push_back( crossE1 );
    calculatedSegmentsE2_.push_back( crossE2 );
    calculatedEnergies_.push_back( energies );
    calculatedAngles_.push_back( angles );
    calculatedAngularDists_.push_back( angularDists );
    calculatedExcitationEnergies_.push_back( excitationEnergies );

  }

  return calculatedSegments_.size( );

}

int AZUREAPI::UpdateSegmentsRWA(vector_r& p) {

  calculatedConv_.clear( );
  calculatedEnergies_.clear( );
  calculatedAngles_.clear( );
  calculatedAngularDists_.clear( );
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
    // Angular-distribution coefficients, self-describing per point:
    // a count followed by that many Legendre coefficients. Points that
    // are not part of an angular-distribution segment contribute a zero.
    std::vector<double> angularDists;

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
        angularDists.push_back( (double)data[k].GetNumAngularDists( ) );
        for( int a = 0; a < data[k].GetNumAngularDists( ); ++a )
          angularDists.push_back( data[k].GetAngularDist( a ) );
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
        angularDists.push_back( (double)data[k].GetNumAngularDists( ) );
        for( int a = 0; a < data[k].GetNumAngularDists( ); ++a )
          angularDists.push_back( data[k].GetAngularDist( a ) );

      }
    }

    calculatedConv_.push_back( conv );
    calculatedSegments_.push_back( cross );
    calculatedSegmentsE1_.push_back( crossE1 );
    calculatedSegmentsE2_.push_back( crossE2 );
    calculatedEnergies_.push_back( energies );
    calculatedAngles_.push_back( angles );
    calculatedAngularDists_.push_back( angularDists );
    calculatedExcitationEnergies_.push_back( excitationEnergies );

  }

  return calculatedSegments_.size( );

}

int AZUREAPI::UpdateSegmentsAllRWA(vector_r& p) {

  calculatedConv_.clear( );
  calculatedEnergies_.clear( );
  calculatedAngles_.clear( );
  calculatedAngularDists_.clear( );
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
    // Angular-distribution coefficients, self-describing per point:
    // a count followed by that many Legendre coefficients. Points that
    // are not part of an angular-distribution segment contribute a zero.
    std::vector<double> angularDists;

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
        angularDists.push_back( (double)data[k].GetNumAngularDists( ) );
        for( int a = 0; a < data[k].GetNumAngularDists( ); ++a )
          angularDists.push_back( data[k].GetAngularDist( a ) );
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
        angularDists.push_back( (double)data[k].GetNumAngularDists( ) );
        for( int a = 0; a < data[k].GetNumAngularDists( ); ++a )
          angularDists.push_back( data[k].GetAngularDist( a ) );

      }
    }

    calculatedConv_.push_back( conv );
    calculatedSegments_.push_back( cross );
    calculatedSegmentsE1_.push_back( crossE1 );
    calculatedSegmentsE2_.push_back( crossE2 );
    calculatedEnergies_.push_back( energies );
    calculatedAngles_.push_back( angles );
    calculatedAngularDists_.push_back( angularDists );
    calculatedExcitationEnergies_.push_back( excitationEnergies );

  }

  return calculatedSegments_.size( );

}

// Transform RWA parameters to physical values
vector_r AZUREAPI::TransformRWAParameters(const vector_r& p) const {

  // p holds the non-fixed parameters only, so the loop has to run over the
  // full parameter array -- running it to p.size() leaves every parameter
  // beyond that index at its .azr value, which silently returns stale widths
  // for the levels at the end of the file.  A shorter p (e.g. only the
  // leading R-matrix block, with the norm/shift tail sliced off) is still
  // accepted: the free parameters it does not cover keep their all_rwa_ value.
  vector_r params = all_rwa_;
  int k = 0;
  for( int i = 0; i < params.size( ); ++i ){
    if( !fixed_[i] ){
      if( k < p.size( ) ) params[i] = p[k];
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

  // p holds every parameter, fixed ones included; a short p updates a prefix.
  vector_r params = all_rwa_;
  for( int i = 0; i < p.size( ) && i < params.size( ); ++i ){
    params[i] = p[i];
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
  for( int i = 0; i < transformedParams.size( ); ++i ){
    if( !fixed_[i] ){
      transformed.push_back( transformedParams[i] );
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

// Set the channel radius of one particle pair (1-based) and rebuild.
//
// The radius enters the penetrabilities, shift functions, boundary conditions,
// Wigner limits, the ANC <-> reduced-width conversion and the lower limit of
// every external-capture integral, so nothing can be reused: the compound
// nucleus and data objects are rebuilt from the .azr with the override, the
// EC-integral cache is bypassed (USE_PREVIOUS_INTEGRALS cleared) so the
// integrals are recomputed on the new radius rather than read back from
// output/intEC*, and the parameter bookkeeping -- values, names, fixed flags
// and the transformed physical parameters -- is refilled to match.
bool AZUREAPI::SetRadius( int idx, double r ) {

  if( compound_ != nullptr ) delete compound_;
  if( data_ != nullptr ) delete data_;

  compound_ = new CNuc;
  data_     = new EData;

  std::pair<int,double> pair = std::make_pair( idx, r );

  if( compound()->Fill( configure( ), pair ) == -1 ) return false;
  if( compound()->NumPairs() == 0 || compound()->NumJGroups() == 0 ) return false;

  // Mirror the branch taken at startup: data mode reads the data segments,
  // extrapolation mode builds its points from <segmentsTest>.
  if( configure().paramMask & Config::CALCULATE_WITH_DATA ) {
    if( data()->Fill( configure( ), compound( ) ) == -1 ) return false;
  } else {
    if( data()->MakePoints( configure( ), compound( ) ) == -1 ) return false;
  }
  if( data()->NumSegments() == 0 ) return false;

  configure().paramMask &= ~Config::USE_PREVIOUS_INTEGRALS;
  try {
    compound( )->Initialize( configure( ) );
  } catch (GSLException e) {
    configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;
    configure().outStream << e.what() << std::endl;
    return false;
  }
  if( data( )->Initialize( compound( ), configure( ) ) == -1 ) {
    configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;
    return false;
  }
  configure().paramMask |= Config::USE_PREVIOUS_INTEGRALS;

  // Without this the client keeps the pre-change parameter vector, names,
  // fixed flags and Wigner limits -- silently mismatched with the new model.
  UpdateParameters( );

  return true;

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

  // One request at a time: operate on the canonical compound/data in place
  // (re-filled below from these parameters) instead of cloning.
  CNuc* localCompound = compound();
  EData* localData = data();

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

  // One request at a time: operate on the canonical compound/data in place
  // (re-filled below from these parameters) instead of cloning.
  CNuc* localCompound = compound();
  EData* localData = data();

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

  return chiSquared;
}

// Map a packed non-fixed RWA parameter vector into the full all_rwa_ layout.
static vector_r MapPackedToFull(const vector_r& packed, const vector_r& all_rwa,
                                const std::vector<bool>& fixed) {
  vector_r full = all_rwa;
  int k = 0;
  for( int i = 0; i < (int)all_rwa.size( ); ++i )
    if( !fixed[i] && k < (int)packed.size() ) full[i] = packed[k++];
  return full;
}

bool AZUREAPI::Chi2GradEGammaNorm(const vector_r& full, vector_r& gradFull,
                                  double& chi2Out) const {
  const bool brune = (configure().paramMask & Config::USE_BRUNE_FORMALISM);
  // Only one API request runs at a time, so operate on the canonical
  // compound/data in place (like UpdateSegments) rather than cloning; every
  // entry point re-fills from its own parameters before use.
  CNuc* lc = compound();
  EData* ld = data();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r* sdp = nullptr;
  if(brune) {
    lc->CalcShiftFunctions(configure());
    shiftDeriv = BuildShiftDerivTable(lc, configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed_);
  GradAccum accum;
  accum.Init(lc);

  // chi2 = sum (fit - data*n)^2/(cmErr*n)^2.  d(chi2)/d(model) = 2 r / err^2;
  // the same model gives the data-term norm gradient, accumulated per segment.
  std::vector<double> normData(ld->NumSegments() + 1, 0.0);
  double chi2 = 0.0;
  FitBarFn fb = [&](ESegment* seg, int i, int pid, double model) -> double {
    EPoint* pt = seg->GetPoint(pid + 1);
    if(!pt) return 0.0;
    double norm = seg->GetNorm();
    double dataval = pt->GetCMCrossSection();
    double cmErr = pt->GetCMCrossSectionError();
    double r = model - dataval * norm;
    double err = cmErr * norm;
    if(err == 0.0) return 0.0;
    // chi2 falls out of the same residual the gradient uses (so value and
    // gradient stay consistent). Runs in the parallel point loop, so guard it.
#pragma omp atomic
    chi2 += (r * r) / (err * err);
    if(seg->IsVaryNorm() && norm != 0.0 && i >= 1 && i < (int)normData.size()) {
      double e2 = cmErr * cmErr;
      // fitBarFn runs inside the parallel point loop of AccumulateEGammaGradient,
      // and all points of a segment share the same normData[i], so guard the
      // accumulation.
      double dNorm = -2.0 * r * dataval / (e2 * norm * norm)
                     - 2.0 * r * r / (e2 * norm * norm * norm);
#pragma omp atomic
      normData[i] += dNorm;
    }
    return 2.0 * r / (err * err);
  };

  bool ok = AccumulateEGammaGradient(lc, ld, configure(), pmap, sdp, fb, accum);
  if(ok) {
    accum.Scatter(pmap, gradFull);
    for(int s = 1; s <= ld->NumSegments(); s++) {
      ESegment* seg = ld->GetSegment(s);
      if(!seg || !seg->IsVaryNorm()) continue;
      int idx = pmap.NormIndex(s);
      if(idx >= 0 && idx < (int)gradFull.size()) gradFull[idx] = normData[s];
    }
    chi2Out = chi2;
  }
  return ok;
}

vector_r AZUREAPI::CalculateChi2GradRWA(const vector_r& params) const {
  vector_r full = MapPackedToFull(params, all_rwa_, fixed_);

  // The analytic gradient pass runs the full forward model at every point, so
  // chi2 comes out as a byproduct -- no separate forward pass. Only fall back to
  // a standalone chi2 evaluation if the analytic path bails.
  vector_r gradFull(all_rwa_.size(), 0.0);
  double chi2 = 0.0;
  bool eg = Chi2GradEGammaNorm(full, gradFull, chi2);
  if(!eg) chi2 = CalculateChi2RWA(params);

  // Finite differences for energy shifts (and the whole block if the analytic
  // path bailed), using the scalar chi-squared.
  ParamIndexMap pmap = BuildParamIndexMap(compound(), data(), fixed_);
  for(int f = 0; f < (int)full.size() && f < pmap.NumFull(); f++) {
    if(fixed_[f]) continue;
    ParamKind kind = pmap.Desc(f).kind;
    if(eg && (kind == ParamKind::LevelEnergy || kind == ParamKind::Gamma ||
              kind == ParamKind::Norm)) continue;
    int packed = pmap.FullToPacked(f);
    if(packed < 0 || packed >= (int)params.size()) continue;
    double x0 = params[packed];
    double h = 1.0e-6 * (std::fabs(x0) + 1.0);
    vector_r pp = params; pp[packed] = x0 + h;
    vector_r pm = params; pm[packed] = x0 - h;
    gradFull[f] = (CalculateChi2RWA(pp) - CalculateChi2RWA(pm)) / (2.0 * h);
  }

  vector_r out;
  out.reserve(1 + params.size());
  out.push_back(chi2);
  for(int i = 0; i < (int)all_rwa_.size(); i++)
    if(!fixed_[i]) out.push_back(gradFull[i]);
  return out;
}

vector_r AZUREAPI::CalculateResidualJacobianRWA(const vector_r& params) const {
  vector_r full = MapPackedToFull(params, all_rwa_, fixed_);

  // One request at a time: operate on the canonical compound/data in place
  // (re-filled here from these parameters) instead of cloning.
  CNuc* lc = compound();
  EData* ld = data();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r* sdp = nullptr;
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) {
    lc->CalcShiftFunctions(configure());
    shiftDeriv = BuildShiftDerivTable(lc, configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed_);
  vector_r residuals, jacobian;
  int nCols = 0;
  bool ok = ComputeResidualJacobian(lc, ld, configure(), pmap, sdp, residuals, jacobian, nCols);

  if(!ok) return vector_r{ -1.0 };

  vector_r out;
  out.reserve(2 + residuals.size() + jacobian.size());
  out.push_back((double)residuals.size());
  out.push_back((double)nCols);
  out.insert(out.end(), residuals.begin(), residuals.end());
  out.insert(out.end(), jacobian.begin(), jacobian.end());
  return out;
}

vector_r AZUREAPI::CalculateModelGradientsRWA(const vector_r& params) const {
  vector_r full = MapPackedToFull(params, all_rwa_, fixed_);

  // Same preparation as CalculateResidualJacobianRWA: one request at a time, so
  // the canonical compound/data are re-filled in place from these parameters.
  CNuc* lc = compound();
  EData* ld = data();
  lc->FillCompoundFromParams(full);
  ld->FillNormsFromParams(full);
  ld->FillEnergyShiftsFromParams(full, ld, lc, &configure());

  vector_matrix_r shiftDeriv;
  const vector_matrix_r* sdp = nullptr;
  if(configure().paramMask & Config::USE_BRUNE_FORMALISM) {
    lc->CalcShiftFunctions(configure());
    shiftDeriv = BuildShiftDerivTable(lc, configure());
    sdp = &shiftDeriv;
  }

  ParamIndexMap pmap = BuildParamIndexMap(lc, ld, fixed_);

  // A band is sensitive only to the R-matrix parameters, and covariance.dat
  // spans exactly those columns -- so reduce each full packed row to them here
  // rather than shipping zero columns for every normalization.
  const std::vector<int> rc = RMatrixPackedColumns(pmap);
  const int nCols = (int)rc.size();

  std::map<EPoint*, vector_r> grad;
  if(!ComputeModelGradients(lc, ld, configure(), pmap, sdp, grad))
    return vector_r{ -1.0 };

  // Walk the segments exactly as UpdateSegments does, so row k of segment s
  // lines up with GET_CALCULATED_SEGMENT's point k: segments sharing a segment
  // key collapse to one calculated segment, and only the first of them is used.
  std::vector<int> counts;
  vector_r rows;
  int prevKey = -1;
  std::vector<ESegment>& segments = ld->GetSegments();
  for( int i = 0; i < (int)segments.size( ); ++i ) {
    const int newKey = segments[i].GetSegmentKey( );
    if( prevKey == newKey ) continue;
    prevKey = newKey;

    std::vector<EPoint>& points = segments[i].GetPoints( );
    int n = 0;
    for( int k = 0; k < (int)points.size( ); ++k ) {
      std::map<EPoint*, vector_r>::const_iterator it = grad.find( &points[k] );
      if( it == grad.end( ) ) continue;
      const vector_r& f = it->second;
      for( int c = 0; c < nCols; ++c )
        rows.push_back( rc[c] < (int)f.size( ) ? f[rc[c]] : 0.0 );
      ++n;
    }
    counts.push_back( n );
  }

  vector_r out;
  out.reserve( 2 + counts.size( ) + rows.size( ) );
  out.push_back( (double)counts.size( ) );
  out.push_back( (double)nCols );
  for( int i = 0; i < (int)counts.size( ); ++i ) out.push_back( (double)counts[i] );
  out.insert( out.end( ), rows.begin( ), rows.end( ) );
  return out;
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

// Structured metadata for every parameter.  The parameters are walked in the
// exact order CNuc::FillMnParams (energies + widths) then EData::FillMnParams
// (norms then energy shifts) emit them, so the records line up one-to-one with
// names_ / all_ / fixed_.  See AZUREAPI.h for the field layout.
vector_r AZUREAPI::GetParameterInfo( ) const {

  vector_r info;

  // Running global parameter index; reads the fixed flag and physical value
  // that UpdateParameters() already cached.
  int gi = 0;
  auto push = [&]( double type, double jgroup, double J, double parity,
                   double level, double levelE, double channel, double L,
                   double S, double pair, double radtype, double segKey,
                   double wignerLimit ) {
    info.push_back( type );
    info.push_back( jgroup );
    info.push_back( J );
    info.push_back( parity );
    info.push_back( level );
    info.push_back( levelE );
    info.push_back( channel );
    info.push_back( L );
    info.push_back( S );
    info.push_back( pair );
    info.push_back( radtype );
    info.push_back( ( gi < (int)fixed_.size() && fixed_[gi] ) ? 1.0 : 0.0 );
    info.push_back( gi < (int)all_.size() ? all_[gi] : 0.0 );
    info.push_back( segKey );
    info.push_back( wignerLimit );
    ++gi;
  };

  // R-matrix parameters: one energy then one width per channel, per level.
  CNuc* nuc = compound();
  for( int j = 1; j <= nuc->NumJGroups(); ++j ) {
    JGroup* jg = nuc->GetJGroup( j );
    double J = jg->GetJ();
    double parity = jg->GetPi();
    for( int la = 1; la <= jg->NumLevels(); ++la ) {
      ALevel* level = jg->GetLevel( la );
      double levelE = level->GetE();
      // energy parameter
      push( 0, j, J, parity, la, levelE, -1, -1, -1, -1, -1, -1, -1 );
      // width parameters (one per channel)
      for( int ch = 1; ch <= jg->NumChannels(); ++ch ) {
        AChannel* chan = jg->GetChannel( ch );
        push( 1, j, J, parity, la, levelE, ch, chan->GetL(), chan->GetS(),
              chan->GetPairNum(), (double)chan->GetRadType(), -1,
              chan->GetWignerLimit() );
      }
    }
  }

  // Normalization parameters: one per segment with IsVaryNorm().
  std::vector<ESegment>& segments = data()->GetSegments();
  for( size_t s = 0; s < segments.size(); ++s ) {
    if( segments[s].IsVaryNorm() )
      push( 2, -1, -1, 0, -1, 0, -1, -1, -1, -1, -1,
            segments[s].GetSegmentKey(), -1 );
  }

  // Energy-shift parameters: one per segment (always emitted).
  for( size_t s = 0; s < segments.size(); ++s ) {
    push( 3, -1, -1, 0, -1, 0, -1, -1, -1, -1, -1,
          segments[s].GetSegmentKey(), -1 );
  }

  return info;
}

vector_r AZUREAPI::GetPairsInfo( ) const {

  vector_r info;
  CNuc* nuc = compound();

  // One record per pair, in 1-based pair-number order so field "pair" of
  // GetParameterInfo() indexes directly into this list.
  for( int p = 1; p <= nuc->NumPairs(); ++p ) {
    PPair* pair = nuc->GetPair( p );
    info.push_back( p );
    info.push_back( pair->GetPairKey() );
    info.push_back( pair->GetPType() );
    info.push_back( pair->IsEntrance() ? 1.0 : 0.0 );
    info.push_back( pair->GetJ( 1 ) );
    info.push_back( pair->GetPi( 1 ) );
    info.push_back( pair->GetZ( 1 ) );
    info.push_back( pair->GetM( 1 ) );
    info.push_back( pair->GetJ( 2 ) );
    info.push_back( pair->GetPi( 2 ) );
    info.push_back( pair->GetZ( 2 ) );
    info.push_back( pair->GetM( 2 ) );
    info.push_back( pair->GetSepE() );
    info.push_back( pair->GetExE() );
    info.push_back( pair->GetChRad() );
    info.push_back( pair->GetI1I2Factor() );
  }

  return info;
}
// ---------------------------------------------------------------------------
//  Diagnostics: the external region, and the caches that make it affordable
// ---------------------------------------------------------------------------

/*!
 * Coulomb wave functions, penetrability, shift function and hard-sphere phase
 * on a requested energy grid.
 *
 * These are ordinary R-matrix quantities that AZURE2 has always computed
 * internally; what is new is being able to ask for them.  The values follow the
 * run's own configuration, so the same call returns the accurate Coulomb
 * routine's answer, GSL's, or the Numerov solution through a nuclear potential,
 * according to how the calculation was set up -- which is how one sees what the
 * hybrid model actually does to the external region.
 */
vector_r AZUREAPI::GetCoulombFunctions(const vector_r& request) const {
  vector_r out;
  if(request.size() < 4) return out;

  int pairKey = static_cast<int>(std::lround(request[0]));
  int lValue  = static_cast<int>(std::lround(request[1]));
  double radius = request[2];
  int nE = static_cast<int>(std::lround(request[3]));
  if(nE < 0 || request.size() < static_cast<size_t>(4 + nE)) return out;

  if(!compound_ || !compound_->IsPairKey(pairKey)) return out;
  PPair* pair = compound_->GetPair(compound_->GetPairNumFromKey(pairKey));
  if(!pair) return out;
  if(radius <= 0.0) radius = pair->GetChRad();   // the channel radius by default

  CoulFunc coul(pair, !!(configure().paramMask & Config::USE_GSL_COULOMB_FUNC));

  out.push_back(static_cast<double>(nE));
  for(int i = 0; i < nE; ++i) {
    double energy = request[4 + i];
    double F = 0., dF = 0., G = 0., dG = 0., P = 0., S = 0., delta = 0.;
    if(energy > 0.0) {
      try {
        CoulWaves w = coul(lValue, radius, energy);
        F = w.F; dF = w.dF; G = w.G; dG = w.dG;
        P = coul.Penetrability(lValue, radius, energy);
        S = coul.PEShift(lValue, radius, energy);
        delta = -std::atan2(w.F, w.G);
      } catch(...) {
        // A failed evaluation returns zeros for that energy rather than
        // aborting the whole grid; the caller can see which points are missing.
        F = dF = G = dG = P = S = delta = 0.;
      }
    }
    out.push_back(F);   out.push_back(dF);
    out.push_back(G);   out.push_back(dG);
    out.push_back(P);   out.push_back(S);
    out.push_back(delta);
  }
  return out;
}

/*!
 * External-capture radial integrals on a requested energy grid.
 *
 * Walks exactly the pathway structure EPoint::CalculateECAmplitudes walks --
 * every EC level whose final pair matches, every KGroup of the entrance pair's
 * decay, every ECMGroup within it -- and evaluates the integral at each
 * requested energy instead of at a data point's energy.  The quantum numbers of
 * each pathway come back with it, so no bookkeeping is needed on the caller's
 * side.
 */
vector_r AZUREAPI::GetECIntegrals(const vector_r& request) const {
  vector_r out;
  if(request.size() < 2) return out;

  int pairKey = static_cast<int>(std::lround(request[0]));
  int nE = static_cast<int>(std::lround(request[1]));
  if(nE <= 0 || request.size() < static_cast<size_t>(2 + nE)) return out;
  std::vector<double> energies(request.begin() + 2, request.begin() + 2 + nE);

  if(!compound_ || !compound_->IsPairKey(pairKey)) return out;
  int aa = compound_->GetPairNumFromKey(pairKey);
  PPair* entrancePair = compound_->GetPair(aa);
  if(!entrancePair || !entrancePair->IsEntrance()) return out;

  // One block per pathway: 6 descriptors then 2*nE numbers.
  std::vector<vector_r> blocks;

  for(int j = 1; j <= compound_->NumJGroups(); j++) {
    for(int la = 1; la <= compound_->GetJGroup(j)->NumLevels(); la++) {
      ALevel* ecLevel = compound_->GetJGroup(j)->GetLevel(la);
      if(!ecLevel->IsECLevel()) continue;
      int ir = ecLevel->GetECPairNum();
      if(ir < 1 || ir > compound_->NumPairs()) continue;

      for(int k = 1; k <= entrancePair->GetDecay(ir)->NumKGroups(); k++) {
        KGroup* theKGroup = entrancePair->GetDecay(ir)->GetKGroup(k);
        for(int ecm = 1; ecm <= theKGroup->NumECMGroups(); ecm++) {
          ECMGroup* g = theKGroup->GetECMGroup(ecm);

          AChannel* finalChannel =
            compound_->GetJGroup(j)->GetChannel(g->GetFinalChannel());
          PPair* finalPair = compound_->GetPair(finalChannel->GetPairNum());

          int liValue;
          double siValue;
          if(g->IsChannelCapture()) {
            MGroup* cc = entrancePair->GetDecay(g->GetChanCapDecay())
                           ->GetKGroup(g->GetChanCapKGroup())
                           ->GetMGroup(g->GetChanCapMGroup());
            liValue = compound_->GetJGroup(cc->GetJNum())
                        ->GetChannel(cc->GetChpNum())->GetL();
            siValue = compound_->GetJGroup(cc->GetJNum())
                        ->GetChannel(cc->GetChpNum())->GetS();
          } else {
            liValue = g->GetL();
            siValue = theKGroup->GetS();
          }

          vector_r block;
          block.push_back(static_cast<double>(liValue));
          block.push_back(static_cast<double>(finalChannel->GetL()));
          block.push_back(2.0 * siValue);
          block.push_back(2.0 * finalChannel->GetS());
          block.push_back(static_cast<double>(g->GetMult()));
          block.push_back(g->GetRadType() == 'E' ? 1.0 : 0.0);

          ECIntegral theECIntegral(finalPair, configure());
          double levelEnergy = ecLevel->GetE();
          for(int i = 0; i < nE; ++i) {
            double inEnergy = energies[i] + entrancePair->GetSepE()
                                          + entrancePair->GetExE();
            complex v(0., 0.);
            if(energies[i] > 0.0) {
              try {
                v = theECIntegral(liValue, finalChannel->GetL(),
                                  siValue, finalChannel->GetS(),
                                  g->GetJ(), compound_->GetJGroup(j)->GetJ(),
                                  g->GetMult(), g->GetRadType(),
                                  inEnergy, levelEnergy,
                                  g->IsChannelCapture());
              } catch(...) {
                v = complex(0., 0.);
              }
            }
            block.push_back(real(v));
            block.push_back(imag(v));
          }
          blocks.push_back(std::move(block));
        }
      }
    }
  }

  out.push_back(static_cast<double>(blocks.size()));
  out.push_back(static_cast<double>(nE));
  for(const vector_r& b : blocks) out.insert(out.end(), b.begin(), b.end());
  return out;
}

/*!
 * Coulomb-function cache counters, aggregated over threads.
 */
vector_r AZUREAPI::GetCacheStats( ) const {
  vector_r out(6, 0.0);
  if(!g_coulFuncCache) return out;
  CoulFuncCache::Stats s = g_coulFuncCache->GetStats();
  out[0] = static_cast<double>(s.queries);
  out[1] = static_cast<double>(s.hits);
  out[2] = static_cast<double>(s.entries);
  out[3] = static_cast<double>(s.keys);
  out[4] = static_cast<double>(s.disabledKeys);
  out[5] = static_cast<double>(s.threads);
  return out;
}
