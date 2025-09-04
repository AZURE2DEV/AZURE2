#include "AMatrixFunc.h"
#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "MatrixInv.h"
#include <assert.h>
#include <iostream>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

/*!
 * The AMatrixFunc object is created with reference to a CNuc object.
 */

AMatrixFunc::AMatrixFunc(CNuc* compound, const Config &configure) :
  compound_(compound), configure_(configure), a_matrices_index_(0), 
  cached_max_levels_(0), cached_max_channels_(0) {}

/*!
 * Returns an A-Matrix element specified by positions in the JGroup and ALevel vectors. 
 */

complex AMatrixFunc::GetAMatrixElement(int jGroupNum, int lambdaNum, int muNum) const {
  return a_matrices_[jGroupNum-1][lambdaNum-1][muNum-1];
}

/*!
 * Returns a pointer to an entire A-Matrix specified by a position in the JGroup vector.
 */

matrix_c *AMatrixFunc::GetJSpecAInvMatrix(int jGroupNum) {
  matrix_c *b=&a_inv_matrices_[jGroupNum-1];
  return b;
}

/*!
 * Clears all matrices associated with the AMatrixFunc object.
 */

void AMatrixFunc::ClearMatrices() {
  a_inv_matrices_.clear();
  a_matrices_.clear();
  tmatrix_.clear();
  ec_tmatrix_.clear();

  const int numJGroups = compound()->NumJGroups();
  a_inv_matrices_.resize(numJGroups);
  a_matrices_.resize(numJGroups);
  level_active_index_.assign(numJGroups, {});   // reset maps

  for (int j = 1; j <= numJGroups; ++j) {
    if (!compound()->GetJGroup(j)->IsInRMatrix()) {
      // keep empty placeholders to preserve external indexing
      level_active_index_[j-1].clear();
      continue;
    }

    JGroup* jg = compound()->GetJGroup(j);
    const int numLevels = jg->NumLevels();

    // Build original->active map (1-based for readability in existing code)
    std::vector<int> map(numLevels + 1, 0);
    int act = 0;
    for (int la = 1; la <= numLevels; ++la) {
      if (jg->GetLevel(la)->IsInRMatrix()) map[la] = ++act;
    }
    level_active_index_[j-1] = std::move(map);

    // Pre-size to dense (act x act) with zeros; no push_backs later
    matrix_c M;
    M.resize(act);
    for (int r = 0; r < act; ++r) M[r].assign(act, complex(0.0, 0.0));
    a_inv_matrices_[j-1] = std::move(M);
  }
}


/*!
 * This function creates the inverted A-Matrix from the parameters in the CNuc object.
 */

void AMatrixFunc::FillMatrices (EPoint *point) {
  double inEnergy;
  if(compound()->
     GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->
     GetPType()==20)
    inEnergy=point->GetCMEnergy()+
      compound()->
      GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->
      GetSepE()+
      compound()->
      GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->
      GetExE();
  else inEnergy = point->GetCMEnergy()+
	 compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetSepE()+
	 compound()->GetPair(compound()->GetPairNumFromKey(point->GetEntranceKey()))->GetExE();
  for(int j=1;j<=compound()->NumJGroups();j++) {
    if(compound()->GetJGroup(j)->IsInRMatrix()) {
      // Cache JGroup pointer to avoid repeated calls
      JGroup *jGroup = compound()->GetJGroup(j);
      int numLevels = jGroup->NumLevels();
      int numChannels = jGroup->NumChannels();
      
      // Ensure pre-allocated buffers are large enough - resize only when needed
      if(numLevels > cached_max_levels_ || numChannels > cached_max_channels_) {
        cached_max_levels_ = std::max(numLevels, cached_max_levels_);
        cached_max_channels_ = std::max(numChannels, cached_max_channels_);
        
        levelGammas_.resize(cached_max_levels_ + 1);
        levelEnergies_.resize(cached_max_levels_ + 1);
        shiftFunctions_.resize(cached_max_levels_ + 1);
        
        for(int la = 0; la <= cached_max_levels_; la++) {
          levelGammas_[la].resize(cached_max_channels_ + 1);
          shiftFunctions_[la].resize(cached_max_channels_ + 1);
        }
      }
      
      // Pre-cache gamma values for all levels and channels using pre-allocated buffers
      for(int la=1; la<=numLevels; la++) {
        if(jGroup->GetLevel(la)->IsInRMatrix()) {
          ALevel *level = jGroup->GetLevel(la);
          levelEnergies_[la] = level->GetFitE();
          
          for(int ch=1; ch<=numChannels; ch++) {
            levelGammas_[la][ch] = level->GetFitGamma(ch);
            shiftFunctions_[la][ch] = level->GetShiftFunction(ch);
          }
        }
      }
      
      for(int la=1;la<=numLevels;la++) {
	if(jGroup->GetLevel(la)->IsInRMatrix()) {
	  ALevel *level=jGroup->GetLevel(la);
	  for(int lap=1;lap<=numLevels;lap++) {
	    if(jGroup->GetLevel(lap)->IsInRMatrix()) {
	      ALevel *levelp=jGroup->GetLevel(lap);
	      complex sum(0.0,0.0);
	      for(int ch=1;ch<=numChannels;ch++) {
		double gammaCh=levelGammas_[la][ch];
		double gammaChp=levelGammas_[lap][ch];
		
		// Early termination for effectively zero gamma values
		if(fabs(gammaCh) < 1.0e-12 || fabs(gammaChp) < 1.0e-12) continue;
		
		complex loElement=point->GetLoElement(j,ch);
		sum+=gammaCh*gammaChp*loElement;
		
		// Cache channel pointer and radiation type to avoid repeated calls
		AChannel *channel = jGroup->GetChannel(ch);
		char radType = channel->GetRadType();
		
		if((radType == 'M' || radType == 'E') && 
		   la==lap &&
		   (configure().paramMask & Config::USE_RMC_FORMALISM)) 
		  sum+=complex(0.0,1.0)*gammaCh*gammaChp;
		if((configure().paramMask & Config::USE_BRUNE_FORMALISM) && radType=='P') {
		  sum+=gammaCh*gammaChp*channel->GetBoundaryCondition();
		  if(la==lap) sum-=gammaCh*gammaChp*shiftFunctions_[la][ch];
		  else sum-=gammaCh*gammaChp*
			 (shiftFunctions_[la][ch]*(inEnergy-levelEnergies_[lap])-shiftFunctions_[lap][ch]*(inEnergy-levelEnergies_[la]))/
			 (levelEnergies_[la]-levelEnergies_[lap]);				
		}
	      }
	      if(la==lap) {
		double resenergy=levelEnergies_[la];
		this->AddAInvMatrixElement(j,la,lap,resenergy-inEnergy-sum);
	      } else this->AddAInvMatrixElement(j,la,lap,-sum);
	    }
	  }
	}     
      }
    }
  }
}

/*!
 * This function inverts the inverse A-Matrix to yeild the A-Matrix.
 */

void AMatrixFunc::InvertMatrices() {
  // Sequential processing to avoid race conditions in matrix data access
  // The matrix inversions themselves are computationally intensive but the data
  // structures are not thread-safe for concurrent access
  for(int j=1;j<=compound()->NumJGroups();j++) {
    if(compound()->GetJGroup(j)->IsInRMatrix()) {
      matrix_c *theAInvMatrix = this->GetJSpecAInvMatrix(j);
      // Add validation to catch corrupted matrices before GSL processing
      if(theAInvMatrix->empty()) {
        continue; // Skip empty matrices
      }
      MatrixInv matrixInv(*theAInvMatrix);
      // Use move semantics to avoid matrix copy and assign directly to correct index
      a_matrices_[j-1] = std::move(matrixInv.inverse());
    }
  }
}

/*!
 * This function calculates the T-Matrix for each reaction pathway based on the A-Matrix.
 */

void AMatrixFunc::CalculateTMatrix(EPoint *point) {
  // Cache frequently accessed values to avoid repeated function calls
  int entranceKey = point->GetEntranceKey();
  int exitKey = point->GetExitKey();
  int aa = compound()->GetPairNumFromKey(entranceKey);
  int exitPairNum = compound()->GetPairNumFromKey(exitKey);
  
  int irEnd;
  int irStart;
  bool isRMC=false;
  if((configure().paramMask & Config::USE_RMC_FORMALISM) && 
     compound()->GetPair(exitPairNum)->GetPType()==10) {
    irStart=1;
    irEnd=compound()->GetPair(aa)->NumDecays();
    isRMC=true;
  } else {
    irStart=0;
    int numDecays = compound()->GetPair(aa)->NumDecays();
    while(irStart<numDecays) {
      irStart++;
      if(compound()->GetPair(aa)->GetDecay(irStart)->GetPairNum()==exitPairNum) break;
    }
    irEnd=irStart;
  }
  for(int ir=irStart;ir<=irEnd;ir++) {
    Decay *theDecay=compound()->GetPair(aa)->GetDecay(ir);
    for(int k=1;k<=theDecay->NumKGroups();k++) {
      for(int m=1;m<=theDecay->GetKGroup(k)->NumMGroups();m++) {
	MGroup *theMGroup=theDecay->GetKGroup(k)->GetMGroup(m);
	// Cache frequently accessed values to avoid repeated function calls
	int jNum = theMGroup->GetJNum();
	int chNum = theMGroup->GetChNum();
	int chpNum = theMGroup->GetChpNum();
	
	JGroup *theJGroup=compound()->GetJGroup(jNum);
	AChannel *entranceChannel=theJGroup->GetChannel(chNum);
	AChannel *exitChannel=theJGroup->GetChannel(chpNum);
	
	// Cache phase calculations 
	complex coulombPhaseEn = point->GetExpCoulombPhase(jNum,chNum);
	complex hardSpherePhaseEn = point->GetExpHardSpherePhase(jNum,chNum);
	complex coulombPhaseEx = point->GetExpCoulombPhase(jNum,chpNum);
	complex hardSpherePhaseEx = point->GetExpHardSpherePhase(jNum,chpNum);
	complex sqrtPenEn = point->GetSqrtPenetrability(jNum,chNum);
	complex sqrtPenEx = point->GetSqrtPenetrability(jNum,chpNum);
	
	complex uphase = coulombPhaseEn * hardSpherePhaseEn * coulombPhaseEx * hardSpherePhaseEx;
	complex umatrix(0.,0.);
	
	int numLevels = theJGroup->NumLevels();
	for(int la=1;la<=numLevels;la++) {
	  if(theJGroup->GetLevel(la)->IsInRMatrix()) {
	    ALevel *level=theJGroup->GetLevel(la);
	    double gammaEn = level->GetFitGamma(chNum);
	    
	    // Skip if entrance gamma is effectively zero
	    if(fabs(gammaEn) < 1.0e-12) continue;
	    
	    for(int lap=1;lap<=numLevels;lap++) {
	      if(theJGroup->GetLevel(lap)->IsInRMatrix()) {
		ALevel *levelp=theJGroup->GetLevel(lap);
		double gammaEx = levelp->GetFitGamma(chpNum);
		
		// Skip if exit gamma is effectively zero
		if(fabs(gammaEx) < 1.0e-12) continue;
		
		umatrix+=2.0*complex(0.0,1.0)*
		  sqrtPenEn * sqrtPenEx *
		  gammaEn * gammaEx *
		  this->GetAMatrixElement(jNum,la,lap);
	      }
	    }
	  }
	}
	complex tphase = coulombPhaseEn * coulombPhaseEn;
	complex tmatrix;
	if(isRMC) this->AddTMatrixElement(k,m,complex(0.0,-1.0)*umatrix,ir);
	else {
	  if(chNum==chpNum) {
	    tmatrix=tphase-uphase*(1.0+umatrix);
	  } else tmatrix=-uphase*umatrix;
	  this->AddTMatrixElement(k,m,tmatrix);
	}
      }
      for(int m=1;m<=theDecay->GetKGroup(k)->NumECMGroups();m++) {
	ECMGroup *theECMGroup=theDecay->GetKGroup(k)->GetECMGroup(m);
	ALevel *finalLevel=compound()->GetJGroup(theECMGroup->GetJGroupNum())
	  ->GetLevel(theECMGroup->GetLevelNum());
	double ecNormParam=finalLevel->GetFitGamma(theECMGroup->GetFinalChannel())*
	  finalLevel->GetSqrtNFFactor()*finalLevel->GetECConversionFactor(theECMGroup->GetFinalChannel());
	// Use energy-shift aware EC amplitude calculation if energy shifts are active
	complex ecAmplitude;
	if(configure().paramMask & Config::USE_EXTERNAL_CAPTURE) {
	  // Use the new method that accounts for energy shifts through interpolation
	  ecAmplitude = point->GetECAmplitudeWithShift(k, m, compound(), configure());
	} else {
	  // Use cached amplitude for backward compatibility
	  ecAmplitude = point->GetECAmplitude(k, m);
	}
	complex tmatrix = ecNormParam * ecAmplitude;
	if(theECMGroup->IsChannelCapture()) {
	  int internalChannel=theECMGroup->GetIntChannelNum();
	  MGroup *chanMGroup=compound()->GetPair(aa)->GetDecay(theECMGroup->GetChanCapDecay())
	    ->GetKGroup(theECMGroup->GetChanCapKGroup())->GetMGroup(theECMGroup->GetChanCapMGroup());
	  AChannel *chanEntranceChannel=compound()->GetJGroup(chanMGroup->GetJNum())
	    ->GetChannel(chanMGroup->GetChNum());
	  AChannel *chanExitChannel=compound()->GetJGroup(chanMGroup->GetJNum())
	    ->GetChannel(chanMGroup->GetChpNum());
	  complex umatrix(0.,0.);
	  for(int la=1;la<=compound()->GetJGroup(chanMGroup->GetJNum())->NumLevels();la++) {
	    if(compound()->GetJGroup(chanMGroup->GetJNum())->GetLevel(la)->IsInRMatrix()) {
	      ALevel *level=compound()->GetJGroup(chanMGroup->GetJNum())->GetLevel(la);
	      if(internalChannel && (configure().paramMask & Config::IGNORE_ZERO_WIDTHS))
		if(fabs(level->GetFitGamma(internalChannel))<1.0e-8) continue;
	      for(int lap=1;lap<=compound()->GetJGroup(chanMGroup->GetJNum())->NumLevels();lap++) {
		if(compound()->GetJGroup(chanMGroup->GetJNum())->GetLevel(lap)->IsInRMatrix()) {
		  ALevel *levelp=compound()->GetJGroup(chanMGroup->GetJNum())->GetLevel(lap);
		  if(internalChannel && (configure().paramMask & Config::IGNORE_ZERO_WIDTHS))
		    if(fabs(levelp->GetFitGamma(internalChannel))<1.0e-8) continue;
		  umatrix+=2.0*complex(0.0,1.0)*
		    point->GetSqrtPenetrability(chanMGroup->GetJNum(),chanMGroup->GetChNum())*
		    level->GetFitGamma(chanMGroup->GetChNum())*
		    levelp->GetFitGamma(chanMGroup->GetChpNum())*
		    this->GetAMatrixElement(chanMGroup->GetJNum(),la,lap);
		}
	      }
	    }
	  }
	  tmatrix=tmatrix*umatrix;
	}
	this->AddECTMatrixElement(k,m,tmatrix);
      }
    }
  }
}

/*!
 * This function adds an inverse A-Matrix element specified by positions in the JGroup and ALevel vectors.
 */

void AMatrixFunc::AddAInvMatrixElement(int jGroupNum, int lambdaNum, int muNum, complex aMatrixElement) {
  // Basic bounds
  if (jGroupNum < 1 || jGroupNum > (int)a_inv_matrices_.size()) return;

  // If this J-group wasn’t included, skip
  if (level_active_index_.empty() || level_active_index_[jGroupNum-1].empty()) return;

  // Translate original level indices to compact active indices
  const auto& map = level_active_index_[jGroupNum-1];
  if (lambdaNum < 1 || lambdaNum >= (int)map.size()) return;
  if (muNum     < 1 || muNum     >= (int)map.size()) return;

  const int row = map[lambdaNum];
  const int col = map[muNum];
  if (row == 0 || col == 0) return; // either level is inactive; logic says: ignore

  // Validate value
  if (!std::isfinite(aMatrixElement.real()) || !std::isfinite(aMatrixElement.imag())) return;

  // Direct assignment into pre-sized dense matrix
  a_inv_matrices_[jGroupNum-1][row-1][col-1] = aMatrixElement;
}

/*!
 * This function adds an entire A-Matrix to a vector.
 */

void AMatrixFunc::AddAMatrix(matrix_c aMatrix) {
  // This method is now primarily used for backward compatibility
  // Direct indexing in InvertMatrices() is preferred for thread safety
  if(a_matrices_index_ < a_matrices_.size()) {
    a_matrices_[a_matrices_index_] = aMatrix;
    a_matrices_index_++;
  }
}

void AMatrixFunc::AddAMatrix(matrix_c&& aMatrix) {
  // Move semantics version - avoids copying large matrices
  if(a_matrices_index_ < a_matrices_.size()) {
    a_matrices_[a_matrices_index_] = std::move(aMatrix);
    a_matrices_index_++;
  }
}
