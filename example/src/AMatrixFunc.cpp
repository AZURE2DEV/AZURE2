#include "AMatrixFunc.h"
#include "CNuc.h"
#include "Config.h"
#include "EPoint.h"
#include "MatrixInv.h"
#include <assert.h>

/*!
 * The AMatrixFunc object is created with reference to a CNuc object.
 */

AMatrixFunc::AMatrixFunc(CNuc* compound, const Config &configure) :
  compound_(compound), configure_(configure) {}

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
      for(int la=1;la<=compound()->GetJGroup(j)->NumLevels();la++) {
	if(compound()->GetJGroup(j)->GetLevel(la)->IsInRMatrix()) {
	  ALevel *level=compound()->GetJGroup(j)->GetLevel(la);
	  for(int lap=1;lap<=compound()->GetJGroup(j)->NumLevels();lap++) {
	    if(compound()->GetJGroup(j)->GetLevel(lap)->IsInRMatrix()) {
	      ALevel *levelp=compound()->GetJGroup(j)->GetLevel(lap);
	      complex sum(0.0,0.0);
	      for(int ch=1;ch<=compound()->GetJGroup(j)->NumChannels();ch++) {
		double gammaCh=level->GetFitGamma(ch);
		double gammaChp=levelp->GetFitGamma(ch);
		complex loElement=point->GetLoElement(j,ch);
		sum+=gammaCh*gammaChp*loElement;
		if((compound()->GetJGroup(j)->GetChannel(ch)->GetRadType() == 'M' || 
		    compound()->GetJGroup(j)->GetChannel(ch)->GetRadType() == 'E' ) && 
		   la==lap &&
		   (configure().paramMask & Config::USE_RMC_FORMALISM)) 
		  sum+=complex(0.0,1.0)*gammaCh*gammaChp;
		if((configure().paramMask & Config::USE_BRUNE_FORMALISM) && compound()->GetJGroup(j)->GetChannel(ch)->GetRadType()=='P') {
		  sum+=gammaCh*gammaChp*compound()->GetJGroup(j)->GetChannel(ch)->GetBoundaryCondition();
		  if(la==lap) sum-=gammaCh*gammaChp*level->GetShiftFunction(ch);
		  else sum-=gammaCh*gammaChp*
			 (level->GetShiftFunction(ch)*(inEnergy-levelp->GetFitE())-levelp->GetShiftFunction(ch)*(inEnergy-level->GetFitE()))/
			 (level->GetFitE()-levelp->GetFitE());				
		}
	      }
	      if(la==lap) {
		double resenergy=level->GetFitE();
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
  for(int j=1;j<=compound()->NumJGroups();j++) {
    if(compound()->GetJGroup(j)->IsInRMatrix()) {
      matrix_c *theAInvMatrix = this->GetJSpecAInvMatrix(j);
      MatrixInv matrixInv(*theAInvMatrix);
      this->AddAMatrix(matrixInv.inverse());
    }
  }
}

/*!
 * This function calculates the T-Matrix for each reaction pathway based on the A-Matrix.
 */

void AMatrixFunc::CalculateTMatrix(EPoint *point) {
  int aa=compound()->GetPairNumFromKey(point->GetEntranceKey());
  int irEnd;
  int irStart;
  bool isRMC=false;
  if((configure().paramMask & Config::USE_RMC_FORMALISM) && 
     compound()->GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->GetPType()==10) {
    irStart=1;
    irEnd=compound()->GetPair(aa)->NumDecays();
    isRMC=true;
  } else {
    irStart=0;
    while(irStart<compound()->GetPair(aa)->NumDecays()) {
      irStart++;
      if(compound()->GetPair(aa)->GetDecay(irStart)->GetPairNum()==
	 compound()->GetPairNumFromKey(point->GetExitKey())) break;
    }
    irEnd=irStart;
  }
  for(int ir=irStart;ir<=irEnd;ir++) {
    Decay *theDecay=compound()->GetPair(aa)->GetDecay(ir);
    for(int k=1;k<=theDecay->NumKGroups();k++) {
      for(int m=1;m<=theDecay->GetKGroup(k)->NumMGroups();m++) {
	MGroup *theMGroup=theDecay->GetKGroup(k)->GetMGroup(m);
	JGroup *theJGroup=compound()->GetJGroup(theMGroup->GetJNum());
	AChannel *entranceChannel=theJGroup->GetChannel(theMGroup->GetChNum());
	AChannel *exitChannel=theJGroup->GetChannel(theMGroup->GetChpNum());
	complex uphase=point->GetExpCoulombPhase(theMGroup->GetJNum(),theMGroup->GetChNum())*
	  point->GetExpHardSpherePhase(theMGroup->GetJNum(),theMGroup->GetChNum())*
	  point->GetExpCoulombPhase(theMGroup->GetJNum(),theMGroup->GetChpNum())*
	  point->GetExpHardSpherePhase(theMGroup->GetJNum(),theMGroup->GetChpNum());
	complex umatrix(0.,0.);
	for(int la=1;la<=compound()->GetJGroup(theMGroup->GetJNum())->NumLevels();la++) {
	  if(compound()->GetJGroup(theMGroup->GetJNum())->GetLevel(la)->IsInRMatrix()) {
	    ALevel *level=compound()->GetJGroup(theMGroup->GetJNum())->GetLevel(la);
	    for(int lap=1;lap<=compound()->GetJGroup(theMGroup->GetJNum())->NumLevels();lap++) {
	      if(compound()->GetJGroup(theMGroup->GetJNum())->GetLevel(lap)->IsInRMatrix()) {
		ALevel *levelp=compound()->GetJGroup(theMGroup->GetJNum())->GetLevel(lap);
		umatrix+=2.0*complex(0.0,1.0)*
		  point->GetSqrtPenetrability(theMGroup->GetJNum(),theMGroup->GetChNum())*
		  point->GetSqrtPenetrability(theMGroup->GetJNum(),theMGroup->GetChpNum())*
		  level->GetFitGamma(theMGroup->GetChNum())*
		  levelp->GetFitGamma(theMGroup->GetChpNum())*
		  this->GetAMatrixElement(theMGroup->GetJNum(),la,lap);
	      }
	    }
	  }
	}
	complex tphase=point->GetExpCoulombPhase(theMGroup->GetJNum(),theMGroup->GetChNum())*
	  point->GetExpCoulombPhase(theMGroup->GetJNum(),theMGroup->GetChNum());
	complex tmatrix;
	if(isRMC) this->AddTMatrixElement(k,m,complex(0.0,-1.0)*umatrix,ir);
	else {
	  if(theMGroup->GetChNum()==theMGroup->GetChpNum()) {
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
	complex tmatrix=ecNormParam*point->GetECAmplitude(k,m);
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
  matrix_c e;
  vector_c f;
  while(jGroupNum>a_inv_matrices_.size()) a_inv_matrices_.push_back(e);
  while(lambdaNum>a_inv_matrices_[jGroupNum-1].size()) a_inv_matrices_[jGroupNum-1].push_back(f);
  a_inv_matrices_[jGroupNum-1][lambdaNum-1].push_back(aMatrixElement);
  assert(muNum=a_inv_matrices_[jGroupNum-1][lambdaNum-1].size());
}

/*!
 * This function adds an entire A-Matrix to a vector.
 */

void AMatrixFunc::AddAMatrix(matrix_c aMatrix) {
  a_matrices_.push_back(aMatrix);
}
