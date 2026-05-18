#include "AngCoeff.h"
#include "CNuc.h"
#include "EPoint.h"
#include "GenMatrixFunc.h"
#include <assert.h>
#include <iostream>

/*!
 * The child classes AMatrixFunc or RMatrixFunc contain functions
 * to calculate the T-Matrix from the fitted R-Matrix parameters.  This function
 * then calculates the cross section from the T-Matrix elements.
 */

void GenMatrixFunc::CalculateCrossSection(EPoint *point) {
  complex sum(0.,0.);
  complex sumE1(0.,0.);
  complex sumE2(0.,0.);
  int aa=compound()->GetPairNumFromKey(point->GetEntranceKey());
  int ir=0;
  while(ir<compound()->GetPair(aa)->NumDecays()) {
    ir++;
    if(compound()->GetPair(aa)->GetDecay(ir)->GetPairNum()==compound()->GetPairNumFromKey(point->GetExitKey())) break;
  }
  Decay *theDecay=compound()->GetPair(aa)->GetDecay(ir);
  if(compound()->GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->GetPType()==10 &&
     (configure().paramMask & Config::USE_RMC_FORMALISM)) {
    int decayNum=0;
    while(decayNum<compound()->GetPair(aa)->NumDecays()) {
      decayNum++;
      if(compound()->GetPair(aa)->GetDecay(decayNum)->GetPairNum()==aa) break;
    }
    for(int k=1;k<=compound()->GetPair(aa)->GetDecay(decayNum)->NumKGroups();k++) {
      for(int m=1;m<=compound()->GetPair(aa)->GetDecay(decayNum)->GetKGroup(k)->NumMGroups();m++) {
	MGroup *theMGroup=compound()->GetPair(aa)->GetDecay(decayNum)->GetKGroup(k)->GetMGroup(m);
	if(theMGroup->GetChNum()==theMGroup->GetChpNum()) {
	  double jValue=compound()->GetJGroup(theMGroup->GetJNum())->GetJ();
	  sum+=2.*point->GetGeometricalFactor()*
	    (2.*jValue+1.)*compound()->GetPair(aa)->GetI1I2Factor()*    
	    imag(this->GetTMatrixElement(k,m,decayNum));
	}
      }
    }
    for(int dp=1;dp<=compound()->GetPair(aa)->NumDecays();dp++) {
      if(compound()->GetPair(compound()->GetPair(aa)->GetDecay(dp)->GetPairNum())->GetPType()==0) {
	for(int k=1;k<=compound()->GetPair(aa)->GetDecay(dp)->NumKGroups();k++) {
	  this->ClearTempTMatrices();
	  for(int m=1;m<=compound()->GetPair(aa)->GetDecay(dp)->GetKGroup(k)->NumMGroups();m++) {
	    MGroup *theMGroup=compound()->GetPair(aa)->GetDecay(dp)->GetKGroup(k)->GetMGroup(m);
	    int lValue=compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChNum())->GetL();
	    int lpValue=compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChpNum())->GetL();
	    double jValue=compound()->GetJGroup(theMGroup->GetJNum())->GetJ();
	    int tempTNum=this->IsTempTMatrix(jValue,lValue,lpValue);
	    if(!tempTNum) {
	      TempTMatrix temptmatrix={jValue,lValue,lpValue,this->GetTMatrixElement(k,m,dp)};
	      this->NewTempTMatrix(temptmatrix);
	    } else this->AddToTempTMatrix(tempTNum,this->GetTMatrixElement(k,m,dp));
	  }
	  for(int temp=1;temp<=this->NumTempTMatrices();temp++) {
	    sum-=point->GetGeometricalFactor()*
	      (2.*this->GetTempTMatrix(temp)->jValue+1.)*
	      compound()->GetPair(aa)->GetI1I2Factor()*
	      (this->GetTempTMatrix(temp)->TMatrix)*conj(this->GetTempTMatrix(temp)->TMatrix);
	  }
	}
      }
    }
    point->SetFitCrossSection(real(sum)/100.);
  } else {
    if(!point->IsPhase()) {
      double angleIntegratedXS=0.;
      double angleIntegratedE1XS=0.;
      double angleIntegratedE2XS=0.;
      if(!point->IsDifferential()) {
	for(int k=1;k<=theDecay->NumKGroups();k++) {
	  // For inelastic particle reactions (aa!=ir, GetPType==0) with UPOS 3-param KGroups,
	  // skip KGroups where sp!=sp2 to avoid over-counting in angle-integrated XS.
	  if(aa!=ir&&compound()->GetPair(compound()->GetPairNumFromKey(point->GetExitKey()))->GetPType()==0) {
	    if(theDecay->GetKGroup(k)->GetSp()!=theDecay->GetKGroup(k)->GetSp2()) continue;
	  }
	  this->ClearTempTMatrices();
          this->ClearTempTMatricesE1();
          this->ClearTempTMatricesE2();
	  for(int m=1;m<=theDecay->GetKGroup(k)->NumMGroups();m++) {
            MGroup *theMGroup=theDecay->GetKGroup(k)->GetMGroup(m);
	    if(compound()->GetPair(aa)->GetPType()==20) {
	      sum+=25.*this->GetTMatrixElement(k,m)*conj(this->GetTMatrixElement(k,m));
            } else { 
	      int lValue=compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChNum())->GetL();
	      int lpValue=compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChpNum())->GetL();
	      double jValue=compound()->GetJGroup(theMGroup->GetJNum())->GetJ();
	      int tempTNum=this->IsTempTMatrix(jValue,lValue,lpValue);
	      if(!tempTNum) {
		TempTMatrix temptmatrix={jValue,lValue,lpValue,this->GetTMatrixElement(k,m)};
		this->NewTempTMatrix(temptmatrix);
	      } else this->AddToTempTMatrix(tempTNum,this->GetTMatrixElement(k,m));
              if (compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChpNum())->GetRadType()=='E' && 
              compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChpNum())->GetL()==1) {
	        int tempTNumE1=this->IsTempTMatrixE1(jValue,lValue,lpValue);
	        if(!tempTNumE1) {
		  TempTMatrix temptmatrixE1={jValue,lValue,lpValue,this->GetTMatrixElement(k,m)};
		  this->NewTempTMatrixE1(temptmatrixE1);
	        } else this->AddToTempTMatrixE1(tempTNumE1,this->GetTMatrixElement(k,m));
              } 
              if (compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChpNum())->GetRadType()=='E' && 
              compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChpNum())->GetL()==2) {
	        int tempTNumE2=this->IsTempTMatrixE2(jValue,lValue,lpValue);
	        if(!tempTNumE2) {
		  TempTMatrix temptmatrixE2={jValue,lValue,lpValue,this->GetTMatrixElement(k,m)};
		  this->NewTempTMatrixE2(temptmatrixE2);
	        } else this->AddToTempTMatrixE2(tempTNumE2,this->GetTMatrixElement(k,m)); 
	      }
            }
	  }
	  if(compound()->GetPair(aa)->GetPType()==20) continue;
	  for(int m=1;m<=theDecay->GetKGroup(k)->NumECMGroups();m++) {
            ECMGroup *theECMGroup=theDecay->GetKGroup(k)->GetECMGroup(m);
            int lValue=theECMGroup->GetL();
	    int lpValue=theECMGroup->GetMult();
	    double jValue=theECMGroup->GetJ();
            int tempTNum=this->IsTempTMatrix(jValue,lValue,lpValue);
	    if(!tempTNum) {
	      TempTMatrix temptmatrix={jValue,lValue,lpValue,this->GetECTMatrixElement(k,m)};
	      this->NewTempTMatrix(temptmatrix);
	    } else this->AddToTempTMatrix(tempTNum,this->GetECTMatrixElement(k,m));
            if (theECMGroup->GetRadType()=='E' && theECMGroup->GetMult()==1) {
	      int tempTNumE1=this->IsTempTMatrixE1(jValue,lValue,lpValue);
	      if(!tempTNumE1) {
	        TempTMatrix temptmatrixE1={jValue,lValue,lpValue,this->GetECTMatrixElement(k,m)};
	        this->NewTempTMatrixE1(temptmatrixE1);
	      } else this->AddToTempTMatrixE1(tempTNumE1,this->GetECTMatrixElement(k,m));
            }
            if (theECMGroup->GetRadType()=='E' && theECMGroup->GetMult()==2) {
	      int tempTNumE2=this->IsTempTMatrixE2(jValue,lValue,lpValue);
	      if(!tempTNumE2) {
	        TempTMatrix temptmatrixE2={jValue,lValue,lpValue,this->GetECTMatrixElement(k,m)};
	        this->NewTempTMatrixE2(temptmatrixE2);
	      } else this->AddToTempTMatrixE2(tempTNumE2,this->GetECTMatrixElement(k,m));     
            }
	  }
	  for(int temp=1;temp<=this->NumTempTMatrices();temp++) {
	    sum+=point->GetGeometricalFactor()*
	      (2.*this->GetTempTMatrix(temp)->jValue+1.)*
	      compound()->GetPair(aa)->GetI1I2Factor()*
	      (this->GetTempTMatrix(temp)->TMatrix)*conj(this->GetTempTMatrix(temp)->TMatrix);
	  }
          for(int temp=1;temp<=this->NumTempTMatricesE1();temp++) {
	    sumE1+=point->GetGeometricalFactor()*
	      (2.*this->GetTempTMatrixE1(temp)->jValue+1.)*
	      compound()->GetPair(aa)->GetI1I2Factor()*
	      (this->GetTempTMatrixE1(temp)->TMatrix)*conj(this->GetTempTMatrixE1(temp)->TMatrix);
	  }
          for(int temp=1;temp<=this->NumTempTMatricesE2();temp++) {
	    sumE2+=point->GetGeometricalFactor()*
	      (2.*this->GetTempTMatrixE2(temp)->jValue+1.)*
	      compound()->GetPair(aa)->GetI1I2Factor()*
	      (this->GetTempTMatrixE2(temp)->TMatrix)*conj(this->GetTempTMatrixE2(temp)->TMatrix);
	  }
	}
	angleIntegratedXS=real(sum)/100.;
        angleIntegratedE1XS=real(sumE1)/100.;
        angleIntegratedE2XS=real(sumE2)/100.;
	if(!point->IsAngularDist()) {
	  point->SetFitCrossSection(angleIntegratedXS);
      point->SetFitE1CrossSection(angleIntegratedE1XS);
      point->SetFitE2CrossSection(angleIntegratedE2XS);
	  return;
	}
      }   
      std::vector<double> angularCoeff(std::min(point->GetMaxLOrder()+1,point->GetMaxAngDistOrder()+1),0.);
      PPair *exitPairDiff=compound()->GetPair(compound()->GetPairNumFromKey(point->GetExitKey()));
      for(int kL=1;kL<=theDecay->NumKLGroups();kL++) {
	for(int inter=1;inter<=theDecay->GetKLGroup(kL)
	      ->NumInterferences();inter++) {
	  Interference *theInterference=theDecay->GetKLGroup(kL)
	    ->GetInterference(inter);
	  complex T1(0.0,0.0),T2(0.0,0.0);
	  std::string interferenceType=theInterference->GetInterferenceType();
	  if(aa!=ir&&exitPairDiff->GetPType()==0) {
	    // Different entrance and exit pairs with particle exit - may be UPOS
	    double sp1=theDecay->GetKGroup(theDecay->GetKLGroup(kL)->GetK())->GetSp();
	    double sp2=theDecay->GetKGroup(theDecay->GetKLGroup(kL)->GetK())->GetSp2();
	    MGroup *theMGroup1=theDecay->GetKGroup(theDecay->GetKLGroup(kL)->GetK())->GetMGroup(theInterference->GetM1());
	    int lp1=compound()->GetJGroup(theMGroup1->GetJNum())->GetChannel(theMGroup1->GetChpNum())->GetL();
	    MGroup *theMGroup2=theDecay->GetKGroup(theDecay->GetKLGroup(kL)->GetK())->GetMGroup(theInterference->GetM2());
	    int lp2=compound()->GetJGroup(theMGroup2->GetJNum())->GetChannel(theMGroup2->GetChpNum())->GetL();
	    if(sp1==sp2&&!point->IsUPOS()) {
	      // Normal angular distribution (filter by sp1==sp2)
	      if(interferenceType=="RR") {
		T1=this->GetTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM1());
		T2=this->GetTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM2());
	      } else if(interferenceType=="ER") {
		T1=this->GetECTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM1());
		T2=this->GetTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM2());
	      } else if(interferenceType=="RE") {
		T1=this->GetTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM1());
		T2=this->GetECTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM2());
	      } else if(interferenceType=="EE") {
		T1=this->GetECTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM1());
		T2=this->GetECTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM2());
	      }
	      int lOrder=theDecay->GetKLGroup(kL)->GetLOrder();
	      sum+=theInterference->GetZ1Z2()*T1*conj(T2)*
		point->GetLegendreP(lOrder);
	      if((lOrder < (int)angularCoeff.size()) && point->IsAngularDist()) {
		double tempCoeff=angularCoeff[lOrder]+
		  real(theInterference->GetZ1Z2()*T1*conj(T2))*point->GetGeometricalFactor()*
		  compound()->GetPair(aa)->GetI1I2Factor()/100.*4./angleIntegratedXS;
		angularCoeff[lOrder]=tempCoeff;
	      }
	    }
	    if(lp1==lp2&&point->IsUPOS()) {
	      // Unobserved primary, observed secondary: only RR interference
	      if(interferenceType=="RR") {
		T1=this->GetTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM1());
		T2=this->GetTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM2());
	      }
	      int lOrder=theDecay->GetKLGroup(kL)->GetLOrder();
	      double finalL=(double)point->GetSecondaryDecayL();
	      double Ic=point->GetIc();
	      double j2f=compound()->GetPair(ir)->GetJ(2);
	      double delta=point->GetDelta();
	      double R_L=0.;
	      if((int)lOrder%2==0) {
		R_L=this->GetRk(j2f,finalL,finalL,Ic,lOrder);
		if(delta!=0.) {
		  double finalLp=finalL+1.;
		  double R_LLp=this->GetRk(j2f,finalL,finalLp,Ic,lOrder);
		  double R_LpLp=this->GetRk(j2f,finalLp,finalLp,Ic,lOrder);
		  if(R_LpLp!=0.) R_L=(R_L+2.*delta*R_LLp+delta*delta*R_LpLp)/(1.+delta*delta);
		}
	      }
	      sum+=theInterference->GetZ1Z2_UPOS()*T1*conj(T2)*
		pow(2.*lOrder+1.,0.5)/(4.)*R_L*point->GetLegendreP(lOrder);
	      if((lOrder < (int)angularCoeff.size()) && point->IsAngularDist()) {
		double tempCoeff=angularCoeff[lOrder]+
		  real(theInterference->GetZ1Z2_UPOS()*T1*conj(T2))*point->GetGeometricalFactor()*
		  compound()->GetPair(aa)->GetI1I2Factor()/100.*pow(2.*lOrder+1.,0.5)*R_L/angleIntegratedXS;
		angularCoeff[lOrder]=tempCoeff;
	      }
	    }
	  } else {
	    // Standard case: entrance==exit or gamma exit
	    if(interferenceType=="RR") {
	      T1=this->GetTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM1());
	      T2=this->GetTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM2());
	    } else if(interferenceType=="ER") {
	      T1=this->GetECTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM1());
	      T2=this->GetTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM2());
	    } else if(interferenceType=="RE") {
	      T1=this->GetTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM1());
	      T2=this->GetECTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM2());
	    } else if(interferenceType=="EE") {
	      T1=this->GetECTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM1());
	      T2=this->GetECTMatrixElement(theDecay->GetKLGroup(kL)->GetK(),theInterference->GetM2());
	    }
	    int lOrder=theDecay->GetKLGroup(kL)->GetLOrder();
	    sum+=theInterference->GetZ1Z2()*T1*conj(T2)*
	      point->GetLegendreP(lOrder);
	    if((lOrder < (int)angularCoeff.size()) && point->IsAngularDist()) {
	      double tempCoeff=angularCoeff[lOrder]+
		real(theInterference->GetZ1Z2()*T1*conj(T2))*point->GetGeometricalFactor()*
		compound()->GetPair(aa)->GetI1I2Factor()/100.*4./angleIntegratedXS;
	      angularCoeff[lOrder]=tempCoeff;
	    }
	  }
	}
      }
      if(point->IsAngularDist()) {
	point->SetAngularDists(angularCoeff);
	return;
      }
      // Identical-particle symmetrization (elastic, aa==ir):
      //   |F|^2 = |F_C|^2 + |F_N|^2 + 2 Re[F_C* F_N]
      // with F_X = f_X(theta) + eps f_X(pi-theta). For two identical 0+
      // bosons the only contributing partial waves have L even, so
      // f_N(pi-theta) = f_N(theta), giving |F_N|^2 = 4|f_N|^2 and
      // 2 Re[F_C* F_N] = 4 Re[F_C* f_N]. The Coulomb piece is already
      // built as F_C in EPoint::CalcCoulombAmplitude. Here we apply the
      // remaining factors 4 to RT and 2 to IT.
      double rtFactor = 1.0;
      double itFactor = 1.0;
      if(aa==ir && compound()->GetPair(aa)->IsIdentical()) {
	rtFactor = 4.0;
	itFactor = 2.0;
      }
      complex RT=sum/pi*point->GetGeometricalFactor()*
	compound()->GetPair(aa)->GetI1I2Factor()*rtFactor;

      complex CT(0.,0.), IT(0.,0.);
      if(aa==ir) {
	complex coulombAmplitude=point->GetCoulombAmplitude();
	CT=coulombAmplitude*conj(coulombAmplitude)*point->GetGeometricalFactor();

	sum=complex(0.,0.);
	for(int k=1;k<=theDecay->NumKGroups();k++) {
	  for(int m=1;m<=theDecay->GetKGroup(k)->NumMGroups();m++) {
	    MGroup *theMGroup=theDecay->GetKGroup(k)->GetMGroup(m);
	    AChannel *entranceChannel=compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChNum());
	    AChannel *exitChannel=compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChpNum());
	    if(entranceChannel==exitChannel)
	      sum+=theMGroup->GetStatSpinFactor()*
		coulombAmplitude*conj(this->GetTMatrixElement(k,m))*
		point->GetLegendreP(compound()->GetJGroup(theMGroup->GetJNum())->
				    GetChannel(theMGroup->GetChNum())->GetL());
	  }
	}
	IT=complex(0.,1.)/sqrt(pi)*sum*point->GetGeometricalFactor()*itFactor;
      }
      point->SetFitCrossSection((real(CT)+real(RT)+real(IT))/100.);
    } else if(aa==ir) {
      double segmentJ=point->GetJ();
      int segmentL=point->GetL();
      PPair *entrancePair=compound()->GetPair(aa);
      this->ClearTempTMatrices();
      for(int k=1;k<=theDecay->NumKGroups();k++) {
	for(int m=1;m<=theDecay->GetKGroup(k)->NumMGroups();m++) {
	  MGroup *theMGroup=theDecay->GetKGroup(k)->GetMGroup(m);
	  double jValue=compound()->GetJGroup(theMGroup->GetJNum())->GetJ();
	  AChannel *entranceChannel=compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChNum());
	  int lValue=entranceChannel->GetL();
	  AChannel *exitChannel=compound()->GetJGroup(theMGroup->GetJNum())->GetChannel(theMGroup->GetChpNum());
	  if(jValue==segmentJ&&lValue==segmentL&&entranceChannel==exitChannel) {
	    complex expCoulPhaseSquared=point->GetExpCoulombPhase(theMGroup->GetJNum(),theMGroup->GetChNum())*
	      point->GetExpCoulombPhase(theMGroup->GetJNum(),theMGroup->GetChNum());
	    complex theUMatrix=(expCoulPhaseSquared-this->GetTMatrixElement(k,m))/expCoulPhaseSquared;
	    int tempTNum=this->IsTempTMatrix(jValue,lValue,lValue);
	    if(!tempTNum) {
	      TempTMatrix temptmatrix={jValue,lValue,lValue,theUMatrix};
	      this->NewTempTMatrix(temptmatrix);
	    } else this->AddToTempTMatrix(tempTNum,theUMatrix);
	  }
	}
      }
      assert(this->NumTempTMatrices()<=1);
      // The U-matrix accumulated above is S_L = exp(2 i delta_L); the nuclear
      // phase shift delta_L is half its argument. For identical particles the
      // physical partial-wave amplitude is symmetrised by [1 + eps (-1)^(L+S)],
      // which doubles allowed waves and removes forbidden ones, but delta_L
      // per allowed wave has the same meaning as in the distinguishable case.
      double phase=0.0;
      if(this->NumTempTMatrices()==1) {
	phase = 180.0/pi/2.0*
	  atan2(imag(this->GetTempTMatrix(1)->TMatrix),
		real(this->GetTempTMatrix(1)->TMatrix));
	// atan2 returns (-pi, pi], so phase lives in (-90, 90]. For
	// identical-particle pairs we wrap to [0, 180), the standard
	// convention in partial-wave phase-shift analyses (e^{2 i delta}
	// is pi-periodic, and delta_L grows monotonically through a
	// resonance in this range). Left untouched for distinguishable
	// pairs to preserve existing fit conventions.
	if(entrancePair->IsIdentical() && phase<0.0) phase+=180.0;
      }
      point->SetFitCrossSection(phase);
    }
  }
  
}


/*!
 * Creates a new temporary T-Matrix element.
 */

void GenMatrixFunc::NewTempTMatrix(TempTMatrix tempTMatrix) {
  temp_t_matrices_.push_back(tempTMatrix);
}

void GenMatrixFunc::NewTempTMatrixE1(TempTMatrix tempTMatrix) {
  temp_t_matrices_E1_.push_back(tempTMatrix);
}

void GenMatrixFunc::NewTempTMatrixE2(TempTMatrix tempTMatrix) {
  temp_t_matrices_E2_.push_back(tempTMatrix);
}

/*!
 * Adds a value to the temporary T-Matrix element specified by its position in the TempTMatrix vector.
 */

void GenMatrixFunc::AddToTempTMatrix(int tempTMatrixNum, complex tempValue) {
  this->GetTempTMatrix(tempTMatrixNum)->TMatrix+=tempValue;
  
}

void GenMatrixFunc::AddToTempTMatrixE1(int tempTMatrixNum, complex tempValue) {
  this->GetTempTMatrixE1(tempTMatrixNum)->TMatrix+=tempValue;
}

void GenMatrixFunc::AddToTempTMatrixE2(int tempTMatrixNum, complex tempValue) {
  this->GetTempTMatrixE2(tempTMatrixNum)->TMatrix+=tempValue;
}

/*!
 * Clears the temporary T-Matrices.
 */

void GenMatrixFunc::ClearTempTMatrices() {
  temp_t_matrices_.clear();
}

void GenMatrixFunc::ClearTempTMatricesE1() {
  temp_t_matrices_E1_.clear();
}

void GenMatrixFunc::ClearTempTMatricesE2() {
  temp_t_matrices_E2_.clear();
}

/*!
 * Adds an internal T-Matrix element to the vector of internal T-matrix elements 
 * corresponding to a specified internal reaction pathway.
 */

 void GenMatrixFunc::AddTMatrixElement(int kGroupNum ,int mGroupNum,complex tMatrixElement, int decayNum) {
  matrix_c d;
  vector_c e;
  while(decayNum>tmatrix_.size()) tmatrix_.push_back(d);
  while(kGroupNum>tmatrix_[decayNum-1].size()) tmatrix_[decayNum-1].push_back(e);
  tmatrix_[decayNum-1][kGroupNum-1].push_back(tMatrixElement);
  assert(kGroupNum==tmatrix_[decayNum-1].size());
  assert(mGroupNum==tmatrix_[decayNum-1][kGroupNum-1].size());
}

/*!
 * Adds an external T-Matrix element to the vector of external T-matrix elements 
 * corresponding to a specified external reaction pathway.
 */

void GenMatrixFunc::AddECTMatrixElement(int kGroupNum ,int mGroupNum,complex tMatrixElement) {
  vector_c d;
  while(kGroupNum>ec_tmatrix_.size()) ec_tmatrix_.push_back(d);
  ec_tmatrix_[kGroupNum-1].push_back(tMatrixElement);
  assert(mGroupNum==ec_tmatrix_[kGroupNum-1].size());
}

/*!
 * Tests if a temporary T-Matrix element already exists for a given \f$ J,l,l' \f$ combination.
 * If the element exists, returns the position in the TempTMatrix vector, otherwise returns 0.
 */

int GenMatrixFunc::IsTempTMatrix(double jValue, int lValue, int lPrimeValue) {
  int d=0;
  bool e=false;
  while(!e&&d<this->NumTempTMatrices()) {
    if(jValue==this->GetTempTMatrix(d+1)->jValue&&
       lValue==this->GetTempTMatrix(d+1)->lValue&&
       lPrimeValue==this->GetTempTMatrix(d+1)->lpValue) e=true;
    d++;
  }
  if(!e) return 0;
  else return d;
}

int GenMatrixFunc::IsTempTMatrixE1(double jValue, int lValue, int lPrimeValue) {
  int d=0;
  bool e=false;
  while(!e&&d<this->NumTempTMatricesE1()) {
    if(jValue==this->GetTempTMatrixE1(d+1)->jValue&&
       lValue==this->GetTempTMatrixE1(d+1)->lValue&&
       lPrimeValue==this->GetTempTMatrixE1(d+1)->lpValue) e=true;
    d++;
  }
  if(!e) return 0;
  else return d;
}

int GenMatrixFunc::IsTempTMatrixE2(double jValue, int lValue, int lPrimeValue) {
  int d=0;
  bool e=false;
  while(!e&&d<this->NumTempTMatricesE2()) {
    if(jValue==this->GetTempTMatrixE2(d+1)->jValue&&
       lValue==this->GetTempTMatrixE2(d+1)->lValue&&
       lPrimeValue==this->GetTempTMatrixE2(d+1)->lpValue) e=true;
    d++;
  }
  if(!e) return 0;
  else return d;
}

/*!
 * Returns the number of temporary T-Matrix elements in the TempTMatrix vector.
 */

int GenMatrixFunc::NumTempTMatrices() const {
  return temp_t_matrices_.size();
}

int GenMatrixFunc::NumTempTMatricesE1() const {
  return temp_t_matrices_E1_.size();
}

int GenMatrixFunc::NumTempTMatricesE2() const {
  return temp_t_matrices_E2_.size();
}

/*!
 * Returns a pointer to the temporary T-Matrix element specified by a position in the TempTMatrix vector.
 */

TempTMatrix *GenMatrixFunc::GetTempTMatrix(int tempTMatrixNum) {
  TempTMatrix *b=&temp_t_matrices_[tempTMatrixNum-1];
  return b;
}

TempTMatrix *GenMatrixFunc::GetTempTMatrixE1(int tempTMatrixNum) {
  TempTMatrix *b=&temp_t_matrices_E1_[tempTMatrixNum-1];
  return b;
}

TempTMatrix *GenMatrixFunc::GetTempTMatrixE2(int tempTMatrixNum) {
  TempTMatrix *b=&temp_t_matrices_E2_[tempTMatrixNum-1];
  return b;
}

/*!
 * Returns the value of the internal T-Matrix element specified by an internal reaction pathway.
 */

 complex GenMatrixFunc::GetTMatrixElement(int kGroupNum, int mGroupNum, int decayNum) const {
  return tmatrix_[decayNum-1][kGroupNum-1][mGroupNum-1];
}

/*!
 * Returns the value of the external T-Matrix element specified by an external reaction pathway.
 */

complex GenMatrixFunc::GetECTMatrixElement(int kGroupNum, int ecMGroupNum) const {
  return ec_tmatrix_[kGroupNum-1][ecMGroupNum-1];
}

/*!
 * Calculates the R_k coefficient for the unobserved primary, observed secondary (UPOS)
 * angular distribution calculation.
 *
 * R_k = sqrt(2j2f+1)*sqrt(2L+1)*sqrt(2L'+1)*(-1)^(j2f-Ic+L-L'+k+1)
 *       * ClebGord(L',L,k,1,-1,0) * Racah(L,L',j2f,j2f,k,Ic)
 */

double GenMatrixFunc::GetRk(double j2f, double finalL, double finalLp, double Ic, int lOrder) {
  AngCoeff angCoeff;
  return sqrt(2.*j2f+1.)*sqrt(2.*finalL+1.)*sqrt(2.*finalLp+1.)*
    pow(-1.,j2f-Ic+finalL-finalLp+lOrder+1)*
    angCoeff.ClebGord(finalLp,finalL,lOrder,1.,-1.,0.)*
    angCoeff.Racah(finalL,finalLp,j2f,j2f,lOrder,Ic);
}
