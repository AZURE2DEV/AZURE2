#ifndef ECINTEGRAL_H
#define ECINTEGRAL_H

#include "CoulFunc.h"
#include "WhitFunc.h"
#include "Config.h"

class EffectiveCharge;

///A function class to calculate external capture integrals

/*!
 * The ECIntegral function class calculates external capture integrals for 
 * both positive and negative energy channels.  The results are returned as 
 * an ECIntResult structure.
 */

class ECIntegral {
 public:
  /*!
   * The ECIntegral object is created with reference to a PPair object. 
   * The PPair object is used to create new instances of the CoulFunc and 
   * WhitFunc objects.  
   */
  ECIntegral(PPair *pPair, const Config& configure) {
    params_.coulFunc = new CoulFunc(pPair,!!(configure.paramMask&Config::USE_GSL_COULOMB_FUNC));
    params_.whitFunc = new WhitFunc(pPair);
    params_.useLongWavelengthApprox = !!(configure.paramMask&Config::USE_LONGWAVELENGTH_APPROX);
    pair_ = pPair;
  };
  /*!
   * The CoulFunc and WhitFunc objects are destroyed with the object.
   */
  ~ECIntegral() {
    delete params_.coulFunc;
    delete params_.whitFunc;
  };
  complex operator()(int,int,double,double,double,double,int,char,double,double,bool);
 private:
  void ResetIntegrals() {FW_=0.;GW_=0.;};
  void Integrate(double);
  static double FWIntegrand(double,void*);
  static double GWIntegrand(double,void*);
  static double WWIntegrand(double,void*);
  CoulFunc *coulfunction() const {return params_.coulFunc;};
  WhitFunc *whitfunction() const {return params_.whitFunc;};
  PPair *pair() const {return pair_;};
  double FW() const {return FW_;};
  double GW() const {return GW_;};  
  typedef struct Params {
    EffectiveCharge* effectiveCharge;
    CoulFunc *coulFunc;
    WhitFunc *whitFunc;
    int liValue;
    int lfValue;
    int multLValue;
    double pairEnergy;
    double bindingEnergy;
    bool useLongWavelengthApprox;
  } Params;
  Params params_;
  PPair *pair_;
  double FW_;
  double GW_;
};




#endif
