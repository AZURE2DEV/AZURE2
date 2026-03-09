#include "ShftFunc.h"
#include <math.h>
#include <gsl/gsl_deriv.h>

double ShftFunc::thisShftFunc(double x, void *p) {
  ShftFunc *shift = (ShftFunc*)p;
  int l= shift->params_.lValue;

  return shift->operator()(l,x);
}

double ShftFunc::theWhitFunc(double x, void * p) {
  Params *params =(Params *)p;
  int l = (params->lValue);
  double bindingenergy = (params->bindingEnergy);
  WhitFunc *whitFunc = (params->whitFunc);
  
  return whitFunc->operator()(l,x,bindingenergy);
}

double ShftFunc::operator()(int l, double energy) {
  params_.bindingEnergy=fabs(energy-totalSepE());
  params_.lValue=l;

  double result;
  double error;
 
  gsl_function F;
  F.function = &theWhitFunc;
  F.params = &params_;
 
  gsl_deriv_central (&F, radius(), 1e-4, &result, &error);

  return radius()*result/theWhitFunc(radius(),&params_);
}

double ShftFunc::EnergyDerivative(int l, double energy) {

  double result;
  double error;
   
  params_.lValue = l;

  gsl_function F;
  F.function = &thisShftFunc;
  F.params = this;

  gsl_deriv_central (&F, energy, 1e-6, &result, &error); 

  return result;
}
