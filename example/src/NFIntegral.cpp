#include "NFIntegral.h"
#include <gsl/gsl_integration.h>

double NFIntegral::Integrand(double x, void * p) {
  Params *params = (Params*)p;
  WhitFunc *whitFunc=(params->whitFunc);
  int lfValue=(params->lfValue);
  double bindingEnergy=(params->bindingEnergy);
  double whitChRadSquaredValue=(params->whitChRadSquaredValue);
  
  double whit=whitFunc->operator()(lfValue,x,bindingEnergy);
  return pow(whit,2.0)/whitChRadSquaredValue;
}

double NFIntegral::operator()(int lf,double levelEnergy) {
  params_.lfValue = lf;
  params_.bindingEnergy = fabs(levelEnergy - totalSepE());
  
  params_.whitChRadSquaredValue = pow(params_.whitFunc->operator()(lf,chanRad(),params_.bindingEnergy),2.0);  
				  
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  
  gsl_function F;
  F.function = &Integrand;
  F.params= &params_;

  double intresult,interror;

  gsl_integration_qagiu(&F,chanRad(),0.0,1e-4,1000,w,&intresult,&interror);
  
  gsl_integration_workspace_free (w);
  
  return intresult;
}
