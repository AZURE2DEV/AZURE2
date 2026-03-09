#include "IntegratedFermiFunc.h"
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>

const double IntegratedFermiFunc::alpha_=1./137.036;
const double IntegratedFermiFunc::pi_=3.14159;
const double IntegratedFermiFunc::electronMass_=0.51099891;
const double IntegratedFermiFunc::hbarc_=197.327;

/*!
 * The constructor for the function takes the charge of the Fermion
 * as the first argument. This is expected to be either -1 or 1 for electrons
 * and positrons, respectively.  The second optional argument is the screening
 * potential \f$V_0\f$, where the actual screening potential, \f$V\f$, will be calculated
 * as \f$V = V_0 \alpha^2 Z^{4/3}\f$.
 */

IntegratedFermiFunc::IntegratedFermiFunc(int charge, double V0) :
  charge_(charge), V0_(V0) {
}

/*!
 * The calling operator for the function class takes the end point energy (in
 * units of \f$m_0c^2\f$), the charge of the daughter nucleus, and the nuclear radius
 * parameter as arguments.  The return value is the integrated Fermi function.
 */

double IntegratedFermiFunc::operator()(double W0, 
				       double Z, 
				       double radius) {
  
  if(W0<=1) return 0.;

  double gamma0 = sqrt(1.-alpha_*alpha_*Z*Z);
  double GammaDenom = gsl_sf_gamma(2.*gamma0+1.);
  Params_ params = {charge_,
  		    gamma0,
		    Z,
		    radius,
		    W0,
		    GammaDenom*GammaDenom,
		    V0_};

  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  gsl_function f;
  f.function = &IntegratedFermiFunc::Integrand;
  f.params = &params;
  
  double result, error;
  gsl_integration_qags (&f, 1., params.W0, 0., 1.e-6, 1000, w, &result, &error);
  
  gsl_integration_workspace_free (w);

  return result;
}

/*!
 * This function returns the actual integrand for use in the GSL dynamic integration
 * routines.  This is a static member function.
 */

double IntegratedFermiFunc::Integrand(double x, void* p) {
  Params_* params = (Params_*)p;

  double W = (params->charge<0) ? x-params->V0*alpha_*alpha_*pow(params->Z,4./3.) : x+params->V0*alpha_*alpha_*pow(params->Z,4./3.);
  double eta = (params->charge<0) ? alpha_*params->Z*W/sqrt(W*W-1.) : 
    -alpha_*params->Z*W/sqrt(W*W-1.);
  gsl_sf_result GammaNumValue;
  gsl_sf_result GammaNumArg;
  gsl_sf_lngamma_complex_e(params->gamma0,eta,&GammaNumValue,&GammaNumArg);
  double GammaNum2 = exp(2.*GammaNumValue.val);
  
  double result = 
    2.*(1.+params->gamma0)*
    exp(pi_*eta-2.*(1.-params->gamma0)*
	log(2.*electronMass_*params->radius*sqrt(W*W-1.)/hbarc_))*
    GammaNum2/params->GammaDenom2*sqrt(W*W-1.)*W*(params->W0-W)*(params->W0-W);

  if(params->V0!=0.) {
    result*=sqrt((W*W-1.)/(x*x-1.))*W/x;
  }

  return result;
}
