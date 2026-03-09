#include "EffectiveCharge.h"
#include "PPair.h"
#include "Constants.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

extern double DoubleFactorial(int);

/*!
 * The constructor takes a PPair object, gamma energy, and multipolarity as arguments.
 */

EffectiveCharge::EffectiveCharge(PPair* pair, double energy, int L) :
  z1_(pair->GetZ(1)),z2_(pair->GetZ(2)),L_(L),m1_(pair->GetM(1)),m2_(pair->GetM(2)),
  energy_(energy) {
}

double EffectiveCharge::Integrand(double x, void* p) {
  int* L = (int*)p;
  return gsl_sf_bessel_jl(*L,x)/x;
}

/*!
 * The operator() function takes the radius as an argument, carries out
 * integration, and returns the rho-dependent effective charge.
 */

double EffectiveCharge::operator()(double r) {
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &Integrand;
  F.params= &L_;

  double intZ1=0.;
  double intZ2=0.;
  double intError;
  gsl_integration_qags (&F, 0., m2_/(m1_+m2_)*energy_/hbarc*r, 0., 1.e-6, 1000, w, &intZ1, &intError);
  gsl_integration_qags (&F, 0., m1_/(m1_+m2_)*energy_/hbarc*r, 0., 1.e-6, 1000, w, &intZ2, &intError);

  gsl_integration_workspace_free (w);

  return sqrt(fstruc*hbarc)*L_*DoubleFactorial(2*L_+1)/
    pow((energy_/hbarc*r),(double)L_)*(z1_*intZ1+pow(-1,L_)*z2_*intZ2);
}
