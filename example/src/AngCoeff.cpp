#include "AngCoeff.h"
#include <gsl/gsl_sf_coupling.h>
#include <math.h>

double AngCoeff::ClebGord(double j1, double j2, double j3, double m1, double m2, double m3) {    
  m3=-m3;
  int j1x2=(int)(2*j1);
  int j2x2=(int)(2*j2);
  int j3x2=(int)(2*j3);
  int m1x2=(int)(2*m1);
  int m2x2=(int)(2*m2);
  int m3x2=(int)(2*m3);

  double w3j=gsl_sf_coupling_3j(j1x2,j2x2,j3x2,m1x2,m2x2,m3x2);

  return pow(-1.0,j1-j2-m3)*sqrt(2.0*j3+1.)*w3j;
}

double AngCoeff::Racah(double j1, double j2, double l2, double l1, double j3, double l3) {

  int j1x2=(int)(2*j1);
  int j2x2=(int)(2*j2);
  int j3x2=(int)(2*j3);
  int l1x2=(int)(2*l1);
  int l2x2=(int)(2*l2);
  int l3x2=(int)(2*l3);
  
  double w6j=gsl_sf_coupling_6j(j1x2,j2x2,j3x2,l1x2,l2x2,l3x2);

  return pow(-1.0,j1+j2+l2+l1)*w6j;
}
