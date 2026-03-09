#ifndef EFFECTIVECHARGE_H
#define EFFECTIVECHARGE_H

class PPair;

///A function class for calculating effective charge without long-wavelength approximation

/*!
 * The EffectiveCharge class calculates the effective charge needed for external capture
 * without using the long wavelength approximation.  The formalism is based on 
 * J.L. Friar and S. Fallieros, Phys. Rev. C 29, 1645 (1984). 
 */

class EffectiveCharge {
 public:
  EffectiveCharge(PPair*,double, int);
  double operator()(double);
 private:
  static double Integrand(double,void*);
  int z1_;
  int z2_;
  int L_;
  double m1_;
  double m2_;
  double energy_;
};

#endif
