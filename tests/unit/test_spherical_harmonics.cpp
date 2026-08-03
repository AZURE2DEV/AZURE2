// Validate Y_l^m against closed forms. Condon-Shortley phase, and the
// Y_l^{-m} = (-1)^m conj(Y_l^m) relation, are the two things that silently
// flip the sign of an analyzing power.
#include <cmath>
#include <cstdio>
#include "AngCoeff.h"
static complex SphericalHarmonic(int l,int m,double t,double p){return AngCoeff::SphericalHarmonic(l,m,t,p);}

static int fails = 0;
static void check(const char* what, complex got, complex want) {
  double e = std::abs(got - want);
  bool ok = e < 1e-10;
  if (!ok) fails++;
  printf("  %-34s got (% .6f,% .6f) want (% .6f,% .6f)  %s\n",
         what, got.real(), got.imag(), want.real(), want.imag(), ok?"ok":"FAIL");
}
int main() {
  const double th = 0.7, ph = 0.4, c = std::cos(th), s = std::sin(th);
  const double pi = M_PI;
  printf("Closed-form spherical harmonics at theta=0.7, phi=0.4\n");
  check("Y_0^0", SphericalHarmonic(0,0,th,ph), complex(0.5/std::sqrt(pi),0));
  check("Y_1^0", SphericalHarmonic(1,0,th,ph), complex(0.5*std::sqrt(3/pi)*c,0));
  check("Y_1^1", SphericalHarmonic(1,1,th,ph),
        -0.5*std::sqrt(1.5/pi)*s*complex(std::cos(ph),std::sin(ph)));
  check("Y_1^-1", SphericalHarmonic(1,-1,th,ph),
        0.5*std::sqrt(1.5/pi)*s*complex(std::cos(-ph),std::sin(-ph)));
  check("Y_2^0", SphericalHarmonic(2,0,th,ph),
        complex(0.25*std::sqrt(5/pi)*(3*c*c-1),0));
  check("Y_2^2", SphericalHarmonic(2,2,th,ph),
        0.25*std::sqrt(7.5/pi)*s*s*complex(std::cos(2*ph),std::sin(2*ph)));
  printf("\nSymmetry Y_l^{-m} = (-1)^m conj(Y_l^m)\n");
  for (int l=1;l<=4;l++) for (int m=1;m<=l;m++) {
    complex a = SphericalHarmonic(l,-m,th,ph);
    complex b = std::conj(SphericalHarmonic(l,m,th,ph)) * ((m%2)? -1.0 : 1.0);
    char buf[64]; snprintf(buf,sizeof buf,"l=%d m=%d",l,m);
    check(buf, a, b);
  }
  printf("\nOut-of-range returns zero\n");
  check("Y_1^2 (|m|>l)", SphericalHarmonic(1,2,th,ph), complex(0,0));
  printf("\n%s (%d failures)\n", fails?"FAILED":"ALL PASSED", fails);
  return fails!=0;
}
