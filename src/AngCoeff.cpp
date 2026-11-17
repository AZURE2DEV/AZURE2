#include "AngCoeff.h"

#include "azure2/physics/AngularCoupling.h"

#include <gsl/gsl_sf_legendre.h>
#include <cmath>
#include <cstdlib>
#include <math.h>

// The two coupling coefficients are those of physics/ (AngularCoupling.h: Racah's closed forms in long double, Condon
// and Shortley's phases), one implementation for the engine and physics/.  Until 2026-11-16 they were GSL's 3j and 6j
// symbols with the phases written here; tests/upgrade/engine_coupling_gsl_test.cpp holds those expressions and
// compares.  An argument is an integer or a half-integer; lround(2 j) is its doubled value.
namespace {
int Twice(double j) { return (int)std::lround(2.0 * j); }
}  // namespace

// (j1 m1 j2 m2 | j3 m3).
double AngCoeff::ClebGord(double j1, double j2, double j3, double m1, double m2, double m3) {
  return azure2::physics::clebsch_gordan(Twice(j1), Twice(m1), Twice(j2), Twice(m2), Twice(j3), Twice(m3));
}

// W(j1 j2 l2 l1; j3 l3) = (-1)^(j1 + j2 + l2 + l1) {j1 j2 j3; l1 l2 l3} (Edmonds 6.2.13).
double AngCoeff::Racah(double j1, double j2, double l2, double l1, double j3, double l3) {
  return azure2::physics::racah_w(Twice(j1), Twice(j2), Twice(l2), Twice(l1), Twice(j3), Twice(l3));
}

double AngCoeff::LegendreP1(int l, double x) {
  if (l < 1) return 0.0;
  // gsl_sf_legendre_Plm includes the Condon-Shortley (-1)^m; Seyler and
  // Weller's P_l^M does not (their Eq. (12) carries the (-1)^M explicitly),
  // so for M = 1 the two differ by a sign.
  return -gsl_sf_legendre_Plm(l, 1, x);
}

complex AngCoeff::SphericalHarmonic(int l, int m, double theta, double phi) {
  if (l < 0 || std::abs(m) > l) return complex(0.0, 0.0);

  // gsl_sf_legendre_sphPlm gives the normalized associated Legendre function
  // for m >= 0 with the Condon-Shortley phase already applied:
  //   sphPlm(l,m,x) = sqrt((2l+1)/(4 pi) (l-m)!/(l+m)!) P_l^m(x)
  // so Y_l^m = sphPlm(l,m,cos theta) exp(i m phi).
  const int am = std::abs(m);
  const double norm = gsl_sf_legendre_sphPlm(l, am, std::cos(theta));
  complex y = norm * complex(std::cos(am * phi), std::sin(am * phi));

  // Y_l^{-m} = (-1)^m conj(Y_l^m).
  if (m < 0) {
    y = std::conj(y);
    if (am % 2) y = -y;
  }
  return y;
}
