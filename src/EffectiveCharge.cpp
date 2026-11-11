#include "EffectiveCharge.h"
#include "PPair.h"
#include "Constants.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_integration.h>

extern double DoubleFactorial(int);

/*!
 * The constructor takes a PPair object, gamma energy, and multipolarity as arguments.
 */

EffectiveCharge::EffectiveCharge(PPair *pair, double energy, int L, int li, int lf) :
  z1_(pair->GetZ(1)),
  z2_(pair->GetZ(2)),
  L_(L),
  m1_(pair->GetM(1)),
  m2_(pair->GetM(2)),
  energy_(energy),
  li_(li),
  lf_(lf) {
}

namespace {
/*!
 * Tangential coefficient of the transverse electric multipole field after its gradient part (the charge term)
 * is removed: (1/k) curl[j_L X_LM] = grad(Phi) + i w_t(x) r grad Y_LM, with
 *   w_t(x) = [(L+1) j_{L-1}(x) - L j_{L+1}(x)] / ((2L+1) sqrt(L(L+1))) - sqrt(L(L+1)) G_L(x)/x,
 * G_L(x) = int_0^x j_L(t)/t dt.  O(x^(L+1)): it vanishes in the long-wavelength limit.
 */
double TangentialW(int L, double x, double G) {
  if (x <= 0.) return 0.;
  double ll = sqrt(L * (L + 1.));
  return ((L + 1) * gsl_sf_bessel_jl(L - 1, x) - L * gsl_sf_bessel_jl(L + 1, x)) / ((2 * L + 1) * ll) - ll * G / x;
}
}  // namespace

double EffectiveCharge::Integrand(double x, void *p) {
  int *L = (int *)p;
  return gsl_sf_bessel_jl(*L, x) / x;
}

/*!
 * G_L(y) = int_0^y j_L(t)/t dt.  Closed forms for the multipolarities external capture uses (E1, E2):
 *   G_1(y) = Si(y)/2 + (y cos y - sin y)/(2 y^2),   G_2(y) = 1/3 + (y cos y - sin y)/y^3,
 * and below y = 0.5, where they cancel, the series sum_k (-1/2)^k y^(L+2k) / (k! (2L+2k+1)!! (L+2k)).
 * The adaptive quadrature used before 2026-10-24 failed ("roundoff error", which AZURE2's GSL handler turns into an
 * exception) at the large y of a near-threshold Whittaker tail; it remains only for L > 2.
 */
double EffectiveCharge::IntJlOverT(int L, double y) {
  if (y <= 0.) return 0.;
  if (y < 0.5) {
    double sum = 0., term = pow(y, L) / DoubleFactorial(2 * L + 1);  // k = 0 term without 1/(L+2k)
    for (int k = 0; k < 8; ++k) {
      sum += term / (L + 2 * k);
      term *= -0.5 * y * y / ((k + 1) * (2 * L + 2 * k + 3));
    }
    return sum;
  }
  if (L == 1) return 0.5 * gsl_sf_Si(y) + (y * cos(y) - sin(y)) / (2. * y * y);
  if (L == 2) return 1. / 3. + (y * cos(y) - sin(y)) / (y * y * y);
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = &Integrand;
  F.params = &L;
  double result = 0., error;
  gsl_integration_qags(&F, 0., y, 0., 1.e-6, 1000, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

/*!
 * The operator() function takes the radius as an argument, carries out
 * integration, and returns the rho-dependent effective charge.
 */

double EffectiveCharge::operator()(double r) {
  double intZ1 = IntJlOverT(L_, m2_ / (m1_ + m2_) * energy_ / hbarc * r);
  double intZ2 = IntJlOverT(L_, m1_ / (m1_ + m2_) * energy_ / hbarc * r);

  double charge = z1_ * intZ1 + pow(-1, L_) * z2_ * intZ2;
  if (li_ >= 0 && lf_ >= 0 && r > 0.) {
    // Convection current of the relative motion, J = sum_i e Z_i p_i/m_i (cluster i at c_i r, c_1 = m2/M,
    // c_2 = -m1/M, p_1 = -p_2 = p), against the tangential field.  For u_f/r Y_lf and u_i/r Y_li,
    //   <f| {p, w(r) r grad Y_LM} |i> = i hbar (lf(lf+1) - li(li+1)) <Y_lf|Y_LM|Y_li> int u_f u_i w/r dr,
    // which is the same angular factor and the same radial functions as the charge term; current conservation
    // fixes its normalization relative to the charge (Siegert) term, checked against exact eigenstates of one
    // Woods-Saxon + Coulomb Hamiltonian to 1e-9 (R-matrix archive 12C+a_onefile/10-23-26_ec_tangential_current).
    double k = energy_ / hbarc;
    double M = m1_ + m2_;
    double x1 = m2_ / M * k * r, x2 = m1_ / M * k * r;
    double dl = lf_ * (lf_ + 1.) - li_ * (li_ + 1.);
    double conv = z1_ * hbarc / (2. * m1_ * uconv) * TangentialW(L_, x1, intZ1) +
        pow(-1, L_) * z2_ * hbarc / (2. * m2_ * uconv) * TangentialW(L_, x2, intZ2);
    charge += -dl / sqrt(L_ * (L_ + 1.)) * conv / r;
  }
  return sqrt(fstruc * hbarc) * L_ * DoubleFactorial(2 * L_ + 1) /
      pow((energy_ / hbarc * r), (double)L_) * charge;
}
