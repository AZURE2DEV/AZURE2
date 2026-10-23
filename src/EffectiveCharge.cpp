#include "EffectiveCharge.h"
#include "PPair.h"
#include "Constants.h"
#include <gsl/gsl_sf_bessel.h>
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
 * The operator() function takes the radius as an argument, carries out
 * integration, and returns the rho-dependent effective charge.
 */

double EffectiveCharge::operator()(double r) {
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  F.function = &Integrand;
  F.params = &L_;

  double intZ1 = 0.;
  double intZ2 = 0.;
  double intError;
  gsl_integration_qags(&F, 0., m2_ / (m1_ + m2_) * energy_ / hbarc * r, 0., 1.e-6, 1000, w, &intZ1, &intError);
  gsl_integration_qags(&F, 0., m1_ / (m1_ + m2_) * energy_ / hbarc * r, 0., 1.e-6, 1000, w, &intZ2, &intError);

  gsl_integration_workspace_free(w);

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
