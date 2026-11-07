// physics/coulomb's default algorithm (coulsm) against GSL's gsl_sf_coulomb_wave_FG_e, an independent algorithm
// (Steed's method with its own continued fractions) that physics/ does not link (docs/upgrade/coulomb.md
// section 5).  Grid: l = 0-10, eta = 0 and 1e-3-50, rho = 1e-3-100.  Compared where GSL is accurate (pycoulomb's
// measurements against mpmath: everywhere for eta = 0, and beyond the outer turning point rho_tp = eta +
// sqrt(eta^2 + l(l+1)) for eta > 0, to 5e-15): G and G' for eta = 0 (coulsm's F there comes from an unstable
// recurrence and is not used alone, coulomb.md), F, F', G, G' relative to the amplitude beyond rho_tp.  coulsm is
// good to 1e-9 there, so a disagreement above 1e-8 means a wrong source, sign or normalisation.
#include <cmath>
#include <cstdio>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coulomb.h>

#include "azure2/physics/coulomb/Coulomb.h"

using namespace azure2::physics::coulomb;

int main() {
  gsl_set_error_handler_off();
  const double etas[] = {0.0, 1e-3, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0};
  const double rhos[] = {1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0};
  double worst = 0.0;
  int compared = 0, failures = 0;
  for (int l = 0; l <= 10; ++l)
    for (double eta : etas)
      for (double rho : rhos) {
        const bool beyond = rho >= eta + std::sqrt(eta * eta + l * (l + 1.0));
        if (eta > 0.0 && !beyond) continue;
        gsl_sf_result F, Fp, G, Gp;
        double eF = 0.0, eG = 0.0;
        if (gsl_sf_coulomb_wave_FG_e(eta, rho, l, 0, &F, &Fp, &G, &Gp, &eF, &eG) != GSL_SUCCESS) continue;
        const double gF = F.val * std::exp(eF), gFp = Fp.val * std::exp(eF), gG = G.val * std::exp(eG), gGp = Gp.val * std::exp(eG);
        const CoulombWaves w = coulomb(l, eta, rho);
        if (w.status != Status::ok) {
          std::printf("FAIL: coulsm failed at l=%d eta=%g rho=%g\n", l, eta, rho);
          ++failures;
          continue;
        }
        double d;
        if (eta == 0.0 && !beyond) {
          d = std::fmax(std::fabs(w.G - gG) / std::fabs(gG), std::fabs(w.Gp - gGp) / std::fabs(gGp));
        } else {
          const double a = std::hypot(gF, gG), ap = std::hypot(gFp, gGp);
          d = std::fmax(std::fmax(std::fabs(w.F - gF), std::fabs(w.G - gG)) / a, std::fmax(std::fabs(w.Fp - gFp), std::fabs(w.Gp - gGp)) / ap);
        }
        ++compared;
        if (!(d <= 1e-8)) {
          std::printf("FAIL: l=%d eta=%g rho=%g: coulsm and GSL differ by %.2e\n", l, eta, rho, d);
          ++failures;
        }
        worst = std::fmax(worst, d);
      }
  std::printf("coulsm against GSL: %d points, largest difference %.2e\n", compared, worst);
  if (compared < 300) {
    std::printf("FAIL: only %d points compared\n", compared);
    ++failures;
  }
  return failures ? 1 : 0;
}
