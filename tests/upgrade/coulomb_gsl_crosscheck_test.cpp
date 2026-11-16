// physics/coulomb's default algorithm (coulsm) against GSL's gsl_sf_coulomb_wave_FG_e, an independent algorithm
// (Steed's method with its own continued fractions) that physics/ does not link (docs/upgrade/coulomb.md
// section 5).  Grid: l = 0-10, eta = 0 and 1e-3-50, rho = 1e-3-100.  Compared where GSL is accurate (pycoulomb's
// measurements against mpmath: everywhere for eta = 0, and beyond the outer turning point rho_tp = eta +
// sqrt(eta^2 + l(l+1)) for eta > 0, to 5e-15): G and G' for eta = 0 (coulsm's F there comes from an unstable
// recurrence and is not used alone, coulomb.md), F, F', G, G' relative to the amplitude beyond rho_tp.  coulsm is
// good to 1e-9 there, so a disagreement above 1e-8 means a wrong source, sign or normalisation.
//
// Inside the outer turning point GSL's routine is not accurate, which is why the engine's option use_gsl_coul was
// retired on 2026-11-16 (docs/upgrade/coulomb.md sections 8.2 and 9.1).  The second part holds that reason as numbers:
// at seven points below the barrier, P = rho/(F^2 + G^2) and S = rho (F F' + G G')/(F^2 + G^2) from mpmath
// (coulombf, coulombg and their derivatives at 40 digits, the Wronskian met to 1e-25); coulsm must give P to 1e-9 and
// S to 1e-9 of max(1, |S|).  GSL's P at the same points is printed, not required (it is off by 0.5 to 5.2 % with GSL
// 2.6; another machine's GSL is not this code's to test): a note says so when a GSL gives the right value.
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
  struct Inside {
    int l;
    double eta, rho, P, S;  // P and S from mpmath
  };
  const Inside inside[] = {{3, 0.6449, 1.2011, 0.0023546071064841554, -2.9657092304755023},
                           {3, 1.0702, 1.2032, 0.0007146450357543889, -3.1374524375945227},
                           {2, 2.288, 1.2627, 0.00035088425288326843, -2.7819978555660874},
                           {1, 0.9, 1.25, 0.18676170227671204, -1.0676045408887106},
                           {3, 25.67, 0.533, 7.460862766382192e-63, -5.9431094761935055},
                           {2, 47.79, 0.2237, 1.5963802885173925e-123, -4.9448933164450665},
                           {2, 14.77, 1.945, 8.448625072439186e-28, -7.446482573868072}};
  double worstOwn = 0.0, gslLeast = 1.0, gslMost = 0.0;
  for (const Inside &p : inside) {
    const CoulombWaves w = coulomb(p.l, p.eta, p.rho);
    const double d = w.F * w.F + w.G * w.G;
    const double eP = std::fabs(p.rho / d / p.P - 1.0);
    const double eS = std::fabs(p.rho * (w.F * w.Fp + w.G * w.Gp) / d - p.S) / std::fmax(1.0, std::fabs(p.S));
    worstOwn = std::fmax(worstOwn, std::fmax(eP, eS));
    if (w.status != Status::ok || !(eP <= 1e-9) || !(eS <= 1e-9)) {
      std::printf("FAIL: inside the turning point, l=%d eta=%g rho=%g: coulsm's P is off by %.2e, S by %.2e (mpmath)\n",
                  p.l, p.eta, p.rho, eP, eS);
      ++failures;
    }
    gsl_sf_result F, Fp, G, Gp;
    double eF = 0.0, eG = 0.0;
    if (gsl_sf_coulomb_wave_FG_e(p.eta, p.rho, p.l, 0, &F, &Fp, &G, &Gp, &eF, &eG) != GSL_SUCCESS || eF != 0.0 || eG != 0.0) {
      std::printf("note: inside the turning point, l=%d eta=%g rho=%g: GSL gave no unscaled value\n", p.l, p.eta, p.rho);
      continue;
    }
    const double g = std::fabs(p.rho / (F.val * F.val + G.val * G.val) / p.P - 1.0);
    gslLeast = std::fmin(gslLeast, g);
    gslMost = std::fmax(gslMost, g);
    if (!(g > 1e-3)) {
      std::printf("note: inside the turning point, l=%d eta=%g rho=%g: GSL's P is off by %.2e only (mpmath); the defect "
                  "that retired use_gsl_coul is not seen with this GSL\n", p.l, p.eta, p.rho, g);
    }
  }
  std::printf("inside the outer turning point, 7 points against mpmath: coulsm within %.2e (P and S); GSL's P off by "
              "%.2e to %.2e\n", worstOwn, gslLeast, gslMost);
  if (compared < 300) {
    std::printf("FAIL: only %d points compared\n", compared);
    ++failures;
  }
  return failures ? 1 : 0;
}
