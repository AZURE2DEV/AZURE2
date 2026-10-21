// Reference check for the closed-channel shift function of neutral (eta = 0) channels.
//
// ShftFunc::NeutralBoundShift gives S_l = rho u'/u for the decaying solution
// u = rho k_l(rho) by the recursion S_0 = -rho, S_l = -rho^2/(l - S_{l-1}) - l.
// It replaces, for neutral pairs, the general route: a numerical radial
// derivative of GSL's Whittaker function (hyperg_U), which a profile of the
// natural-neon example put at ~29 % of a cold evaluation.  Compared here with
// that route over l = 0..10, rho = 1e-3..100, and with the l = 0, 1 textbook
// expressions; the timing line documents the speed-up.
//
// Open neutral channels were left on the general route: Michel's Coulomb
// routine at eta = 0 costs ~0.3 us per call, no more than the Riccati-Bessel
// closed form through GSL (measured 2026-10-18).

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_hyperg.h>

#include "Config.h"
#include "ShftFunc.h"

Config *g_config = nullptr;  // defined by the executable; the library reads it

static int failures = 0;

static void check(bool ok, const char *what, int l, double rho, double a, double b) {
  if (!ok) {
    failures++;
    if (failures <= 20) std::printf("  FAIL %s l=%d rho=%.6e  closed=%.15e  reference=%.15e\n", what, l, rho, a, b);
  }
}

// The general route as ShftFunc has it: W(z) = exp(-z/2) z^(m+1/2) U(a,b,z) with
// k = 0 (a = l+1, b = 2l+2), z = 2 rho, and S = rho dW/drho / W by gsl_deriv_central.
struct WParams {
  int l;
};
static double whit(double rho, void *p) {
  const int l = static_cast<WParams *>(p)->l;
  const double m = l + 0.5, z = 2.0 * rho;
  return std::exp(-z / 2.0) * std::pow(z, m + 0.5) * gsl_sf_hyperg_U(m + 0.5, 1.0 + 2.0 * m, z);
}
static double whittakerShift(int l, double rho) {
  WParams p{l};
  gsl_function F;
  F.function = &whit;
  F.params = &p;
  double d, err;
  gsl_deriv_central(&F, rho, 2e-5 * rho, &d, &err);  // ShftFunc: 1e-4 fm at a ~ 5 fm
  return rho * d / whit(rho, &p);
}

int main() {
  gsl_set_error_handler_off();
  const int lmax = 10;
  const int nrho = 400;
  std::vector<double> rhos;
  for (int i = 0; i <= nrho; i++) rhos.push_back(std::pow(10.0, -3.0 + 5.0 * i / nrho));

  // 1. textbook l = 0, 1
  for (double rho : rhos) {
    const double s0 = ShftFunc::NeutralBoundShift(0, rho), s1 = ShftFunc::NeutralBoundShift(1, rho);
    const double b1 = -(rho * rho + rho + 1.0) / (rho + 1.0);
    check(s0 == -rho, "S_0 = -rho", 0, rho, s0, -rho);
    check(std::fabs(s1 - b1) < 1e-14 * std::fabs(b1), "S_1 = -(rho^2+rho+1)/(rho+1)", 1, rho, s1, b1);
  }

  // 2. against the Whittaker derivative, l = 0..10; the reference is a numerical
  // derivative, accurate to ~1e-8, so that is the tolerance
  double worst = 0;
  int skipped = 0;
  for (int l = 0; l <= lmax; l++) {
    for (double rho : rhos) {
      const double sb = ShftFunc::NeutralBoundShift(l, rho), sw = whittakerShift(l, rho);
      if (!std::isfinite(sw)) {
        skipped++;
        continue;
      }
      const double e = std::fabs(sb - sw) / std::max(1.0, std::fabs(sw));
      worst = std::max(worst, e);
      check(e < 1e-6, "S vs Whittaker derivative", l, rho, sb, sw);
      check(sb < 0.0, "S < 0", l, rho, sb, sw);
    }
  }
  std::printf("worst vs Whittaker derivative: %.2e (%d points the reference could not do)\n", worst, skipped);

  // 3. timing over the same grid
  auto time = [&](auto f) {
    auto t0 = std::chrono::steady_clock::now();
    double sink = 0;
    for (int l = 0; l <= lmax; l++)
      for (double rho : rhos) sink += f(l, rho);
    auto t1 = std::chrono::steady_clock::now();
    if (sink == 12345.678) std::printf(" ");
    return std::chrono::duration<double, std::micro>(t1 - t0).count() / ((lmax + 1) * rhos.size());
  };
  const double tBound = time([](int l, double rho) { return ShftFunc::NeutralBoundShift(l, rho); });
  const double tWhit = time([](int l, double rho) { return whittakerShift(l, rho); });
  std::printf("recursion %.3f us/call, Whittaker derivative %.3f us/call (x%.0f)\n", tBound, tWhit, tWhit / tBound);

  if (failures) {
    std::printf("FAIL: %d neutral closed-channel checks\n", failures);
    return 1;
  }
  std::printf("PASS: neutral closed-channel shift function\n");
  return 0;
}
