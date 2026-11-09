// physics/coulomb's Whittaker function against the engine's GSL path (docs/upgrade/coulomb.md section 1.4).
//
// GSL's U (gsl_sf_hyperg_U_e10_e, as in WhitFunc::Scaled) and physics' whittaker_scaled on the mpmath grid of
// physics/tests/data/whittaker.csv: both errors against mpmath are reported, and wherever GSL gives a value the two
// must agree within GSL's own error plus physics' tolerance -- so a physics value that drifted from the function the
// engine has used would fail, while GSL's inaccuracy does not.  The time per call of both is printed on an
// engine-like range (l = 0-3, eta_b = 0-20, rho = 0.1-20), for the engine part of step 3b (i); it is not a pass/fail.
// Usage: whittaker_gsl_test <path to physics/tests/data/whittaker.csv>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_hyperg.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "azure2/physics/coulomb/Whittaker.h"

using namespace azure2::physics::coulomb;

namespace {

struct Row {
  int l;
  double eta, rho, log10W;
};

/// log10 W by GSL, the engine's formula (WhitFunc::Scaled); false where GSL reports an error.
bool GslLog10W(int l, double eta, double rho, double &log10W) {
  gsl_sf_result_e10 u;
  const double a = l + 1.0 + eta, b = 2.0 * l + 2.0, z = 2.0 * rho;
  if (gsl_sf_hyperg_U_e10_e(a, b, z, &u) != GSL_SUCCESS || !(u.val > 0.0)) return false;
  log10W = -rho / std::log(10.0) + (l + 1.0) * std::log10(z) + std::log10(u.val) + (double)u.e10;
  return std::isfinite(log10W);
}

double PhysicsLog10W(int l, double eta, double rho) {
  const WhittakerScaled s = whittaker_scaled(l, eta, rho);
  return std::log10(s.mantissa) + (double)s.log10_scale;
}

double RelFromLog10(double d) { return std::fabs(std::expm1(d * std::log(10.0))); }

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s whittaker.csv\n", argv[0]);
    return 2;
  }
  gsl_set_error_handler_off();
  std::vector<Row> rows;
  {
    std::ifstream in(argv[1]);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      for (char &c : line)
        if (c == ',') c = ' ';
      std::istringstream ss(line);
      Row r{};
      std::string t[3];
      ss >> r.l >> t[0] >> t[1] >> t[2];
      r.eta = std::strtod(t[0].c_str(), nullptr);
      r.rho = std::strtod(t[1].c_str(), nullptr);
      r.log10W = std::strtod(t[2].c_str(), nullptr);
      rows.push_back(r);
    }
  }
  int failures = 0, gslMissing = 0, compared = 0;
  double worstGsl = 0.0, worstPhys = 0.0, worstExcess = 0.0;
  const double kPhysTol = 2e-11;  // physics_whittaker's tolerance for W
  for (const Row &r : rows) {
    const double p = PhysicsLog10W(r.l, r.eta, r.rho);
    worstPhys = std::fmax(worstPhys, RelFromLog10(p - r.log10W));
    double g;
    if (!GslLog10W(r.l, r.eta, r.rho, g)) {
      gslMissing++;
      continue;
    }
    compared++;
    const double eG = RelFromLog10(g - r.log10W);
    worstGsl = std::fmax(worstGsl, eG);
    // physics within GSL's own error (plus physics' tolerance) of GSL
    worstExcess = std::fmax(worstExcess, RelFromLog10(p - g) - (eG + kPhysTol));
  }
  std::printf("  %d points: GSL without a value at %d\n", (int)rows.size(), gslMissing);
  std::printf("  W against mpmath: GSL %.2e, physics %.2e\n", worstGsl, worstPhys);
  const bool agree = compared > 0 && worstExcess <= 0.0;
  std::printf("  %s  physics within GSL's error of GSL on the %d points where GSL gives a value (excess %.2e)\n",
              agree ? "ok  " : "FAIL", compared, worstExcess);
  if (!agree) failures++;

  // Timing on an engine-like range (not a pass/fail).
  std::vector<Row> grid;
  for (int l = 0; l <= 3; l++)
    for (double eta : {0.0, 0.5, 2.0, 5.0, 20.0})
      for (double rho : {0.1, 0.5, 1.0, 3.0, 10.0, 20.0}) grid.push_back({l, eta, rho, 0.0});
  const int reps = 200;
  volatile double sink = 0.0;
  auto t0 = std::chrono::steady_clock::now();
  for (int k = 0; k < reps; k++)
    for (const Row &r : grid) {
      double g;
      if (GslLog10W(r.l, r.eta, r.rho, g)) sink = sink + g;
    }
  auto t1 = std::chrono::steady_clock::now();
  for (int k = 0; k < reps; k++)
    for (const Row &r : grid) sink = sink + PhysicsLog10W(r.l, r.eta, r.rho);
  auto t2 = std::chrono::steady_clock::now();
  const double n = (double)reps * grid.size();
  const double usG = std::chrono::duration<double, std::micro>(t1 - t0).count() / n;
  const double usP = std::chrono::duration<double, std::micro>(t2 - t1).count() / n;
  std::printf("  time per call, engine-like range: GSL U (W only) %.2f us, physics whittaker_scaled (W and rho W'/W) %.2f us\n",
              usG, usP);
  std::printf("\n%s\n", failures ? "FAIL" : "all checks passed");
  return failures ? 1 : 0;
}
