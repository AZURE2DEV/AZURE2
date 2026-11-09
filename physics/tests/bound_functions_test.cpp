// The bound channel functions and the derivatives of physics/ (ChannelFunctions.h, Derivatives.h; docs/upgrade/
// coulomb.md section 1.4, "Bound channel functions in physics/"), checked against things that do not use them:
//  1. bessel_k_scaled against mpmath (tests/data/bessel_k.csv, make_bessel_k.py there), n = 0-21, x = 1e-3-500;
//  2. zero_energy_shift against the limit of whittaker_scaled's rho W'/W at threshold: S is linear in the binding
//     near 0, so 2 S(B) - S(2B) at B = 1e-8 MeV extrapolates to S(0) (p+12C, p+p, a+12C, l = 0-3);
//  3. neutral_bound_shift against rho u'/u, u = rho k_l(rho), from the closed forms of k_0 ... k_2;
//  4. coulomb_bound_shift: whittaker_scaled itself for eta_b <= 1000, linear in B below, continuous at B_min;
//  5. central_derivative and second_derivative on polynomials, where the rules are exact up to rounding, and on
//     exp and sin.
// Argument: the path of bessel_k.csv.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <string>

#include "azure2/physics/ChannelFunctions.h"
#include "azure2/physics/Constants.h"
#include "azure2/physics/Derivatives.h"
#include "azure2/physics/coulomb/Whittaker.h"

using namespace azure2::physics;

namespace {

int failures = 0;
long checks = 0;

void close(const std::string &what, double got, double want, double rel) {
  checks++;
  const double d = std::fabs(got - want);
  if (!(d <= rel * std::fabs(want))) {
    std::printf("  FAIL  %s: got %.17g, want %.17g (rel %.3g)\n", what.c_str(), got, want, d / std::fabs(want));
    failures++;
  }
}

struct Pair {
  const char *name;
  double z1z2, mu_u, a_fm;
};
const Pair kPairs[] = {{"p+12C", 6, 12.0 / 13.0, 4.2}, {"p+p", 1, 0.5, 4.0}, {"a+12C", 12, 3.0, 5.5}};

BoundShiftArgs Args(const Pair &p) {
  BoundShiftArgs a;
  a.mu_MeV = p.mu_u * kAtomicMassUnit_MeV;
  a.c_fm = p.z1z2 * kFineStructure * a.mu_MeV / kHbarC_MeV_fm;
  a.a_fm = p.a_fm;
  const double mu = a.mu_MeV, z = p.z1z2, r = p.a_fm;
  a.eta_rho = [mu, z, r](double B, double &eta, double &rho) {
    eta = z * kFineStructure * std::sqrt(mu / (2.0 * B));
    rho = std::sqrt(2.0 * mu * B) / kHbarC_MeV_fm * r;
  };
  return a;
}

double RhoDlog(int l, const BoundShiftArgs &a, double B) {
  double eta, rho;
  a.eta_rho(B, eta, rho);
  return coulomb::whittaker_scaled(l, eta, rho).rho_dlog;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::printf("usage: %s bessel_k.csv\n", argv[0]);
    return 2;
  }
  std::printf("bessel_k_scaled vs mpmath\n");
  {
    std::ifstream in(argv[1]);
    std::string line;
    long n = 0;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream s(line);
      int order;
      double x, ref;
      s >> order >> x >> ref;
      close("K_" + std::to_string(order) + "(" + std::to_string(x) + ")", bessel_k_scaled(order, x), ref, 2e-15);
      n++;
    }
    if (n < 300) {
      std::printf("  FAIL  only %ld reference values read\n", n);
      failures++;
    }
  }

  std::printf("zero_energy_shift vs the threshold limit of rho W'/W\n");
  for (const Pair &p : kPairs) {
    const BoundShiftArgs a = Args(p);
    for (int l = 0; l <= 3; l++) {
      const double B = 1e-8;
      const double limit = 2.0 * RhoDlog(l, a, B) - RhoDlog(l, a, 2.0 * B);
      close(std::string(p.name) + " S(0) l=" + std::to_string(l), zero_energy_shift(a.c_fm, a.a_fm, l), limit, 1e-8);
    }
  }
  close("S_2(0) without Coulomb", zero_energy_shift(0.0, 4.0, 2), -2.0, 0.0);

  std::printf("neutral_bound_shift vs closed forms\n");
  for (double rho : {0.01, 0.3, 1.0, 2.5, 7.0}) {
    // u = rho k_l: e^-r, e^-r (1 + r) / r, e^-r (r^2 + 3 r + 3) / r^2; S = r u'/u.
    close("l=0", neutral_bound_shift(0, rho), -rho, 1e-15);
    close("l=1", neutral_bound_shift(1, rho), -1.0 - rho * rho / (1.0 + rho), 1e-14);
    const double q = rho * rho + 3.0 * rho + 3.0;
    close("l=2", neutral_bound_shift(2, rho), -rho - 2.0 + rho * (2.0 * rho + 3.0) / q, 1e-14);
  }

  std::printf("coulomb_bound_shift\n");
  for (const Pair &p : kPairs) {
    const BoundShiftArgs a = Args(p);
    const double Bmin = std::pow(kHbarC_MeV_fm * (a.c_fm / 1000.0), 2.0) / (2.0 * a.mu_MeV);  // as the function
    for (int l = 0; l <= 3; l++) {
      const std::string t = std::string(p.name) + " l=" + std::to_string(l);
      for (double B : {Bmin, 3.0 * Bmin, 1e-3, 0.1, 2.0}) close(t + " exact", coulomb_bound_shift(l, B, a), RhoDlog(l, a, B), 0.0);
      const double s0 = zero_energy_shift(a.c_fm, a.a_fm, l), sMin = RhoDlog(l, a, Bmin);
      for (double f : {0.0, 0.25, 0.5, 0.999999}) {
        if (f == 0.0) continue;
        close(t + " band", coulomb_bound_shift(l, f * Bmin, a), s0 + (sMin - s0) * f, 1e-14);
      }
      close(t + " continuity at B_min", coulomb_bound_shift(l, Bmin * (1.0 - 1e-12), a), sMin, 1e-11);
    }
  }

  std::printf("derivatives\n");
  auto cubic = [](double x) { return 2.0 * x * x * x - x * x + 3.0 * x - 1.0; };
  for (double x : {-1.5, 0.0, 0.7, 4.0}) {
    close("d/dx cubic", central_derivative(cubic, x, 1e-3), 6.0 * x * x - 2.0 * x + 3.0, 1e-9);
    close("d2/dx2 cubic", second_derivative(cubic, x, 1e-2, 1e-3, false, Side::Below), 12.0 * x - 2.0, 1e-8);
    close("d2/dx2 quadratic one-sided below", second_derivative([](double y) { return 3.0 * y * y - y; }, x, 1e-2, 1e-3, true, Side::Below), 6.0, 1e-6);
    close("d2/dx2 quadratic one-sided above", second_derivative([](double y) { return 3.0 * y * y - y; }, x, 1e-2, 1e-3, true, Side::Above), 6.0, 1e-6);
  }
  for (double x : {0.1, 1.0, 3.0}) {
    close("d/dx exp", central_derivative([](double y) { return std::exp(y); }, x, 1e-3), std::exp(x), 1e-10);
    close("d/dx sin", central_derivative([](double y) { return std::sin(y); }, x, 1e-3), std::cos(x), 1e-10);
    close("d2/dx2 exp", second_derivative([](double y) { return std::exp(y); }, x, 1e-2, 1e-3, false, Side::Above), std::exp(x), 1e-9);
  }

  std::printf("%s: %ld checks, %d failure(s)\n", failures ? "FAIL" : "PASS", checks, failures);
  return failures ? 1 : 0;
}
