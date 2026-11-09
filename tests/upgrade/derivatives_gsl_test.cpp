// physics/ central_derivative (Derivatives.h) against GSL's gsl_deriv_central, which it ports (docs/upgrade/coulomb.md
// section 1.4, "Bound channel functions in physics/"): the same bits for the same function on 10^4 cases -- smooth
// functions (exp, sin, a shift-function-like ratio), a noisy one (rounding noise of 1e-10 relative, which makes GSL
// take its second, optimal step) and a kinked one, at steps 1e-6 to 1e-2.  The error estimate is compared too.
#include <gsl/gsl_deriv.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>

#include "azure2/physics/Derivatives.h"

namespace {

double Adaptor(double x, void *p) { return (*static_cast<std::function<double(double)> *>(p))(x); }

// A deterministic pseudo-noise of relative size 1e-10, a function of x alone.
double Noise(double x) {
  std::uint64_t u;
  std::memcpy(&u, &x, sizeof u);
  u ^= u >> 33;
  u *= 0xff51afd7ed558ccdULL;
  u ^= u >> 33;
  return ((double)(u >> 11) / 9007199254740992.0 - 0.5) * 2e-10;
}

}  // namespace

int main() {
  const std::function<double(double)> functions[] = {
      [](double x) { return std::exp(0.7 * x); },
      [](double x) { return std::sin(3.0 * x) + 0.1 * x; },
      [](double x) { return -1.0 - x * x / (1.0 + std::fabs(x)); },
      [](double x) { return std::exp(0.7 * x) * (1.0 + Noise(x)); },
      [](double x) { return std::fabs(x - 0.3) + 0.5 * x; },
  };
  const double steps[] = {1e-6, 1e-5, 1e-4, 1e-3, 1e-2};
  long cases = 0, differ = 0, optimal = 0;
  for (const auto &fn : functions)
    for (double h : steps)
      for (int i = 0; i < 400; i++) {
        const double x = -2.0 + 4.0 * i / 399.0;
        std::function<double(double)> f = fn;
        gsl_function F;
        F.function = &Adaptor;
        F.params = &f;
        double gr, ge, pe;
        gsl_deriv_central(&F, x, h, &gr, &ge);
        const double pr = azure2::physics::central_derivative(f, x, h, &pe);
        // GSL takes the optimal step when rounding < truncation: count those cases to show the branch is exercised.
        int evaluations = 0;
        std::function<double(double)> counting = [&](double y) { evaluations++; return f(y); };
        azure2::physics::central_derivative(counting, x, h);
        if (evaluations > 4) optimal++;
        cases++;
        if (std::memcmp(&gr, &pr, sizeof gr) != 0 || std::memcmp(&ge, &pe, sizeof ge) != 0) {
          if (differ < 5) std::printf("  FAIL  x = %.17g h = %g: GSL %.17g (%.3g), physics %.17g (%.3g)\n", x, h, gr, ge, pr, pe);
          differ++;
        }
      }
  std::printf("  %ld cases (%ld with the second, optimal step), %ld differ\n", cases, optimal, differ);
  const bool ok = differ == 0 && optimal > 100;
  std::printf("%s\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
