// Reference check for Equation::Evaluate (compiled postfix program) against
// Equation::EvaluateTokens (the original string-rewriting evaluator).
//
// The two must agree for every operator, function, parameter and negation the
// parser supports.  EvaluateTokens rounds each intermediate result to 15
// significant digits, so agreement is required to 1e-12 relative, not bitwise.
// The timing line documents why Evaluate was rewritten: energy-dependent
// convolution kernels are evaluated ~1e7 times by the target-effect Jacobian.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "Config.h"
#include "Equation.h"

static int failures = 0;

int main() {
  Config configure(std::cout);
  struct Case {
    std::string expr;
    std::vector<double> params;
    std::vector<double> xs;
  };
  const std::vector<double> xsPos = {0.05, 0.3, 0.7, 1.0, 1.4142, 2.5, 7.0};
  const std::vector<Case> cases = {
      // the energy-dependent Gaussian sigma of the 22Ne+a gas-target data
      {"a0+a1*x+a2*x^2+a3*x^3+a4*x^4", {5.087936e-03, 9.953659e-03, -1.052552e-02, 4.064479e-03, -5.446409e-04}, xsPos},
      {"a0+a1*x+a2*x^2+a3*x^3+a4*x^4", {2.821e-20, 4.669e-20, -4.177e-20, 1.3646e-20, -1.547e-21}, xsPos},
      {"x", {}, xsPos},
      {"3.5e-3", {}, xsPos},
      {"2.5E+2*x-1.0e-1", {}, xsPos},
      {"-x+2", {}, xsPos},
      {"-(x*x)/(1+x)", {}, xsPos},
      {"2^x^2", {}, xsPos},             // right associative
      {"(2^x)^2", {}, xsPos},
      {"x-1-2-3", {}, xsPos},           // left associative
      {"x/2/3", {}, xsPos},
      {"exp(-x/3)*sqrt(x)", {}, xsPos},
      {"e^(x)-exp(x)", {}, xsPos},
      {"ln(x)+log(x)", {}, xsPos},
      {"sin(x)*cos(x)-tan(x/10)", {}, xsPos},
      {"asin(x/10)+acos(x/10)+atan(x)", {}, xsPos},
      {"a0*exp(-a1*(x-a2)^2)+a3", {1.7, 0.4, 1.1, -0.2}, xsPos},
      {"sqrt(a0+a1*ln(x))", {3.0, 0.5}, xsPos},
  };
  for (const Case &c : cases) {
    Equation eq(c.expr, c.params, configure);
    for (double x : c.xs) {
      double fast = eq.Evaluate(configure, x);
      double ref = eq.EvaluateTokens(configure, x);
      double denom = std::fabs(ref) > 1e-300 ? std::fabs(ref) : 1.0;
      bool ok = std::fabs(fast - ref) / denom <= 1e-12 || (std::isnan(fast) && std::isnan(ref));
      if (!ok) {
        failures++;
        std::printf("FAIL  %-34s x=%-7g compiled %.17g  tokens %.17g\n", c.expr.c_str(), x, fast, ref);
      }
    }
    std::printf("%-34s %s\n", c.expr.c_str(), "checked");
  }
  // a copied Equation (TargetEffect stores them by value) must keep its program
  Equation original("a0+a1*x", std::vector<double>{1.0, 2.0}, configure);
  Equation copy = original;
  if (std::fabs(copy.Evaluate(configure, 3.0) - 7.0) > 1e-15) {
    failures++;
    std::printf("FAIL  copied equation\n");
  }
  // SetParameter after construction must be seen by the compiled program
  Equation withParam("a0*x", 1, configure);
  withParam.SetParameter(0, 4.0, configure);
  if (std::fabs(withParam.Evaluate(configure, 2.0) - 8.0) > 1e-15) {
    failures++;
    std::printf("FAIL  SetParameter\n");
  }

  Equation poly(cases[0].expr, cases[0].params, configure);
  const int n = 200000;
  double sink = 0.0;
  auto t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < n; i++) sink += poly.Evaluate(configure, 1.0 + 1e-6 * i);
  auto t1 = std::chrono::steady_clock::now();
  for (int i = 0; i < n; i++) sink -= poly.EvaluateTokens(configure, 1.0 + 1e-6 * i);
  auto t2 = std::chrono::steady_clock::now();
  double fastNs = std::chrono::duration<double, std::nano>(t1 - t0).count() / n;
  double refNs = std::chrono::duration<double, std::nano>(t2 - t1).count() / n;
  std::printf("timing: compiled %.0f ns, token rewriting %.0f ns per evaluation (x%.0f); residual %.1e\n",
              fastNs, refNs, refNs / fastNs, sink);

  std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
