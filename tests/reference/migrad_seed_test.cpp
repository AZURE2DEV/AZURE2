// MigradSeed (MigradSeed.h): the bounded curvature of MIGRAD's seed against hand-computed values, and the watch of
// the first step as a state machine (no model, no fit).
//
// BoundCurvature, per free parameter with gradient g, Gauss-Newton curvature G, initial step e and limits [l, u]:
// with a = (dx/du)^2 and d2 = d2x/du2 of Minuit2's transformation, the internal curvature is
//   max( a G + max(0, g d2),  a |g| / (10 e),  |g| sqrt(a) / 0.5 [two limits only] )
// and the curvature handed back is that over a, never below G.
#include "MigradSeed.h"

#include <atomic>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

int fails = 0;

void Check(bool ok, const char *what) {
  if (!ok) {
    std::printf("FAIL  %s\n", what);
    fails++;
  }
}

bool Near(double a, double b, double rel = 1e-12) { return std::fabs(a - b) <= rel * std::fabs(b); }

MigradSeed::Parameter Free(double step) {
  MigradSeed::Parameter p;
  p.free = true;
  p.step = step;
  return p;
}

MigradSeed::Parameter Limited(double step, double lower, double upper) {
  MigradSeed::Parameter p = Free(step);
  p.hasLower = p.hasUpper = true;
  p.lower = lower;
  p.upper = upper;
  return p;
}

}  // namespace

int main() {
  // ---- the curvature -------------------------------------------------------------------------------------------
  {
    MigradSeed::Parameter fixed = Free(0.1);
    fixed.free = false;
    MigradSeed::Parameter lower = Free(0.5);
    lower.hasLower = true;
    lower.lower = 0.0;
    const std::vector<MigradSeed::Parameter> parameters = {Free(0.1), Free(1.0), Limited(0.1, 0.0, 1.0), Free(0.1),
                                                           fixed, lower};
    const std::vector<double> p = {3.0, 3.0, 0.99, 3.0, 3.0, 1.0};
    const std::vector<double> grad = {100.0, 1.0, -5.0, 0.0, 100.0, 4.0};
    std::vector<double> g2 = {1e-6, 10.0, 2.0, 1e-9, 1e-6, 0.0};
    const MigradSeed::Bounded done = MigradSeed::BoundCurvature(parameters, p, grad, g2);
    // 0: hardly constrained: |g| / (10 e) = 100 / 1 = 100; its unbounded step was 100 / 1e-6 / 0.1 = 1e9 steps
    Check(Near(g2[0], 100.0), "a near-zero curvature is raised to |g| / (10 step)");
    // 1: a step of 0.1 initial steps: left alone
    Check(g2[1] == 10.0, "a well-conditioned parameter is left alone");
    // 2: two limits, pushed against the upper one: a = 0.99 * 0.01, d2 = -0.49, g d2 = 2.45;
    //    a G + g d2 = 2.4698 > the step bound 0.0495 and the angle bound 5 sqrt(a) / 0.5 = 0.99499
    Check(Near(g2[2], (0.0099 * 2.0 + 2.45) / 0.0099), "the limit term g d2x/du2 is included");
    // 3: zero gradient: nothing to bound
    Check(g2[3] == 1e-9, "a zero gradient changes nothing");
    // 4: fixed
    Check(g2[4] == 1e-6, "a fixed parameter is left alone");
    // 5: lower limit only, x - l = 1: w = 2, a = 3/4, d2 = 1/8, g d2 = 0.5; step bound a |g| / (10 e) = 0.6
    Check(Near(g2[5], 0.6 / 0.75), "one limit: the step bound over (dx/du)^2");
    Check(done.free == 5 && done.raised == 3, "counts: 5 free, 3 raised");
    Check(done.worst == 5 && std::isinf(done.worstSteps), "the longest unbounded step: the parameter without curvature");
    const double edm = 0.5 * (1e4 / g2[0] + 1.0 / 10.0 + 25.0 / g2[2] + 16.0 / g2[5]);
    Check(Near(done.edm, edm), "the estimated distance to the minimum uses the bounded curvatures");

    std::vector<double> atMinimum = {1e-6, 10.0, 2.0, 1e-9, 1e-6, 0.0};
    const std::vector<double> before = atMinimum;
    const MigradSeed::Bounded none =
        MigradSeed::BoundCurvature(parameters, p, std::vector<double>(6, 0.0), atMinimum);
    Check(atMinimum == before && none.raised == 0 && none.edm == 0.0, "at a minimum (zero gradient) nothing is raised");

    std::vector<double> shortG2 = {1.0};
    const MigradSeed::Bounded mismatch = MigradSeed::BoundCurvature(parameters, p, grad, shortG2);
    Check(shortG2[0] == 1.0 && mismatch.free == 0, "vectors of another length are left alone");
  }

  // ---- the watch of the first step -----------------------------------------------------------------------------
  const std::vector<MigradSeed::Parameter> one = {Free(100.0)};
  const std::vector<double> at = {1.0};
  auto seeded = [&](MigradSeed &seed, double g) {
    std::vector<double> g2 = {1.0};
    seed.Bound(at, std::vector<double>{g}, g2);  // step 100: no bound; edm = g^2 / 2
    seed.Evaluated(10.0, false);                 // the value at the seed point
    seed.GradientCalled();
  };
  {
    // a line search without improvement, then Minuit2's parallel numerical gradient: stalled, stop requested
    MigradSeed seed(one, 1.0);
    seeded(seed, 100.0);
    for (int i = 0; i < 12; i++) seed.Evaluated(11.0 + i, false);
    Check(!seed.Stalled() && !seed.Stop()->load(), "no verdict while the line search runs");
    seed.Evaluated(10.0, true);
    Check(seed.Stalled() && seed.Stop()->load(), "no improvement, then a parallel evaluation: stalled and stopped");
  }
  {
    // the same on one thread: recognised by the length of the run
    MigradSeed seed(one, 1.0);
    seeded(seed, 100.0);
    bool notice = false;
    for (int i = 0; i < MigradSeed::kHesseRun - 1; i++) notice = seed.Evaluated(11.0, false) || notice;
    Check(!seed.Stalled() && !notice, "39 serial evaluations: still a line search");
    notice = seed.Evaluated(10.0, false);
    Check(seed.Stalled() && notice, "the 40th serial evaluation: stalled, and the HESSE notice once");
    Check(!seed.Evaluated(10.0, false), "the notice is given once");
  }
  {
    // an improvement in the line search: Minuit2 goes on; a later HESSE is its own business
    MigradSeed seed(one, 1.0);
    seeded(seed, 100.0);
    seed.Evaluated(12.0, false);
    seed.Evaluated(9.5, false);
    seed.Evaluated(10.0, true);
    Check(!seed.Stalled() && !seed.Stop()->load(), "an improvement in the line search: not stalled");
  }
  {
    // the start is a minimum (edm = 1e-8 / 2 < 1e-3 up): HESSE is the verification, not a stall
    MigradSeed seed(one, 1.0);
    seeded(seed, 1e-4);
    seed.Evaluated(10.0, false);
    seed.Evaluated(10.0, true);
    Check(!seed.Stalled(), "a start at a minimum is not a stall");
  }
  {
    // after the second gradient call the watch is over
    MigradSeed seed(one, 1.0);
    seeded(seed, 100.0);
    seed.Evaluated(9.0, false);
    seed.GradientCalled();
    for (int i = 0; i < 100; i++) seed.Evaluated(20.0, i > 50);
    Check(!seed.Stalled(), "after the second gradient call nothing is a stall");
  }
  {
    // the caller's stop request is forwarded
    MigradSeed seed(one, 1.0);
    std::atomic<bool> caller{false};
    seed.Forward(&caller);
    seed.Forward(nullptr);
    Check(!seed.Stop()->load(), "no stop without a request");
    caller.store(true);
    seed.Forward(&caller);
    Check(seed.Stop()->load() && !seed.Stalled(), "the caller's stop request is forwarded, and is not a stall");
  }
  {
    // the first Bound is the seed's
    MigradSeed seed(one, 1.0);
    std::vector<double> g2 = {1.0};
    seed.Bound(at, std::vector<double>{3.0}, g2);
    seed.Bound(at, std::vector<double>{30.0}, g2);
    Check(Near(seed.Seed().edm, 4.5), "Seed() is the first Bound's record");
  }

  if (fails == 0) std::printf("ok  migrad_seed\n");
  return fails == 0 ? 0 : 1;
}
