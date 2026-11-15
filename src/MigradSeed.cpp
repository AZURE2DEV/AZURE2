// The seed of a MIGRAD fit with the analytic gradient: the bounded curvature and the watch of the first step
// (MigradSeed.h; docs/upgrade/fit.md, "The seed of MIGRAD").
#include "MigradSeed.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include <Minuit2/MnUserParameters.h>

namespace {

std::vector<MigradSeed::Parameter> FromMinuit(const ROOT::Minuit2::MnUserParameters &start) {
  const auto &rows = start.Parameters();
  std::vector<MigradSeed::Parameter> parameters(rows.size());
  for (size_t i = 0; i < rows.size(); i++) {
    MigradSeed::Parameter &s = parameters[i];
    s.free = !rows[i].IsFixed() && !rows[i].IsConst();
    s.hasLower = rows[i].HasLowerLimit();
    s.hasUpper = rows[i].HasUpperLimit();
    s.lower = rows[i].LowerLimit();
    s.upper = rows[i].UpperLimit();
    s.step = std::fabs(rows[i].Error());
  }
  return parameters;
}

}  // namespace

MigradSeed::Bounded MigradSeed::BoundCurvature(const std::vector<Parameter> &parameters, const std::vector<double> &p,
                                               const std::vector<double> &grad, std::vector<double> &g2) {
  Bounded done;
  if (parameters.size() != p.size() || grad.size() != p.size() || g2.size() != p.size()) return done;
  for (size_t i = 0; i < p.size(); i++) {
    const Parameter &s = parameters[i];
    if (!s.free) continue;
    done.free++;
    const double g = grad[i];
    if (!std::isfinite(g) || !std::isfinite(g2[i]) || g == 0.0) continue;
    // a = (dx/du)^2 and d2 = d2x/du2 of Minuit2's transformations, written with the external value x:
    // two limits, x = l + (u - l)(sin t + 1)/2; one limit, x = l - 1 + sqrt(t^2 + 1) (and its mirror image).
    double a = 1.0, d2 = 0.0;
    if (s.hasLower && s.hasUpper) {
      a = std::max(0.0, (p[i] - s.lower) * (s.upper - p[i]));
      d2 = -(p[i] - 0.5 * (s.lower + s.upper));
      if (a <= 1e-24 * (s.upper - s.lower) * (s.upper - s.lower)) continue;  // on the limit: the internal gradient is 0
    } else if (s.hasLower || s.hasUpper) {
      const double w = s.hasLower ? p[i] - s.lower + 1.0 : s.upper - p[i] + 1.0;
      if (w <= 1.0 + 1e-12) continue;
      a = (w * w - 1.0) / (w * w);
      d2 = (s.hasLower ? 1.0 : -1.0) / (w * w * w);
    }
    double internal = a * std::max(0.0, g2[i]) + std::max(0.0, g * d2);
    if (s.step > 0.0) internal = std::max(internal, a * std::fabs(g) / (kSteps * s.step));
    if (s.hasLower && s.hasUpper) internal = std::max(internal, std::fabs(g) * std::sqrt(a) / kInternal);
    const double bounded = internal / a;
    if (bounded > g2[i] * (1.0 + 1e-12)) {
      done.raised++;
      const double steps = (g2[i] > 0.0 && s.step > 0.0) ? std::fabs(g) / g2[i] / s.step
                                                         : std::numeric_limits<double>::infinity();
      if (steps > done.worstSteps) {
        done.worstSteps = steps;
        done.worst = (int)i;
      }
      g2[i] = bounded;
    }
    if (g2[i] > 0.0) done.edm += 0.5 * g * g / g2[i];
  }
  return done;
}

MigradSeed::MigradSeed(const ROOT::Minuit2::MnUserParameters &start, double up) : parameters_(FromMinuit(start)), up_(up) {}

MigradSeed::MigradSeed(const std::vector<Parameter> &parameters, double up) : parameters_(parameters), up_(up) {}

void MigradSeed::Bound(const std::vector<double> &p, const std::vector<double> &grad, std::vector<double> &g2) {
  const Bounded done = BoundCurvature(parameters_, p, grad, g2);
  std::lock_guard<std::mutex> lock(mutex_);
  if (!seeded_) {
    seeded_ = true;
    seed_ = done;
  }
}

void MigradSeed::GradientCalled() {
  gradients_.fetch_add(1);
  sinceGradient_.store(0);
}

bool MigradSeed::Evaluated(double objective, bool inParallel) {
  const long gradients = gradients_.load();
  if (gradients == 0) {
    // before the first gradient call: the value at the seed point is the last of these
    std::lock_guard<std::mutex> lock(mutex_);
    startObjective_ = objective;
    return false;
  }
  const long run = sinceGradient_.fetch_add(1) + 1;
  const bool hesse = inParallel || run >= kHesseRun;
  if (gradients == 1) {
    double start, edm;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      start = startObjective_;
      edm = seed_.edm;
    }
    if (!hesse && objective < start) improved_.store(true);
    // The line search is over and found nothing; a start that is already a minimum (estimated distance below
    // 1e-3 up) is left to Minuit2's HESSE, which is then the verification of the minimum.
    if (hesse && !improved_.load() && edm > 1e-3 * up_) {
      stalled_.store(true);
      stop_.store(true);
    }
  }
  return hesse && run >= kHesseRun && !hesseSeen_.exchange(true);
}

void MigradSeed::Forward(const std::atomic<bool> *request) {
  if (request && request->load(std::memory_order_relaxed)) stop_.store(true);
}

MigradSeed::Bounded MigradSeed::Seed() {
  std::lock_guard<std::mutex> lock(mutex_);
  return seed_;
}
