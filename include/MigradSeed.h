#ifndef MIGRADSEED_H
#define MIGRADSEED_H

// The start of a MIGRAD fit with the analytic gradient (docs/upgrade/fit.md, "The seed of MIGRAD").
//
// MIGRAD starts from the inverse Hessian diag(1 / G2) in Minuit2's internal coordinates and takes the step -g / G2.
// With the Gauss-Newton G2 of the residual Jacobian that first step has no bound for a parameter the data hardly
// constrain (G2 near zero) and for a limited parameter near its limit (Minuit2 forms the internal G2 as
// (dx/du)^2 G2 and leaves out g d2x/du2).  When the line search that follows finds no improvement, Minuit2 starts a
// numerical HESSE of n(n-1)/2 evaluations with no gradient call.  A MigradSeed, owned by the fitter for one MIGRAD,
// (1) raises the curvature the objective hands to Minuit2 so that the first step is bounded, and (2) watches the
// evaluations until the second gradient call: a first line search without improvement sets its stop request, MIGRAD
// returns, and the fitter starts again with Minuit2's numerical seed.

#include <atomic>
#include <mutex>
#include <vector>

namespace ROOT {
namespace Minuit2 {
class MnUserParameters;
}
}  // namespace ROOT

class MigradSeed {
 public:
  /// A parameter as MIGRAD sees it: free or not, its limits, and its initial step (Minuit's "error").
  struct Parameter {
    bool free = false, hasLower = false, hasUpper = false;
    double lower = 0.0, upper = 0.0, step = 0.0;
  };
  /// What BoundCurvature did: free parameters, how many curvatures were raised, the parameter whose unbounded step
  /// was the longest in units of its initial step (-1: none; `worstSteps` infinite for a zero curvature), and the
  /// estimated distance to the minimum (1/2) sum g^2 / G2 with the bounded curvatures.
  struct Bounded {
    int free = 0, raised = 0, worst = -1;
    double worstSteps = 0.0, edm = 0.0;
  };
  /// No parameter moves by more than kSteps of its initial step in the unit step from the seed.
  static constexpr double kSteps = 10.0;
  /// A parameter with two limits moves by at most this in Minuit2's internal angle.
  static constexpr double kInternal = 0.5;
  /// A run of this many evaluations without a gradient call is Minuit2's numerical HESSE (a line search takes at
  /// most about 15).
  static constexpr long kHesseRun = 40;

  /// Raises `g2` (never lowers it) where the diagonal Newton step grad / g2 has no bound: the term g d2x/du2 of the
  /// limit transformation where it is positive, kSteps, kInternal.  Every correction is proportional to the
  /// gradient, so nothing changes at a minimum.  Vectors of other lengths than `parameters` are left alone.
  static Bounded BoundCurvature(const std::vector<Parameter> &parameters, const std::vector<double> &p,
                                const std::vector<double> &grad, std::vector<double> &g2);

  /// `up`: the objective's error definition (the scale of "already at a minimum").
  MigradSeed(const ROOT::Minuit2::MnUserParameters &start, double up);
  MigradSeed(const std::vector<Parameter> &parameters, double up);

  /// The objective, from the residual Jacobian's gradient and G2 at `p`: bounds g2; the first call is the seed's.
  void Bound(const std::vector<double> &p, const std::vector<double> &grad, std::vector<double> &g2);
  /// The objective: a gradient was asked for (Minuit2's, not an evaluation inside one).
  void GradientCalled();
  /// The fitter's observer: an evaluation gave `objective`; `inParallel`: it was made inside a parallel region of
  /// the caller (Minuit2's numerical gradient).  Returns true once, when it recognizes a numerical HESSE.
  bool Evaluated(double objective, bool inParallel);
  /// Set the stop request as well when `request` is set (the caller's own stop).
  void Forward(const std::atomic<bool> *request);

  /// The stop request the objective obeys while this seed watches (EvaluationObserver::stop).
  const std::atomic<bool> *Stop() const { return &stop_; }
  /// The first line search found no improvement and the start is not a minimum: the stop request was set.
  bool Stalled() const { return stalled_.load(); }
  /// What the first Bound did (free = 0 before it).
  Bounded Seed();

 private:
  const std::vector<Parameter> parameters_;
  const double up_;
  std::mutex mutex_;  ///< seed_, startObjective_: evaluations arrive from several threads
  Bounded seed_;
  bool seeded_ = false;
  double startObjective_ = 0.0;
  std::atomic<long> gradients_{0}, sinceGradient_{0};
  std::atomic<bool> improved_{false}, stalled_{false}, stop_{false}, hesseSeen_{false};
};

#endif
