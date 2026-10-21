#ifndef AZURE_PROFILE_H
#define AZURE_PROFILE_H

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iosfwd>
#include <string>

/// Per-stage wall-clock accounting of the calculation, switched on with AZURE_PROFILE=1.

/*!
 * Each stage accumulates the time spent inside it, summed over threads, until
 * the next Report().  With AZURE_PROFILE unset a Timer costs one test of a
 * cached flag; nothing is read from the clock.
 *
 * Stages nest: the level-matrix stages (fill, factor, T-matrix) run inside
 * the point stage, and the channel solves inside the T-matrix stage, so the
 * report lists inclusive times.
 */
namespace AZUREProfile {

enum Stage {
  kEvaluation,       ///< one chi-squared evaluation (AZURECalc::operator())
  kPoints,           ///< the parallel pass over the points of an evaluation
  kEDependent,       ///< Coulomb functions, penetrabilities, shifts (EPoint::CalcEDependentValues)
  kFill,             ///< level-matrix (or channel-matrix) fill
  kFactor,           ///< factorization / channel-space solve
  kTMatrix,          ///< T-matrix from the solved matrices
  kChannelSolve,     ///< level-space solves A g_c (inside kTMatrix and the gradient)
  kCrossSection,     ///< angular and cross-section sums from the T-matrix
  kGradient,         ///< one gradient or residual-Jacobian pass
  kGradForward,      ///< forward solves inside the gradient
  kGradTBar,         ///< observable -> T-matrix cotangents (PointTBar)
  kGradBackprop,     ///< T-matrix cotangents -> parameters (TBarAdjoint, bilinear basis)
  kGradLevelMatrix,  ///< structural level-matrix adjoint (inside kGradBackprop)
  kGradRows,         ///< per-member rows: basis combination, scatter into the Jacobian
  kNumStages
};

/// Events counted alongside the stages: the route each J-group's level matrix takes at an energy.
enum Counter {
  kChannelRoute,       ///< solved in channel space
  kLevelRouteSize,     ///< level space: no fewer coupled channels than levels
  kLevelRouteSetup,    ///< level space: no formal representation (eigen-decomposition, RMC with Brune/Park)
  kLevelRoutePole,     ///< level space: an energy too close to a formal pole
  kLevelRoutePivot,    ///< level space: 1 - R Lambda ill-conditioned
  kLevelRouteForced,   ///< level space: AZURE_LEVEL_SPACE or SetLevelSpaceOnly
  kDeflated,           ///< channel space with levels near a pole deflated
  kNumCounters
};

/// True when AZURE_PROFILE is set to anything but "0" (read once).
inline bool Enabled() {
  static const bool on = [] {
    const char *e = std::getenv("AZURE_PROFILE");
    return e && *e && std::string(e) != "0";
  }();
  return on;
}

/// Adds `ns` nanoseconds to a stage (and one call).
void Add(Stage stage, long long ns);

/// Counts one event (only when profiling).
void CountEvent(Counter counter);
inline void Count(Counter counter) {
  if (Enabled()) CountEvent(counter);
}

/// Writes the stages accumulated since the last report under `label`, then clears them.
void Report(std::ostream &out, const char *label);

/// Times the enclosing scope into one stage.
class Timer {
 public:
  explicit Timer(Stage stage) : stage_(stage), on_(Enabled()) {
    if (on_) start_ = std::chrono::steady_clock::now();
  }
  ~Timer() { Stop(); }
  /// Ends the timing before the end of the scope.
  void Stop() {
    if (on_)
      Add(stage_, std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start_)
                      .count());
    on_ = false;
  }
  Timer(const Timer &) = delete;
  Timer &operator=(const Timer &) = delete;

 private:
  Stage stage_;
  bool on_;
  std::chrono::steady_clock::time_point start_;
};

/*!
 * Reports (to std::cerr) when the enclosing scope ends.  Declared before the
 * scope's Timer, it runs after the Timer has added its time.  A report inside
 * another one's scope (an evaluation inside a gradient) is left to the outer one.
 */
class ScopedReport {
 public:
  explicit ScopedReport(const char *label);
  ~ScopedReport();
  ScopedReport(const ScopedReport &) = delete;
  ScopedReport &operator=(const ScopedReport &) = delete;

 private:
  const char *label_;
  bool outer_;
};

}  // namespace AZUREProfile

#endif
