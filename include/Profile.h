#ifndef AZURE_PROFILE_H
#define AZURE_PROFILE_H

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iosfwd>
#include <string>

#include "azure2/support/Environment.h"

/// Per-stage wall-clock accounting of the calculation, switched on with AZURE_PROFILE=1.

/*!
 * Each stage accumulates the time spent inside it, summed over threads, until
 * the next Report(), in the run's Account.  With AZURE_PROFILE unset a Timer
 * costs one test of a null pointer; nothing is read from the clock.
 *
 * Stages nest: the level-matrix stages (fill, factor, T-matrix) run inside
 * the point stage, and the channel solves inside the T-matrix stage, so the
 * report lists inclusive times.
 */
namespace AZUREProfile {

enum Stage {
  kEvaluation,       ///< one chi-squared evaluation (AZURECalc::operator())
  kPoints,           ///< the parallel pass over the points of an evaluation
  kEDependent,       ///< Coulomb functions, penetrabilities, shifts (Unit::CalcEDependentValues)
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
  kShiftLocal,       ///< energy-shift Jacobian columns differenced on their own segments' rows
  kShiftFull,        ///< energy-shift / density Jacobian columns differenced on every row (forward only)
  kEnergyTaskSetup,  ///< the energy tasks of a pass: checked against the units, rebuilt when they changed (serial)
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
  kWhitMemoHit,        ///< final-state Whittaker of ECIntegral's data path from the memo (WhitFunc::ScaledMemo)
  kWhitMemoMiss,       ///< ... computed and stored
  kEnergyTasks,        ///< energy tasks of the point pass: one level-matrix fill and factorisation each (CalculationPlan.h)
  kEnergyUnits,        ///< units (points and sub-points) calculated by those tasks
  kEnergyTaskBuilds,   ///< passes whose energy tasks were rebuilt (a unit, its E_x, an amplitude source or the work changed)
  kNumCounters
};

/// True when AZURE_PROFILE is set to anything but "0" (the process's environment record, read once).
inline bool Enabled() {
  static const bool on = azure2::support::ProcessEnvironment().On(azure2::support::Switch::kProfile);
  return on;
}

/*!
 * The accounting of one run: stage times and calls, event counts, and the nesting depth of the reports.  Owned by the
 * run (RunDiagnostics::profile); Config::profile points to it when AZURE_PROFILE is on and is null otherwise, so a
 * Timer of a run without profiling reads no clock and writes nothing.
 */
class Account {
 public:
  /// Adds `ns` nanoseconds to a stage (and one call).
  void Add(Stage stage, long long ns) {
    ns_[stage].fetch_add(ns, std::memory_order_relaxed);
    calls_[stage].fetch_add(1, std::memory_order_relaxed);
  }
  /// Counts `n` events.
  void CountEvents(Counter counter, long long n) { count_[counter].fetch_add(n, std::memory_order_relaxed); }
  /// Writes the stages accumulated since the last report under `label`, then clears them.
  void Report(std::ostream &out, const char *label);
  /// ScopedReport's nesting: true for the outermost report.
  bool Enter() { return depth_.fetch_add(1) == 0; }
  void Leave() { depth_.fetch_sub(1); }

 private:
  std::atomic<long long> ns_[kNumStages] = {};
  std::atomic<long long> calls_[kNumStages] = {};
  std::atomic<long long> count_[kNumCounters] = {};
  std::atomic<int> depth_{0};
};

/// Counts one event in `account` (null: not profiling).
inline void Count(Account *account, Counter counter) {
  if (account) account->CountEvents(counter, 1);
}

/// Counts `n` events at once in `account` (null: not profiling).
inline void CountMany(Account *account, Counter counter, long long n) {
  if (account) account->CountEvents(counter, n);
}

/// Writes `account`'s stages under `label` (nothing when null).
inline void Report(Account *account, std::ostream &out, const char *label) {
  if (account) account->Report(out, label);
}

/// Times the enclosing scope into one stage of `account` (null: nothing).
class Timer {
 public:
  Timer(Account *account, Stage stage) : account_(account), stage_(stage) {
    if (account_) start_ = std::chrono::steady_clock::now();
  }
  ~Timer() { Stop(); }
  /// Ends the timing before the end of the scope.
  void Stop() {
    if (account_)
      account_->Add(stage_, std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() -
                                                                                  start_).count());
    account_ = nullptr;
  }
  Timer(const Timer &) = delete;
  Timer &operator=(const Timer &) = delete;

 private:
  Account *account_;
  Stage stage_;
  std::chrono::steady_clock::time_point start_;
};

/*!
 * Reports (to std::cerr) when the enclosing scope ends.  Declared before the
 * scope's Timer, it runs after the Timer has added its time.  A report inside
 * another one's scope (an evaluation inside a gradient) is left to the outer one.
 */
class ScopedReport {
 public:
  ScopedReport(Account *account, const char *label);
  ~ScopedReport();
  ScopedReport(const ScopedReport &) = delete;
  ScopedReport &operator=(const ScopedReport &) = delete;

 private:
  Account *account_;
  const char *label_;
  bool outer_;
};

}  // namespace AZUREProfile

#endif
