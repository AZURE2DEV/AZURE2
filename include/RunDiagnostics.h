#ifndef AZURE2_RUNDIAGNOSTICS_H
#define AZURE2_RUNDIAGNOSTICS_H

#include <atomic>
#include <mutex>

#include "CoulombTable.h"
#include "Profile.h"

namespace azure2 {

/*!
 * What the run's diagnostics count and accumulate (docs/upgrade/evaluation.md section 8.4): the profile accounting
 * (AZURE_PROFILE), the Coulomb table's lookups, and the call counters and worst values of the printing and checking
 * switches (azure2/support/Environment.h).  Owned by the run's Config (Config::diagnostics), as its scratch is, and
 * shared by the Config's copies, so two runs in one process -- two pyazr sessions, the web UI's server -- count apart;
 * they were process-wide before 2026-11-13.  Nothing a calculation returns reads it.
 *
 * Add a counter: a member here, written where the diagnostic prints (CLAUDE.md, recipe "Add a diagnostic switch").
 */
struct RunDiagnostics {
  RunDiagnostics() = default;
  /// Prints the Coulomb table's use (when it was used) on standard error: the run ends.  The CLI's run ends with its
  /// process, where this line was printed before.
  ~RunDiagnostics();
  RunDiagnostics(const RunDiagnostics &) = delete;
  RunDiagnostics &operator=(const RunDiagnostics &) = delete;

  /// AZURE_PROFILE: stage times and counters (Profile.h); Config::profile points here when the switch is on.
  AZUREProfile::Account profile;
  /// The Coulomb table's lookups by the run's channel functions (PairFunctionSettings::tableUse points here).
  CoulombTableUse coulombTable;

  // The call counters of the printing switches: each prints its first calls and then every 50th.
  std::atomic<int> evaluationTimings{0};   ///< AZURE_GRAD_TIMING, one chi-squared evaluation (AZURECalc)
  std::atomic<int> jacobianPathTimings{0}; ///< AZURE_GRAD_TIMING, a gradient by the Jacobian path
  std::atomic<int> gradientTimings{0};     ///< AZURE_GRAD_TIMING or AZURE_GRAD_CHECK, an adjoint gradient
  std::atomic<int> g2Timings{0};           ///< AZURE_GRAD_TIMING, Minuit's G2
  std::atomic<int> residualJacobianTimings{0};  ///< AZURE_GRAD_TIMING, a residual Jacobian (AZURECalc)
  std::atomic<int> jacobianRowTimings{0};  ///< AZURE_GRAD_TIMING, the Jacobian rows (AZUREGrad)
  std::atomic<long long> foldWeightNs{0};  ///< AZURE_GRAD_TIMING, the target-effect combiner weights (ns, all threads)
  std::atomic<int> gridPlacements{0};      ///< AZURE_GRID_TIMING, a re-placement of the sub-point grids
  std::atomic<int> nonFiniteDumps{0};      ///< AZURE_DUMP_NONFINITE, dumps written (at most 20)
  std::atomic<int> poolChecks{0};          ///< AZURE_POOL_CHECK, comparisons made (the first five)
  std::atomic<int> foldWeightChecks{0};    ///< AZURE_TE_WEIGHTS_CHECK, points checked
  /// AZURE_GRAD_CHECK: the comparison of the Jacobian path with the adjoint, and the check against differences, run
  /// once per run.
  std::atomic<bool> jacobianPathChecked{false};
  std::atomic<bool> gradientChecked{false};
  /// AZURE_CHANNEL_VERIFY: the largest relative difference reported so far.
  std::mutex channelVerifyMutex;
  double channelVerifyWorst = 0.0;
};

}  // namespace azure2

#endif
