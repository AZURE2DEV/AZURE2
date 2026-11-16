#ifndef AZURE2_RUNSCRATCH_H
#define AZURE2_RUNSCRATCH_H

#include <map>
#include <memory>
#include <mutex>
#include <thread>
#include <tuple>
#include <vector>

#include "ECAmplitudeCache.h"
#include "ECIntegral.h"
#include "ECTable.h"
#include "LineWindow.h"
#include "ThmVertexBoundary.h"

namespace azure2 {

/*!
 * Park's numerical shift-function derivatives at level energies (CompoundSystem::CalcShiftFunctions) by their exact arguments:
 * the pair's parameters, l, 2J and the level energy, compared bit for bit, and the settings that change the shift
 * function at those arguments (GSL or own Coulomb functions, a nuclear potential, the Coulomb table, relativistic
 * kinematics).  The pooled copies of a fit (evaluation, gradient, residual Jacobian) need the same values at the
 * same parameters.  An entry is what the computation returns, so a hit changes nothing.  The Coulomb table is fixed
 * within a run (PPair::FunctionSettings), so the key holds whether it is used, not its file.  Cleared when it grows
 * past kMax.
 */
struct ShiftDerivativeMemo {
  struct Key {
    int z1, z2, l, two_j;  ///< two_j: PPair::SpinOrbitTwoJ
    bool gsl;
    std::vector<double> hybrid;  ///< azure2::physics::PotentialKey of the pair's potential
    bool table, relativistic;
    double redmass, mass, radius, threshold, energy;  ///< u, u, fm, MeV, MeV
    bool operator<(const Key &o) const {
      return std::tie(z1, z2, l, two_j, gsl, hybrid, table, relativistic, redmass, mass, radius, threshold, energy) <
          std::tie(o.z1, o.z2, o.l, o.two_j, o.gsl, o.hybrid, o.table, o.relativistic, o.redmass, o.mass, o.radius,
                   o.threshold, o.energy);
    }
  };
  struct Value {
    bool has1 = false, has2 = false;
    double d1 = 0.0, d2 = 0.0;  ///< dS/dE (1/MeV), d2S/dE2 (1/MeV^2)
  };
  static constexpr std::size_t kMax = 1 << 16;
  std::mutex mutex;
  std::map<Key, Value> values;
};

/// What one thread of a run keeps between calls that have no energy task's scratch at hand (set-up, hand-over,
/// external-capture integrals, the THM vertex boundary): values keyed by everything they depend on.
struct ThreadScratch {
  ECScratch ec;
  ThmShiftMemo thmShift;
};

/*!
 * The run's scratch (docs/upgrade/evaluation.md: no hidden state): per thread, and tables shared by the threads.
 * Owned by the run's Config (Config::scratch), so two runs in one process -- two pyazr sessions, the web UI's
 * server -- never share a value, and released with the run.  A thread's slot is made on its first request and used
 * by that thread only; the lookup takes a lock, so callers take the slot once per object (an integral, a vertex
 * boundary), not per value.
 */
class RunScratch {
 public:
  ThreadScratch &ForThisThread();

  /// Shared by the run's threads, each behind its own lock: the external-capture integrals' energy tables
  /// (ECTable.h) and Park's shift-function derivatives at level energies.
  ECTermStore ecTerms;
  ShiftDerivativeMemo shiftDerivatives;
  /// The quadrature nodes of the line-shape integration of a broad capture final state (LineWindow.h; its own lock).
  LineNodeStore lineNodes;
  /// The external-capture amplitudes at the energies of each segment's points, interpolated at a shifted energy with
  /// --no-ec-table (ECAmplitudeCache.h; its own lock).  Null until the run's set-up starts it
  /// (InitializeECAmplitudeCache, single-threaded); a site then returns its own amplitudes.
  std::unique_ptr<ECAmplitudeCache> ecAmplitudes;

 private:
  std::mutex mutex_;
  std::map<std::thread::id, std::unique_ptr<ThreadScratch>> threads_;
};

}  // namespace azure2

#endif
