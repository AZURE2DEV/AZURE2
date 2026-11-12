#ifndef AZURE2_RUNSCRATCH_H
#define AZURE2_RUNSCRATCH_H

#include <map>
#include <memory>
#include <mutex>
#include <thread>

#include "ECIntegral.h"
#include "ThmVertexBoundary.h"

namespace azure2 {

/// What one thread of a run keeps between calls that have no energy task's scratch at hand (set-up, hand-over,
/// external-capture integrals, the THM vertex boundary): values keyed by everything they depend on.
struct ThreadScratch {
  ECScratch ec;
  ThmShiftMemo thmShift;
};

/*!
 * The run's per-thread scratch (docs/upgrade/evaluation.md: no hidden state).  Owned by the run's Config
 * (Config::scratch), so two runs in one process -- two pyazr sessions, the web UI's server -- never share a value,
 * and released with the run.  A thread's slot is made on its first request and used by that thread only; the lookup
 * takes a lock, so callers take the slot once per object (an integral, a vertex boundary), not per value.
 */
class RunScratch {
 public:
  ThreadScratch &ForThisThread();

 private:
  std::mutex mutex_;
  std::map<std::thread::id, std::unique_ptr<ThreadScratch>> threads_;
};

}  // namespace azure2

#endif
