#include "RunScratch.h"

// Here, not in ECIntegral.cpp: Config makes the run's scratch, and a program that links Config alone must not pull
// in the integrals.
ECScratch::ECScratch() : workspace(gsl_integration_workspace_alloc(1000)) {}
ECScratch::~ECScratch() {
  if (workspace) gsl_integration_workspace_free(workspace);
}

namespace azure2 {

ThreadScratch &RunScratch::ForThisThread() {
  std::lock_guard<std::mutex> lock(mutex_);
  std::unique_ptr<ThreadScratch> &slot = threads_[std::this_thread::get_id()];
  if (!slot) slot.reset(new ThreadScratch());
  return *slot;
}

}  // namespace azure2
