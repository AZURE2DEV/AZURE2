#include "TaskScratch.h"

#include <stdexcept>
#include <string>

#include "AMatrixFunc.h"
#include "THMMatrixFunc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace azure2 {

TaskScratch::TaskScratch() = default;
TaskScratch::~TaskScratch() = default;

AMatrixFunc &TaskScratch::PointSolver(CNuc *compound, const LevelState &state, const Config &configure) {
  if (!point_) {
    point_.reset(new AMatrixFunc(compound, state, configure));
    point_->ShareChannelMemo(&channel_);
  }
  point_->Reset(compound, state, configure);
  return *point_;
}

THMMatrixFunc &TaskScratch::PointThmSolver(CNuc *compound, const LevelState &state, const Config &configure) {
  if (!pointThm_) {
    pointThm_.reset(new THMMatrixFunc(compound, state, configure));
    pointThm_->ShareChannelMemo(&channel_);
  }
  pointThm_->Reset(compound, state, configure);
  return *pointThm_;
}

AMatrixFunc &TaskScratch::GradientSolver(CNuc *compound, const LevelState &state, const Config &configure) {
  if (!gradient_) {
    gradient_.reset(new AMatrixFunc(compound, state, configure));
    gradient_->ShareChannelMemo(&channel_);
  }
  gradient_->Reset(compound, state, configure);
  return *gradient_;
}

THMMatrixFunc &TaskScratch::GradientThmSolver(CNuc *compound, const LevelState &state, const Config &configure) {
  if (!gradientThm_) {
    gradientThm_.reset(new THMMatrixFunc(compound, state, configure));
    gradientThm_->ShareChannelMemo(&channel_);
  }
  return *gradientThm_;
}

TaskScratchSet::TaskScratchSet() : slots_(kMaxThreads) {}
TaskScratchSet::~TaskScratchSet() = default;
TaskScratchSet::TaskScratchSet(const TaskScratchSet &) : slots_(kMaxThreads) {}
TaskScratchSet &TaskScratchSet::operator=(const TaskScratchSet &) { return *this; }

TaskScratch &TaskScratchSet::ForThisThread() {
#ifdef _OPENMP
  const int t = omp_get_thread_num();
#else
  const int t = 0;
#endif
  if (t < 0 || t >= kMaxThreads)
    throw std::logic_error("TaskScratchSet: thread " + std::to_string(t) + " beyond " + std::to_string(kMaxThreads));
  std::unique_ptr<TaskScratch> &slot = slots_[t];
  if (!slot) slot.reset(new TaskScratch());
  return *slot;
}

}  // namespace azure2
