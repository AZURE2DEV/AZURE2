// AZURECalc, the idle states: the level states and working copies (each with its buffer) that evaluations running
// at the same time take from and give back to the caller's FitterStates (AZURECalc.h; docs/upgrade/evaluation.md
// section 7.7).  One unit because this is the only place a fit's evaluation gets its own copy: a copy of the
// caller's state is made when none is idle, so two evaluations never share one.
#include "AZURECalc.h"
#include <memory>
#include "WorkingCopy.h"

FitterStates::FitterStates() = default;
FitterStates::~FitterStates() = default;

std::unique_ptr<LevelState> AZURECalc::GetPooledState() const {
  // None idle: one more evaluation runs at the same time as the others (a thread of a parallel gradient), on a copy of
  // the caller's state as it is now.
  return states_.idle.Take(state_);
}

/*!
 * Get a working copy of the data set and its buffer from the pool, creating new if pool is empty
 */
WorkingData *AZURECalc::GetPooledData() const {
  // None idle: as GetPooledState, on a copy of the caller's working copy.
  return states_.idleData
      .TakeOrMake([this] {
        return std::unique_ptr<WorkingData>(new WorkingData{std::unique_ptr<WorkingCopy>(data_->Clone()), {}});
      })
      .release();
}

void AZURECalc::ReturnPooledState(std::unique_ptr<LevelState> state) const { states_.idle.Give(std::move(state)); }

/*!
 * Return a working copy and its buffer to the pool for reuse
 */
void AZURECalc::ReturnPooledData(WorkingData *obj) const { states_.idleData.Give(std::unique_ptr<WorkingData>(obj)); }
