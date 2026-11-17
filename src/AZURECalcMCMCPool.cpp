// AZURECalcMCMC, the idle states: the walkers' level states and working copies (each with its buffer) that
// evaluations running at the same time take from and give back to the caller's SamplerStates (AZURECalcMCMC.h;
// docs/upgrade/evaluation.md section 7.8).  One unit because this is the only place a walker's evaluation gets
// its own copy: a copy of the sampler's is made when none is idle.
#include "AZURECalcMCMC.h"
#include <memory>
#include "WorkingCopy.h"

/*!
 * Get a CompoundSystem object from the pool, creating new if pool is empty
 */
WalkerLevels *AZURECalcMCMC::GetPooledLevels() const {
  // None idle: a copy of the sampler's levels as they are now (one per walker evaluated at the same time).
  return states_.levels.Take(levels_).release();
}

/*!
 * Get a working copy of the data set and its buffer from the pool, creating new if pool is empty
 */
WorkingData *AZURECalcMCMC::GetPooledData() const {
  // None idle: one more walker is evaluated at the same time as the others, on a copy of the sampler's working copy.
  return states_.data
      .TakeOrMake([this] {
        return std::unique_ptr<WorkingData>(new WorkingData{std::unique_ptr<WorkingCopy>(data_->Clone()), {}});
      })
      .release();
}

/*!
 * Return a CompoundSystem object to the pool for reuse
 */
void AZURECalcMCMC::ReturnPooledLevels(WalkerLevels *walker) const {
  states_.levels.Give(std::unique_ptr<WalkerLevels>(walker));
}

/*!
 * Return a working copy and its buffer to the pool for reuse
 */
void AZURECalcMCMC::ReturnPooledData(WorkingData *obj) const { states_.data.Give(std::unique_ptr<WorkingData>(obj)); }
