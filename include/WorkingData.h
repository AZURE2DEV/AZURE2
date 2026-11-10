#ifndef AZURE2_WORKINGDATA_H
#define AZURE2_WORKINGDATA_H

#include <memory>

#include "Results.h"

class EData;

/*!
 * A working copy of the data set, as the pools of AZURECalc and AZURECalcMCMC hand them to an evaluation, with the
 * buffer that evaluation writes (Results.h, docs/upgrade/data_model.md S3).  The two travel together, so a copy never
 * writes into another owner's values.  The buffer starts empty and is sized by the copy's first hand-over
 * (EData::BeginEvaluation).  Goes with the pools (Phase 6).
 */
struct WorkingData {
  std::unique_ptr<EData> data;
  azure2::Results results;
};

#endif
