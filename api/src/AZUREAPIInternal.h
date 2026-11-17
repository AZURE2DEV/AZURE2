#ifndef AZUREAPIINTERNAL_H
#define AZUREAPIINTERNAL_H

// Two helpers the files of AZUREAPI share (AZUREAPIReports.cpp, AZUREAPIChi2.cpp, AZUREAPIGradient.cpp): a private
// header of api/src, not installed.  They were file-local functions of AZUREAPI.cpp while it was one file.

#include <vector>
#include "CompoundSystem.h"
#include "Config.h"
#include "Constants.h"
#include "LevelState.h"
#include "Segment.h"
#include "WorkingCopy.h"

// The points of an API calculation in one parallel pass, as the CLI's
// AZURECalc::PrecalculatePoints (before 2026-10-30 the API calculated them one
// at a time).  Used only when every point calculated: after a failure the
// caller repeats the calculation point by point, so a failed point keeps the
// serial path's handling (a NaN slot in UpdateSegments*, 0 or a skipped
// component in a chi-squared) and the results equal the serial path's.
inline bool PrecalculateAll(WorkingCopy &data, const std::vector<const Segment *> &segments, const CompoundSystem *lc,
                            const LevelState &state, const Config &cfg, azure2::Results &results) {
  return data.CalculatePointsInParallel(segments, lc, state, cfg, results) == 0;
}

// Map a packed non-fixed RWA parameter vector into the full all_rwa_ layout.
inline vector_r MapPackedToFull(const vector_r &packed, const vector_r &all_rwa,
                                const std::vector<bool> &fixed) {
  vector_r full = all_rwa;
  int k = 0;
  for (int i = 0; i < (int)all_rwa.size(); ++i)
    if (!fixed[i] && k < (int)packed.size()) full[i] = packed[k++];
  return full;
}

#endif
