// The capture analyzing-power requests of a run, read from its data sets before the pathway table is built.
#include "CaptureAyRequests.h"

#include "CompoundSystem.h"
#include "CompoundSystemSetUp.h"
#include "WorkingCopy.h"
#include "Segment.h"
#include "PPair.h"

void RequestCaptureAyTerms(const CompoundSystem &compound, const WorkingCopy &data, CompoundSetUpRequests &requests) {
  for (int i = 1; i <= data.Data().NumSegments(); i++) {
    const Segment &segment = *data.Data().GetSegment(i);
    if (!segment.IsAnalyzingPower()) continue;
    for (int t = 0; t < segment.NumTerms(); t++) {
      const int entranceKey = segment.TermEntranceKey(t), exitKey = segment.TermExitKey(t);
      if (!compound.IsPairKey(entranceKey) || !compound.IsPairKey(exitKey)) continue;
      const int exitPair = compound.GetPairNumFromKey(exitKey);
      if (compound.GetPair(exitPair)->GetPType() != 10) continue;
      requests.AddCaptureAy(compound.GetPairNumFromKey(entranceKey), exitPair);
    }
  }
}
