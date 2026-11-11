// The capture analyzing-power requests of a run, read from its data sets before the pathway table is built.
#include "CaptureAyRequests.h"

#include "CNuc.h"
#include "EData.h"
#include "ESegment.h"
#include "PPair.h"

void RequestCaptureAyTerms(CNuc &compound, EData &data) {
  for (int i = 1; i <= data.NumSegments(); i++) {
    const ESegment &segment = *data.GetSegment(i);
    if (!segment.IsAnalyzingPower()) continue;
    for (int t = 0; t < segment.NumTerms(); t++) {
      const int entranceKey = segment.TermEntranceKey(t), exitKey = segment.TermExitKey(t);
      if (!compound.IsPairKey(entranceKey) || !compound.IsPairKey(exitKey)) continue;
      const int exitPair = compound.GetPairNumFromKey(exitKey);
      if (compound.GetPair(exitPair)->GetPType() != 10) continue;
      compound.RequestCaptureAy(compound.GetPairNumFromKey(entranceKey), exitPair);
    }
  }
}
