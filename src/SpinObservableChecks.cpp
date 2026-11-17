// Segment::SpinObservableAllowed (Segment.h): whether the general spin observable of a data set can be computed for
// its pairs, asked by the segment readers of WorkingCopy.cpp beside A_y's own check.
#include <iostream>
#include <string>

#include "CompoundSystem.h"
#include "Config.h"
#include "PPair.h"
#include "Segment.h"

bool Segment::SpinObservableAllowed(const CompoundSystem *compound, int number, bool isTest,
                                    const Config &configure) const {
  const azure2::physics::SpinObservable *observable = GetSpinObservable();
  if (!observable) return true;
  // A general spin observable: particle pairs on both sides, and every particle its operator names with the spin
  // the operator needs (physics/SpinObservable.h, Check).
  const char *what = isTest ? "test segment #" : "segment #";
  const PPair *entrancePair = compound->GetPair(compound->GetPairNumFromKey(GetEntranceKey()));
  const PPair *exitPair =
      compound->IsPairKey(GetExitKey()) ? compound->GetPair(compound->GetPairNumFromKey(GetExitKey())) : nullptr;
  std::string why;
  if (!exitPair)
    why = "needs one entrance and one exit pair";
  else if (entrancePair->GetPType() != 0 || exitPair->GetPType() != 0)
    why = "needs particle pairs on both sides (for capture only A_y, quantity analyzing_power, is computed)";
  else
    observable->Check(entrancePair->GetJ(1), entrancePair->GetJ(2), exitPair->GetJ(1), exitPair->GetJ(2), why);
  if (why.empty()) return true;
  configure.outStream << "ERROR: " << what << number << ": spin observable \"" << observable->token << "\": " << why
                      << "." << std::endl;
  return false;
}
