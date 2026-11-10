#ifndef EXPERIMENT_ADAPTER_H
#define EXPERIMENT_ADAPTER_H

#include "PPair.h"
#include "azure2/experiment/Frames.h"

/// A pair as experiment/'s frame functions read it: masses (u), separation and excitation energy (MeV).
inline azure2::experiment::Partition PartitionOf(const PPair *p) {
  azure2::experiment::Partition r;
  r.m1 = p->GetM(1);
  r.m2 = p->GetM(2);
  r.sepE = p->GetSepE();
  r.exE = p->GetExE();
  return r;
}

#endif
