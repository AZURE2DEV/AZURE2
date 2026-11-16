#ifndef AZURECALCMCMCRUN_H
#define AZURECALCMCMCRUN_H

// The state of one run of the sampler and the one reader of a number in its files: a private header of
// src/AZURECalcMCMCRun.cpp (which makes the Run and drives the chain) and src/AZURECalcMCMCChain.cpp (which
// writes each step from it).  It was in AZURECalcMCMC.cpp while these were one file.

#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include "AZURECalcMCMC.h"

// The state of one run of RunMCMCSampling, owned by that call and reached by the callbacks it hands the ensemble
// sampler (lambdas capturing it): the sample file, the step numbering and the best walker of the last step.
struct AZURECalcMCMC::Run {
  std::ofstream *sampleFile = nullptr;
  int nwalkers = 0;
  int currentStep = 0;  // absolute: steps of earlier runs of a resumed chain included
  int stepOffset = 0;   // the step this run started from, so a resumed run keeps numbering the chain
  bool useRWA = true;
  // Best walker of the most recent completed step, for the progress callback and the snapshot files
  double lastLikelihood = 0.0;
  double lastPrior = 0.0;
  double lastLogprob = -std::numeric_limits<double>::infinity();
  std::vector<double> lastParams;
};

// Custom robust string-to-double parser for scientific notation (src/AZURECalcMCMCChain.cpp)
double robust_stod(const std::string &str);

#endif
