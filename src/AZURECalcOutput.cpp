// AZURECalc, the files written during a fit: the intermediate output of a parameter vector (WriteIterationOutput:
// the model on copies of its own in the caller's FitterStates, then the data files, param.fit and parameters.out)
// and param.fit itself (WriteParameters).  One unit because these are the fitter's only writes to disk, and they
// must not disturb the objects the minimizer is stepping.
#include "AZURECalc.h"
#include "ModelPass.h"
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include "Config.h"
#include "CompoundSystem.h"
#include "LevelReport.h"
#include "WorkingCopy.h"
#include "AZUREParams.h"

void AZURECalc::WriteIterationOutput(const vector_r &p) const {
  // Work on private copies: this runs in the middle of a fit and must not
  // disturb the objects the minimizer is stepping.  The copies are made once
  // and reused (every snapshot refills them from p below, exactly as the pooled
  // objects are refilled in operator()), instead of being cloned and deleted at
  // every snapshot -- see FitterStates in AZURECalc.h.  The lock also keeps two
  // threads from writing the snapshot files at the same time.
  std::lock_guard<std::mutex> outputLock(states_.outputMutex);
  if (!states_.output) states_.output = std::make_unique<LevelState>(state_);
  if (!states_.outputData) states_.outputData.reset(data()->Clone());
  const CompoundSystem *lc = compound();
  LevelState &ls = *states_.output;
  WorkingCopy *ld = states_.outputData.get();
  azure2::Results &results = states_.outputResults;

  try {
    lc->FillCompoundFromParams(p, ls);
    ld->FillNormsFromParams(p);
    ld->FillEnergyShiftsFromParams(p, results, lc, &ls, &inputs_, &configure());
    if (configure().paramMask & Config::USE_BRUNE_FORMALISM) lc->CalcShiftFunctions(configure(), ls);

    azure2::EvaluateModel(*ld, lc, ls, configure(), results, azure2::LeafValues::Calculate);

    AZUREParams params;
    lc->FillMnParams(params, ls, inputs_, &configure());
    ld->FillMnParams(params);
    // FillMnParams reads the levels' input values (GetE/GetGamma), not the
    // point being evaluated; copy p in so param.fit matches this snapshot.
    ROOT::Minuit2::MnUserParameters &mp = params.GetMinuitParams();
    if (mp.Params().size() == p.size())
      for (unsigned int i = 0; i < p.size(); i++) mp.SetValue(i, p[i]);
    ld->SyncShiftGroupParams(mp);
    ld->SyncNormGroupParams(mp);
    WriteParameters(params, configure());
    ld->WriteOutputFiles(configure(), results, true);
    WriteParametersOut(*lc, ls, lc->TransformOut(configure(), ls, inputs_), configure());
  } catch (...) {
    // An intermediate snapshot is a convenience, never a reason to abort a fit.
    // A failure may leave the copies half-filled: drop them so the next
    // snapshot starts from fresh clones.
    states_.output.reset();
    states_.outputData.reset();
    states_.outputResults.slot.clear();
  }
}

/*!
 * Write parameters to file
 */
void AZURECalc::WriteParameters(AZUREParams &params, const Config &configure) const {
  char filename[256];
  snprintf(filename, sizeof(filename), "%sparam.fit", configure.outputdir.c_str());
  std::ofstream out;
  out.open(configure.OpenPath(filename));
  if (out) {
    out.precision(7);
    // The same first line as param.par/param.sav: which amplitudes the file
    // holds.  Before 2026-11-04 param.fit had none, so a --use-park fit
    // resumed from it took its Park amplitudes for Brune's and converted
    // them again (every amplitude times sqrt(J), fixed ones included;
    // 11B+a 11-04-26_park_gpSM, chi2 10733 -> 49532 after three resumes).
    out << std::setw(20) << "#parametrization"
        << std::scientific << std::setw(20) << (double)AZUREParams::BasisForConfig(configure)
        << std::scientific << std::setw(20) << 0.0 << std::endl;
    for (int i = 0; i < params.GetMinuitParams().Params().size(); i++) {
      out << std::setw(20) << params.GetMinuitParams().GetName(i)
          << std::scientific << std::setw(20) << params.GetMinuitParams().Value(i)
          << std::scientific << std::setw(20) << params.GetMinuitParams().Error(i) << std::endl;
    }
    out.flush();
    out.close();
  } else
    configure.outStream << "Could not save param.fit file." << std::endl;
}
