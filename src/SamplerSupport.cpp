// The sampler's side of the support of the posterior (docs/upgrade/objective.md section 6, "The support";
// PosteriorSupport.h, the rule): the fit's vector in both sampling bases, and the walkers' start drawn inside it.
#include <limits>
#include <random>
#include <vector>

#include "AZURECalcMCMC.h"
#include "AZUREParams.h"
#include "CNuc.h"
#include "Config.h"
#include "WorkingCopy.h"
#include "WorkingData.h"

bool AZURECalcMCMC::InsideFitLimits(const vector_r &fitVector) const { return support_.Contains(fitVector); }

vector_r AZURECalcMCMC::FitVectorOfObserved(CNuc *transformed, WorkingCopy *copy, const vector_r &observedFull) const {
  AZUREParams params;
  transformed->FillMnParams(params, &configure());
  copy->FillMnParams(params);
  // The tail of the observed vector (normalizations, energy shifts, every data parameter) is the fit's.
  vector_r fitVector = params.GetMinuitParams().Params();
  for (size_t i = (size_t)copy->Data().GetNormParamOffset(); i < fitVector.size() && i < observedFull.size(); i++)
    fitVector[i] = observedFull[i];
  return fitVector;
}

bool AZURECalcMCMC::StartInsideFitLimits(const std::vector<double> &sampled, bool useRWA) const {
  if (useRWA) return InsideFitLimits(ReconstructFullParameters(sampled));
  CNuc *c = GetPooledCNuc();
  WorkingData *working = GetPooledData();
  bool inside = true;  // a conversion that fails: outside no limit (its posterior is zero, as before)
  try {
    const vector_r observed = ReconstructFullParametersPhysical(sampled);
    c->FillCompoundFromParamsPhysical(observed);
    if (c->TransformIn(configure())) inside = InsideFitLimits(FitVectorOfObserved(c, working->data.get(), observed));
  } catch (...) {
  }
  ReturnPooledCNuc(c);
  ReturnPooledData(working);
  return inside;
}

std::vector<std::vector<double>> AZURECalcMCMC::DrawInitialEnsemble(int nwalkers, const std::vector<double> &centre,
                                                                    const std::vector<double> &spread,
                                                                    std::default_random_engine &generator,
                                                                    bool useRWA) const {
  // A walker drawn outside the fit's limits is drawn again, up to kMaxDraws times; a draw inside the limits uses the
  // random stream as a ball without limits does, so an ensemble none of whose walkers fell outside is that ball.
  const int kMaxDraws = 1000;
  int nRedrawn = 0, nOutside = 0;
  std::vector<std::vector<double>> positions;
  for (int i = 0; i < nwalkers; i++) {
    std::vector<double> pos;
    int draws = 0;
    do {
      pos = centre;
      for (size_t j = 0; j < pos.size(); j++) {
        std::normal_distribution<double> init_distribution(0.0, spread[j]);
        pos[j] += init_distribution(generator);
      }
    } while (++draws < kMaxDraws && !StartInsideFitLimits(pos, useRWA));
    if (draws > 1) nRedrawn++;
    if (draws == kMaxDraws && !StartInsideFitLimits(pos, useRWA)) nOutside++;
    positions.push_back(pos);
  }
  if (nRedrawn > 0)
    configure().outStream << "Initial spread: " << nRedrawn << " walker" << (nRedrawn == 1 ? "" : "s")
                          << " drawn again inside the fit's limits\n";
  if (nOutside > 0)
    configure().outStream << "Warning: " << nOutside << " walker" << (nOutside == 1 ? "" : "s") << " still outside "
                          << "the fit's limits after " << kMaxDraws << " draws (zero posterior probability); "
                          << "reduce the initial spread or start inside the limits.\n";
  return positions;
}
