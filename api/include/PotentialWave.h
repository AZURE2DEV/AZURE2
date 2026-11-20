#ifndef AZURE2_API_POTENTIALWAVE_H
#define AZURE2_API_POTENTIALWAVE_H

#include "azure2/physics/HybridPotential.h"

class PPair;

/// The arguments of the physics/ functions of one partial wave (l, 2J) of a pair's potential: the pair, the wave's
/// potential (ForPartialWave), the 2J of its keys, <2 l.i> and the reduced mass (MeV).  One place for the session's calls
/// (AZUREAPI::GetPotentialBoxStates, GetPotentialBoundState, GetPotentialEffectiveRange) and the Qt program's (api/front
/// FrontPotential.cpp), so that both ask the same wave.
struct PotentialWaveArgs {
  const PPair *pair = nullptr;
  azure2::physics::HybridPotential potential;
  int twoJ = azure2::physics::kNoSpinOrbit;
  double ldoti = 0.0, mu = 0.0;
};

/// The wave (l, requestedTwoJ; < 0: none) of the potential of `pair`; throws std::invalid_argument, its message headed by
/// `call`, for a pair without potential, l < 0, a J the potential needs (it depends on J) or that l and the spin cannot
/// make.
PotentialWaveArgs PotentialWaveOf(const PPair &pair, int l, int requestedTwoJ, const char *call);

#endif
