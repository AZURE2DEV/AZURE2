// AZUREAPI: the channel functions of a pair on a requested energy grid (pyazr's coulomb_functions).
#include <cmath>

#include "AZUREAPI.h"
#include "CNuc.h"
#include "ChannelFunc.h"
#include "Config.h"
#include "CoulFunc.h"
#include "PPair.h"

/*!
 * Coulomb wave functions, penetrability, shift function and hard-sphere phase
 * on a requested energy grid; below threshold the shift function alone.
 *
 * These are ordinary R-matrix quantities that AZURE2 has always computed
 * internally; what is new is being able to ask for them.  The values follow the
 * run's own configuration, so the same call returns the accurate Coulomb
 * routine's answer, GSL's, or the Numerov solution through a nuclear potential,
 * according to how the calculation was set up -- which is how one sees what the
 * hybrid model actually does to the external region.
 */
vector_r AZUREAPI::GetCoulombFunctions(const vector_r &request) const {
  vector_r out;
  if (request.size() < 4) return out;

  int pairKey = static_cast<int>(std::lround(request[0]));
  int lValue = static_cast<int>(std::lround(request[1]));
  double radius = request[2];
  int nE = static_cast<int>(std::lround(request[3]));
  if (nE < 0 || request.size() < static_cast<size_t>(4 + nE)) return out;

  if (!compound_ || !compound_->IsPairKey(pairKey)) return out;
  PPair *pair = compound_->GetPair(compound_->GetPairNumFromKey(pairKey));
  if (!pair) return out;
  if (radius <= 0.0) radius = pair->GetChRad();  // the channel radius by default

  CoulFunc coul(pair, !!(configure().paramMask & Config::USE_GSL_COULOMB_FUNC));

  out.push_back(static_cast<double>(nE));
  for (int i = 0; i < nE; ++i) {
    double energy = request[4 + i];
    double F = 0., dF = 0., G = 0., dG = 0., P = 0., S = 0., delta = 0.;
    if (energy > 0.0 && !ChannelFunc::IsCoulombThreshold(pair, energy)) {
      try {
        CoulWaves w = coul(lValue, azure2::physics::kNoSpinOrbit, radius, energy);
        F = w.F, dF = w.dF, G = w.G, dG = w.dG;
        P = coul.Penetrability(lValue, azure2::physics::kNoSpinOrbit, radius, energy);
        S = coul.PEShift(lValue, azure2::physics::kNoSpinOrbit, radius, energy);
        delta = -std::atan2(w.F, w.G);
      } catch (...) {
        // A failed evaluation returns zeros for that energy rather than aborting the grid.
        F = dF = G = dG = P = S = delta = 0.;
      }
    } else if (radius == pair->GetChRad()) {  // closed (E <= 0, eta > 100): S through threshold (ChannelFunc)
      S = ChannelFunc(pair, !!(configure().paramMask & Config::USE_GSL_COULOMB_FUNC)).Shift(lValue, azure2::physics::kNoSpinOrbit, energy);
    }
    out.insert(out.end(), {F, dF, G, dG, P, S, delta});
  }
  return out;
}
