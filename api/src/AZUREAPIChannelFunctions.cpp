// AZUREAPI: the channel functions of a pair on a requested energy grid (pyazr's channel_functions, through the
// binding coulomb_functions).
#include <cmath>
#include <stdexcept>
#include <string>

#include "AZUREAPI.h"
#include "CompoundSystem.h"
#include "ChannelFunc.h"
#include "Config.h"
#include "CoulFunc.h"
#include "PPair.h"
#include "azure2/experiment/SpinOrbitChecks.h"

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
 *
 * request = {pair key, l, 2J (< 0: none), radius (fm, <= 0: the channel radius), n, E_1 ... E_n (MeV, c.m.)}.  A pair
 * whose potential has a spin-orbit term needs J (its functions depend on it: E_SO_NEEDS_J); any other pair refuses
 * one (E_SO_NO_TERM).  Throws std::invalid_argument for those and for a J that l and the spin cannot make.
 */
vector_r AZUREAPI::GetCoulombFunctions(const vector_r &request) const {
  vector_r out;
  if (request.size() < 5) return out;

  int pairKey = static_cast<int>(std::lround(request[0]));
  int lValue = static_cast<int>(std::lround(request[1]));
  const int requestedTwoJ = static_cast<int>(std::lround(request[2]));
  double radius = request[3];
  int nE = static_cast<int>(std::lround(request[4]));
  if (nE < 0 || request.size() < static_cast<size_t>(5 + nE)) return out;

  if (!compound_ || !compound_->IsPairKey(pairKey)) return out;
  const PPair *pair = compound_->GetPair(compound_->GetPairNumFromKey(pairKey));
  if (!pair) return out;
  if (radius <= 0.0) radius = pair->GetPairRadius();  // the pair's radius by default

  const bool spinOrbit = pair->GetHybridPotential() && pair->GetHybridPotential()->spin_orbit;
  const azure2::experiment::SpinOrbitRefusal refusal =
      azure2::experiment::ChannelFunctionsJRefusal(spinOrbit, requestedTwoJ >= 0, std::to_string(pairKey));
  if (refusal) throw std::invalid_argument("ERROR (" + refusal.code + "): " + refusal.message);
  const int twoJ = spinOrbit ? pair->SpinOrbitTwoJ(requestedTwoJ / 2.0) : azure2::physics::kNoSpinOrbit;
  pair->SpinOrbitFactor(lValue, twoJ);  // throws std::invalid_argument for a J that l and the spin cannot make

  CoulFunc coul(pair, !!(configure().paramMask & Config::USE_GSL_COULOMB_FUNC));

  out.push_back(static_cast<double>(nE));
  for (int i = 0; i < nE; ++i) {
    double energy = request[5 + i];
    double F = 0., dF = 0., G = 0., dG = 0., P = 0., S = 0., delta = 0.;
    if (energy > 0.0 && !ChannelFunc::IsCoulombThreshold(pair, energy)) {
      try {
        CoulWaves w = coul(lValue, twoJ, radius, energy);
        F = w.F, dF = w.dF, G = w.G, dG = w.dG;
        P = coul.Penetrability(lValue, twoJ, radius, energy);
        S = coul.PEShift(lValue, twoJ, radius, energy);
        delta = -std::atan2(w.F, w.G);
      } catch (...) {
        // A failed evaluation returns zeros for that energy rather than aborting the grid.
        F = dF = G = dG = P = S = delta = 0.;
      }
    } else if (radius == pair->GetPairRadius()) {  // closed (E <= 0, eta > 100): S through threshold (ChannelFunc)
      S = ChannelFunc(pair, !!(configure().paramMask & Config::USE_GSL_COULOMB_FUNC), radius).Shift(lValue, twoJ, energy);
    }
    out.insert(out.end(), {F, dF, G, dG, P, S, delta});
  }
  return out;
}
