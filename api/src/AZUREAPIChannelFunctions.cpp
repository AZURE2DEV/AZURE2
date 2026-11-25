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
#include "PotentialWave.h"
#include "azure2/experiment/SpinOrbitChecks.h"
#include "azure2/physics/HybridPotential.h"

#include <algorithm>
#include <vector>

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
 * request = {pair key, l, 2J (< 0: none), radius (fm, <= 0: the radius of the pair's channels of this l), n, E_1 ... E_n (MeV, c.m.)}.  A pair
 * whose functions depend on J -- a spin-orbit term, or a `wave` of its potential that names a J -- needs J
 * (E_SO_NEEDS_J); any other pair refuses one (E_SO_NO_TERM).  Throws std::invalid_argument for those and for a J that
 * l and the spin cannot make.
 *
 * A partial wave with an interior background (potential_background.md section 2): F, F', G, G', P, S and the phase are
 * all the dressed ones of Eq. (9), the channel functions the levels on top of the background see (delta = -phi~,
 * P~ = rho/|H~|^2, S~); without a background the plain (or hybrid) ones.  Before 2026-11-20 F, G and the phase were
 * the undressed waves while P and S were dressed.  A wave with `only bound_state`: the point-Coulomb ones.
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
  // The radius by default: the one the model's channels of this pair and l have (a channel_radius record gives a
  // channel its own); the pair's when the model has no such channel.  Channels of the l with more than one radius
  // (records with s or J) have no one default: the caller names the radius.
  double defaultRadius = pair->GetPairRadius();
  {
    const int pairNum = compound_->GetPairNumFromKey(pairKey);
    bool found = false;
    for (int j = 1; j <= compound_->NumJGroups(); j++) {
      const JGroup *group = compound_->GetJGroup(j);
      for (int ch = 1; ch <= group->NumChannels(); ch++) {
        const AChannel *channel = group->GetChannel(ch);
        if (channel->GetPairNum() != pairNum || channel->GetL() != lValue || channel->GetRadType() != 'P') continue;
        if (found && channel->GetRadius() != defaultRadius && radius <= 0.0)
          throw std::invalid_argument("ERROR: the channels of pair " + std::to_string(pairKey) + " with l = " +
                                      std::to_string(lValue) + " have more than one radius (channel_radius records "
                                      "with s or J): give the radius");
        if (!found) defaultRadius = channel->GetRadius();
        found = true;
      }
    }
  }
  if (radius <= 0.0) radius = defaultRadius;

  // J matters for a spin-orbit term and for a `wave` that names a J (its own depth or background): PPair::SpinOrbitTwoJ
  const bool spinOrbit = pair->GetHybridPotential() && pair->GetHybridPotential()->spin_orbit;
  const bool dependsOnJ = pair->GetHybridPotential() && azure2::physics::DependsOnJ(*pair->GetHybridPotential());
  azure2::experiment::SpinOrbitRefusal refusal =
      azure2::experiment::ChannelFunctionsJRefusal(dependsOnJ, requestedTwoJ >= 0, std::to_string(pairKey));
  if (refusal && dependsOnJ && !spinOrbit)
    refusal.message = "the potential of pair " + std::to_string(pairKey) + " has a wave with its own J: its channel "
                      "functions depend on J; give J";
  if (refusal) throw std::invalid_argument("ERROR (" + refusal.code + "): " + refusal.message);
  const int twoJ = dependsOnJ ? pair->SpinOrbitTwoJ(requestedTwoJ / 2.0) : azure2::physics::kNoSpinOrbit;
  pair->SpinOrbitFactor(lValue, twoJ);  // throws std::invalid_argument for a J that l and the spin cannot make

  CoulFunc coul(pair);

  out.push_back(static_cast<double>(nE));
  for (int i = 0; i < nE; ++i) {
    double energy = request[5 + i];
    double F = 0., dF = 0., G = 0., dG = 0., P = 0., S = 0., delta = 0.;
    if (energy > 0.0 && !ChannelFunc::IsCoulombThreshold(pair, energy)) {
      try {
        CoulWaves w = coul.ChannelWaves(lValue, twoJ, radius, energy);  // dressed where a background acts
        F = w.F, dF = w.dF, G = w.G, dG = w.dG;
        P = coul.Penetrability(lValue, twoJ, radius, energy);
        S = coul.PEShift(lValue, twoJ, radius, energy);
        delta = -std::atan2(w.F, w.G);
      } catch (...) {
        // A failed evaluation returns zeros for that energy rather than aborting the grid.
        F = dF = G = dG = P = S = delta = 0.;
      }
    } else if (radius == defaultRadius) {  // closed (E <= 0, eta > 100): S through threshold (ChannelFunc)
      S = ChannelFunc(pair, radius).Shift(lValue, twoJ, energy);
    }
    out.insert(out.end(), {F, dF, G, dG, P, S, delta});
  }
  return out;
}

PotentialWaveArgs PotentialWaveOf(const PPair &pair, int l, int requestedTwoJ, const char *call) {
  const std::string key = std::to_string(pair.GetPairKey());
  PotentialWaveArgs a;
  a.pair = &pair;
  if (!pair.GetHybridPotential())
    throw std::invalid_argument(std::string(call) + ": pair " + key + " has no potential record");
  if (l < 0) throw std::invalid_argument(std::string(call) + ": l must be >= 0");
  const bool dependsOnJ = azure2::physics::DependsOnJ(*pair.GetHybridPotential());
  if (dependsOnJ && requestedTwoJ < 0)
    throw std::invalid_argument(std::string(call) + ": the potential of pair " + key +
                                " depends on J (a spin-orbit term or a wave with its own J): give J");
  a.twoJ = dependsOnJ ? pair.SpinOrbitTwoJ(requestedTwoJ / 2.0) : azure2::physics::kNoSpinOrbit;
  a.ldoti = pair.SpinOrbitFactor(l, a.twoJ);  // throws for a J that l and the spin cannot make
  a.potential = azure2::physics::ForPartialWave(*pair.GetHybridPotential(), l, a.twoJ);
  a.mu = pair.GetRedMass() * uconv;
  return a;
}

namespace {

/// PotentialWaveOf for the pair of key `pairKey` of the model; throws std::invalid_argument for no such pair.
PotentialWaveArgs WaveArgs(const CompoundSystem *compound, int pairKey, int l, int requestedTwoJ, const char *call) {
  if (!compound || !compound->IsPairKey(pairKey))
    throw std::invalid_argument(std::string(call) + ": no pair " + std::to_string(pairKey) + " in this model");
  return PotentialWaveOf(*compound->GetPair(compound->GetPairNumFromKey(pairKey)), l, requestedTwoJ, call);
}

}  // namespace

std::vector<azure2::physics::BoxState> AZUREAPI::GetPotentialBoxStates(int pairKey, int l, int twoJ, double radius,
                                                                       double boundary, double eMin, double eMax) const {
  const PotentialWaveArgs a = WaveArgs(compound_.get(), pairKey, l, twoJ, "potential_box_states");
  if (!(radius > 0.0)) radius = a.pair->GetPairRadius();
  if (!(eMax > eMin)) throw std::invalid_argument("potential_box_states: e_max must exceed e_min");
  return azure2::physics::InteriorBoxStates(a.potential, a.pair->GetZ(1), a.pair->GetZ(2), a.mu, l, a.ldoti, radius,
                                            boundary, eMin, eMax);
}

std::vector<double> AZUREAPI::GetPotentialBoundState(int pairKey, int l, int twoJ, double energy, double window) const {
  const PotentialWaveArgs a = WaveArgs(compound_.get(), pairKey, l, twoJ, "potential_bound_state");
  if (!(energy < 0.0)) throw std::invalid_argument("potential_bound_state: the energy must be below threshold (< 0 MeV)");
  if (!(window > 0.0)) throw std::invalid_argument("potential_bound_state: the window must be positive");
  // The matching radius: inside it the potential, beyond it its Whittaker tail (the eigenvalue does not depend on it).
  const double a_fm = std::min(a.pair->GetPairRadius(), a.potential.matching_radius_fm);
  const azure2::physics::PotentialEigenvalue ev = azure2::physics::NearestEigenvalue(
      a.potential, a.pair->GetZ(1), a.pair->GetZ(2), a.mu, l, a.ldoti, a_fm, energy, std::min(window, -energy));
  if (!ev.found) return {0.0, 0.0, 0.0};
  const double anc = azure2::physics::PotentialBoundStateANC(a.potential, a.pair->GetZ(1), a.pair->GetZ(2), a.mu, l,
                                                             a.ldoti, a_fm, ev.E_MeV);
  return {1.0, ev.E_MeV, anc};
}

azure2::physics::EffectiveRange AZUREAPI::GetPotentialEffectiveRange(int pairKey, int l, int twoJ) const {
  const PotentialWaveArgs a = WaveArgs(compound_.get(), pairKey, l, twoJ, "potential_scattering_length");
  return azure2::physics::PotentialEffectiveRange(a.potential, a.pair->GetZ(1), a.pair->GetZ(2), a.mu, l, a.ldoti);
}

azure2::physics::DepthWidthForRange AZUREAPI::GetPotentialDepthFor(int pairKey, int l, int twoJ, double a, double r) const {
  const PotentialWaveArgs w = WaveArgs(compound_.get(), pairKey, l, twoJ, "potential_depth_for");
  const int z1 = w.pair->GetZ(1), z2 = w.pair->GetZ(2);
  if (std::isfinite(r)) return azure2::physics::DepthWidthForEffectiveRange(w.potential, z1, z2, w.mu, l, w.ldoti, a, r);
  const azure2::physics::DepthForLength d = azure2::physics::DepthForScatteringLength(w.potential, z1, z2, w.mu, l, w.ldoti, a);
  azure2::physics::DepthWidthForRange out;
  out.found = d.found;
  out.why = d.why;
  out.depth_MeV = d.depth_MeV;
  out.width_fm = w.potential.central.width_fm;
  out.a = d.a;
  out.r = d.r;
  out.boundStates = d.boundStates;
  out.nodes = d.nodes;
  return out;
}
