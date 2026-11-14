#ifndef AZURE2_EXPERIMENT_SPINORBITCHECKS_H
#define AZURE2_EXPERIMENT_SPINORBITCHECKS_H

#include <string>

#include "azure2/physics/HybridPotential.h"

namespace azure2 {
namespace experiment {

/*!
 * The refusals of a spin-orbit term (docs/upgrade/spin_orbit_potential.md section 3.4), written once for the three
 * places that need them: the format-3 build of a `potential` record (BuildPotentials), the pair's set-up (PPair)
 * and the engine's API for channel functions (api/).  A refusal has a code of DiagnosticCodes() and a message; an
 * empty code means allowed.
 */
struct SpinOrbitRefusal {
  std::string code;  ///< "" when allowed
  std::string message;
  explicit operator bool() const { return !code.empty(); }
};

/// A spin-orbit term acting on `on` of a pair whose particles have spins 2j = twoJLight, twoJHeavy: E_SO_TWO_SPINS
/// when both spins are non-zero (the term couples the channel spins of one (l, J): option (b), not implemented),
/// E_SO_ZERO_SPIN when the particle it acts on has spin 0 (the term would vanish).  pairText names the pair.
SpinOrbitRefusal SpinOrbitPairRefusal(int twoJLight, int twoJHeavy, physics::SpinOrbitParticle on, const std::string &pairText);

/// A request for the channel functions of a pair, with (jGiven) or without J: E_SO_NEEDS_J when the pair's
/// potential has a spin-orbit term and no J is given (its functions depend on J), E_SO_NO_TERM when a J is given
/// for a pair without one (its functions do not depend on J; a J there is a mistake about the model).
SpinOrbitRefusal ChannelFunctionsJRefusal(bool pairHasSpinOrbit, bool jGiven, const std::string &pairText);

}  // namespace experiment
}  // namespace azure2

#endif
