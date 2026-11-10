// The refusals of a spin-orbit term (SpinOrbitChecks.h).
#include "azure2/experiment/SpinOrbitChecks.h"

#include "azure2/experiment/Diagnostics.h"

namespace azure2 {
namespace experiment {

namespace {

std::string SpinText(int twoJ) { return twoJ % 2 ? std::to_string(twoJ) + "/2" : std::to_string(twoJ / 2); }

}  // namespace

SpinOrbitRefusal SpinOrbitPairRefusal(int twoJLight, int twoJHeavy, physics::SpinOrbitParticle on, const std::string &pairText) {
  SpinOrbitRefusal r;
  if (twoJLight != 0 && twoJHeavy != 0) {
    r.code = codes::kSoTwoSpins;
    r.message = "pair " + pairText + " has two spinning particles (spins " + SpinText(twoJLight) + " and " +
                SpinText(twoJHeavy) + "): a spin-orbit term on one of them couples the channel spins of one (l, J), "
                "which is not implemented (docs/upgrade/spin_orbit_potential.md section 2.3); only pairs with one "
                "spinning particle can have a spin-orbit term";
    return r;
  }
  const bool light = on == physics::SpinOrbitParticle::Light;
  if ((light ? twoJLight : twoJHeavy) == 0) {
    r.code = codes::kSoZeroSpin;
    r.message = std::string("the spin-orbit term of pair ") + pairText + " acts on the " + (light ? "light" : "heavy") +
                " particle, which has spin 0 (the term would vanish); the spinning particle is the " +
                (light ? "heavy" : "light") + " one" + ((light ? twoJHeavy : twoJLight) == 0 ? " -- neither has spin" : "");
  }
  return r;
}

SpinOrbitRefusal ChannelFunctionsJRefusal(bool pairHasSpinOrbit, bool jGiven, const std::string &pairText) {
  SpinOrbitRefusal r;
  if (pairHasSpinOrbit && !jGiven) {
    r.code = codes::kSoNeedsJ;
    r.message = "pair " + pairText + " has a spin-orbit term: its channel functions depend on J; give J";
  } else if (!pairHasSpinOrbit && jGiven) {
    r.code = codes::kSoNoTerm;
    r.message = "pair " + pairText + " has no spin-orbit term: its channel functions do not depend on J; leave J out";
  }
  return r;
}

}  // namespace experiment
}  // namespace azure2
