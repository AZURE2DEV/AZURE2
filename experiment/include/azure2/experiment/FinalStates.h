#ifndef AZURE2_EXPERIMENT_FINALSTATES_H
#define AZURE2_EXPERIMENT_FINALSTATES_H

#include <vector>

#include "azure2/experiment/LevelTables.h"

namespace azure2 {
namespace experiment {

/*!
 * The final state of the external capture of a photon pair (user guide, "Gamma-Ray Calculations"): a level of the
 * project with the J^pi of the pair's heavy particle, in its compound nucleus, within 1 keV of its E_x; it carries the
 * ANC and the gamma widths.  The one rule: the engine's set-up (CompoundSystem::ParseExternalCapture) and the validation
 * (Schema::Build) both take their refusals from ExternalCaptureFinalStates.
 */
struct FinalStateOfPair {
  int pairKey = 0;
  int twoJ = 0, parity = 1;  ///< the final state's 2J and parity (the pair's heavy particle)
  double energy = 0.0;       ///< MeV, its E_x
  bool groupExists = false;  ///< a J^pi group of the compound nucleus with that J^pi exists
  /// The levels within 1 keV of `energy` in that group (indices into LevelTables::levels), in level order.  None: no
  /// final state; more than one: an ambiguous one; one: the final state.
  std::vector<int> levels;
  /// When the final state is one level that an earlier pair of the list already has: that pair's key, else 0.
  int sharedWithPairKey = 0;
};

/// The final states of the photon pairs `pairKeys` (the pairs with external capture, in the order the set-up takes
/// them: an earlier pair keeps a level a later one shares), one per key.  A key that is no photon pair is left out.
std::vector<FinalStateOfPair> ExternalCaptureFinalStates(const LevelTables &tables, const std::vector<int> &pairKeys);
/// The entry of pair `key` in `states`, or nullptr.
const FinalStateOfPair *FinalStateOf(const std::vector<FinalStateOfPair> &states, int key);

/// A capture asked of the model: from the particle pair `entranceKey` into the photon pair `exitKey` (-1: into every
/// photon pair of the entrance's compound nucleus, a total capture).  `request` is the caller's index (its record).
struct CaptureRequest {
  int entranceKey = 0;
  int exitKey = 0;
  int request = 0;
};

/*!
 * A capture into a final state that is not bound where the engine takes it as bound (E_FINAL_STATE_UNBOUND):
 *
 *  - Bound: the photon pair has external capture and its final state (a level of <levels>) lies at or above the
 *    threshold of a particle pair with a channel in the final state's J^pi group.  The N_f integral, the
 *    ANC-to-amplitude factor and the radial integral take |E_f - S| of that pair as a binding energy
 *    (src/NFIntegral.cpp, src/ExternalCaptureFinalStates.cpp, src/ECIntegral.cpp).
 *  - Entrance: the final state lies above the threshold of the entrance pair of a capture asked into it, internal or
 *    external.  Between the two energies the photon energy E_gamma = E_x - E_f is negative and E_gamma^(L+1/2) is NaN
 *    (the photon channel's penetrability, src/ChannelTables.cpp; the external amplitude, src/UnitExternalCapture.cpp).
 *
 * Energies are E_x of the compound nucleus (MeV): a pair's threshold is S + e1 + e2 (PairInfo::threshold); the final
 * state's energy is the larger of the photon pair's threshold and its level's energy (they agree within 1 keV).
 */
struct UnboundFinalState {
  enum class Kind { Bound, Entrance };
  Kind kind = Kind::Bound;
  int photonPairKey = 0;
  int particlePairKey = 0;   ///< the pair whose threshold lies at or below the final state (the lowest one)
  double finalEnergy = 0.0;  ///< MeV, E_x of the final state
  double threshold = 0.0;    ///< MeV, E_x at which the particle pair opens
  int level = -1;            ///< Bound: the final state's index into LevelTables::levels
  /// CaptureRequest::request of the capture asked (Bound: the first one into the photon pair; -1 when none is)
  int request = -1;
};

/// The captures of `tables` into final states above a particle threshold that the engine would compute as bound or
/// at a negative photon energy: one Bound finding per photon pair with external capture (the lowest threshold), then
/// one Entrance finding per request and photon pair that has no Bound finding.  Empty: none.
std::vector<UnboundFinalState> UnboundFinalStates(const LevelTables &tables, const std::vector<CaptureRequest> &requests);

}  // namespace experiment
}  // namespace azure2

#endif
