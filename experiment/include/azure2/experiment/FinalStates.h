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
 * An external-capture final state the engine cannot compute as unbound (E_FINAL_STATE_UNBOUND).
 *
 * A final state above the threshold of a particle pair with a channel in its J^pi group is unbound in that channel:
 * the engine treats it as a narrow resonance, its exterior the outgoing wave O_l(k_f r)/O_l(k_f a) at
 * E_f - threshold (src/ECIntegral.cpp, IntegrateUnbound; the N_f integral continued as (hbar^2/2 mu a) dS/dE,
 * src/ExternalCaptureFinalStates.cpp).  That wave cannot be evaluated at the threshold itself, nor just above the
 * threshold of a charged pair, where its Sommerfeld parameter eta exceeds kMaxOutgoingWaveEta (within about a keV for
 * light nuclei): such a final state is refused.  Everything else is computed: an unbound final state farther above
 * the threshold; a capture whose entrance pair opens below the final state (no photon below E_f: the pathway is
 * closed and its cross section 0 there); internal capture into any final state.
 *
 * Energies are E_x of the compound nucleus (MeV): a pair's threshold is S + e1 + e2 (PairInfo::threshold); the final
 * state's energy is the larger of the photon pair's threshold and its level's energy (they agree within 1 keV).  eta
 * is the non-relativistic one (the engine's check at set-up, with the run's kinematics, is the same to 1e-3).
 */
struct UnboundFinalState {
  int photonPairKey = 0;
  int particlePairKey = 0;   ///< the pair too close above whose threshold the final state lies (the largest eta)
  double finalEnergy = 0.0;  ///< MeV, E_x of the final state
  double threshold = 0.0;    ///< MeV, E_x at which the particle pair opens
  double eta = 0.0;          ///< the Sommerfeld parameter of the pair at finalEnergy - threshold (0 at the threshold)
  int level = -1;            ///< the final state's index into LevelTables::levels
  /// CaptureRequest::request of the first capture asked into the photon pair; -1 when none is
  int request = -1;
};

/// Above this Sommerfeld parameter the outgoing wave of an unbound final state is not evaluated (the engine's
/// open-channel limit, ECIntegral::kMaxOpenChannelEta).
constexpr double kMaxOutgoingWaveEta = 100.0;

/// The final states of `tables` with external capture that lie at a particle threshold, or too close above the
/// threshold of a charged pair, with a channel in their J^pi group: one finding per photon pair.  Empty: none.
std::vector<UnboundFinalState> UnboundFinalStates(const LevelTables &tables, const std::vector<CaptureRequest> &requests);

}  // namespace experiment
}  // namespace azure2

#endif
