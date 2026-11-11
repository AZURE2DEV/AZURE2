#ifndef AZURE2_EXPERIMENT_FINALSTATES_H
#define AZURE2_EXPERIMENT_FINALSTATES_H

#include <vector>

#include "azure2/experiment/LevelTables.h"

namespace azure2 {
namespace experiment {

/*!
 * The final state of the external capture of a photon pair (user guide, "Gamma-Ray Calculations"): a level of the
 * project with the J^pi of the pair's heavy particle, in its compound nucleus, within 1 keV of its E_x; it carries the
 * ANC and the gamma widths.  The one rule: the engine's set-up (CNuc::ParseExternalCapture) and the validation
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

}  // namespace experiment
}  // namespace azure2

#endif
