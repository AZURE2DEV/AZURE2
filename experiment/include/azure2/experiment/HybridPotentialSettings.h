#ifndef AZURE2_EXPERIMENT_HYBRIDPOTENTIALSETTINGS_H
#define AZURE2_EXPERIMENT_HYBRIDPOTENTIALSETTINGS_H

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "azure2/physics/HybridPotential.h"

namespace azure2 {
namespace experiment {

/*!
 * The hybrid potential of one pair as a project states it: a `potential` record of format 3 (Schema::Build), a
 * section of the old <potential> block (convert/ ReadPotentialBlock); pyazr edits the record in a Project.  An input
 * record; the pair's potential is the value Potential() makes from it.  The old block has no Coulomb radius, matching
 * radius or step: it keeps the defaults below, which are the physics value's (3.5 fm, 15 fm, 0.0275 fm).
 */
struct HybridPotentialSetting {
  bool enabled = false;
  physics::FormShape shape = physics::FormShape::WoodsSaxon;
  double V0 = 150.0;  ///< depth (MeV), both shapes
  double R = 3.6;     ///< Woods-Saxon radius (fm)
  double a = 0.6;     ///< Woods-Saxon diffuseness (fm)
  double r0 = 5.0;    ///< Gaussian width (fm)
  double coulomb_radius_fm = physics::HybridPotential().coulomb_radius_fm;    ///< R_C of the charged sphere (fm)
  double matching_radius_fm = physics::HybridPotential().matching_radius_fm;  ///< r_m (fm)
  double step_fm = physics::HybridPotential().step_fm;                        ///< Numerov step (fm)
  /// The spin-orbit term (a format-3 record's spin_orbit); the old block has none.
  std::optional<physics::SpinOrbitTerm> spin_orbit;
  /// The boundary constant b of the potential's interior R-function (a format-3 record's interior_boundary;
  /// docs/upgrade/interior_potential.md); none: the potential acts outside the channel radius only.
  std::optional<double> interior_boundary;
  /// The potential of this setting; throws std::invalid_argument when a parameter of its shape, R_C, r_m, the
  /// step or of the spin-orbit term is not positive.
  physics::HybridPotential Potential() const;
  /// "WoodsSaxon" or "Gaussian", the names pyazr uses.
  std::string ShapeName() const;
};

/*!
 * The hybrid potentials of a project, per pair key: a default for every pair without a section of its own (the old
 * <potential> block and a pyazr session only; a format-3 project has none), and the pairs' own settings.  Owned by the
 * Config of a run or a pyazr session; CompoundSystem::Fill gives each pair the potential that applies to it
 * (PPair::GetHybridPotential).
 */
struct HybridPotentialSettings {
  HybridPotentialSetting defaults;
  std::map<int, HybridPotentialSetting> pairs;

  /// The setting that applies to a pair: its own, else the default.
  const HybridPotentialSetting &For(int pairKey) const;
  /// Is the hybrid model on for any pair?
  bool AnyEnabled() const;
  /// The potential of a pair, or none when the model is off for it.
  std::optional<physics::HybridPotential> PotentialFor(int pairKey) const;
  /// The keys of the pairs with a section of their own, ascending.
  std::vector<int> ConfiguredPairs() const;
};

}  // namespace experiment
}  // namespace azure2

#endif
