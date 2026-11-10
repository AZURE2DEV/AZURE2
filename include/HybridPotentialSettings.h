#ifndef HYBRIDPOTENTIALSETTINGS_H
#define HYBRIDPOTENTIALSETTINGS_H

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "azure2/physics/HybridPotential.h"

/*!
 * One section of the <potential> block (Config::ReadPotentialBlock) or of pyazr's set_nuclear_potential: whether the
 * hybrid model is on and the potential's parameters.  An input record; the pair's potential is the value
 * Potential() makes from it.
 */
struct HybridPotentialSetting {
  bool enabled = false;
  azure2::physics::FormShape shape = azure2::physics::FormShape::WoodsSaxon;
  double V0 = 150.0;  ///< depth (MeV), both shapes
  double R = 3.6;     ///< Woods-Saxon radius (fm)
  double a = 0.6;     ///< Woods-Saxon diffuseness (fm)
  double r0 = 5.0;    ///< Gaussian width (fm)
  /// The potential of this setting; throws std::invalid_argument when a parameter of its shape is not positive.
  azure2::physics::HybridPotential Potential() const;
  /// "WoodsSaxon" or "Gaussian", the names pyazr uses.
  std::string ShapeName() const;
};

/*!
 * The hybrid potentials of a project, per pair key: a default for every pair without a section of its own, and
 * the pairs' sections.  Owned by the Config of a run or a pyazr session; CNuc::Fill gives each pair the potential
 * that applies to it (PPair::GetHybridPotential).
 */
struct HybridPotentialSettings {
  HybridPotentialSetting defaults;
  std::map<int, HybridPotentialSetting> pairs;

  /// The setting that applies to a pair: its own, else the default.
  const HybridPotentialSetting &For(int pairKey) const;
  /// Is the hybrid model on for any pair?
  bool AnyEnabled() const;
  /// The potential of a pair, or none when the model is off for it.
  std::optional<azure2::physics::HybridPotential> PotentialFor(int pairKey) const;
  /// The keys of the pairs with a section of their own, ascending.
  std::vector<int> ConfiguredPairs() const;
};

#endif
