#include "HybridPotentialSettings.h"

azure2::physics::HybridPotential HybridPotentialSetting::Potential() const {
  if (shape == azure2::physics::FormShape::Gaussian) return azure2::physics::GaussianPotential(V0, r0);
  return azure2::physics::WoodsSaxonPotential(V0, R, a);
}

std::string HybridPotentialSetting::ShapeName() const {
  return shape == azure2::physics::FormShape::Gaussian ? "Gaussian" : "WoodsSaxon";
}

const HybridPotentialSetting &HybridPotentialSettings::For(int pairKey) const {
  std::map<int, HybridPotentialSetting>::const_iterator it = pairs.find(pairKey);
  return it != pairs.end() ? it->second : defaults;
}

bool HybridPotentialSettings::AnyEnabled() const {
  if (defaults.enabled) return true;
  for (const auto &p : pairs)
    if (p.second.enabled) return true;
  return false;
}

std::optional<azure2::physics::HybridPotential> HybridPotentialSettings::PotentialFor(int pairKey) const {
  const HybridPotentialSetting &s = For(pairKey);
  if (!s.enabled) return std::nullopt;
  return s.Potential();
}

std::vector<int> HybridPotentialSettings::ConfiguredPairs() const {
  std::vector<int> keys;
  for (const auto &p : pairs) keys.push_back(p.first);
  return keys;
}
