#include "azure2/experiment/HybridPotentialSettings.h"

#include <stdexcept>

namespace azure2 {
namespace experiment {

physics::HybridPotential HybridPotentialSetting::Potential() const {
  physics::HybridPotential p =
      shape == physics::FormShape::Gaussian ? physics::GaussianPotential(V0, r0) : physics::WoodsSaxonPotential(V0, R, a);
  if (!(coulomb_radius_fm > 0.0) || !(matching_radius_fm > 0.0) || !(step_fm > 0.0))
    throw std::invalid_argument("the Coulomb radius, the matching radius and the step must be positive");
  p.coulomb_radius_fm = coulomb_radius_fm;
  p.matching_radius_fm = matching_radius_fm;
  p.step_fm = step_fm;
  if (spin_orbit)
    p.spin_orbit = physics::WoodsSaxonSpinOrbit(spin_orbit->depth_MeV, spin_orbit->form.radius_fm,
                                                spin_orbit->form.diffuseness_fm, spin_orbit->acts_on);
  p.interior_boundary = interior_boundary;
  return p;
}

std::string HybridPotentialSetting::ShapeName() const {
  return shape == physics::FormShape::Gaussian ? "Gaussian" : "WoodsSaxon";
}

const HybridPotentialSetting &HybridPotentialSettings::For(int pairKey) const {
  std::map<int, HybridPotentialSetting>::const_iterator it = pairs.find(pairKey);
  return it != pairs.end() ? it->second : defaults;
}

bool HybridPotentialSettings::AnyEnabled() const {
  if (defaults.enabled) return true;
  for (const auto &p : pairs)
    if (p.second.enabled) return true;
  for (const auto &p : byL)
    if (p.second.enabled) return true;
  return false;
}

std::map<int, physics::HybridPotential> HybridPotentialSettings::PotentialsByL(int pairKey) const {
  std::map<int, physics::HybridPotential> out;
  for (const auto &p : byL)
    if (p.first.first == pairKey && p.second.enabled) out.emplace(p.first.second, p.second.Potential());
  return out;
}

std::optional<physics::HybridPotential> HybridPotentialSettings::PotentialFor(int pairKey) const {
  const HybridPotentialSetting &s = For(pairKey);
  if (!s.enabled) return std::nullopt;
  return s.Potential();
}

std::vector<int> HybridPotentialSettings::ConfiguredPairs() const {
  std::vector<int> keys;
  for (const auto &p : pairs) keys.push_back(p.first);
  return keys;
}

}  // namespace experiment
}  // namespace azure2
