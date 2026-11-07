#include <algorithm>
#include <stdexcept>
#include "NuclearPotentialManager.h"

NuclearPotentialManager &NuclearPotentialManager::instance() {
  static NuclearPotentialManager manager;
  return manager;
}

NuclearPotentialManager::NuclearPotentialManager() :
  defaultTag_(1),
  tagCounter_(1) {
  // Initialize with default Woods-Saxon potential (V=150 MeV, R=3.6 fm, a=0.6 fm)
  potential_ = build(defaultSetting_);
}

std::shared_ptr<NuclearPotential>
NuclearPotentialManager::build(const NuclearPotentialSetting &s) {
  if (s.type == "Gaussian")
    return std::make_shared<GaussianPotential>(s.V0, s.r0);
  return std::make_shared<WoodsSaxonPotential>(s.V0, s.R, s.a);
}

// -- default -----------------------------------------------------------------

std::shared_ptr<NuclearPotential> NuclearPotentialManager::getPotential() const {
  return potential_;
}

void NuclearPotentialManager::setWoodsSaxonPotential(double V0, double R, double a) {
  try {
    potential_ = std::make_shared<WoodsSaxonPotential>(V0, R, a);
    defaultSetting_.type = "WoodsSaxon";
    defaultSetting_.V0 = V0;
    defaultSetting_.R = R;
    defaultSetting_.a = a;
    defaultTag_ = ++tagCounter_;
  } catch (const std::exception &e) {
    // Log error but maintain previous valid state
    throw std::runtime_error(std::string("Failed to set Woods-Saxon potential: ") + e.what());
  }
}

void NuclearPotentialManager::setGaussianPotential(double V0, double r0) {
  try {
    potential_ = std::make_shared<GaussianPotential>(V0, r0);
    defaultSetting_.type = "Gaussian";
    defaultSetting_.V0 = V0;
    defaultSetting_.r0 = r0;
    defaultTag_ = ++tagCounter_;
  } catch (const std::exception &e) {
    throw std::runtime_error(std::string("Failed to set Gaussian potential: ") + e.what());
  }
}

NuclearPotentialManager::State NuclearPotentialManager::exportState() const {
  State s;
  s.defaultSetting = defaultSetting_;
  s.pairSettings = pairSettings_;
  s.pairTags = pairTags_;
  s.defaultTag = defaultTag_;
  return s;
}

void NuclearPotentialManager::importState(const State &state) {
  defaultSetting_ = state.defaultSetting;
  potential_ = build(defaultSetting_);
  pairSettings_ = state.pairSettings;
  pairPotentials_.clear();
  for (std::map<int, NuclearPotentialSetting>::const_iterator it = pairSettings_.begin();
       it != pairSettings_.end(); ++it)
    pairPotentials_[it->first] = build(it->second);
  pairTags_ = state.pairTags;
  for (std::map<int, NuclearPotentialSetting>::const_iterator it = pairSettings_.begin();
       it != pairSettings_.end(); ++it)
    if (pairTags_.find(it->first) == pairTags_.end()) pairTags_[it->first] = ++tagCounter_;
  defaultTag_ = state.defaultTag ? state.defaultTag : ++tagCounter_;
}

void NuclearPotentialManager::resetToDefault() {
  defaultSetting_ = NuclearPotentialSetting();
  potential_ = build(defaultSetting_);
  pairSettings_.clear();
  pairPotentials_.clear();
  pairTags_.clear();
  defaultTag_ = ++tagCounter_;
}

// -- per pair ----------------------------------------------------------------

void NuclearPotentialManager::setWoodsSaxonPotential(int pairKey, double V0,
                                                     double R, double a) {
  NuclearPotentialSetting s = getSetting(pairKey);
  s.type = "WoodsSaxon";
  s.V0 = V0;
  s.R = R;
  s.a = a;
  setSetting(pairKey, s);
}

void NuclearPotentialManager::setGaussianPotential(int pairKey, double V0, double r0) {
  NuclearPotentialSetting s = getSetting(pairKey);
  s.type = "Gaussian";
  s.V0 = V0;
  s.r0 = r0;
  setSetting(pairKey, s);
}

std::shared_ptr<NuclearPotential>
NuclearPotentialManager::getPotential(int pairKey) const {
  std::map<int, std::shared_ptr<NuclearPotential>>::const_iterator it =
      pairPotentials_.find(pairKey);
  if (it != pairPotentials_.end()) return it->second;
  return potential_;
}

// -- activation --------------------------------------------------------------

void NuclearPotentialManager::setDefaultEnabled(bool enabled) {
  defaultSetting_.enabled = enabled;
}

void NuclearPotentialManager::setPairEnabled(int pairKey, bool enabled) {
  NuclearPotentialSetting s = getSetting(pairKey);
  s.enabled = enabled;
  setSetting(pairKey, s);
}

bool NuclearPotentialManager::isPairEnabled(int pairKey) const {
  std::map<int, NuclearPotentialSetting>::const_iterator it = pairSettings_.find(pairKey);
  if (it != pairSettings_.end()) return it->second.enabled;
  return defaultSetting_.enabled;
}

bool NuclearPotentialManager::isAnyEnabled() const {
  if (defaultSetting_.enabled) return true;
  for (std::map<int, NuclearPotentialSetting>::const_iterator it = pairSettings_.begin();
       it != pairSettings_.end(); ++it)
    if (it->second.enabled) return true;
  return false;
}

// -- inspection and bulk edit ------------------------------------------------

NuclearPotentialSetting NuclearPotentialManager::getSetting(int pairKey) const {
  std::map<int, NuclearPotentialSetting>::const_iterator it = pairSettings_.find(pairKey);
  if (it != pairSettings_.end()) return it->second;
  return defaultSetting_;
}

NuclearPotentialSetting NuclearPotentialManager::getDefaultSetting() const {
  return defaultSetting_;
}

void NuclearPotentialManager::setSetting(int pairKey,
                                         const NuclearPotentialSetting &setting) {
  std::shared_ptr<NuclearPotential> built;
  try {
    built = build(setting);
  } catch (const std::exception &e) {
    throw std::runtime_error(std::string("Failed to set potential for pair: ") + e.what());
  }
  pairSettings_[pairKey] = setting;
  pairPotentials_[pairKey] = built;
  pairTags_[pairKey] = ++tagCounter_;
}

void NuclearPotentialManager::setDefaultSetting(const NuclearPotentialSetting &setting) {
  std::shared_ptr<NuclearPotential> built;
  try {
    built = build(setting);
  } catch (const std::exception &e) {
    throw std::runtime_error(std::string("Failed to set default potential: ") + e.what());
  }
  defaultSetting_ = setting;
  potential_ = built;
  defaultTag_ = ++tagCounter_;
}

bool NuclearPotentialManager::hasPairSetting(int pairKey) const {
  return pairSettings_.find(pairKey) != pairSettings_.end();
}

void NuclearPotentialManager::clearPairSetting(int pairKey) {
  pairSettings_.erase(pairKey);
  pairPotentials_.erase(pairKey);
  pairTags_.erase(pairKey);
}

long NuclearPotentialManager::tagFor(int pairKey) const {
  if (!isPairEnabled(pairKey)) return 0;
  std::map<int, long>::const_iterator it = pairTags_.find(pairKey);
  if (it != pairTags_.end()) return it->second;
  return defaultTag_;
}

std::vector<int> NuclearPotentialManager::configuredPairs() const {
  std::vector<int> keys;
  for (std::map<int, NuclearPotentialSetting>::const_iterator it = pairSettings_.begin();
       it != pairSettings_.end(); ++it)
    keys.push_back(it->first);
  std::sort(keys.begin(), keys.end());
  return keys;
}
