#include "NuclearPotentialManager.h"

NuclearPotentialManager& NuclearPotentialManager::instance() {
  static NuclearPotentialManager manager;
  return manager;
}

NuclearPotentialManager::NuclearPotentialManager()
    : currentType_("WoodsSaxon") {
  // Initialize with default Woods-Saxon potential (V=150 MeV, R=3.6 fm, a=0.6 fm)
  potential_ = std::make_shared<WoodsSaxonPotential>(150.0, 3.6, 0.6);
}

std::shared_ptr<NuclearPotential> NuclearPotentialManager::getPotential() const {
  return potential_;
}

std::string NuclearPotentialManager::getCurrentPotentialType() const {
  return currentType_;
}

void NuclearPotentialManager::setWoodsSaxonPotential(double V0, double R, double a) {
  try {
    potential_ = std::make_shared<WoodsSaxonPotential>(V0, R, a);
    currentType_ = "WoodsSaxon";
  } catch(const std::exception& e) {
    // Log error but maintain previous valid state
    throw std::runtime_error(std::string("Failed to set Woods-Saxon potential: ") + e.what());
  }
}

bool NuclearPotentialManager::getWoodsSaxonParameters(double& V0, double& R, double& a) const {
  auto ws = std::dynamic_pointer_cast<WoodsSaxonPotential>(potential_);
  if(!ws) {
    return false;
  }
  V0 = ws->get_V0();
  R = ws->get_R();
  a = ws->get_a();
  return true;
}

void NuclearPotentialManager::setGaussianPotential(double V0, double r0) {
  try {
    potential_ = std::make_shared<GaussianPotential>(V0, r0);
    currentType_ = "Gaussian";
  } catch(const std::exception& e) {
    throw std::runtime_error(std::string("Failed to set Gaussian potential: ") + e.what());
  }
}

bool NuclearPotentialManager::getGaussianParameters(double& V0, double& r0) const {
  auto gauss = std::dynamic_pointer_cast<GaussianPotential>(potential_);
  if(!gauss) {
    return false;
  }
  V0 = gauss->get_V0();
  r0 = gauss->get_r0();
  return true;
}

void NuclearPotentialManager::resetToDefault() {
  potential_ = std::make_shared<WoodsSaxonPotential>(150.0, 3.6, 0.6);
  currentType_ = "WoodsSaxon";
}
