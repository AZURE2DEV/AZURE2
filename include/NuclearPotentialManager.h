#ifndef NUCLEARPOTENTIALMANAGER_H
#define NUCLEARPOTENTIALMANAGER_H

#include <memory>
#include "NuclearPotential.h"

/*!
 * @brief Singleton manager for the shared nuclear potential used throughout the application.
 *
 * This class maintains a single, globally-accessible nuclear potential instance that is
 * used by all CoulFunc objects and other calculations. It ensures consistency across
 * all components of the application.
 *
 * Usage:
 *   auto& manager = NuclearPotentialManager::instance();
 *   manager.setWoodsSaxonPotential(50.0, 5.5, 0.65);
 *   auto potential = manager.getPotential();
 */
class NuclearPotentialManager {
 public:
  /*!
   * @brief Get the singleton instance of the nuclear potential manager
   */
  static NuclearPotentialManager &instance();

  /*!
   * @brief Prevent copying
   */
  NuclearPotentialManager(const NuclearPotentialManager &) = delete;
  NuclearPotentialManager &operator=(const NuclearPotentialManager &) = delete;

  /*!
   * @brief Get the current nuclear potential
   */
  std::shared_ptr<NuclearPotential> getPotential() const;

  /*!
   * @brief Get the current potential type name
   */
  std::string getCurrentPotentialType() const;

  // Woods-Saxon Potential Methods
  /*!
   * @brief Set Woods-Saxon potential parameters
   * @param V0 Potential depth in MeV
   * @param R Nuclear radius in fm
   * @param a Surface diffuseness in fm
   */
  void setWoodsSaxonPotential(double V0, double R, double a);

  /*!
   * @brief Get Woods-Saxon potential parameters
   */
  bool getWoodsSaxonParameters(double &V0, double &R, double &a) const;

  // Gaussian Potential Methods
  /*!
   * @brief Set Gaussian potential parameters
   * @param V0 Potential depth in MeV
   * @param r0 Gaussian width parameter in fm
   */
  void setGaussianPotential(double V0, double r0);

  /*!
   * @brief Get Gaussian potential parameters
   */
  bool getGaussianParameters(double &V0, double &r0) const;

  /*!
   * @brief Reset to default Woods-Saxon potential (V=150, R=3.6, a=0.6)
   */
  void resetToDefault();

 private:
  /*!
   * @brief Private constructor for singleton pattern
   */
  NuclearPotentialManager();

  std::shared_ptr<NuclearPotential> potential_;
  std::string currentType_;
};

#endif  // NUCLEARPOTENTIALMANAGER_H
