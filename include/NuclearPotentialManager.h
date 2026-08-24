#ifndef NUCLEARPOTENTIALMANAGER_H
#define NUCLEARPOTENTIALMANAGER_H

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "NuclearPotential.h"

/*!
 * @brief One pair's hybrid-model setting: whether it is on, and which
 * potential it uses.
 *
 * Held by value so a caller can read a setting out, edit it and hand it back
 * without touching the manager's own state.
 */
struct NuclearPotentialSetting {
  /// Is the hybrid model active for this pair?
  bool enabled = false;
  ///"WoodsSaxon" or "Gaussian".
  std::string type = "WoodsSaxon";
  /// Depth (MeV), shared by both shapes.
  double V0 = 150.0;
  /// Woods-Saxon radius (fm).
  double R = 3.6;
  /// Woods-Saxon surface diffuseness (fm).
  double a = 0.6;
  /// Gaussian width (fm).
  double r0 = 5.0;
};

/*!
 * @brief The nuclear potentials the hybrid Coulomb method uses, per particle
 * pair.
 *
 * A potential is a property of the pair -- it modifies the radial wave
 * functions of that channel and nothing else -- so the manager keeps one
 * setting per pair key, plus a default that stands in for every pair without
 * one of its own.  A project that names no pair therefore behaves exactly as
 * it did when the model was global: the default is what every CoulFunc sees.
 *
 * The pair-less overloads read and write the default, so existing callers
 * (Config::ReadPotentialBlock's flat form, the GUI's single-potential path)
 * keep working unchanged.
 *
 * Usage:
 *   auto& manager = NuclearPotentialManager::instance();
 *   manager.setWoodsSaxonPotential(50.0, 5.5, 0.65);   // every pair
 *   manager.setGaussianPotential(2, 80.0, 4.0);        // pair 2 only
 *   auto potential = manager.getPotential(2);
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

  // -- default (every pair without an override) -----------------------------

  /*!
   * @brief Get the default nuclear potential
   */
  std::shared_ptr<NuclearPotential> getPotential() const;

  /*!
   * @brief Get the default potential type name
   */
  std::string getCurrentPotentialType() const;

  /*!
   * @brief Set the default Woods-Saxon potential
   * @param V0 Potential depth in MeV
   * @param R Nuclear radius in fm
   * @param a Surface diffuseness in fm
   */
  void setWoodsSaxonPotential(double V0, double R, double a);

  /*!
   * @brief Get the default Woods-Saxon parameters
   */
  bool getWoodsSaxonParameters(double &V0, double &R, double &a) const;

  /*!
   * @brief Set the default Gaussian potential
   * @param V0 Potential depth in MeV
   * @param r0 Gaussian width parameter in fm
   */
  void setGaussianPotential(double V0, double r0);

  /*!
   * @brief Get the default Gaussian parameters
   */
  bool getGaussianParameters(double &V0, double &r0) const;

  /*!
   * @brief Reset to default Woods-Saxon potential (V=150, R=3.6, a=0.6) and
   * drop every per-pair override
   */
  void resetToDefault();

  // -- per pair -------------------------------------------------------------

  /*!
   * @brief Set a Woods-Saxon potential for one pair only
   * @param pairKey The pair's key, as PPair::GetPairKey() reports it
   */
  void setWoodsSaxonPotential(int pairKey, double V0, double R, double a);

  /*!
   * @brief Set a Gaussian potential for one pair only
   */
  void setGaussianPotential(int pairKey, double V0, double r0);

  /*!
   * @brief Get the potential that applies to a pair -- its own if it has one,
   * otherwise the default
   */
  std::shared_ptr<NuclearPotential> getPotential(int pairKey) const;

  /*!
   * @brief Get the type name that applies to a pair
   */
  std::string getCurrentPotentialType(int pairKey) const;

  /*!
   * @brief Get the Woods-Saxon parameters that apply to a pair.  False if the
   * potential in force for it is not a Woods-Saxon.
   */
  bool getWoodsSaxonParameters(int pairKey, double &V0, double &R, double &a) const;

  /*!
   * @brief Get the Gaussian parameters that apply to a pair.  False if the
   * potential in force for it is not a Gaussian.
   */
  bool getGaussianParameters(int pairKey, double &V0, double &r0) const;

  // -- activation -----------------------------------------------------------

  /*!
   * @brief Turn the hybrid model on or off for the pairs that have no setting
   * of their own
   */
  void setDefaultEnabled(bool enabled);

  /*!
   * @brief Is the hybrid model on by default?
   */
  bool getDefaultEnabled() const;

  /*!
   * @brief Turn the hybrid model on or off for one pair.  Giving a pair an
   * activation state also gives it a setting, seeded from the default.
   */
  void setPairEnabled(int pairKey, bool enabled);

  /*!
   * @brief Is the hybrid model on for this pair?  Falls back to the default
   * for a pair with no setting of its own.
   */
  bool isPairEnabled(int pairKey) const;

  /*!
   * @brief Is the hybrid model on for any pair at all?
   */
  bool isAnyEnabled() const;

  // -- inspection and bulk edit ---------------------------------------------

  /*!
   * @brief The setting that applies to a pair, resolved against the default
   */
  NuclearPotentialSetting getSetting(int pairKey) const;

  /*!
   * @brief The default setting
   */
  NuclearPotentialSetting getDefaultSetting() const;

  /*!
   * @brief Install a setting for one pair
   */
  void setSetting(int pairKey, const NuclearPotentialSetting &setting);

  /*!
   * @brief Install the default setting
   */
  void setDefaultSetting(const NuclearPotentialSetting &setting);

  /*!
   * @brief Does this pair carry a setting of its own?
   */
  bool hasPairSetting(int pairKey) const;

  /*!
   * @brief Drop a pair's setting; it falls back to the default again
   */
  void clearPairSetting(int pairKey);

  /*!
   * @brief The keys of every pair carrying a setting of its own, ascending
   */
  std::vector<int> configuredPairs() const;

  /*!
   * @brief An identifier for the potential in force for a pair.
   *
   * 0 when the hybrid model is off for that pair, otherwise a value that
   * changes whenever the setting behind it does.  The Coulomb function memo
   * keys on it, so that switching a potential on, off, or to different
   * parameters cannot serve back waves computed under the previous one.
   */
  long tagFor(int pairKey) const;

 private:
  /*!
   * @brief Private constructor for singleton pattern
   */
  NuclearPotentialManager();

  /// Build the potential object a setting describes.
  static std::shared_ptr<NuclearPotential> build(const NuclearPotentialSetting &);

  NuclearPotentialSetting defaultSetting_;
  std::shared_ptr<NuclearPotential> potential_;
  std::string currentType_;
  std::map<int, NuclearPotentialSetting> pairSettings_;
  std::map<int, std::shared_ptr<NuclearPotential>> pairPotentials_;
  std::map<int, long> pairTags_;
  long defaultTag_;
  long tagCounter_;
};

#endif  // NUCLEARPOTENTIALMANAGER_H
